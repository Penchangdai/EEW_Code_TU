"""
synthetic.py
สร้างคลื่นแผ่นดินไหวจำลอง 3 แกน (Z, N, E) หน่วย gal สำหรับทดสอบระบบ
(เมื่อมีข้อมูลจริง เช่นไฟล์ miniSEED ให้ใช้ ObsPy อ่านแทนไฟล์นี้)
"""
import numpy as np

import config as C
from signal_utils import haversine_km, azimuth_deg


def _wavelet(t, t_on, amp, freq, rise=0.1, decay=3.0):
    """การกระจัดรูป sine ที่มี envelope ขึ้นเร็ว-ลดลงช้า เริ่มที่เวลา t_on"""
    tt = t - t_on
    out = np.zeros_like(t)
    m = tt > 0
    env = (1 - np.exp(-tt[m] / rise)) * np.exp(-tt[m] / decay)
    out[m] = amp * env * np.sin(2 * np.pi * freq * tt[m])
    return out / max(np.max(np.abs(out)), 1e-30) * amp


def make_station_record(sta_lat, sta_lon, eq, fs=100.0, duration=120.0, pre=20.0,
                        noise_gal=0.01, seed=None):
    """
    eq : dict {lat, lon, depth, origin_time, M}
    คืนค่า dict {t0 (เวลา sample แรก), z, n, e (gal), tp_true, ts_true}
    """
    rng = np.random.default_rng(seed)
    t_start = eq["origin_time"] - pre
    t = t_start + np.arange(int(duration * fs)) / fs

    delta = haversine_km(eq["lat"], eq["lon"], sta_lat, sta_lon)
    R = np.sqrt(delta ** 2 + eq["depth"] ** 2)
    tp = eq["origin_time"] + R / C.VP_KM_S
    ts = eq["origin_time"] + R / C.VS_KM_S
    az = np.radians(azimuth_deg(eq["lat"], eq["lon"], sta_lat, sta_lon))   # ทิศ source->station
    M = eq["M"]

    # P: แอมพลิจูดตาม Pd (Wu & Zhao) และความถี่ตาม τc (Wu & Kanamori)
    pd = 10 ** (-3.463 + 0.729 * M - 1.374 * np.log10(R))                # cm
    fp = 1.0 / 10 ** ((M - 5.787) / 3.373)
    inc = np.arctan2(delta, eq["depth"]) * 0.4                           # มุมตกกระทบ (ประมาณ)
    dP = _wavelet(t, tp, pd, fp, decay=4.0)
    # S: แอมพลิจูดตาม Tsuboi (JMA) แนวราบ, ขั้วตั้งฉากกับแนวรัศมี (SH)
    A_s = 10 ** (M - 1.73 * np.log10(max(delta, 1.0)) + 0.83) * 1e-4     # µm -> cm
    dS = _wavelet(t, ts, A_s, max(fp * 1.5, 0.5), decay=2.0 + 2.0 * max(M - 5, 0))

    radial_n, radial_e = np.cos(az), np.sin(az)
    trans_n, trans_e = -np.sin(az), np.cos(az)
    d_z = dP * np.cos(inc) + 0.15 * dS
    d_n = dP * np.sin(inc) * radial_n + dS * (0.9 * trans_n + 0.3 * radial_n)
    d_e = dP * np.sin(inc) * radial_e + dS * (0.9 * trans_e + 0.3 * radial_e)

    def to_acc(d):
        return np.gradient(np.gradient(d, 1 / fs), 1 / fs) + rng.normal(0, noise_gal, len(d))

    return {"t0": t_start, "fs": fs, "z": to_acc(d_z), "n": to_acc(d_n), "e": to_acc(d_e),
            "tp_true": tp, "ts_true": ts, "epi_dist_true": delta}
