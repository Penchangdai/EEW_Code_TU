"""
run_magnitude.py
คำนวณขนาดแผ่นดินไหว (Magnitude) จากสถานีที่เลือก 1-5 สถานี ด้วย 3 วิธี
  1) M_JMA  (สูตร Tsuboi, 1954)        ใช้การกระจัดสูงสุดแนวราบ (หลังคลื่น P ทั้งหมด) + ระยะ Δ
  2) M_Pd   (Wu & Zhao, 2006)          ใช้การกระจัดสูงสุดแนวดิ่งในคลื่น P 3 วินาทีแรก + ระยะ R
  3) M_τc   (Wu & Kanamori, 2005)      ใช้คาบเด่นของคลื่น P 3 วินาทีแรก (ไม่ใช้ระยะทาง)
แล้วรวมหลายสถานีด้วยค่ามัธยฐาน
และรายงานค่าที่วัดได้ที่สถานี: PGA, PGV, ความรุนแรงตามเครื่องมือ JMA

ระยะทางที่ใช้ในสูตร:
  - ถ้าใส่ EPICENTER (จุดศูนย์กลางที่รู้แล้ว) จะคำนวณระยะจากพิกัดสถานีใน stations.csv
  - ถ้าไม่ใส่ จะใช้ระยะจากเวลา S-P ของแต่ละสถานี (โลก 2 ชั้น)

รัน:  ใส่รหัสสถานีใน SELECTED_STATIONS ด้านล่าง แล้วกด Run ใน VS Code
"""
import os
import sys

import numpy as np

import config as C
from data_loader import load_station_txt, resolve_data_path, find_station_file
from p_wave import detect_p_wave
from s_wave import detect_s_wave
from epicenter import delta_from_sp
from magnitude import (pd_value, magnitude_from_pd, tau_c, magnitude_from_tauc,
                       max_horizontal_displacement_um, jma_magnitude_from_amplitude,
                       network_magnitude)
from intensity import jma_instrumental_intensity, jma_intensity_class
from signal_utils import haversine_km, acc_to_vel, acc_to_disp, remove_offset
from station_db import get_coords
from run_p_wave import to_utc

# ======================================================================
#  ตั้งค่าตรงนี้
# ======================================================================
# รหัสสถานีที่ต้องการใช้ 1-5 สถานี (ต้องมีไฟล์ในโฟลเดอร์ station)
SELECTED_STATIONS = ["BKKA", "CMMT", "PHRA"]
# ถ้า SELECTED_STATIONS = [] จะใช้ไฟล์นี้ไฟล์เดียวแทน
DATA_FILE = "BKKA_20250328_061959.txt"

PROFILE = "regional"      # "regional" = สถานีไกล (> ~300 km), "local" = สถานีใกล้

# จุดศูนย์กลางที่รู้แล้ว (lat, lon, ความลึก km) เช่น จาก run_network.py หรือ USGS
# None = ใช้ระยะจาก S-P ของแต่ละสถานี
EPICENTER = None
DEPTH_KM = 10.0           # ความลึกสมมติ เมื่อ EPICENTER = None
MAX_STATIONS = 5


def analyze(path, profile, plot=True):
    rec = load_station_txt(path)
    code, fs = rec["station"], rec["fs"]
    z, n, e = rec["z"], rec["n"], rec["e"]
    regional = profile == "regional"
    print(f"\n========== สถานี {code} ==========")

    p = detect_p_wave(z, fs, start_time=rec["t0"], **(C.P_REGIONAL if regional else {}))
    if not p["detected"]:
        print("ไม่พบคลื่น P"); return None
    ip = p["index"]
    n_pre = max(ip - int(1.0 * fs), 10)          # ใช้ช่วง noise ก่อนคลื่น P ลบ offset
    print(f"P-wave arrival : {to_utc(p['time'])} UTC")

    # ---------- ระยะทาง ----------
    if EPICENTER is not None:
        coords = get_coords(code)
        if coords is None:
            print(f"ไม่พบพิกัดของ {code} ใน stations.csv"); return None
        delta = haversine_km(EPICENTER[0], EPICENTER[1], coords[0], coords[1])
        depth = EPICENTER[2]
        src = "จากพิกัดจุดศูนย์กลาง"
    else:
        s = detect_s_wave(z, n, e, fs, ip, start_time=rec["t0"], **(C.S_REGIONAL if regional else {}))
        if not s["detected"]:
            print("ไม่พบคลื่น S จึงหาระยะทางไม่ได้ (ใส่ EPICENTER แทนได้)"); return None
        delta = delta_from_sp(s["sp_time"], DEPTH_KM)
        depth = DEPTH_KM
        src = f"จาก S-P = {s['sp_time']:.2f} s"
    R = float(np.sqrt(delta ** 2 + depth ** 2))
    print(f"ระยะทาง        : Δ = {delta:.0f} km, R = {R:.0f} km  ({src})")

    # ---------- 1) M_JMA (Tsuboi) ----------
    A_um = max_horizontal_displacement_um(n, e, fs, i0=ip, n_pre=n_pre)
    m_jma = float(jma_magnitude_from_amplitude(A_um, delta))
    # ---------- 2) M_Pd ----------
    pd = pd_value(z, fs, ip, n_pre=n_pre)
    m_pd = float(magnitude_from_pd(pd, R))
    # ---------- 3) M_τc ----------
    tc = tau_c(z, fs, ip, n_pre=n_pre)
    m_tc = float(magnitude_from_tauc(tc))

    # ---------- ค่าที่วัดได้ที่สถานี ----------
    za, na, ea = (remove_offset(x, n_pre) for x in (z, n, e))
    pga = float(np.max(np.sqrt(za ** 2 + na ** 2 + ea ** 2)))
    vn, ve = acc_to_vel(n, fs, C.DISP_HIGHPASS_HZ, n_pre), acc_to_vel(e, fs, C.DISP_HIGHPASS_HZ, n_pre)
    pgv = float(np.max(np.sqrt(vn ** 2 + ve ** 2)))
    inst_i = jma_instrumental_intensity(z, n, e, fs)

    print(f"M_JMA (Tsuboi) : {m_jma:4.1f}   A = {A_um:,.0f} µm (การกระจัดสูงสุดแนวราบ)")
    print(f"M_Pd           : {m_pd:4.1f}   Pd = {pd:.4f} cm (คลื่น P {C.PD_TAUC_WINDOW_SEC:.0f} s แรก)")
    print(f"M_τc           : {m_tc:4.1f}   τc = {tc:.2f} s")
    print(f"ที่สถานี       : PGA = {pga:.2f} gal, PGV = {pgv:.3f} cm/s, "
          f"ความรุนแรง JMA = {inst_i:.1f} (ระดับ {jma_intensity_class(inst_i)})")

    res = {"code": code, "delta": delta, "R": R, "A_um": A_um, "pd": pd, "tc": tc,
           "M_jma": m_jma, "M_pd": m_pd, "M_tc": m_tc, "pga": pga, "pgv": pgv, "I": inst_i}
    if plot:
        plot_station(rec, ip, n_pre, res)
    return res


def plot_station(rec, ip, n_pre, res):
    """กราฟการกระจัด: แนวดิ่ง (ช่วง Pd/τc) และแนวราบ (ค่าที่ใช้ใน M_JMA)"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fs = rec["fs"]
    t = np.arange(len(rec["z"])) / fs
    uz = acc_to_disp(rec["z"], fs, C.DISP_HIGHPASS_HZ, n_pre) * 1e4      # cm -> µm
    un = acc_to_disp(rec["n"], fs, C.DISP_HIGHPASS_HZ, n_pre) * 1e4
    ue = acc_to_disp(rec["e"], fs, C.DISP_HIGHPASS_HZ, n_pre) * 1e4
    tp, w = ip / fs, C.PD_TAUC_WINDOW_SEC

    fig, ax = plt.subplots(3, 1, figsize=(12, 8))
    ax[0].plot(t, uz, lw=0.5, color="0.3")
    ax[0].axvspan(tp, tp + w, color="orange", alpha=0.3, label=f"P {w:.0f} s window (Pd, τc)")
    win = (t > tp - 10) & (t < tp + 30)
    lim = 1.2 * max(np.max(np.abs(uz[win])), 1e-6)
    ax[0].set_xlim(tp - 10, tp + 30); ax[0].set_ylim(-lim, lim)
    ax[0].set_ylabel("Z disp (µm)"); ax[0].legend(loc="upper right")
    ax[0].set_title(f"{res['code']}   M_JMA {res['M_jma']:.1f} | M_Pd {res['M_pd']:.1f} | "
                    f"M_τc {res['M_tc']:.1f}   (Δ = {res['delta']:.0f} km)")
    for a, u, name in ((ax[1], un, "N"), (ax[2], ue, "E")):
        a.plot(t, u, lw=0.5, color="0.3")
        a.axvline(tp, color="r", lw=0.8, label="P")
        k = ip + int(np.argmax(np.abs(u[ip:])))
        a.plot(k / fs, u[k], "bo", label=f"max |{name}| = {abs(u[k]):,.0f} µm")
        a.set_ylabel(f"{name} disp (µm)"); a.legend(loc="upper right")
    ax[2].set_xlabel("Time from start (s)")
    fig.tight_layout()
    out = f"magnitude_{res['code']}.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"บันทึกกราฟ: {out}")


def main(stations=SELECTED_STATIONS, profile=PROFILE, plot=True):
    if stations:
        paths = []
        for code in stations[:MAX_STATIONS]:
            path = find_station_file(code)
            if path is None:
                print(f"⚠ ไม่พบไฟล์ข้อมูลของสถานี {code} (ต้องมี '# Station: {code}' ในโฟลเดอร์ station)")
            else:
                paths.append(path)
    else:
        paths = [resolve_data_path(DATA_FILE)]
    print(f"Profile : {profile}" + (f"   จุดศูนย์กลาง: {EPICENTER}" if EPICENTER else ""))

    results = [r for r in (analyze(p, profile, plot) for p in paths) if r is not None]
    if not results:
        print("\nไม่มีสถานีที่คำนวณได้"); return None

    print(f"\n================ สรุปขนาดจาก {len(results)} สถานี ================")
    print(f"{'Station':<8}{'Δ (km)':>8}{'M_JMA':>7}{'M_Pd':>7}{'M_τc':>7}{'PGA(gal)':>10}{'I_JMA':>7}")
    for r in results:
        print(f"{r['code']:<8}{r['delta']:8.0f}{r['M_jma']:7.1f}{r['M_pd']:7.1f}{r['M_tc']:7.1f}"
              f"{r['pga']:10.2f}{r['I']:7.1f}")
    summary = {k: network_magnitude([r[k] for r in results]) for k in ("M_jma", "M_pd", "M_tc")}
    print(f"{'มัธยฐาน':<8}{'':8}{summary['M_jma']:7.1f}{summary['M_pd']:7.1f}{summary['M_tc']:7.1f}")

    # ---------- คำเตือนเรื่องช่วงที่สูตรใช้ได้ ----------
    far = [r["code"] for r in results if r["delta"] > 500]
    if far:
        print(f"\n⚠ {', '.join(far)} อยู่ไกลเกิน 500 km: สูตร Tsuboi ปรับเทียบไว้สำหรับ Δ ≤ 500 km")
    if any(r["R"] > 200 for r in results):
        print("⚠ M_Pd ปรับเทียบจากสถานีใกล้ (ส่วนใหญ่ R < ~100-200 km) สถานีไกลอาจคลาดมาก")
    if max(summary["M_pd"], summary["M_tc"]) >= 6.5:
        print("⚠ M_Pd และ M_τc ใช้คลื่นแค่ 3 วินาที จะอิ่มตัว (ต่ำกว่าจริง) เมื่อ M > ~7")
    return summary


if __name__ == "__main__":
    if len(sys.argv) > 1:              # python run_magnitude.py BKKA BSPA [local]
        main([a for a in sys.argv[1:] if a not in ("local", "regional")],
             "local" if "local" in sys.argv else PROFILE)
    else:
        main()
