"""
magnitude.py
คำนวณขนาดแผ่นดินไหว 2 ชนิด แล้วรวมเป็นขนาดของเครือข่าย

  Mw_P   ขนาดช่วงแรก จากคลื่น P 3 วินาทีแรก (ได้เร็ว แต่อิ่มตัวที่ราว M 6-7)
         M_Pd = (log10(Pd) + 3.463 + 1.374 log10(R)) / 0.729         (Wu & Zhao 2006)
         M_τc = 3.373 log10(τc) + 5.787,  τc = 2π / sqrt(∫v²/∫u²)    (Wu & Kanamori 2005)
         Mw_P = (M_Pd + M_τc) / 2
  Ms_BB  ขนาดสุดท้าย จากคลื่นผิวพื้น (มาช้ากว่า แต่แม่นกว่าเมื่อแผ่นดินไหวใหญ่)
         Ms_BB = log10(Vmax / 2π) + 1.66 log10(Δ) + 0.3   2° - 160°  (IASPEI 2013)
         Vmax (nm/s) = ความเร็วแนวดิ่งสูงสุด คาบ 3-60 s ในหน้าต่างคลื่น Rayleigh 4.5-2.5 km/s
         แปลงเป็น Mw: Mw = 0.67 Ms + 2.07 (Ms < 6.15), Mw = 0.99 Ms + 0.08 (Scordilis 2006) แล้วบวกค่าปรับเครือข่ายไทย

การเลือกขนาดหลัก:  มี Ms_BB ใช้ Ms_BB  ถ้ายังไม่มีใช้ Mw_P
ขนาดของเครือข่าย = มัธยฐานถ่วงน้ำหนักของทุกสถานี   น้ำหนัก = SNR x ระยะอยู่ในช่วงของสูตร x คุณภาพ pick P

หน่วยข้อมูลเข้า: ความเร่ง (gal = cm/s²) เปลี่ยนได้ด้วย config.ACC_UNIT_M_S2
"""
import numpy as np

import config as C
from signal_utils import acc_to_disp, acc_to_vel


def _cfg(name, default):
    return getattr(C, name, default)


ACC_UNIT_M_S2 = _cfg("ACC_UNIT_M_S2", 0.01)         # 1 gal = 0.01 m/s²
KM_PER_DEG = 111.195
MIN_SNR = _cfg("MAG_MIN_SNR", 3.0)                  # SNR ขั้นต่ำ (สัญญาณ / noise ก่อนคลื่น P)
SNR_GOOD = _cfg("MAG_SNR_GOOD", 10.0)               # SNR ตั้งแต่ค่านี้ = น้ำหนักเต็ม

MSBB_BAND_HZ = (1.0 / 60.0, 1.0 / 3.0)              # Ms_BB: ความเร็วคาบ 3-60 s
RAYLEIGH_VMAX, RAYLEIGH_VMIN = 4.5, 2.5             # หน้าต่างความเร็วกลุ่มของคลื่น Rayleigh (km/s)
MW_P_SATURATION = 7.0                               # Mw_P อิ่มตัวที่ราวค่านี้ จึงแสดงเป็น "ขึ้นไป" เมื่อเข้าใกล้

TYPES = ("Ms_BB", "Mw_P")                           # เรียงตามลำดับที่เลือกเป็นขนาดหลัก
TYPE_NOTE = {
    "Ms_BB": "คลื่นผิวพื้นคาบ 3-60 s (IASPEI 2013)",
    "Mw_P": "คลื่น P 3 วินาทีแรก = ค่าเฉลี่ย M_Pd และ M_τc (Wu & Zhao 2006; Wu & Kanamori 2005)",
}


# ======================================================================
#  ฟังก์ชันช่วย
# ======================================================================
def _butter_mag(f, fc, order, kind):
    """ขนาดของตัวกรอง Butterworth (zero phase)"""
    r = (fc / f) if kind == "high" else (f / fc)
    return 1.0 / np.sqrt(1.0 + r ** (2 * order))


def broadband_velocity(acc, fs, band=MSBB_BAND_HZ, n_pre=None, order=2, taper_frac=0.05):
    """ความเร็ว (m/s) ในช่วง band จากความเร่ง  ทำในโดเมนความถี่"""
    x = np.asarray(acc, dtype=float) * ACC_UNIT_M_S2
    x = x - (x[:n_pre].mean() if n_pre else x.mean())
    n = len(x)
    m = max(int(taper_frac * n), 1)
    w = np.ones(n)
    w[:m] = 0.5 * (1 - np.cos(np.pi * np.arange(m) / m))
    w[n - m:] = w[:m][::-1]
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfftfreq(nfft, 1.0 / fs)
    s = 2j * np.pi * f
    H = np.zeros_like(s)
    H[1:] = (1.0 / s[1:]) * _butter_mag(f[1:], band[0], order, "high") * _butter_mag(f[1:], band[1], order, "low")
    return np.fft.irfft(np.fft.rfft(x * w, nfft) * H, nfft)[:n]


def _snr(x, p_index, k, fs, noise_sec=30.0):
    """SNR = |x[k]| / ค่าสูงสุดของ |x| ในช่วง noise ก่อนคลื่น P (None ถ้า noise สั้นกว่า 2 s)"""
    i0 = max(p_index - int(noise_sec * fs), 0)
    if p_index - i0 < int(2 * fs):
        return None
    return float(abs(x[k]) / max(np.max(np.abs(x[i0:p_index])), 1e-30))


def _result(mtype, mag=None, reason=None, **kw):
    d = {"type": mtype, "mag": None if mag is None else float(mag), "valid": mag is not None}
    if reason:
        d["reason"] = reason
    d.update(kw)
    return d


def _clip(i, n):
    return int(min(max(i, 0), n))


# ======================================================================
#  1) Ms_BB  (IASPEI 2013)
# ======================================================================
def msbb_from_velocity(Vmax_nm_s, epi_dist_km):
    d = epi_dist_km / KM_PER_DEG
    return np.log10(Vmax_nm_s / (2 * np.pi)) + 1.66 * np.log10(d) + 0.3


def msbb_magnitude(z_acc, fs, epi_dist_km, depth_km, t0_rel=None, p_index=0, n_pre=None):
    """Vmax = ความเร็วสูงสุดแนวดิ่ง (คาบ 3-60 s) ในหน้าต่างคลื่น Rayleigh (4.5 - 2.5 km/s)"""
    d = epi_dist_km / KM_PER_DEG
    if not 2.0 <= d <= 160.0:
        return _result("Ms_BB", reason=f"Δ = {d:.1f}° อยู่นอกช่วง 2° - 160°")
    if depth_km >= 60:
        return _result("Ms_BB", reason=f"ความลึก {depth_km:.0f} km >= 60 km")
    v = broadband_velocity(z_acc, fs, n_pre=n_pre)
    if t0_rel is not None:
        i0 = _clip(round((t0_rel + epi_dist_km / RAYLEIGH_VMAX) * fs), len(v))
        i1 = _clip(round((t0_rel + epi_dist_km / RAYLEIGH_VMIN) * fs), len(v))
    else:
        i0, i1 = max(p_index, 0), len(v)
    if i1 - i0 < int(10 * fs):
        return _result("Ms_BB", reason="ข้อมูลไม่ครอบคลุมหน้าต่างคลื่นผิวพื้น")
    seg = v[i0:i1]
    k = int(np.argmax(np.abs(seg)))
    Vmax = abs(seg[k]) * 1e9
    snr = _snr(v, p_index, i0 + k, fs) if p_index and p_index > 0 else None
    info = dict(vmax_nm_s=Vmax, index=i0 + k, snr=snr)
    if snr is not None and snr < MIN_SNR:
        return _result("Ms_BB", reason=f"SNR = {snr:.1f} < {MIN_SNR:g}", **info)
    return _result("Ms_BB", msbb_from_velocity(Vmax, epi_dist_km), **info)


# ======================================================================
#  2) Mw_P  (คลื่น P 3 วินาทีแรก)
# ======================================================================
def pd_value(z_acc, fs, p_index, window=C.PD_TAUC_WINDOW_SEC, n_pre=None):
    """Pd = ค่าสูงสุดของการกระจัดแนวดิ่งใน window วินาทีแรกของคลื่น P (cm)"""
    d = acc_to_disp(z_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)
    return float(np.max(np.abs(d[p_index:p_index + int(window * fs)])))


def magnitude_from_pd(pd_cm, hypo_dist_km):
    """(Wu & Zhao 2006)  M = (log10(Pd) + 3.463 + 1.374 log10(R)) / 0.729"""
    return (np.log10(pd_cm) + 3.463 + 1.374 * np.log10(max(hypo_dist_km, 1.0))) / 0.729


def tau_c(z_acc, fs, p_index, window=C.PD_TAUC_WINDOW_SEC, n_pre=None):
    """τc = 2π / sqrt(∫v² dt / ∫u² dt)   (Kanamori 2005)"""
    v = acc_to_vel(z_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)
    u = acc_to_disp(z_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)
    sl = slice(p_index, p_index + int(window * fs))
    return float(2 * np.pi / np.sqrt(np.sum(v[sl] ** 2) / max(np.sum(u[sl] ** 2), 1e-30)))


def magnitude_from_tauc(tc):
    """(Wu & Kanamori 2005)  M = 3.373 log10(τc) + 5.787"""
    return 3.373 * np.log10(tc) + 5.787


def mw_p_magnitude(z_acc, fs, p_index, epi_dist_km, depth_km=10.0, n_pre=None):
    """Mw_P = (M_Pd + M_τc) / 2   Pd บอกแอมพลิจูด τc บอกคาบเด่น จึงให้ข้อมูลต่างกัน"""
    if p_index is None or p_index < 0 or p_index >= len(z_acc) - int(2 * fs):
        return _result("Mw_P", reason="ข้อมูลไม่ครอบคลุมช่วงคลื่น P")
    R = float(np.hypot(epi_dist_km, depth_km))
    pd = pd_value(z_acc, fs, p_index, n_pre=n_pre)
    tc = tau_c(z_acc, fs, p_index, n_pre=n_pre)
    m_pd, m_tc = float(magnitude_from_pd(pd, R)), float(magnitude_from_tauc(tc))
    return _result("Mw_P", 0.5 * (m_pd + m_tc), pd_cm=pd, tau_c=tc, m_pd=m_pd, m_tc=m_tc)


def station_magnitudes(z_acc, fs, p_index, epi_dist_km, depth_km=10.0, t0_rel=None, n_pre=None):
    """t0_rel : เวลาเกิด (s นับจากต้นไฟล์)   return : {"Ms_BB": {...}, "Mw_P": {...}}"""
    return {
        "Ms_BB": msbb_magnitude(z_acc, fs, epi_dist_km, depth_km, t0_rel, p_index, n_pre=n_pre),
        "Mw_P": mw_p_magnitude(z_acc, fs, p_index, epi_dist_km, depth_km, n_pre=n_pre),
    }


# ======================================================================
#  น้ำหนักระดับสถานี
# ======================================================================
def _snr_factor(snr):
    """SNR: ไม่ทราบ -> 0.7, MIN_SNR -> 0.5, ตั้งแต่ SNR_GOOD -> 1.0"""
    if snr is None:
        return 0.7
    if snr < MIN_SNR:
        return 0.0
    return float(np.clip(0.5 + 0.5 * (snr - MIN_SNR) / (SNR_GOOD - MIN_SNR), 0.5, 1.0))


def _distance_factor(mtype, delta_km, R_km):
    """1 = ระยะอยู่ในช่วงที่สูตรปรับเทียบไว้ดี, น้อยกว่า 1 = ใช้ได้แต่ไม่ดีนัก"""
    if mtype == "Ms_BB":              # ช่วงมาตรฐาน Ms_20 คือ 20°-160°
        return 1.0 if delta_km / KM_PER_DEG >= 20 else 0.7
    return 1.0 if R_km <= 200 else (0.3 if R_km <= 500 else 0.1)      # Mw_P ปรับเทียบจากสถานีใกล้


def _pick_factor(p_weight):
    """คุณภาพ pick คลื่น P (weight 0-4) -> ตัวคูณ"""
    return {0: 1.0, 1: 1.0, 2: 0.9, 3: 0.7}.get(int(p_weight if p_weight is not None else 2), 0.3)


def station_weight(mtype, mag, delta_km, depth_km=10.0, p_weight=None):
    """น้ำหนัก 0-1 ของค่าขนาดหนึ่งค่าจากสถานีหนึ่ง"""
    if mag is None or not mag.get("valid"):
        return 0.0
    R = float(np.hypot(delta_km, depth_km))
    return float(_snr_factor(mag.get("snr")) * _distance_factor(mtype, delta_km, R) * _pick_factor(p_weight))


# ======================================================================
#  ขนาดของเครือข่าย
# ======================================================================
def weighted_median(vals, wts):
    v, w = np.asarray(vals, float), np.asarray(wts, float)
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w)
    return float(v[np.searchsorted(c, 0.5 * c[-1])])


def network_magnitude(values, weights=None):
    """ขนาดของเครือข่าย = มัธยฐาน (ถ่วงน้ำหนักถ้าให้ weights) ของ station magnitude  -> dict(mag, n, mad)"""
    pairs = [(v, 1.0 if weights is None else w) for v, w in zip(values, weights or [None] * len(values))
             if v is not None and np.isfinite(v) and (weights is None or w > 0)]
    if not pairs:
        return {"mag": None, "n": 0, "mad": None}
    vals, wts = [p[0] for p in pairs], [p[1] for p in pairs]
    med = weighted_median(vals, wts)
    mad = weighted_median(np.abs(np.array(vals) - med), wts) if len(vals) > 1 else None
    return {"mag": med, "n": len(vals), "mad": mad}


# ======================================================================
#  แปลงเป็น Mw
# ======================================================================
#   ค่าปรับเฉพาะพื้นที่ (ลบความเอนเอียงของเครือข่ายไทยเทียบกับ Mw ของ USGS) วัดจากเหตุการณ์ทดสอบด้วย run_validation.py
#   ยังเป็นค่าจากข้อมูลน้อย ควรคำนวณใหม่เมื่อมีเหตุการณ์มากขึ้น   {} = ไม่ปรับ
MW_REGIONAL_CORR = getattr(C, "MW_REGIONAL_CORR", {"Ms_BB": -0.3})


def to_mw(mtype, mag, regional=True):
    """แปลงขนาดชนิด mtype เป็น Mw  -> (Mw, คำอธิบายสูตร, อยู่ในช่วงที่สูตรใช้ได้หรือไม่)"""
    if mag is None:
        return None, "", False
    if mtype != "Ms_BB":
        return float(mag), "ปรับเทียบกับ Mw (Wu & Zhao 2006; Wu & Kanamori 2005)", True
    mw = float(0.67 * mag + 2.07 if mag < 6.15 else 0.99 * mag + 0.08)      # Scordilis (2006)
    note = "Mw = 0.67 Ms + 2.07 / 0.99 Ms + 0.08 (Scordilis 2006)"
    corr = MW_REGIONAL_CORR.get(mtype, 0.0) if regional else 0.0
    if corr:
        mw += corr
        note = f"{note} {corr:+.2f} (ค่าปรับเครือข่ายไทย)"
    return mw, note, bool(3.0 <= mag <= 8.2)


def network_by_type(stations, types=TYPES):
    """
    stations: list ของ dict ที่มี "mags" (จาก station_magnitudes), "delta", "depth", "p_weight"
    คำนวณน้ำหนัก (เก็บใน station["w"][ชนิด]) แล้วรวมเป็นขนาดของเครือข่ายแต่ละชนิด
    """
    net = {}
    for t in types:
        vals, wts = [], []
        for st in stations:
            m = st["mags"].get(t)
            w = station_weight(t, m, st["delta"], st.get("depth", 10.0), st.get("p_weight"))
            st.setdefault("w", {})[t] = w
            if w > 0:
                vals.append(m["mag"]); wts.append(w)
        if not vals:
            reasons = sorted({st["mags"][t].get("reason", "") for st in stations if st["mags"].get(t)})
            net[t] = {"mag": None, "n": 0, "reason": "; ".join(r for r in reasons if r) or "ไม่มีสถานีที่ใช้ได้"}
            continue
        r = network_magnitude(vals, weights=wts)
        mw, note, ok = to_mw(t, r["mag"])
        r.update(mag_raw=r["mag"], mag=mw, mw_note=note, mw_in_range=ok)
        net[t] = r
    return net


# ======================================================================
#  ขนาดหลัก
# ======================================================================
def choose_preferred(net):
    """มี Ms_BB ใช้ Ms_BB  ถ้ายังไม่มีใช้ Mw_P   (Mw_P ที่เข้าใกล้ช่วงอิ่มตัว = ค่าต่ำสุด ขนาดจริงอาจใหญ่กว่า)"""
    for t in TYPES:
        v = net.get(t)
        if v and v.get("mag") is not None:
            lower_bound = t == "Mw_P" and v["mag_raw"] >= MW_P_SATURATION - 0.3
            return {"type": t, "mag": v["mag"], "mag_raw": v["mag_raw"], "mw_note": v.get("mw_note", ""),
                    "lower_bound": bool(lower_bound)}
    return None


def size_class(mag):
    """คำบอกขนาดตามการแบ่งกลุ่มของ USGS (Great / Major / Strong / Moderate / Light / Minor / Micro)"""
    for lim, w in ((8.0, "ใหญ่มาก"), (7.0, "ใหญ่"), (6.0, "รุนแรง"), (5.0, "ปานกลาง"), (4.0, "เล็ก"), (3.0, "เล็กมาก")):
        if mag >= lim:
            return w
    return "จิ๋ว (ส่วนใหญ่ไม่รู้สึก)"


def reliability(pick, net, stations, gap=None, located=True):
    """
    ระดับความน่าเชื่อถือ ("สูง" / "ปานกลาง" / "ต่ำ"), เหตุผล, ช่วงคลาดเคลื่อน ±
      เริ่มจาก "สูง" แล้วลดลงทีละขั้นเมื่อพบข้อจำกัด   ช่วงคลาด = max(0.3, MAD) + เผื่อกรณีที่รู้ว่ามี bias
    """
    levels = ("ต่ำ", "ปานกลาง", "สูง")
    lv = 2
    why = []
    n = net[pick["type"]]["n"]
    if pick["type"] == "Mw_P":
        lv -= 1; why.append("เป็นค่าเบื้องต้นจากคลื่น P 3 วินาทีแรก")
    if n < 3:
        lv -= 1; why.append(f"ใช้ได้เพียง {n} สถานี")
    if gap is not None and gap > 180:
        lv -= 1; why.append("สถานีอยู่ด้านเดียวของจุดเกิด")
    if not located:
        lv -= 1; why.append("ยังระบุตำแหน่งไม่ได้")
    unc = max(0.3, net[pick["type"]].get("mad") or 0.0)
    if gap is not None and gap > 180:
        unc += 0.2                                    # directivity ทำให้แอมพลิจูดเพี้ยนได้
    if pick["type"] == "Ms_BB" and all(st["delta"] / KM_PER_DEG < 20 for st in stations):
        lv -= 1; why.append("ใช้ Ms จากระยะน้อยกว่า 20°")
        unc += 0.2                                    # นอกช่วงมาตรฐาน Ms_20
    return levels[max(lv, 0)], why, unc
