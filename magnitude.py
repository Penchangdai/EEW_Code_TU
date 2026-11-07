"""
magnitude.py
คำนวณขนาดแผ่นดินไหว (Magnitude) ตามสูตรมาตรฐาน IASPEI (2013) ที่ USGS / NEIC ใช้
แล้วให้น้ำหนัก (weight) ความน่าเชื่อถือ และเลือกขนาดหลัก (preferred magnitude)

------------------------------------------------------------------------------------------------
ชนิดขนาด  (A = แอมพลิจูดการเคลื่อนที่ของพื้นจริง หลังจำลองเครื่องวัดมาตรฐาน)          ระยะที่ใช้ได้
------------------------------------------------------------------------------------------------
มาตรฐาน USGS / IASPEI
  ML     ML = log10(A) + 1.11 log10(R) + 0.00189 R - 2.09                           R <= 1000 km
         A (nm) บน Wood-Anderson จำลอง (T0 0.8 s, h 0.7, กำลังขยาย 1) แนวราบ, R = ระยะถึงจุดกำเนิด
         (IASPEI 2013 = Hutton & Boore 1987)
  mb_Lg  mbLg = 3.75 + 0.90 log10(Δ) + log10(A/T)     0.5° <= Δ < 4°               0.5° - 30°
         mbLg = 3.30 + 1.66 log10(Δ) + log10(A/T)     4°  <= Δ <= 30°   (Nuttli 1973)
         A (µm) = ยอดที่ 3 ของคลื่น Lg แนวดิ่งบน WWSSN-SP จำลอง ในหน้าต่าง 3.6-3.2 km/s, T 0.7-1.3 s
  Ms_BB  Ms_BB = log10(Vmax / 2π) + 1.66 log10(Δ) + 0.3                             2° - 160°
         Vmax (nm/s) = ความเร็วแนวดิ่งสูงสุดของคลื่น Rayleigh คาบ 3-60 s, ความลึก < 60 km (IASPEI 2013)
  Mwp    Mo = max|∫uz dt| · 4πρα³r / Fp ,  Mwp = (log10 Mo - 9.1) / 1.5            5° - 105°
         แก้ค่าแบบ Whitmore et al. (2002): Mwp = (Mwp - 1.03) / 0.843
         ρ = 3400 kg/m³, α = 7900 m/s, Fp = 0.52 (Tsuboi et al. 1995)
เตือนภัยล่วงหน้า / เดิม (ไม่ใช่ขนาดมาตรฐาน USGS)
  M_JMA  M = log10(A) + 1.73 log10(Δ) - 0.83   A (µm) การกระจัดสูงสุดแนวราบ (Tsuboi 1954)  Δ <= 500 km
  M_Pd   M = (log10(Pd) + 3.463 + 1.374 log10(R)) / 0.729                (Wu & Zhao 2006)
  M_tc   M = 3.373 log10(τc) + 5.787,  τc = 2π / sqrt(∫v²/∫u²)            (Wu & Kanamori 2005)
------------------------------------------------------------------------------------------------
ไม่ได้ทำ: mb (20°-100°), Ms_20 (20°-160°), Mww / Mwr (ต้อง invert moment tensor), Md (ขึ้นกับเครือข่าย)

แปลงเป็น Mw ก่อนเลือก (ดู MW_CONVERSION): ขนาดของเครือข่ายทุกชนิดถูกแปลงเป็น Mw แล้วจึงให้คะแนนและเลือก
  ผลที่แสดงเป็น Mw เสมอ (ชนิดเดิมและค่าก่อนแปลงเก็บไว้ใน "mag_raw")
  ช่วงแรกของการเตือนภัยมีแค่คลื่น P 3 s -> Mw_P = ค่าเฉลี่ย M_Pd และ M_τc  (อิ่มตัวที่ ~M 7 จึงแสดงเป็น "อย่างน้อย")

การเลือกขนาดหลัก
  1) น้ำหนักสถานี  w = SNR x ระยะอยู่ในช่วงที่สูตรปรับเทียบไว้ x คุณภาพ pick P (x ความเสถียรของ Mwp)
  2) ขนาดของเครือข่ายแต่ละชนิด = มัธยฐานถ่วงน้ำหนัก (USGS ใช้มัธยฐานของ station magnitude)
  3) คะแนนของชนิด = พื้นฐาน x ไม่อิ่มตัว x ช่วงขนาด x จำนวนสถานี x ความกระจาย x ความสอดคล้อง (x หน้าต่าง Mwp)
     พื้นฐานเรียงตามลำดับความสำคัญของ USGS: Mw (Mwp) > Ms > mb_Lg, ML > ขนาดสำหรับเตือนภัยเร็ว
  4) ขนาดหลัก = ชนิดที่คะแนนสูงสุด

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

# เครื่องวัดมาตรฐาน
WA_T0, WA_H = 0.8, 0.7                               # Wood-Anderson (IASPEI 2013, กำลังขยาย 1)
WWSSN_TS, WWSSN_HS = 1.0, 1.0                        # WWSSN-SP seismometer
WWSSN_TG, WWSSN_HG = 0.75, 1.0                       # WWSSN-SP galvanometer
MSBB_BAND_HZ = (1.0 / 60.0, 1.0 / 3.0)               # Ms_BB: ความเร็วคาบ 3-60 s
# หน้าต่างความเร็วกลุ่ม (km/s)
LG_VMAX, LG_VMIN = 3.6, 3.2
RAYLEIGH_VMAX, RAYLEIGH_VMIN = 4.5, 2.5
# Mwp
MWP_RHO, MWP_ALPHA, MWP_FP = 3400.0, 7900.0, 0.52
MWP_HIGHPASS_HZ = _cfg("MWP_HIGHPASS_HZ", None)      # None = อินทิเกรตในโดเมนเวลา ไม่กรอง
MWP_WHITMORE = _cfg("MWP_WHITMORE", True)
MWP_HP_TEST = _cfg("MWP_HP_TEST", (None, 0.005, 0.01, 0.02))   # ตัวกรองที่ใช้ตรวจความเสถียรของ Mwp

# ----------------------------------------------------------------------
#  ค่าที่ใช้ให้น้ำหนัก
# ----------------------------------------------------------------------
TYPES = ("Mwp", "Ms_BB", "mb_Lg", "ML", "M_JMA", "Mw_P", "M_Pd", "M_tc")
DISPLAY_ONLY = ("M_Pd", "M_tc")      # แสดงไว้เทียบ ไม่ใช้เลือกขนาดหลัก (รวมเป็น Mw_P แล้ว)
USGS_TYPES = ("Mwp", "Ms_BB", "mb_Lg", "ML")
#   ชนิด: (ความน่าเชื่อถือพื้นฐาน, ขนาดที่อิ่มตัว, ขนาดต่ำสุดที่สูตรเหมาะ, คำอธิบาย)
#   ขนาดที่อิ่มตัวโดยประมาณจาก Kanamori (1983) และ Bormann (NMSOP-2, บทที่ 3)
TYPE_INFO = {
    "Mwp":   (1.00, 8.0, 5.0, "Mw จากคลื่น P (Tsuboi 1995 + Whitmore 2002)"),
    "Ms_BB": (0.90, 8.0, 4.5, "คลื่นผิวพื้นคาบ 3-60 s (IASPEI 2013)"),
    "mb_Lg": (0.80, 6.5, 3.0, "คลื่น Lg คาบ ~1 s (Nuttli 1973)"),
    "ML":    (0.80, 6.5, 0.0, "Wood-Anderson คาบ 0.8 s (IASPEI 2013)"),
    "M_JMA": (0.55, 8.0, 3.0, "การกระจัดสูงสุดแนวราบ (Tsuboi 1954) ปรับเทียบในญี่ปุ่น"),
    "Mw_P":  (0.50, 7.0, 3.0, "Mw จากคลื่น P 3 s แรก = ค่าเฉลี่ย M_Pd และ M_τc (เตือนภัยล่วงหน้า)"),
    "M_Pd":  (0.40, 6.5, 3.0, "การกระจัด P 3 s แรก (Wu & Zhao 2006)"),
    "M_tc":  (0.30, 6.5, 3.0, "คาบเด่น P 3 s แรก (Wu & Kanamori 2005)"),
}
SNR_GOOD = _cfg("MAG_SNR_GOOD", 10.0)       # SNR ตั้งแต่ค่านี้ = น้ำหนักเต็ม
SAT_HALF_WIDTH = 0.5                        # เริ่มลดคะแนนที่ (ขนาดอิ่มตัว - 0.5) จนเหลือ 0.1 ที่ (+0.5)
SCATTER_REF = 0.3                           # MAD 0.3 = คะแนนความกระจายเหลือครึ่งหนึ่ง
CONSIST_REF = 0.7                           # ต่างจากค่ากลางของทุกชนิด 0.7 = คะแนนความสอดคล้องเหลือครึ่งหนึ่ง
MWP_SPREAD_REF = 0.3                        # Mwp ต่างกันตามตัวกรอง 0.3 = น้ำหนักเหลือครึ่งหนึ่ง
MWP_DRIFT_MAX = 0.2                         # baseline ลอย / แอมพลิจูด P เกินนี้ = ลดน้ำหนัก


# ======================================================================
#  เครื่องมือ: จำลองการตอบสนองของเครื่องวัดในโดเมนความถี่
# ======================================================================
def _apply_response(acc, fs, resp, n_pre=None, taper_frac=0.05):
    """
    acc  : ความเร่ง (หน่วยตาม ACC_UNIT_M_S2)
    resp : ฟังก์ชัน resp(s, f) -> H(s) ของ "ความเร่ง -> ผลลัพธ์"  (s = 2πif)
    return: สัญญาณผลลัพธ์ในหน่วย SI (m, m/s, ...)
    """
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
    H[1:] = resp(s[1:], f[1:])                        # ตัด DC
    return np.fft.irfft(np.fft.rfft(x * w, nfft) * H, nfft)[:n]


def _osc(s, T, h):
    """ตัวส่วนของ oscillator อันดับ 2:  s² + 2hω0 s + ω0²"""
    w0 = 2 * np.pi / T
    return s ** 2 + 2 * h * w0 * s + w0 ** 2


def _butter_mag(f, fc, order, kind):
    """ขนาดของตัวกรอง Butterworth (zero phase)"""
    r = (fc / f) if kind == "high" else (f / fc)
    return 1.0 / np.sqrt(1.0 + r ** (2 * order))


def wood_anderson(acc, fs, n_pre=None):
    """การกระจัดบน Wood-Anderson จำลอง กำลังขยาย 1 (m)"""
    return _apply_response(acc, fs, lambda s, f: 1.0 / _osc(s, WA_T0, WA_H), n_pre)


def _wwssn_sp_disp_resp(s):
    """WWSSN-SP (การกระจัด -> เอาต์พุต) ปรับให้ |H| = 1 ที่ 1 Hz"""
    raw = lambda ss: ss ** 3 / (_osc(ss, WWSSN_TS, WWSSN_HS) * _osc(ss, WWSSN_TG, WWSSN_HG))
    return raw(s) / abs(raw(2j * np.pi * 1.0))


def wwssn_sp(acc, fs, n_pre=None):
    """เอาต์พุต WWSSN-SP จำลอง (หน่วยเท่าการกระจัด m ที่ 1 Hz)"""
    return _apply_response(acc, fs, lambda s, f: _wwssn_sp_disp_resp(s) / s ** 2, n_pre)


def wwssn_sp_gain(period):
    """|H| ของ WWSSN-SP ที่คาบ T (แปลงแอมพลิจูดบนเครื่องกลับเป็นการกระจัดของพื้น)"""
    return float(abs(_wwssn_sp_disp_resp(np.array([2j * np.pi / period]))[0]))


def broadband_velocity(acc, fs, band=MSBB_BAND_HZ, n_pre=None, order=2):
    """ความเร็ว (m/s) ราบในช่วง band"""
    return _apply_response(
        acc, fs,
        lambda s, f: (1.0 / s) * _butter_mag(f, band[0], order, "high") * _butter_mag(f, band[1], order, "low"),
        n_pre)


def displacement(acc, fs, highpass=MWP_HIGHPASS_HZ, n_pre=None, order=4):
    """
    การกระจัด (m) จากความเร่ง
      highpass = None : ลบค่าเฉลี่ยช่วงก่อนเหตุการณ์ แล้วอินทิเกรต 2 ครั้งในโดเมนเวลา
                        (เก็บพลังงานคาบยาวไว้ครบ แต่ baseline ลอยง่าย โดยเฉพาะเครื่องวัดความเร่ง)
      highpass = fc   : high-pass Butterworth แบบ zero-phase แล้วอินทิเกรตในโดเมนความถี่
    """
    if highpass is None:
        x = np.asarray(acc, dtype=float) * ACC_UNIT_M_S2
        x = x - (x[:n_pre].mean() if n_pre else x[:max(int(fs), 1)].mean())
        v = np.cumsum(x) / fs
        return np.cumsum(v) / fs
    return _apply_response(acc, fs, lambda s, f: _butter_mag(f, highpass, order, "high") / s ** 2, n_pre)


# ======================================================================
#  ฟังก์ชันช่วย
# ======================================================================
def _snr(x, p_index, k, fs, noise_sec=30.0):
    """SNR = |x[k]| / ค่าสูงสุดของ |x| ในช่วง noise ก่อนคลื่น P (None ถ้า noise สั้นกว่า 2 s)"""
    i0 = max(p_index - int(noise_sec * fs), 0)
    if p_index - i0 < int(2 * fs):
        return None
    return float(abs(x[k]) / max(np.max(np.abs(x[i0:p_index])), 1e-30))


def _period_at(x, k, fs):
    """คาบของยอดที่ index k = 2 × ระยะระหว่างจุดตัดศูนย์ที่ล้อมยอดนั้น"""
    sgn = np.sign(x)
    a = k
    while a > 0 and sgn[a - 1] == sgn[k]:
        a -= 1
    b = k
    while b < len(x) - 1 and sgn[b + 1] == sgn[k]:
        b += 1
    return 2.0 * (b - a + 1) / fs


def _nth_peak(x, nth=1):
    """index ของยอดที่ใหญ่เป็นอันดับ nth (นับยอดละครึ่งรอบ)"""
    if len(x) < 3:
        return None
    sgn = np.sign(x)
    cut = np.flatnonzero(np.diff(sgn) != 0) + 1
    bounds = np.concatenate(([0], cut, [len(x)]))
    peaks = [lo + int(np.argmax(np.abs(x[lo:hi]))) for lo, hi in zip(bounds[:-1], bounds[1:]) if hi > lo]
    peaks.sort(key=lambda k: -abs(x[k]))
    return peaks[nth - 1] if len(peaks) >= nth else (peaks[-1] if peaks else None)


def _result(mtype, mag=None, reason=None, **kw):
    d = {"type": mtype, "mag": None if mag is None else float(mag), "valid": mag is not None}
    if reason:
        d["reason"] = reason
    d.update(kw)
    return d


def _clip(i, n):
    return int(min(max(i, 0), n))


def _check_index(p_index, n, fs, mtype):
    """กัน index ติดลบ / เกินท้ายข้อมูล (ไฟล์ไม่ครอบคลุมคลื่น P)"""
    if p_index is None or p_index < 0 or p_index >= n - int(2 * fs):
        return _result(mtype, reason="ข้อมูลไม่ครอบคลุมช่วงคลื่น P")
    return None


# ======================================================================
#  1) ML  (IASPEI 2013 / Hutton & Boore 1987)
# ======================================================================
def ml_from_amplitude(A_nm, hypo_dist_km):
    """ML = log10(A) + 1.11 log10(R) + 0.00189 R - 2.09   (A nm บน WA กำลังขยาย 1, R km)"""
    R = max(hypo_dist_km, 1.0)
    return np.log10(A_nm) + 1.11 * np.log10(R) + 0.00189 * R - 2.09


def ml_magnitude(n_acc, e_acc, fs, hypo_dist_km, p_index, i1=None, n_pre=None):
    """ML ของสถานี = ค่าเฉลี่ย ML ของแกน N และ E (วัดแยกแต่ละแกน แบบ USGS / CISN)"""
    if hypo_dist_km > 1000:
        return _result("ML", reason=f"R = {hypo_dist_km:.0f} km > 1000 km")
    bad = _check_index(p_index, len(n_acc), fs, "ML")
    if bad:
        return bad
    comp = {}
    for name, x in (("N", n_acc), ("E", e_acc)):
        wa = wood_anderson(x, fs, n_pre)
        k = p_index + int(np.argmax(np.abs(wa[p_index:i1])))
        A_nm = float(abs(wa[k]) * 1e9)
        snr = _snr(wa, p_index, k, fs)
        comp[name] = {"amp_nm": A_nm, "snr": snr, "index": k,
                      "mag": float(ml_from_amplitude(A_nm, hypo_dist_km)),
                      "used": snr is None or snr >= MIN_SNR}
    used = [c for c in comp.values() if c["used"]]
    if not used:
        return _result("ML", reason=f"SNR < {MIN_SNR:g} ทั้งสองแกน", components=comp)
    snrs = [c["snr"] for c in used if c["snr"] is not None]
    return _result("ML", np.mean([c["mag"] for c in used]), components=comp,
                   amp_nm=max(c["amp_nm"] for c in comp.values()), snr=min(snrs) if snrs else None)


# ======================================================================
#  2) mb_Lg  (Nuttli 1973)
# ======================================================================
def mblg_from_amplitude(A_um, T, epi_dist_km):
    d = epi_dist_km / KM_PER_DEG
    if d < 4.0:
        return 3.75 + 0.90 * np.log10(d) + np.log10(A_um / T)
    return 3.30 + 1.66 * np.log10(d) + np.log10(A_um / T)


def mblg_magnitude(z_acc, fs, epi_dist_km, t0_rel, p_index=None, n_pre=None, nth=3):
    """t0_rel = เวลาเกิด (s นับจากต้นไฟล์)  ยอดที่ 3 ในหน้าต่างคลื่น Lg บน WWSSN-SP จำลอง"""
    d = epi_dist_km / KM_PER_DEG
    if not 0.5 <= d <= 30.0:
        return _result("mb_Lg", reason=f"Δ = {d:.1f}° อยู่นอกช่วง 0.5° - 30°")
    sp = wwssn_sp(z_acc, fs, n_pre)
    i0 = _clip(round((t0_rel + epi_dist_km / LG_VMAX) * fs), len(sp))
    i1 = _clip(round((t0_rel + epi_dist_km / LG_VMIN) * fs), len(sp))
    if i1 - i0 < int(2 * fs):
        return _result("mb_Lg", reason="ข้อมูลไม่ครอบคลุมหน้าต่างคลื่น Lg")
    seg = sp[i0:i1]
    k = _nth_peak(seg, nth)
    T = float(np.clip(_period_at(seg, k, fs), 0.7, 1.3))
    A_um = abs(seg[k]) / wwssn_sp_gain(T) * 1e6
    snr = _snr(sp, p_index, i0 + k, fs) if p_index is not None and p_index > 0 else None
    info = dict(amp_um=A_um, period=T, index=i0 + k, snr=snr)
    if snr is not None and snr < MIN_SNR:
        return _result("mb_Lg", reason=f"SNR = {snr:.1f} < {MIN_SNR:g}", **info)
    return _result("mb_Lg", mblg_from_amplitude(A_um, T, epi_dist_km), **info)


# ======================================================================
#  3) Ms_BB  (IASPEI 2013)
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
    T = _period_at(seg, k, fs)
    Vmax = abs(seg[k]) * 1e9
    snr = _snr(v, p_index, i0 + k, fs) if p_index and p_index > 0 else None
    info = dict(vmax_nm_s=Vmax, period=T, index=i0 + k, snr=snr)
    if snr is not None and snr < MIN_SNR:
        return _result("Ms_BB", reason=f"SNR = {snr:.1f} < {MIN_SNR:g}", **info)
    return _result("Ms_BB", msbb_from_velocity(Vmax, epi_dist_km), **info)


# ======================================================================
#  4) Mwp  (Tsuboi et al. 1995; Whitmore et al. 2002)
# ======================================================================
def mwp_magnitude(z_acc, fs, epi_dist_km, p_index, s_index=None, n_pre=None, whitmore=MWP_WHITMORE,
                  highpass=MWP_HIGHPASS_HZ):
    """
    1) ความเร่งแนวดิ่ง -> การกระจัด uz   2) ∫uz dt ตั้งแต่คลื่น P ถึงก่อนคลื่น S
    3) Mo = max|∫uz dt| · 4πρα³r / Fp ,  Mwp = (log10 Mo - 9.1) / 1.5
    highpass   : ตัวกรอง high-pass ของการกระจัด (None = ไม่กรอง)
    drift_ratio: การลอยของการกระจัดใน 30 s ก่อนคลื่น P / แอมพลิจูดคลื่น P  (> ~0.2 = baseline ลอย)
    """
    d = epi_dist_km / KM_PER_DEG
    if not 5.0 <= d <= 105.0:
        return _result("Mwp", reason=f"Δ = {d:.1f}° อยู่นอกช่วง 5° - 105°")
    bad = _check_index(p_index, len(z_acc), fs, "Mwp")
    if bad:
        return bad
    u = displacement(z_acc, fs, highpass=highpass, n_pre=n_pre)
    if s_index is None:                                  # ประมาณเวลา S-P จาก Pn / Sn
        s_index = p_index + int(epi_dist_km * (1 / 4.6 - 1 / 8.0) * fs)
    i1 = _clip(s_index, len(u))
    if i1 - p_index < int(5 * fs):
        return _result("Mwp", reason="ช่วงคลื่น P ก่อนคลื่น S สั้นเกินไป")
    kp = p_index + int(np.argmax(np.abs(u[p_index:i1] - u[p_index])))
    snr = _snr(u - u[p_index], p_index, kp, fs)
    if snr is not None and snr < MIN_SNR:
        return _result("Mwp", reason=f"SNR การกระจัด = {snr:.1f} < {MIN_SNR:g}", snr=snr)
    seg = u[p_index:i1] - u[p_index]
    j0 = max(p_index - int(30 * fs), 0)
    drift_ratio = float(abs(u[p_index] - u[j0]) / max(np.max(np.abs(seg)), 1e-30))
    iu = np.cumsum(seg) / fs                             # ∫uz dt (m·s)
    k = int(np.argmax(np.abs(iu)))
    Mo = abs(iu[k]) * 4 * np.pi * MWP_RHO * MWP_ALPHA ** 3 * (epi_dist_km * 1e3) / MWP_FP
    m = raw = (np.log10(Mo) - 9.1) / 1.5
    if whitmore:
        m = (m - 1.03) / 0.843
    return _result("Mwp", m, mo_nm=float(Mo), mwp_raw=float(raw), index=p_index + k,
                   window_sec=(i1 - p_index) / fs, snr=snr, drift_ratio=drift_ratio, highpass=highpass)


def mwp_stable(z_acc, fs, epi_dist_km, p_index, s_index=None, n_pre=None, filters=MWP_HP_TEST):
    """
    Mwp จากหลายตัวกรอง (filters): ค่าของสถานี = มัธยฐาน, spread = สูงสุด - ต่ำสุด
    ถ้าข้อมูลคาบยาวดีพอ ผลแทบไม่ขึ้นกับตัวกรอง (spread เล็ก)  ถ้า spread ใหญ่ = เชื่อไม่ได้
    """
    runs = [mwp_magnitude(z_acc, fs, epi_dist_km, p_index, s_index, n_pre, highpass=hp) for hp in filters]
    ok = [r for r in runs if r.get("valid")]
    if not ok:
        return runs[0] if runs else _result("Mwp", reason="คำนวณไม่ได้")
    vals = [r["mag"] for r in ok]
    out = dict(ok[0])
    nofilt = [r for r in ok if r.get("highpass") is None]
    out.update(mag=float(np.median(vals)), variants=vals, spread=float(max(vals) - min(vals)),
               drift_ratio=(nofilt[0] if nofilt else ok[0]).get("drift_ratio"))
    return out


# ======================================================================
#  5) ขนาดแบบเดิม / เตือนภัยล่วงหน้า
# ======================================================================
def jma_magnitude_from_amplitude(A_um, epi_dist_km):
    """M_JMA (Tsuboi 1954) = log10(A) + 1.73 log10(Δ) - 0.83   (A µm, Δ <= 500 km)"""
    return np.log10(A_um) + 1.73 * np.log10(max(epi_dist_km, 1.0)) - 0.83


def max_horizontal_displacement_um(n_acc, e_acc, fs, i0=0, i1=None, n_pre=None):
    dn = acc_to_disp(n_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)[i0:i1]
    de = acc_to_disp(e_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)[i0:i1]
    return float(np.sqrt(np.max(np.abs(dn)) ** 2 + np.max(np.abs(de)) ** 2) * 1e4)     # cm -> µm


def jma_magnitude(n_acc, e_acc, fs, epi_dist_km, p_index, n_pre=None):
    A = max_horizontal_displacement_um(n_acc, e_acc, fs, i0=p_index, n_pre=n_pre)
    return float(jma_magnitude_from_amplitude(A, epi_dist_km)), A


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


def legacy_magnitudes(z_acc, n_acc, e_acc, fs, p_index, epi_dist_km, depth_km=10.0, n_pre=None):
    """M_JMA, M_Pd, M_tc พร้อมตรวจช่วงระยะ"""
    bad = _check_index(p_index, len(z_acc), fs, "M_JMA")
    if bad:
        return {t: _result(t, reason=bad["reason"]) for t in ("M_JMA", "M_Pd", "M_tc")}
    R = float(np.hypot(epi_dist_km, depth_km))
    A = max_horizontal_displacement_um(n_acc, e_acc, fs, i0=p_index, n_pre=n_pre)
    ok = epi_dist_km <= 500
    pd = pd_value(z_acc, fs, p_index, n_pre=n_pre)
    tc = tau_c(z_acc, fs, p_index, n_pre=n_pre)
    m_pd, m_tc = float(magnitude_from_pd(pd, R)), float(magnitude_from_tauc(tc))
    return {
        "M_JMA": _result("M_JMA", jma_magnitude_from_amplitude(A, epi_dist_km) if ok else None,
                         None if ok else f"Δ = {epi_dist_km:.0f} km > 500 km", amp_um=A),
        "M_Pd": _result("M_Pd", m_pd, pd_cm=pd),
        "M_tc": _result("M_tc", m_tc, tau_c=tc),
        # Mw จากคลื่น P 3 s แรก: Pd (แอมพลิจูด) กับ τc (คาบเด่น) ให้ข้อมูลต่างกัน จึงใช้ค่าเฉลี่ยของทั้งสอง
        # (แนวทางของ Kanamori 2005; Wu & Kanamori 2005)  สมการทั้งสองปรับเทียบกับ Mw ของเหตุการณ์ใหญ่
        "Mw_P": _result("Mw_P", 0.5 * (m_pd + m_tc), pd_cm=pd, tau_c=tc, m_pd=m_pd, m_tc=m_tc),
    }


# ======================================================================
#  ทุกชนิดของสถานีเดียว
# ======================================================================
def station_magnitudes(z_acc, n_acc, e_acc, fs, p_index, epi_dist_km, depth_km=10.0,
                       t0_rel=None, s_index=None, n_pre=None, legacy=False, mwp_check=False):
    """
    t0_rel   : เวลาเกิด (s นับจากต้นไฟล์)  ถ้าไม่รู้ ประมาณจากเวลา P (Pn 8.0 / Pg 6.0 km/s)
    legacy   : True = คำนวณ M_JMA, M_Pd, M_tc ด้วย
    mwp_check: True = Mwp จากหลายตัวกรอง (mwp_stable) แทนตัวกรองเดียว
    return   : {"ML": {...}, "mb_Lg": {...}, "Ms_BB": {...}, "Mwp": {...}, (+ M_JMA, M_Pd, M_tc)}
    """
    hypo = float(np.hypot(epi_dist_km, depth_km))
    if t0_rel is None:
        v = 8.0 if epi_dist_km > 200 else 6.0
        t0_rel = p_index / fs - hypo / v
    mwp = (mwp_stable if mwp_check else mwp_magnitude)(z_acc, fs, epi_dist_km, p_index, s_index, n_pre=n_pre)
    out = {
        "ML": ml_magnitude(n_acc, e_acc, fs, hypo, p_index, n_pre=n_pre),
        "mb_Lg": mblg_magnitude(z_acc, fs, epi_dist_km, t0_rel, p_index, n_pre=n_pre),
        "Ms_BB": msbb_magnitude(z_acc, fs, epi_dist_km, depth_km, t0_rel, p_index, n_pre=n_pre),
        "Mwp": mwp,
    }
    if legacy:
        out.update(legacy_magnitudes(z_acc, n_acc, e_acc, fs, p_index, epi_dist_km, depth_km, n_pre))
    return out


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
    deg = delta_km / KM_PER_DEG
    if mtype == "ML":                 # Hutton & Boore (1987) ปรับเทียบด้วยข้อมูลถึง ~600 km
        return 1.0 if R_km <= 600 else 0.5
    if mtype == "Ms_BB":              # ช่วงมาตรฐาน Ms_20 คือ 20°-160°
        return 1.0 if deg >= 20 else 0.7
    if mtype in ("M_Pd", "M_tc", "Mw_P"):     # ปรับเทียบจากสถานีใกล้ (R < ~100-200 km)
        return 1.0 if R_km <= 200 else (0.3 if R_km <= 500 else 0.1)
    return 1.0


def _pick_factor(p_weight):
    """คุณภาพ pick คลื่น P (weight 0-4) -> ตัวคูณ"""
    return {0: 1.0, 1: 1.0, 2: 0.9, 3: 0.7}.get(int(p_weight if p_weight is not None else 2), 0.3)


def station_weight(mtype, mag, delta_km, depth_km=10.0, p_weight=None):
    """น้ำหนัก 0-1 ของค่าขนาดหนึ่งค่าจากสถานีหนึ่ง"""
    if mag is None or not mag.get("valid"):
        return 0.0
    R = float(np.hypot(delta_km, depth_km))
    w = _snr_factor(mag.get("snr")) * _distance_factor(mtype, delta_km, R) * _pick_factor(p_weight)
    if mtype == "Mwp":
        w /= 1 + (mag.get("spread", 0.0) / MWP_SPREAD_REF) ** 2           # ไม่เสถียรต่อตัวกรอง
        d = mag.get("drift_ratio")
        if d is not None and d > MWP_DRIFT_MAX:
            w *= MWP_DRIFT_MAX / d                                         # baseline ลอย
    return float(w)


# ======================================================================
#  ขนาดของเครือข่าย
# ======================================================================
def weighted_median(vals, wts):
    v, w = np.asarray(vals, float), np.asarray(wts, float)
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w)
    return float(v[np.searchsorted(c, 0.5 * c[-1])])


def network_magnitude(values, min_stations=1, return_stats=False, weights=None):
    """
    ขนาดของเครือข่าย = มัธยฐาน (ถ่วงน้ำหนักถ้าให้ weights) ของ station magnitude
    return_stats=True -> dict(mag, n, mad)
    """
    pairs = [(v, 1.0 if weights is None else w) for v, w in zip(values, weights or [None] * len(values))
             if v is not None and np.isfinite(v) and (weights is None or w > 0)]
    if len(pairs) < max(min_stations, 1):
        return None if not return_stats else {"mag": None, "n": len(pairs), "mad": None}
    vals, wts = [p[0] for p in pairs], [p[1] for p in pairs]
    med = weighted_median(vals, wts)
    if not return_stats:
        return med
    mad = weighted_median(np.abs(np.array(vals) - med), wts) if len(vals) > 1 else None
    return {"mag": med, "n": len(vals), "mad": mad, "w_sum": float(np.sum(wts))}


# ======================================================================
#  แปลงขนาดทุกชนิดเป็น Mw (moment magnitude) ก่อนเลือกและแสดงผล
# ======================================================================
#   Mwp, Mw_P (Pd + τc)  : เป็น Mw อยู่แล้ว (สมการปรับเทียบกับ Mw)
#   Ms_BB -> Mw : Scordilis (2006)   Mw = 0.67 Ms + 2.07  (3.0 ≤ Ms ≤ 6.1)
#                                    Mw = 0.99 Ms + 0.08  (6.2 ≤ Ms ≤ 8.2)
#   mb_Lg -> Mw : Scordilis (2006)   Mw = 0.85 mb + 1.03  (3.5 ≤ mb ≤ 6.2)   (mb_Lg ปรับให้เทียบเท่า mb ตาม Nuttli 1973)
#   ML    -> Mw : Akkar et al. (2010) Mw = 0.953 ML + 0.422 (3.9 ≤ ML ≤ 6.8)  (ปรับเทียบในตุรกี ควรปรับเทียบใหม่ในไทย)
#   M_JMA       : ไม่มีสูตรแปลงสากล ใช้ค่าเดิม (M_JMA ≈ Mw ที่ M 5-7)
MW_CONVERSION = {
    "Mwp": ((lambda m: m), (-9, 99), "เป็น Mw อยู่แล้ว (Tsuboi 1995)"),
    "Mw_P": ((lambda m: m), (-9, 99), "ปรับเทียบกับ Mw (Wu & Zhao 2006; Wu & Kanamori 2005)"),
    "M_Pd": ((lambda m: m), (-9, 99), "ปรับเทียบกับ Mw (Wu & Zhao 2006)"),
    "M_tc": ((lambda m: m), (-9, 99), "ปรับเทียบกับ Mw (Wu & Kanamori 2005)"),
    "Ms_BB": ((lambda m: 0.67 * m + 2.07 if m < 6.15 else 0.99 * m + 0.08), (3.0, 8.2),
              "Mw = 0.67 Ms + 2.07 / 0.99 Ms + 0.08 (Scordilis 2006)"),
    "mb_Lg": ((lambda m: 0.85 * m + 1.03), (3.5, 6.2), "Mw = 0.85 mb + 1.03 (Scordilis 2006)"),
    "ML": ((lambda m: 0.953 * m + 0.422), (3.9, 6.8), "Mw = 0.953 ML + 0.422 (Akkar et al. 2010)"),
    "M_JMA": ((lambda m: m), (5.0, 7.0), "ใช้ค่าเดิม (M_JMA ≈ Mw ที่ M 5-7)"),
}


# ค่าปรับเฉพาะพื้นที่ (Mw ของระบบ − Mw จริง ที่วัดจากเหตุการณ์ทดสอบ ด้วย run_validation.py) ลบออกหลังแปลงเป็น Mw
#   สูตร Scordilis / Akkar ปรับเทียบจากข้อมูลทั่วโลก/ยุโรป  เครือข่ายไทยได้สูงกว่า Mw ของ USGS อย่างเป็นระบบ
#   ค่าตั้งต้นมาจาก 2 เหตุการณ์ (Mw 5.2 18 พ.ค. 2569, mb 4.6 21 ส.ค. 2568):
#     Ms_BB +0.34 / +0.29,  mb_Lg +0.28 / +0.19,  ML +0.68 / +0.50,  M_JMA +0.42 / +0.10
#   ยังเป็นค่าจากข้อมูลน้อย  ควรคำนวณใหม่เมื่อมีเหตุการณ์มากขึ้น (run_validation.py พิมพ์ค่าที่แนะนำให้)
#   {} = ไม่ปรับ
MW_REGIONAL_CORR = getattr(C, "MW_REGIONAL_CORR", {"Ms_BB": -0.3, "mb_Lg": -0.25, "ML": -0.6, "M_JMA": -0.25})


def to_mw(mtype, mag, regional=True):
    """แปลงขนาดชนิด mtype เป็น Mw  -> (Mw, คำอธิบายสูตร, อยู่ในช่วงที่สูตรใช้ได้หรือไม่)"""
    if mag is None:
        return None, "", False
    f, (lo, hi), note = MW_CONVERSION.get(mtype, ((lambda m: m), (-9, 99), "ใช้ค่าเดิม"))
    mw = float(f(mag))
    corr = MW_REGIONAL_CORR.get(mtype, 0.0) if regional else 0.0
    if corr:
        mw += corr
        note = f"{note} {corr:+.2f} (ค่าปรับเครือข่ายไทย)"
    return mw, note, bool(lo <= mag <= hi)


def network_by_type(stations, types=TYPES):
    """
    stations: list ของ dict ที่มี "mags" (จาก station_magnitudes), "delta", "depth", "p_weight"
    คำนวณน้ำหนัก (เก็บใน station["w"][ชนิด]) แล้วรวมเป็นขนาดของเครือข่ายแต่ละชนิด
    """
    net = {}
    for t in types:
        vals, wts, used = [], [], []
        for st in stations:
            m = st["mags"].get(t)
            w = station_weight(t, m, st["delta"], st.get("depth", 10.0), st.get("p_weight"))
            st.setdefault("w", {})[t] = w
            if w > 0:
                vals.append(m["mag"]); wts.append(w); used.append(st)
        if not vals:
            reasons = sorted({st["mags"][t].get("reason", "") for st in stations if st["mags"].get(t)})
            net[t] = {"mag": None, "n": 0, "reason": "; ".join(r for r in reasons if r) or "ไม่มีสถานีที่ใช้ได้"}
            continue
        r = network_magnitude(vals, return_stats=True, weights=wts)
        mw, note, ok = to_mw(t, r["mag"])               # ค่ากลางของชนิดนั้น -> Mw
        r.update(mag_raw=r["mag"], mag=mw, mw_note=note, mw_in_range=ok)
        windows = [st["mags"][t].get("window_sec") for st in used if st["mags"][t].get("window_sec")]
        spreads = [st["mags"][t].get("spread") for st in used if st["mags"][t].get("spread") is not None]
        r.update(min_window=min(windows) if windows else None,
                 spread=float(np.median(spreads)) if spreads else None)
        net[t] = r
    return net


# ======================================================================
#  คะแนนของแต่ละชนิด และขนาดหลัก
# ======================================================================
def rupture_duration(mag):
    """เวลาที่รอยเลื่อนแตกโดยประมาณ (s):  log10 T ≈ 0.5 M - 2.0  (M 7 ~ 30 s, M 8 ~ 100 s)"""
    return 10 ** (0.5 * mag - 2.0)


def type_scores(net):
    """
    score = base x saturation x range x data x scatter x consistency (x หน้าต่างเวลาของ Mwp)
    M0 = ค่ากลางของทุกชนิด (มัธยฐานถ่วงด้วย base x data x scatter)
         ใช้ค่ากลางแทนค่าสูงสุด เพราะบางชนิดอาจสูงเกินจริงจากข้อมูลไม่ดี (เช่น Mwp ที่ baseline ลอย)
    return (M0, {ชนิด: {score, base, sat, range, data, scatter, consist, window}})
    """
    have = {t: v for t, v in net.items() if v.get("mag") is not None and t not in DISPLAY_ONLY}
    if not have:
        return None, {}
    pre = {}
    for t, v in have.items():
        data = 1 - np.exp(-v["w_sum"] / 1.5)            # 1 สถานีน้ำหนักเต็ม ~0.5, 3 สถานี ~0.86
        scat = 1.0 if v["mad"] is None else 1.0 / (1 + (v["mad"] / SCATTER_REF) ** 2)
        pre[t] = (data, scat)
    ts = list(have)
    M0 = weighted_median([have[t]["mag"] for t in ts], [TYPE_INFO[t][0] * pre[t][0] * pre[t][1] for t in ts])
    scores = {}
    for t, v in have.items():
        base, msat, mmin, _ = TYPE_INFO[t]
        data, scat = pre[t]
        sat = float(np.clip(1 - 0.9 * (M0 - (msat - SAT_HALF_WIDTH)) / (2 * SAT_HALF_WIDTH), 0.1, 1.0))
        rng = 1.0 if M0 >= mmin else 0.3
        cons = 1.0 / (1 + ((v["mag"] - M0) / CONSIST_REF) ** 2)
        win = 1.0
        if t == "Mwp" and v.get("min_window"):
            win = float(np.clip(v["min_window"] / rupture_duration(M0), 0.2, 1.0))
        scores[t] = {"score": base * sat * rng * data * scat * cons * win, "base": base, "sat": sat,
                     "range": rng, "data": data, "scatter": scat, "consist": cons, "window": win}
    return M0, scores


def choose_preferred(net, scores):
    """
    ขนาดหลัก = ชนิดที่ score สูงสุด
    combined = ค่าเฉลี่ยถ่วง score ของชนิดที่ score >= ครึ่งหนึ่งของอันดับ 1
    """
    if not scores:
        return None
    best = max(scores, key=lambda t: scores[t]["score"])
    top = [t for t in scores if scores[t]["score"] >= 0.5 * scores[best]["score"]]
    wsum = sum(scores[t]["score"] for t in top)
    combo = sum(scores[t]["score"] * net[t]["mag"] for t in top) / wsum
    spread = max(net[t]["mag"] for t in top) - min(net[t]["mag"] for t in top)
    raw = net[best].get("mag_raw", net[best]["mag"])
    # ขนาดจากคลื่นคาบสั้นที่ถึงช่วงอิ่มตัวแล้ว = ค่าต่ำสุดที่เป็นไปได้ (ขนาดจริงอาจใหญ่กว่า)
    #   ใช้กับชนิดที่อิ่มตัวเร็ว (≤ M 7: Mw_P, ML, mb_Lg) ซึ่งให้ค่าต่ำกว่าจริงเมื่ออิ่มตัว
    lower_bound = TYPE_INFO[best][1] <= 7.0 and raw >= TYPE_INFO[best][1] - 0.3
    return {"type": best, "mag": net[best]["mag"], "mag_raw": raw, "mw_note": net[best].get("mw_note", ""),
            "lower_bound": bool(lower_bound), "score": scores[best]["score"],
            "combined": combo, "combined_types": top, "spread": spread}


PREFERRED_ORDER = USGS_TYPES


def preferred_magnitude(net):
    """
    (ใช้แทนฟังก์ชันเดิมได้)  net : {"ML": 4.1, "mb_Lg": 4.3, ...} (ค่าตัวเลข) หรือผลจาก network_by_type
    ถ้าเป็นค่าตัวเลขอย่างเดียว ใช้ลำดับแบบ USGS:  M >= 5.5 -> Mwp, Ms_BB, mb_Lg, ML
                                              3.5-5.5 -> mb_Lg, ML      M < 3.5 -> ML
    return (ชนิด, ค่า) หรือ (None, None)
    """
    if net and isinstance(next(iter(net.values())), dict):
        _, sc = type_scores(net)
        p = choose_preferred(net, sc)
        return (p["type"], p["mag"]) if p else (None, None)
    have = {k: v for k, v in net.items() if v is not None and np.isfinite(v)}
    if not have:
        return None, None
    big = max(have.values())
    order = (("Mwp", "Ms_BB", "mb_Lg", "ML") if big >= 5.5 else
             ("mb_Lg", "ML", "Mwp", "Ms_BB") if big >= 3.5 else ("ML", "mb_Lg", "Mwp", "Ms_BB"))
    for t in order:
        if t in have:
            return t, float(have[t])
    t = next(iter(have))
    return t, float(have[t])


def size_class(mag):
    """คำบอกขนาดตามการแบ่งกลุ่มของ USGS (Great / Major / Strong / Moderate / Light / Minor / Micro)"""
    for lim, w in ((8.0, "ใหญ่มาก"), (7.0, "ใหญ่"), (6.0, "รุนแรง"), (5.0, "ปานกลาง"), (4.0, "เล็ก"), (3.0, "เล็กมาก")):
        if mag >= lim:
            return w
    return "จิ๋ว (ส่วนใหญ่ไม่รู้สึก)"


def reliability(pick, net, stations, gap=None, located=True):
    """
    ระดับความน่าเชื่อถือ ("สูง" / "ปานกลาง" / "ต่ำ"), เหตุผล, ช่วงคลาดเคลื่อน ±
      ระดับเริ่มจาก score ของขนาดหลัก แล้วลดลงเมื่อ สถานีน้อย / สถานีอยู่ด้านเดียว / ไม่รู้ตำแหน่ง
      ช่วงคลาด = max(0.3, MAD, ครึ่งของผลต่างระหว่างชนิด) + เผื่อกรณีที่รู้ว่ามี bias
    """
    levels = ("ต่ำ", "ปานกลาง", "สูง")
    lv = 2 if pick["score"] >= 0.5 else (1 if pick["score"] >= 0.2 else 0)
    why = []
    n = net[pick["type"]]["n"]
    if n < 3:
        lv -= 1; why.append(f"ใช้ได้เพียง {n} สถานี")
    if gap is not None and gap > 180:
        lv -= 1; why.append("สถานีอยู่ด้านเดียวของจุดเกิด")
    if not located:
        lv -= 1; why.append("ยังระบุตำแหน่งไม่ได้")
    unc = max(0.3, net[pick["type"]].get("mad") or 0.0, 0.5 * pick["spread"])
    if gap is not None and gap > 180:
        unc += 0.2                                    # directivity ทำให้แอมพลิจูดเพี้ยนได้
    if pick["type"] == "Ms_BB" and all(st["delta"] / KM_PER_DEG < 20 for st in stations):
        unc += 0.2                                    # Ms จากระยะ < 20° (นอกช่วงมาตรฐาน Ms_20)
    if not located:
        unc += 0.5                                    # ระยะจาก S-P สถานีเดียว คลาดได้มาก
    return levels[max(lv, 0)], why, unc
