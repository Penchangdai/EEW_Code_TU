"""
intensity.py
ความรุนแรงแผ่นดินไหว (Seismic Intensity)  มาตรา MMI แบบ USGS  และมาตรา JMA (0 - 7)

  1) ความรุนแรงที่วัดได้ที่สถานี
       - MMI จาก PGV / PGA ด้วย GMICE ของ Worden et al. (2012)   <- ที่ USGS ShakeMap / ShakeAlert ใช้
       - JMA instrumental intensity (FFT -> F1·F2·F3 -> a0 ที่เกินรวม 0.3 s -> I = 2 log a0 + 0.94)
  2) ความรุนแรงที่ทำนายได้ (ก่อนคลื่นมาถึง) รวม 2 วิธีที่เป็นอิสระต่อกัน
       ก) GMPE -> GMICE : PGV บนหิน (Si & Midorikawa 1999) x อัตราขยายดิน (Vs30) -> MMI (Worden 2012)
       ข) IPE          : MMI จากขนาดและระยะโดยตรง (Allen, Wald & Worden 2012 ที่ USGS PAGER ใช้) + ผลของดิน
       MMI ทำนาย = w ก) + (1-w) ข)   w = 1 ที่ระยะ <= 200 km ลดเป็น 0 ที่ 300 km (นอกช่วงของ GMPE ใช้ IPE)
       ความไม่แน่นอน = ของสมการ + ของขนาด + ความต่างระหว่าง 2 วิธี
       ถ้ารู้แนวรอยแตก (finite fault) ใช้ระยะถึงรอยแตกแทนระยะถึงจุดศูนย์กลาง (rjb_km)
  3) ระดับการแจ้งเตือนตามเกณฑ์ ShakeAlert (USGS)
       M >= 5.0 และ MMI >= IV  -> แจ้งเตือนสาธารณะ (แบบ Wireless Emergency Alert)
       M >= 4.5 และ MMI >= III -> แจ้งเตือนผ่านแอป
     (MMI IV = ค่าที่ปัดแล้วเป็น 4 คือ MMI >= 3.5,  MMI III คือ MMI >= 2.5)

ฟังก์ชันเดิม (jma_instrumental_intensity, jma_intensity_detail, jma_intensity_class, predict_intensity,
mmi_from_pgv, ...) ยังใช้ได้เหมือนเดิม
"""
import numpy as np

import config as C


def _cfg(name, default):
    return getattr(C, name, default)


# ======================================================================
#  1) JMA instrumental intensity
# ======================================================================
def _jma_filter(f):
    f = np.maximum(f, 1e-6)
    F1 = np.sqrt(1.0 / f)                                    # period-effect filter
    X = f / 10.0
    F2 = 1.0 / np.sqrt(1 + 0.694 * X**2 + 0.241 * X**4 + 0.0557 * X**6
                       + 0.009664 * X**8 + 0.00134 * X**10 + 0.000155 * X**12)  # high-cut
    F3 = np.sqrt(1.0 - np.exp(-(f / 0.5) ** 3))              # low-cut
    return F1 * F2 * F3


def jma_instrumental_intensity(z, n, e, fs):
    """ความรุนแรงตามเครื่องมือของ JMA (คืนค่าเฉพาะ I) รายละเอียดดู jma_intensity_detail"""
    return jma_intensity_detail(z, n, e, fs)["I"]


def jma_intensity_detail(z, n, e, fs):
    """
    z, n, e : ความเร่ง 3 แกน (gal)
    1) FFT แต่ละแกน แล้วคูณ filter F1*F2*F3, inverse FFT
    2) a(t) = sqrt(az^2 + an^2 + ae^2)
    3) หา a0 ที่ทำให้ เวลารวมที่ a(t) >= a0 เท่ากับ 0.3 วินาที
    4) I = 2 log10(a0) + 0.94
    return: I, a0, a (สัญญาณหลังกรอง), i_peak (index ที่ a สูงสุด)
    """
    nfft = 1 << int(np.ceil(np.log2(2 * len(z))))           # x2 กันการวนรอบของ FFT
    f = np.fft.rfftfreq(nfft, 1.0 / fs)
    W = _jma_filter(f)
    W[0] = 0.0
    comps = [np.fft.irfft(np.fft.rfft(np.asarray(x, float) - np.mean(x), nfft) * W, nfft)[:len(z)]
             for x in (z, n, e)]
    a = np.sqrt(sum(c ** 2 for c in comps))
    k = int(round(0.3 * fs))                     # จำนวน sample ใน 0.3 s
    a0 = np.sort(a)[::-1][max(k - 1, 0)]
    I = 2.0 * np.log10(max(a0, 1e-10)) + 0.94
    I = np.floor(round(I, 2) * 10) / 10.0        # ปัดเศษหลักที่ 3 แล้วตัดเหลือ 1 ตำแหน่ง (ตาม JMA)
    return {"I": float(I), "a0": float(a0), "a": a, "i_peak": int(np.argmax(a))}


def jma_intensity_class(I):
    """แปลงค่าความรุนแรงเป็นระดับ JMA (Shindo)"""
    table = [(0.5, "0"), (1.5, "1"), (2.5, "2"), (3.5, "3"), (4.5, "4"),
             (5.0, "5-"), (5.5, "5+"), (6.0, "6-"), (6.5, "6+")]
    for upper, cls in table:
        if I < upper:
            return cls
    return "7"


# ======================================================================
#  2) MMI จากการเคลื่อนที่ของพื้น  (GMICE, Worden et al. 2012 = WGRW12 ใน USGS ShakeMap)
# ======================================================================
#   log Y < t1b :  MMI = c1b + c2b log Y           (ช่วงสั่นเบามาก)
#   t1b <= log Y < t1 :  MMI = c1 + c2 log Y
#   log Y >= t1 :  MMI = c3 + c4 log Y
#   Y = PGA (cm/s² = gal) หรือ PGV (cm/s) ของแกนราบที่ใหญ่กว่า
WGRW12 = {
    "PGA": dict(c1=1.78, c2=1.55, c3=-1.60, c4=3.70, t1=1.57, c1b=1.71, c2b=2.08, t1b=0.14, sigma=0.66),
    "PGV": dict(c1=3.78, c2=1.47, c3=2.89, c4=3.16, t1=0.53, c1b=4.62, c2b=2.17, t1b=-1.21, sigma=0.63),
}


def mmi_from_gm(value, imt="PGV"):
    """MMI จาก PGA (gal) หรือ PGV (cm/s)  ตาม Worden et al. (2012)  จำกัดในช่วง 1-10"""
    c = WGRW12[imt]
    ly = np.log10(max(float(value), 1e-6))
    if ly < c["t1b"]:
        mmi = c["c1b"] + c["c2b"] * ly
    elif ly < c["t1"]:
        mmi = c["c1"] + c["c2"] * ly
    else:
        mmi = c["c3"] + c["c4"] * ly
    return float(np.clip(mmi, 1.0, 10.0))


def mmi_slope(value, imt="PGV"):
    """dMMI / dlog10(Y) ณ ค่านั้น (ใช้แปลงอัตราขยายของดินเป็น MMI ที่เพิ่ม)"""
    c = WGRW12[imt]
    ly = np.log10(max(float(value), 1e-6))
    return c["c2b"] if ly < c["t1b"] else (c["c2"] if ly < c["t1"] else c["c4"])


def mmi_from_pgv(pgv_cm_s):
    """(เดิม) MMI = 3.47 log10(PGV) + 2.35  (Wald et al., 1999) ใช้ได้เฉพาะ MMI >= V  ปัจจุบัน USGS ใช้ mmi_from_gm"""
    return float(np.clip(3.47 * np.log10(max(pgv_cm_s, 1e-6)) + 2.35, 1, 12))


def observed_mmi(pga_gal, pgv_cm_s):
    """
    MMI ที่วัดได้ที่สถานี
      ใช้ PGV เป็นหลัก (สัมพันธ์กับ MMI ดีที่สุดในช่วงที่สร้างความเสียหาย) และแสดงค่าจาก PGA ไว้เทียบ
    """
    m_v, m_a = mmi_from_gm(pgv_cm_s, "PGV"), mmi_from_gm(pga_gal, "PGA")
    return {"mmi": m_v, "mmi_pgv": m_v, "mmi_pga": m_a, "sigma": WGRW12["PGV"]["sigma"]}


# ======================================================================
#  3) MMI ที่ทำนายได้
# ======================================================================
# ---- ก) GMPE: PGV บนชั้นหิน Vs 600 m/s (Si & Midorikawa 1999) ----
SM99_MAX_MW = 8.3            # ขนาดสูงสุดที่สมการปรับเทียบไว้
SM99_MAX_X = 300.0           # ระยะสูงสุดที่สมการปรับเทียบไว้ (km)


def pgv_bedrock_si_midorikawa(Mw, X_km, depth_km, eq_type="crustal"):
    """
    log10(PGV600) = 0.58 Mw + 0.0038 D + d - 1.29 - log10(X + 0.0028*10^(0.50 Mw)) - 0.002 X
      d = 0 (crustal), -0.02 (interplate), 0.12 (intraplate)   PGV หน่วย cm/s
    """
    d = {"crustal": 0.0, "interplate": -0.02, "intraplate": 0.12}[eq_type]
    D = min(depth_km, 200.0)
    logv = (0.58 * Mw + 0.0038 * D + d - 1.29
            - np.log10(X_km + 0.0028 * 10 ** (0.50 * Mw)) - 0.002 * X_km)
    return 10 ** logv


def site_amplification(vs30=C.DEFAULT_VS30):
    """อัตราขยาย PGV เทียบกับหิน Vs 600 m/s:  log10(ARV) = 1.83 - 0.66 log10(Vs30)   (Midorikawa et al., 1994)"""
    return 10 ** (1.83 - 0.66 * np.log10(vs30))


def intensity_from_pgv(pgv_cm_s):
    """I_JMA = 2.68 + 1.72 log10(PGV)   (Midorikawa et al., 1999)"""
    return 2.68 + 1.72 * np.log10(max(pgv_cm_s, 1e-6))


def fault_distance_from_hypo(R_km, Mw):
    """
    ประมาณระยะถึงรอยเลื่อนจากระยะ hypocentral โดยถือว่ารอยเลื่อนเป็นทรงกลมรัศมี L/2
      log10(L) = 0.5 Mw - 1.85  (km)  ->  X = max(R - L/2, 3)
    """
    L = 10 ** (0.5 * Mw - 1.85)
    return max(R_km - L / 2.0, 3.0)


# ---- ข) IPE: Allen, Wald & Worden (2012) แบบระยะ hypocentral (USGS PAGER) ----
AWW12 = dict(c0=2.085, c1=1.428, c2=-1.402, c4=0.078, m1=-0.209, m2=2.042, s1=0.82, s2=0.37, s3=22.9)
IPE_REF_VS30 = _cfg("IPE_REF_VS30", 400.0)   # สภาพดิน "เฉลี่ย" ที่ IPE แทน (ข้อมูล Did You Feel It?)  ใช้แปลงผลของดิน


def ipe_allen2012(M, R_hypo_km):
    """
    MMI = c0 + c1 M + c2 ln( √(R² + rm²) ) + c4 ln(R/50) [R > 50 km],   rm = m1 + m2 exp(M - 5)
    σ   = s1 + s2 / (1 + (R/s3)²)
    """
    c = AWW12
    R = max(float(R_hypo_km), 1.0)
    rm = c["m1"] + c["m2"] * np.exp(M - 5.0)
    mmi = c["c0"] + c["c1"] * M + c["c2"] * np.log(np.sqrt(R ** 2 + rm ** 2))
    if R > 50.0:
        mmi += c["c4"] * np.log(R / 50.0)
    sigma = c["s1"] + c["s2"] / (1 + (R / c["s3"]) ** 2)
    return float(np.clip(mmi, 1.0, 10.0)), float(sigma)


def site_mmi_term(vs30, pgv_level_cm_s):
    """MMI ที่เพิ่ม/ลดเพราะดิน (เทียบกับดินเฉลี่ย IPE_REF_VS30):  ΔMMI = (dMMI/dlogPGV) x log10(ARV(Vs30)/ARV(ref))"""
    ratio = site_amplification(vs30) / site_amplification(IPE_REF_VS30)
    return float(mmi_slope(pgv_level_cm_s, "PGV") * np.log10(ratio))


# น้ำหนักของ GMPE ในค่าเฉลี่ย: 1 เมื่อระยะ <= GMPE_TAPER_KM[0], ลดลงเป็นเส้นตรงจนเป็น 0 ที่ GMPE_TAPER_KM[1]
#   Si & Midorikawa (1999) ปรับเทียบจากข้อมูลญี่ปุ่นที่ระยะไม่เกิน ~300 km จึงไม่ใช้นอกช่วงนี้ (ใช้ IPE แทน)
#   ทดสอบกับเหตุการณ์ Mw 7.7 พม่า 2025 (75 สถานี 320-1,400 km): GMPE ต่ำกว่าค่าวัด +2.0 MMI, IPE +0.3 MMI
GMPE_TAPER_KM = _cfg("GMPE_TAPER_KM", (200.0, 300.0))
IPE_FAULT_DISTANCE = _cfg("IPE_FAULT_DISTANCE", True)   # IPE ใช้ R − L/2 แทน R (เมื่อไม่มี FAULT_TRACE)


def gmpe_weight(X_km, taper=None):
    x1, x2 = taper or GMPE_TAPER_KM
    if x2 <= x1:
        return 1.0 if X_km <= x1 else 0.0
    return float(np.clip((x2 - X_km) / (x2 - x1), 0.0, 1.0))


def distance_to_trace_km(lat, lon, trace, n_sub=50):
    """
    ระยะในแนวราบ (Rjb) จากตำแหน่ง (lat, lon) ถึงแนวรอยแตก trace = [(lat, lon), ...] (km)
    ใช้เมื่อรู้แนวรอยเลื่อนที่แตก เช่น จากแบบจำลอง finite fault ของ USGS หรือรอยเลื่อนที่รู้จัก
    """
    tr = np.asarray(trace, float)
    if tr.ndim != 2 or len(tr) < 2:
        raise ValueError("trace ต้องมีอย่างน้อย 2 จุด [(lat, lon), ...]")
    f = np.linspace(0.0, 1.0, n_sub)[:, None]
    pts = np.vstack([a + f * (b - a) for a, b in zip(tr[:-1], tr[1:])])
    la1, lo1 = np.radians(lat), np.radians(lon)
    la2, lo2 = np.radians(pts[:, 0]), np.radians(pts[:, 1])
    h = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return float(np.min(2 * 6371.0 * np.arcsin(np.sqrt(h))))


def _predict_core(M, delta_km, depth_km, vs30, eq_type, rjb_km=None):
    Mg = min(M, SM99_MAX_MW)
    if rjb_km is None:                                       # จุดกำเนิดเป็นจุด: ประมาณระยะถึงรอยเลื่อนจากขนาด
        R = float(np.hypot(delta_km, depth_km))
        X = fault_distance_from_hypo(R, Mg)
        # IPE ก็ใช้ระยะถึงรอยเลื่อนโดยประมาณ (R − L/2) เหมือน GMPE: แผ่นดินไหวใหญ่รอยแตกยาว (Mw 7.9 L ≈ 130 km)
        # ทดสอบ Mw 7.7 2025: bias +0.41 -> +0.29 MMI  ส่วนแผ่นดินไหว M < 6 ต่างจากเดิมไม่ถึง 3 km
        R_ipe = X if IPE_FAULT_DISTANCE else R
    else:                                                    # รู้แนวรอยแตก: ใช้ระยะถึงรอยแตกทั้งสองสมการ
        X = R = R_ipe = float(np.hypot(rjb_km, depth_km))
    pgv_b = float(pgv_bedrock_si_midorikawa(Mg, X, depth_km, eq_type))
    pgv_s = pgv_b * float(site_amplification(vs30))
    mmi_a = mmi_from_gm(pgv_s, "PGV")
    ipe, s_ipe = ipe_allen2012(M, R_ipe)
    mmi_b = float(np.clip(ipe + site_mmi_term(vs30, pgv_s), 1.0, 10.0))
    w = gmpe_weight(X)
    return dict(R=R, X=X, pgv_bedrock=pgv_b, pgv_surface=pgv_s, mmi_gmpe=mmi_a, mmi_ipe=mmi_b,
                sigma_ipe=s_ipe, w_gmpe=w, mmi=w * mmi_a + (1.0 - w) * mmi_b)


def predict_mmi(M, delta_km, depth_km=10.0, vs30=C.DEFAULT_VS30, eq_type="crustal", m_unc=0.3, rjb_km=None):
    """
    ทำนาย MMI ที่ตำแหน่งหนึ่ง (ก่อนคลื่นมาถึง)
      M       : ขนาดหลัก (จาก magnitude.py)     m_unc : ความไม่แน่นอนของขนาด (±)
      delta_km: ระยะจากจุดเหนือศูนย์กลาง,  vs30 : m/s
      rjb_km  : (ถ้ามี) ระยะในแนวราบถึงแนวรอยแตก (distance_to_trace_km)  None = ถือเป็นจุดกำเนิด
    return: mmi (เฉลี่ยถ่วงน้ำหนัก 2 วิธี), mmi_gmpe, mmi_ipe, w_gmpe, sigma, mmi_low / mmi_high (±1σ),
            pgv_surface, jma (I ทำนาย), jma_class, flags (ข้อจำกัด)
    """
    p = _predict_core(M, delta_km, depth_km, vs30, eq_type, rjb_km)
    up = _predict_core(M + 0.1, delta_km, depth_km, vs30, eq_type, rjb_km)["mmi"]
    dm = (up - p["mmi"]) / 0.1                               # dMMI/dM
    w = p["w_gmpe"]
    s_model = w * WGRW12["PGV"]["sigma"] + (1 - w) * p["sigma_ipe"]
    s_methods = w * (1 - w) * abs(p["mmi_gmpe"] - p["mmi_ipe"])
    sigma = float(np.sqrt(s_model ** 2 + (dm * m_unc) ** 2 + s_methods ** 2))
    flags = []
    if p["X"] > SM99_MAX_X:
        flags.append(f"ระยะ {p['X']:.0f} km เกินช่วงของ GMPE (~{SM99_MAX_X:.0f} km) จึงใช้ IPE")
    if M > SM99_MAX_MW:
        flags.append(f"M {M:.1f} เกินช่วงของ Si & Midorikawa (ใช้ {SM99_MAX_MW})")
    I = float(intensity_from_pgv(p["pgv_surface"]))
    return dict(p, sigma=sigma, mmi_low=max(p["mmi"] - sigma, 1.0), mmi_high=min(p["mmi"] + sigma, 10.0),
                dmmi_dm=dm, jma=I, jma_class=jma_intensity_class(I), flags=flags)


def predict_intensity(Mw, hypo_dist_km, depth_km, vs30=C.DEFAULT_VS30, eq_type="crustal"):
    """(เดิม) ทำนายความรุนแรง JMA จากระยะ hypocentral  -> dict เหมือนเดิม + mmi แบบใหม่"""
    delta = float(np.sqrt(max(hypo_dist_km ** 2 - depth_km ** 2, 0.0)))
    p = predict_mmi(Mw, delta, depth_km, vs30, eq_type)
    return {"fault_dist_km": p["X"], "pgv_bedrock": p["pgv_bedrock"], "pgv_surface": p["pgv_surface"],
            "arv": float(site_amplification(vs30)), "intensity": p["jma"], "jma_class": p["jma_class"],
            "mmi": p["mmi"]}


# ======================================================================
#  4) คำอธิบาย MMI และระดับการแจ้งเตือน (USGS)
# ======================================================================
#   (เลขโรมัน, การสั่น, ความเสียหาย)  ตามตารางของ USGS ShakeMap
MMI_TABLE = {
    1: ("I", "ไม่รู้สึก", "ไม่มี"),
    2: ("II", "อ่อน", "ไม่มี"),
    3: ("III", "อ่อน", "ไม่มี"),
    4: ("IV", "เบา", "ไม่มี"),
    5: ("V", "ปานกลาง", "น้อยมาก"),
    6: ("VI", "แรง", "เล็กน้อย"),
    7: ("VII", "แรงมาก", "ปานกลาง"),
    8: ("VIII", "รุนแรง", "ปานกลางถึงหนัก"),
    9: ("IX", "รุนแรงมาก", "หนัก"),
    10: ("X+", "รุนแรงที่สุด", "หนักมาก"),
}
MMI_EN = {1: "Not felt", 2: "Weak", 3: "Weak", 4: "Light", 5: "Moderate", 6: "Strong", 7: "Very strong",
          8: "Severe", 9: "Violent", 10: "Extreme"}


def mmi_class(mmi):
    """MMI (ทศนิยม) -> (ระดับจำนวนเต็ม, เลขโรมัน, การสั่น, ความเสียหาย)   ระดับ = ค่าที่ปัดแล้ว"""
    k = int(np.clip(np.floor(mmi + 0.5), 1, 10))
    return (k,) + MMI_TABLE[k]


# เกณฑ์ ShakeAlert (USGS):  (ชื่อ, M ขั้นต่ำ, MMI ขั้นต่ำ)
ALERT_LEVELS = (
    ("แจ้งเตือนสาธารณะ", _cfg("ALERT_PUBLIC_M", 5.0), _cfg("ALERT_PUBLIC_MMI", 3.5)),   # WEA: M5+, MMI IV
    ("แจ้งเตือนผ่านแอป", _cfg("ALERT_APP_M", 4.5), _cfg("ALERT_APP_MMI", 2.5)),          # แอป: M4.5+, MMI III
)


def alert_level(M, mmi):
    """ระดับการแจ้งเตือน (0 = ไม่แจ้ง, 1 = แอป, 2 = สาธารณะ) และชื่อ"""
    for k, (name, m_min, mmi_min) in enumerate(ALERT_LEVELS):
        if M >= m_min and mmi >= mmi_min:
            return 2 - k, name
    return 0, "ไม่แจ้งเตือน"


def action_text(mmi):
    """ข้อความสิ่งที่ควรทำ"""
    k = mmi_class(mmi)[0]
    if k >= 6:
        return "หมอบลง ป้องศีรษะและคอ เกาะสิ่งที่มั่นคง  อาจมีสิ่งของตกหล่นและอาคารเสียหาย"
    if k >= 4:
        return "หมอบลง ป้องศีรษะ เกาะสิ่งที่มั่นคง  อยู่ห่างจากกระจกและของที่อาจล้ม"
    return "อาจรู้สึกสั่นเล็กน้อย โดยเฉพาะชั้นสูงของอาคาร ไม่ต้องตื่นตระหนก"
