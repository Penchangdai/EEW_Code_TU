"""
intensity.py
ความรุนแรงแผ่นดินไหว (Seismic Intensity)  มาตรา MMI แบบ USGS

  1) MMI ที่วัดได้ที่สถานี : จาก PGV ด้วย GMICE ของ Worden et al. (2012)
  2) MMI ที่ทำนายได้ (ก่อนคลื่นมาถึง)
       IPE ของ Allen, Wald & Worden (2012) : MMI จากขนาดและระยะโดยตรง
       + ผลของดินจาก Vs30 (Midorikawa et al., 1994)
       ถ้ารู้แนวรอยแตก (finite fault) ใช้ระยะถึงรอยแตกแทนระยะถึงจุดกำเนิด (rjb_km)
  3) ระดับการแจ้งเตือนตามเกณฑ์ ShakeAlert (USGS)
       M >= 5.0 และ MMI >= IV  -> แจ้งเตือนสาธารณะ (แบบ Wireless Emergency Alert)
       M >= 4.5 และ MMI >= III -> แจ้งเตือนผ่านแอป
     (MMI IV = ค่าที่ปัดแล้วเป็น 4 คือ MMI >= 3.5,  MMI III คือ MMI >= 2.5)
"""
import numpy as np

import config as C


def _cfg(name, default):
    return getattr(C, name, default)


# ======================================================================
#  1) MMI จาก PGV  (GMICE, Worden et al. 2012)
# ======================================================================
#   log PGV < -1.21        :  MMI = 4.62 + 2.17 log PGV        (ช่วงสั่นเบามาก)
#   -1.21 <= log PGV < 0.53 :  MMI = 3.78 + 1.47 log PGV
#   log PGV >= 0.53        :  MMI = 2.89 + 3.16 log PGV
#   PGV (cm/s) ของแกนราบที่ใหญ่กว่า
WGRW12_SIGMA = 0.63


def mmi_from_pgv(pgv_cm_s):
    """MMI จาก PGV (cm/s)  จำกัดในช่วง 1-10"""
    ly = np.log10(max(float(pgv_cm_s), 1e-6))
    if ly < -1.21:
        mmi = 4.62 + 2.17 * ly
    elif ly < 0.53:
        mmi = 3.78 + 1.47 * ly
    else:
        mmi = 2.89 + 3.16 * ly
    return float(np.clip(mmi, 1.0, 10.0))


def mmi_slope(mmi):
    """dMMI / dlog10(PGV) ของสูตร Worden ณ ระดับ MMI นั้น (จุดเปลี่ยนสูตรอยู่ที่ MMI 2.0 และ 4.56)"""
    return 2.17 if mmi < 2.0 else (1.47 if mmi < 4.56 else 3.16)


# ======================================================================
#  2) MMI ที่ทำนายได้  (IPE: Allen, Wald & Worden 2012)
# ======================================================================
AWW12 = dict(c0=2.085, c1=1.428, c2=-1.402, c4=0.078, m1=-0.209, m2=2.042, s1=0.82, s2=0.37, s3=22.9)
IPE_REF_VS30 = _cfg("IPE_REF_VS30", 400.0)   # สภาพดิน "เฉลี่ย" ที่ IPE แทน (ข้อมูล Did You Feel It?)
FAULT_MAX_M = 8.3                            # ใช้คำนวณความยาวรอยแตก ไม่เกินขนาดนี้


def ipe_allen2012(M, R_km):
    """
    MMI = c0 + c1 M + c2 ln( √(R² + rm²) ) + c4 ln(R/50) [R > 50 km],   rm = m1 + m2 exp(M - 5)
    σ   = s1 + s2 / (1 + (R/s3)²)
    """
    c = AWW12
    R = max(float(R_km), 1.0)
    rm = c["m1"] + c["m2"] * np.exp(M - 5.0)
    mmi = c["c0"] + c["c1"] * M + c["c2"] * np.log(np.sqrt(R ** 2 + rm ** 2))
    if R > 50.0:
        mmi += c["c4"] * np.log(R / 50.0)
    sigma = c["s1"] + c["s2"] / (1 + (R / c["s3"]) ** 2)
    return float(np.clip(mmi, 1.0, 10.0)), float(sigma)


def site_amplification(vs30=C.DEFAULT_VS30):
    """อัตราขยาย PGV ตาม Vs30:  log10(ARV) = 1.83 - 0.66 log10(Vs30)   (Midorikawa et al., 1994)"""
    return 10 ** (1.83 - 0.66 * np.log10(vs30))


def site_mmi_term(vs30, mmi):
    """MMI ที่เพิ่ม/ลดเพราะดิน (เทียบกับดินเฉลี่ย IPE_REF_VS30):  ΔMMI = (dMMI/dlogPGV) x log10(ARV(Vs30)/ARV(ref))"""
    ratio = site_amplification(vs30) / site_amplification(IPE_REF_VS30)
    return float(mmi_slope(mmi) * np.log10(ratio))


def fault_distance_from_hypo(R_km, Mw):
    """
    ประมาณระยะถึงรอยเลื่อนจากระยะ hypocentral โดยถือว่ารอยเลื่อนเป็นทรงกลมรัศมี L/2
      log10(L) = 0.5 Mw - 1.85  (km)  ->  X = max(R - L/2, 3)
    """
    L = 10 ** (0.5 * Mw - 1.85)
    return max(R_km - L / 2.0, 3.0)


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


def _predict_core(M, delta_km, depth_km, vs30, rjb_km=None):
    if rjb_km is None:                                       # จุดกำเนิดเป็นจุด: ประมาณระยะถึงรอยเลื่อนจากขนาด
        R = float(np.hypot(delta_km, depth_km))
        X = fault_distance_from_hypo(R, min(M, FAULT_MAX_M))
    else:                                                    # รู้แนวรอยแตก: ใช้ระยะถึงรอยแตก
        X = R = float(np.hypot(rjb_km, depth_km))
    ipe, s_ipe = ipe_allen2012(M, X)
    return dict(R=R, X=X, sigma_ipe=s_ipe, mmi=float(np.clip(ipe + site_mmi_term(vs30, ipe), 1.0, 10.0)))


def predict_mmi(M, delta_km, depth_km=10.0, vs30=C.DEFAULT_VS30, m_unc=0.3, rjb_km=None):
    """
    ทำนาย MMI ที่ตำแหน่งหนึ่ง (ก่อนคลื่นมาถึง)
      M       : ขนาดหลัก (จาก magnitude.py)     m_unc : ความไม่แน่นอนของขนาด (±)
      delta_km: ระยะจากจุดเหนือศูนย์กลาง,  vs30 : m/s
      rjb_km  : (ถ้ามี) ระยะในแนวราบถึงแนวรอยแตก (distance_to_trace_km)  None = ถือเป็นจุดกำเนิด
    return: mmi, sigma, mmi_low / mmi_high (±1σ), R, X
    """
    p = _predict_core(M, delta_km, depth_km, vs30, rjb_km)
    dm = (_predict_core(M + 0.1, delta_km, depth_km, vs30, rjb_km)["mmi"] - p["mmi"]) / 0.1      # dMMI/dM
    sigma = float(np.sqrt(p["sigma_ipe"] ** 2 + (dm * m_unc) ** 2))
    return dict(p, sigma=sigma, mmi_low=max(p["mmi"] - sigma, 1.0), mmi_high=min(p["mmi"] + sigma, 10.0), dmmi_dm=dm)


# ======================================================================
#  3) คำอธิบาย MMI และระดับการแจ้งเตือน (USGS)
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
