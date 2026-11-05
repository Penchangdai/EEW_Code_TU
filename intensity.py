"""
intensity.py
ความรุนแรงแผ่นดินไหว (Seismic Intensity) ตามมาตรา JMA (0 - 7)

  1) Instrumental intensity ของ JMA (วัดจากคลื่นความเร่งจริงที่สถานี)
  2) Predicted intensity (ทำนายที่ตำแหน่งของเรา) ตามวิธี JMA EEW:
       PGV (บนชั้นหินฐาน Vs=600 m/s) จากสมการ Si & Midorikawa (1999)
       -> ขยายตามสภาพดิน (ARV จาก Vs30)
       -> แปลง PGV เป็น JMA Intensity (Midorikawa et al., 1999)
"""
import numpy as np

import config as C


# ============================ 1) JMA instrumental intensity ============================
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
    """
    nfft = 1 << int(np.ceil(np.log2(len(z))))
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
    return {"I": float(I), "a0": float(a0), "a": a}


def jma_intensity_class(I):
    """แปลงค่าความรุนแรงเป็นระดับ JMA (Shindo)"""
    table = [(0.5, "0"), (1.5, "1"), (2.5, "2"), (3.5, "3"), (4.5, "4"),
             (5.0, "5-"), (5.5, "5+"), (6.0, "6-"), (6.5, "6+")]
    for upper, cls in table:
        if I < upper:
            return cls
    return "7"


# ============================ 2) Predicted intensity ============================
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
    """log10(ARV) = 1.83 - 0.66 log10(Vs30)   (Midorikawa et al., 1994)"""
    return 10 ** (1.83 - 0.66 * np.log10(vs30))


def intensity_from_pgv(pgv_cm_s):
    """I_JMA = 2.68 + 1.72 log10(PGV)   (Midorikawa et al., 1999)"""
    return 2.68 + 1.72 * np.log10(max(pgv_cm_s, 1e-6))


def mmi_from_pgv(pgv_cm_s):
    """MMI = 3.47 log10(PGV) + 2.35  (Wald et al., 1999) - สำหรับเทียบกับมาตรา Mercalli"""
    return float(np.clip(3.47 * np.log10(max(pgv_cm_s, 1e-6)) + 2.35, 1, 12))


def fault_distance_from_hypo(R_km, Mw):
    """
    ประมาณระยะถึงรอยเลื่อนจากระยะ hypocentral โดยถือว่ารอยเลื่อนเป็นทรงกลมรัศมี L/2
      log10(L) = 0.5 Mw - 1.85  (km)  ->  X = max(R - L/2, 3)
    """
    L = 10 ** (0.5 * Mw - 1.85)
    return max(R_km - L / 2.0, 3.0)


def predict_intensity(Mw, hypo_dist_km, depth_km, vs30=C.DEFAULT_VS30, eq_type="crustal"):
    X = fault_distance_from_hypo(hypo_dist_km, Mw)
    pgv_b = pgv_bedrock_si_midorikawa(Mw, X, depth_km, eq_type)
    pgv_s = pgv_b * site_amplification(vs30)
    I = intensity_from_pgv(pgv_s)
    return {"fault_dist_km": X, "pgv_bedrock": pgv_b, "pgv_surface": pgv_s, "arv": site_amplification(vs30),
            "intensity": float(I), "jma_class": jma_intensity_class(I),
            "mmi": mmi_from_pgv(pgv_s)}
