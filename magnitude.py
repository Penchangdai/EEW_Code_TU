"""
magnitude.py
คำนวณขนาดแผ่นดินไหว (Magnitude)

  1) JMA displacement magnitude (สูตร Tsuboi, 1954) - วิธีที่ JMA ใช้
       M = log10(A) + 1.73 log10(Δ) - 0.83
       A = แอมพลิจูดการกระจัดสูงสุดแนวราบ (ไมครอน, µm), Δ = ระยะทางจุดเหนือศูนย์ (km)
  2) Pd magnitude (Wu & Zhao, 2006) - ใช้คลื่น P 3 วินาทีแรก (เร็วที่สุด)
       log10(Pd) = -3.463 + 0.729 M - 1.374 log10(R)
  3) τc magnitude (Wu & Kanamori, 2005) - ใช้คาบเด่นของคลื่น P 3 วินาทีแรก
       M = 3.373 log10(τc) + 5.787
"""
import numpy as np

import config as C
from signal_utils import acc_to_disp, acc_to_vel


# ---------------- 1) JMA (Tsuboi) ----------------
def jma_magnitude_from_amplitude(A_um, epi_dist_km):
    """M_JMA = log10(A) + 1.73 log10(Δ) - 0.83   (A ใน µm, ใช้ได้ Δ <= 500 km, ลึก <= 60 km)"""
    return np.log10(A_um) + 1.73 * np.log10(max(epi_dist_km, 1.0)) - 0.83


def max_horizontal_displacement_um(n_acc, e_acc, fs, i0=0, i1=None, n_pre=None):
    """A = sqrt(An^2 + Ae^2) ; An, Ae = ค่าสูงสุดของการกระจัดแต่ละแกน (cm -> µm)"""
    dn = acc_to_disp(n_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)[i0:i1]
    de = acc_to_disp(e_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)[i0:i1]
    An, Ae = np.max(np.abs(dn)), np.max(np.abs(de))
    return float(np.sqrt(An ** 2 + Ae ** 2) * 1e4)          # 1 cm = 10^4 µm


def jma_magnitude(n_acc, e_acc, fs, epi_dist_km, p_index, n_pre=None):
    A = max_horizontal_displacement_um(n_acc, e_acc, fs, i0=p_index, n_pre=n_pre)
    return float(jma_magnitude_from_amplitude(A, epi_dist_km)), A


# ---------------- 2) Pd ----------------
def pd_value(z_acc, fs, p_index, window=C.PD_TAUC_WINDOW_SEC, n_pre=None):
    """Pd = ค่าสูงสุดของการกระจัดแนวดิ่งในช่วง window วินาทีแรกของคลื่น P (cm)"""
    d = acc_to_disp(z_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)
    seg = d[p_index:p_index + int(window * fs)]
    return float(np.max(np.abs(seg)))


def magnitude_from_pd(pd_cm, hypo_dist_km):
    """M = (log10(Pd) + 3.463 + 1.374 log10(R)) / 0.729"""
    return (np.log10(pd_cm) + 3.463 + 1.374 * np.log10(max(hypo_dist_km, 1.0))) / 0.729


# ---------------- 3) τc ----------------
def tau_c(z_acc, fs, p_index, window=C.PD_TAUC_WINDOW_SEC, n_pre=None):
    """
    r  = ∫ v(t)^2 dt / ∫ u(t)^2 dt   (v = ความเร็ว, u = การกระจัด, ช่วง 0..τ0)
    τc = 2π / sqrt(r)
    """
    v = acc_to_vel(z_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)
    u = acc_to_disp(z_acc, fs, C.DISP_HIGHPASS_HZ, n_pre)
    sl = slice(p_index, p_index + int(window * fs))
    r = np.sum(v[sl] ** 2) / max(np.sum(u[sl] ** 2), 1e-30)
    return float(2 * np.pi / np.sqrt(r))


def magnitude_from_tauc(tc):
    """M = 3.373 log10(τc) + 5.787"""
    return 3.373 * np.log10(tc) + 5.787


# ---------------- รวมผลจากหลายสถานี ----------------
def network_magnitude(values):
    """ใช้ค่ามัธยฐาน (median) เพื่อลดผลจากสถานีที่ผิดปกติ"""
    vals = [v for v in values if v is not None and np.isfinite(v)]
    return float(np.median(vals)) if vals else None
