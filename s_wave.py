"""
s_wave.py
Pick คลื่น S และ Lg ในหน้าต่างรอบเวลาที่คาด (หลังหาตำแหน่งคร่าว ๆ จากคลื่น P แล้ว)
  1) กรองสัญญาณ แล้วหมุนแกน N, E เป็น R (radial), T (transverse) จากทิศไปยังจุดเกิด
  2) AIC บนแกน R และ T รวมกัน หาจุดเริ่มคลื่น
  3) ให้ weight ตามอัตราส่วนพลังงานหลังและก่อน pick
"""
import numpy as np

import config as C
from signal_utils import bandpass
from p_wave import aic_pick


# ======================================================================
#  วิธี "usgs"
# ======================================================================
def _cfg(name, default):
    return getattr(C, name, default)




def rotate_ne_rt(n, e, baz):
    """หมุนแกน N, E -> R, T (สูตรเดียวกับ ObsPy)"""
    b = np.radians(baz)
    r = -e * np.sin(b) - n * np.cos(b)
    t = -e * np.cos(b) + n * np.sin(b)
    return r, t


# ======================================================================
#  Pick คลื่น S / Lg ในช่วงเวลาที่คาด (ใช้หลังหาตำแหน่งจากคลื่น P ได้แล้ว)
# ======================================================================
S_WIN_SNR = _cfg("S_WIN_SNR", (4.0, 2.5, 1.7))     # อัตราส่วนแอมพลิจูด หลัง/ก่อน ขั้นต่ำของ weight 1, 2, 3
S_PRE_SEC = _cfg("S_PRE_SEC", 5.0)                 # ช่วงก่อน pick ที่ใช้เทียบ (คือ P coda) (s)
S_POST_SEC = _cfg("S_POST_SEC", 3.0)               # ช่วงหลัง pick ที่ใช้เทียบ (s)
S_WIN_BAND = _cfg("S_WIN_BAND", (0.5, 5.0))        # ช่วงความถี่สำหรับ pick คลื่น S (s) ที่สถานีไกล


def s_quality(rf, tf, zf, pick, fs, pre_sec=S_PRE_SEC, post_sec=S_POST_SEC, phase="S"):
    """
    คุณภาพของ pick คลื่น S ที่อยู่ใน P coda (ไม่ได้เทียบกับ noise ก่อนเหตุการณ์ แต่เทียบกับช่วงก่อนหน้า)
      ratio  = RMS แนวราบ หลัง pick / RMS แนวราบ ก่อน pick   (พลังงานกระโดดขึ้นกี่เท่า)
      hv_gain = (H/V หลัง) / (H/V ก่อน)                        (S เด่นแนวราบ จึงควร > 1)
    weight: ratio >= 4 -> 1, >= 2.5 -> 2, >= 1.7 -> 3, น้อยกว่า -> 4   (S ไม่ให้ weight 0)
    """
    nb, na = max(int(pre_sec * fs), 2), max(int(post_sec * fs), 2)
    i0, i1 = max(pick - nb, 0), min(pick + na, len(rf))
    hb = np.mean(rf[i0:pick] ** 2 + tf[i0:pick] ** 2) if pick - i0 >= 2 else 1e-30
    ha = np.mean(rf[pick:i1] ** 2 + tf[pick:i1] ** 2) if i1 - pick >= 2 else 0.0
    vb = np.mean(zf[i0:pick] ** 2) if pick - i0 >= 2 else 1e-30
    va = np.mean(zf[pick:i1] ** 2) if i1 - pick >= 2 else 1e-30
    ratio = float(np.sqrt(ha / max(hb, 1e-30)))
    hv_gain = float((ha / max(va, 1e-30)) / max(hb / max(vb, 1e-30), 1e-30))
    n_early = max(int(0.5 * fs), 2)
    early = float(np.sqrt(np.mean(rf[pick:pick + n_early] ** 2 + tf[pick:pick + n_early] ** 2) / max(hb, 1e-30)))
    weight = 4
    for w, r_min in zip((1, 2, 3), S_WIN_SNR):
        if ratio >= r_min:
            weight = w
            break
    if phase != "Lg" and hv_gain < 1.0 and weight < 4:
        weight += 1                     # S: พลังงานเพิ่มแต่ไม่เด่นแนวราบ -> อาจเป็น P ระลอกหลัง
                                        # (Lg มีพลังงานแนวดิ่งมากด้วย จึงไม่ใช้เงื่อนไขนี้)
    onset = "i" if early >= 3.0 else "e"
    return {"ratio": ratio, "hv_gain": hv_gain, "weight": weight, "onset": onset}


def pick_s_window(z, n, e, fs, rec_t0, t_pred, baz, half_win=10.0, band=S_WIN_BAND, phase="S",
                  max_weight=3, not_before=None):
    """
    pick คลื่น S (หรือ Lg) ในช่วง t_pred ± half_win โดยรู้ back-azimuth จริงจากตำแหน่งที่หาได้
      1) กรองสัญญาณ แล้วหมุนแกน N, E -> R, T ด้วย baz (แม่นกว่า baz จากคลื่น P ที่ค่อย ๆ โผล่)
      2) AIC บน R + T ในหน้าต่าง -> จุดที่พลังงานแนวราบเปลี่ยนชัดที่สุด
         (ทดสอบกับข้อมูลจริงแล้ว AIC ตรงกับจุดที่ตาเห็นว่าพลังงานเพิ่ม ส่วน STA/LTA มัก trigger
          บนความผันผวนของ P coda ก่อนหน้า)
      3) คุณภาพจากพลังงานที่กระโดดขึ้นเทียบช่วงก่อนหน้า (P coda) และ H/V ที่เพิ่มขึ้น
    not_before : เวลา (epoch s) ที่ pick ต้องไม่ก่อนหน้านี้ เช่น เวลา P + 1 s
    return dict (index, time, weight, ratio, hv_gain, label) หรือ None
    """
    z, n, e = (np.asarray(x, dtype=float) for x in (z, n, e))
    i_pred = int(round((t_pred - rec_t0) * fs))
    w = int(half_win * fs)
    if not (0 <= i_pred < len(z)):
        return None
    s0 = max(i_pred - 2 * w - int((S_PRE_SEC + 15) * fs), 0)      # +15 s ให้ตัวกรองนิ่งก่อนถึงหน้าต่าง
    s1 = min(i_pred + 2 * w + int((S_POST_SEC + 5) * fs), len(z))
    zf, nf, ef = (bandpass(x[s0:s1], fs, *band) if band else x[s0:s1] for x in (z, n, e))
    rf, tf = rotate_ne_rt(nf, ef, baz)
    i0, i1 = max(i_pred - w - s0, 0), min(i_pred + w - s0, len(rf))
    if not_before is not None:
        i0 = max(i0, int(round((not_before - rec_t0) * fs)) - s0)
    if i1 - i0 < int(4 * fs):
        return None
    # AIC ใช้ข้อมูลเลยขอบหน้าต่างออกไปอีก half_win (ให้เห็นทั้งก่อนและหลังจุดเปลี่ยนเต็มที่)
    # แต่เลือกจุด pick เฉพาะในหน้าต่าง  -> ผลไม่เปลี่ยนมากเมื่อเวลาที่คาดขยับเล็กน้อย
    a0, a1 = max(i0 - w, 0), min(i1 + w, len(rf))
    _, a_r = aic_pick(rf[a0:a1])
    _, a_t = aic_pick(tf[a0:a1])
    aic = (a_r + a_t)[i0 - a0:i1 - a0]
    pick = i0 + int(np.argmin(aic))
    q = s_quality(rf, tf, zf, pick, fs, phase=phase)
    if q["weight"] > max_weight:
        return None
    return {"index": s0 + pick, "time": rec_t0 + (s0 + pick) / fs, "weight": q["weight"],
            "ratio": q["ratio"], "hv_gain": q["hv_gain"], "onset": q["onset"], "phase": phase,
            "label": f"{q['onset']}{phase}{q['weight']}"}
