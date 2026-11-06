"""
s_wave.py
ตรวจจับคลื่น S (S-wave picking) จากสัญญาณแนวราบ N, E
  1) CF แนวราบ  H(i) = N(i)^2 + E(i)^2
  2) STA/LTA โดย LTA เริ่มนับตั้งแต่คลื่น P มาถึง (post-P LTA)
  3) เงื่อนไขอัตราส่วนพลังงานแนวราบ/แนวดิ่ง (S เด่นในแนวราบ)
  4) ปรับเวลาด้วย AIC
และฟังก์ชันทำนายเวลาคลื่น S เมื่อยังไม่เห็นคลื่น S จริง

มี 2 วิธี (เลือกด้วย method=)
  "first" : จุดแรกที่ STA/LTA แนวราบเกิน threshold (เหมาะกับสถานีใกล้, เดิม)
  "jump"  : หาจุดที่ "พลังงานแนวราบกระโดดขึ้น" และ "H/V เพิ่มขึ้น" ชัดที่สุด (เหมาะกับสถานีไกล)
            แก้ปัญหาแผ่นดินไหวใหญ่ที่รอยเลื่อนแตกยาว: คลื่น P ระลอกหลัง (P coda) ทำให้ STA/LTA
            เกินเกณฑ์ก่อนคลื่น S จริง แต่คลื่น P ระลอกหลังเพิ่มพลังงานทั้งแนวดิ่งและแนวราบ (H/V ไม่เพิ่ม)
            ส่วนคลื่น S เพิ่มพลังงานแนวราบเป็นหลัก (H/V เพิ่ม)
            score(k) = log10(Ha/Hb) + w·log10[(Ha/Va)/(Hb/Vb)]
              Ha, Va = พลังงานเฉลี่ยแนวราบ/แนวดิ่ง ช่วง jump_after_sec หลังจุด k
              Hb, Vb = พลังงานเฉลี่ยแนวราบ/แนวดิ่ง ช่วง jump_before_sec ก่อนจุด k
            เลือก "ยอดแรก" ที่ score >= rel_score × score สูงสุด (กันไปจับคลื่น Lg ที่ใหญ่กว่าแต่มาทีหลัง Sn)
  "usgs"  : ตามขั้นตอนแบบเครือข่าย USGS/ANSS
            1) หา back-azimuth จากการสั่นของคลื่น P (covariance ของ Z กับ N, E)
               แล้วหมุนแกน N, E -> R (radial), T (transverse)
            2) Polarization: ในหน้าต่างเลื่อน หาทิศสั่นหลัก (eigenvector ของ covariance 3 แกน)
               P สั่นใกล้แนวดิ่ง/แนว R, S สั่นแนวราบ -> ค่า horizontality = sqrt(vN^2 + vE^2)
               CF = (R^2 + T^2) x horizontality^2  (ลดน้ำหนักช่วงที่ยังเป็น P coda)
            3) หน้าต่างค้นหา: ถ้ารู้ origin time ใช้ Ts = T0 + 1.73 (Tp - T0) ± s_tol_sec
               ถ้าไม่รู้ ใช้ช่วง min_sp_sec ... max_sp_sec หลัง P
            4) STA/LTA (post-P LTA) บน CF -> trigger, แล้วปรับเวลาด้วย AIC บน R + T
            5) คุณภาพ pick: onset i/e, weight 0-4, ความไม่แน่นอน (label เช่น "eS2")
"""
import numpy as np

import config as C
from signal_utils import bandpass
from p_wave import aic_pick, first_sustained, pick_quality


def detect_s_wave(z, n, e, fs, p_index, start_time=0.0,
                  sta_sec=C.S_STA_SEC, threshold=C.S_TRIGGER_ON, band=C.P_BANDPASS,
                  min_lta_sec=C.S_MIN_LTA_SEC, min_sp_sec=C.S_MIN_SP_SEC,
                  max_sp_sec=C.S_MAX_SP_SEC, hv_min=C.S_HV_RATIO_MIN,
                  min_trigger_sec=0.0, aic_pre_sec=1.0, aic_post_sec=0.5,
                  method="first", jump_after_sec=5.0, jump_before_sec=15.0, hv_weight=1.0,
                  rel_score=0.5, min_score=0.5, origin_time=None):
    """
    z, n, e : ความเร่ง 3 แกน, p_index : index ของคลื่น P (จาก p_wave.detect_p_wave)
    method  : "first" (STA/LTA จุดแรก), "jump" (จุดกระโดดของพลังงานแนวราบ + H/V)
              หรือ "usgs" (R/T + polarization + AIC + น้ำหนัก pick)
    origin_time : เวลาเกิดแผ่นดินไหว (epoch s) ถ้ารู้ ใช้จำกัดหน้าต่างค้นหาของวิธี "usgs"
    return dict: detected, index, time, sp_time (+ candidates เมื่อใช้ "jump",
                 + baz, onset, weight, uncertainty, label เมื่อใช้ "usgs")
    """
    if method == "usgs":
        return _detect_s_usgs(z, n, e, fs, p_index, start_time, sta_sec, threshold, band,
                              min_lta_sec, min_sp_sec, max_sp_sec, hv_min, min_trigger_sec,
                              aic_pre_sec, aic_post_sec, origin_time)
    if method == "jump":
        return _detect_s_jump(z, n, e, fs, p_index, start_time, band, min_sp_sec, max_sp_sec, hv_min,
                              aic_pre_sec, aic_post_sec, jump_after_sec, jump_before_sec, hv_weight,
                              rel_score, min_score)
    zf = bandpass(z, fs, *band)
    nf = bandpass(n, fs, *band)
    ef = bandpass(e, fs, *band)
    h = nf ** 2 + ef ** 2
    v = zf ** 2

    nsta = max(int(sta_sec * fs), 1)
    i_start = p_index + max(int(min_sp_sec * fs), int(min_lta_sec * fs) + nsta)
    i_end = min(p_index + int(max_sp_sec * fs), len(h))
    if i_start >= i_end:
        return {"detected": False}

    ch = np.concatenate(([0.0], np.cumsum(h)))
    cv = np.concatenate(([0.0], np.cumsum(v)))
    i = np.arange(i_start, i_end)
    sta_h = (ch[i + 1] - ch[i + 1 - nsta]) / nsta
    lta_h = (ch[i + 1 - nsta] - ch[p_index]) / (i + 1 - nsta - p_index)   # post-P LTA
    sta_v = (cv[i + 1] - cv[i + 1 - nsta]) / nsta
    ratio = sta_h / np.maximum(lta_h, 1e-20)
    hv = sta_h / np.maximum(sta_v, 1e-20)

    # เงื่อนไข: STA/LTA >= threshold และ H/V >= hv_min ต่อเนื่องอย่างน้อย min_trigger_sec
    score = np.where(hv >= hv_min, ratio, 0.0)
    k = first_sustained(score, threshold, int(min_trigger_sec * fs))
    if k is None:
        return {"detected": False, "max_ratio": float(ratio.max())}

    trig = int(i[k])
    # AIC บนสัญญาณแนวราบทั้ง 2 แกน แล้วรวมกัน
    i0 = max(trig - int(aic_pre_sec * fs), p_index + int(min_sp_sec * fs))
    i1 = min(trig + int(aic_post_sec * fs), len(h))
    _, aic_n = aic_pick(nf[i0:i1])
    _, aic_e = aic_pick(ef[i0:i1])
    pick = i0 + int(np.argmin(aic_n + aic_e))

    return {
        "detected": True,
        "index": pick,
        "time": start_time + pick / fs,
        "trigger_index": trig,
        "sp_time": (pick - p_index) / fs,
        "max_ratio": float(ratio.max()),
    }


def _detect_s_jump(z, n, e, fs, p_index, start_time, band, min_sp_sec, max_sp_sec, hv_min,
                   aic_pre_sec, aic_post_sec, after_sec, before_sec, hv_weight, rel_score, min_score):
    """วิธี "jump": จุดที่พลังงานแนวราบกระโดดขึ้นพร้อม H/V เพิ่มขึ้น (ดูคำอธิบายด้านบนของไฟล์)"""
    zf, nf, ef = (bandpass(x, fs, *band) for x in (z, n, e))
    h = nf ** 2 + ef ** 2
    v = zf ** 2
    na, nb = max(int(after_sec * fs), 2), max(int(before_sec * fs), 2)
    i_start = p_index + max(int(min_sp_sec * fs), 2)
    i_end = min(p_index + int(max_sp_sec * fs), len(h) - na)
    if i_start >= i_end:
        return {"detected": False}

    ch = np.concatenate(([0.0], np.cumsum(h)))
    cv = np.concatenate(([0.0], np.cumsum(v)))
    i = np.arange(i_start, i_end)
    b0 = np.maximum(i - nb, p_index)                      # ช่วงก่อนหน้าไม่ย้อนไปก่อนคลื่น P
    ha = (ch[i + na] - ch[i]) / na
    va = (cv[i + na] - cv[i]) / na
    hb = (ch[i] - ch[b0]) / (i - b0)
    vb = (cv[i] - cv[b0]) / (i - b0)
    eps = 1e-20
    hv_a = ha / np.maximum(va, eps)
    hv_b = hb / np.maximum(vb, eps)
    score = np.log10(np.maximum(ha, eps) / np.maximum(hb, eps)) \
        + hv_weight * np.log10(np.maximum(hv_a, eps) / np.maximum(hv_b, eps))
    score = np.where(hv_a >= hv_min, score, -np.inf)       # หลังจุด S ต้องเด่นแนวราบ
    if not np.isfinite(score).any() or np.nanmax(score) < min_score:
        return {"detected": False, "max_score": float(np.nanmax(np.where(np.isfinite(score), score, -9)))}

    # ยอดของ score (ห่างกันอย่างน้อย jump_after_sec) เรียงตามเวลา
    peaks, k, sep = [], 0, na
    order = np.argsort(score)[::-1]
    taken = np.zeros(len(score), bool)
    for k in order:
        if not np.isfinite(score[k]) or score[k] < min_score or len(peaks) >= 8:
            break
        if taken[max(0, k - sep):k + sep].any():
            continue
        taken[k] = True
        peaks.append(int(k))
    smax = score[peaks[0]]
    good = sorted(p for p in peaks if score[p] >= rel_score * smax)
    k = good[0]                                            # ยอดแรกที่ใหญ่พอ (Sn มาก่อน Lg)
    trig = int(i[k])

    # ปรับเวลาด้วย AIC บนแนวราบ
    i0 = max(trig - int(aic_pre_sec * fs), i_start)
    i1 = min(trig + int(aic_post_sec * fs), len(h))
    _, aic_n = aic_pick(nf[i0:i1])
    _, aic_e = aic_pick(ef[i0:i1])
    pick = i0 + int(np.argmin(aic_n + aic_e))

    cands = [{"sp_time": (int(i[p]) - p_index) / fs, "score": float(score[p]),
              "hv": float(hv_a[p]), "chosen": p == k} for p in sorted(peaks)]
    return {
        "detected": True,
        "index": pick,
        "time": start_time + pick / fs,
        "trigger_index": trig,
        "sp_time": (pick - p_index) / fs,
        "max_ratio": float(10 ** smax),
        "score": float(score[k]),
        "candidates": cands,
        "method": "jump",
    }


# ======================================================================
#  วิธี "usgs"
# ======================================================================
def _cfg(name, default):
    return getattr(C, name, default)


S_VPVS = _cfg("S_VPVS", 1.73)                 # Vp/Vs สำหรับทำนายเวลา S
S_TOL_SEC = _cfg("S_TOL_SEC", 10.0)           # ± หน้าต่างรอบเวลา S ที่ทำนาย (s)
S_BAZ_SEC = _cfg("S_BAZ_SEC", 2.0)            # ช่วงหลัง P ที่ใช้หา back-azimuth (s)
S_POL_WIN_SEC = _cfg("S_POL_WIN_SEC", 1.0)    # หน้าต่างเลื่อนของ polarization (s)


def estimate_backazimuth(z, n, e, p_index, fs, win_sec=S_BAZ_SEC):
    """
    back-azimuth (องศาจากทิศเหนือ ตามเข็ม) จากการสั่นของคลื่น P
    คลื่น P ที่ขึ้น (Z+) จะผลักพื้นออกจากแหล่งกำเนิด
    -> ทิศของ (cov(N,Z), cov(E,Z)) ชี้ "ออกจาก" แหล่ง ดังนั้น baz = ทิศนั้น + 180
    return: baz (องศา), ความเชื่อมั่น 0-1 (สัดส่วนพลังงานแนวราบที่อยู่ในแนว R)
    """
    i1 = min(p_index + max(int(win_sec * fs), 4), len(z))
    zz, nn, ee = (np.asarray(x[p_index:i1], dtype=float) for x in (z, n, e))
    zz, nn, ee = zz - zz.mean(), nn - nn.mean(), ee - ee.mean()
    baz = (np.degrees(np.arctan2(np.sum(ee * zz), np.sum(nn * zz))) + 180.0) % 360.0
    r, t = rotate_ne_rt(nn, ee, baz)
    conf = float(np.sum(r ** 2) / max(np.sum(r ** 2) + np.sum(t ** 2), 1e-20))
    return float(baz), conf


def rotate_ne_rt(n, e, baz):
    """หมุนแกน N, E -> R, T (สูตรเดียวกับ ObsPy)"""
    b = np.radians(baz)
    r = -e * np.sin(b) - n * np.cos(b)
    t = -e * np.cos(b) + n * np.sin(b)
    return r, t


def horizontality(z, n, e, fs, win_sec=S_POL_WIN_SEC):
    """
    polarization แบบหน้าต่างเลื่อน: ทิศสั่นหลัก = eigenvector ที่ค่า eigenvalue สูงสุด
    return: horizontality = sqrt(vN^2 + vE^2) ในช่วง 0 (สั่นแนวดิ่ง) - 1 (สั่นแนวราบ)
    ค่าที่ index i ใช้ข้อมูลช่วง [i, i + win) (มองไปข้างหน้า เพื่อให้ขึ้นตอน S เริ่ม)
    """
    x = np.vstack([z, n, e]).astype(float)
    w = max(int(win_sec * fs), 4)
    N = x.shape[1]
    c = {}
    for a in range(3):
        for b in range(a, 3):
            c[a, b] = np.concatenate(([0.0], np.cumsum(x[a] * x[b])))
    m = N - w
    hor = np.zeros(N)
    if m <= 0:
        return hor
    cov = np.empty((m, 3, 3))
    for (a, b), cs in c.items():
        v = (cs[w:w + m] - cs[:m]) / w
        cov[:, a, b] = v
        cov[:, b, a] = v
    _, vec = np.linalg.eigh(cov)          # eigenvalue เรียงน้อย -> มาก
    v1 = vec[:, :, -1]                    # ทิศสั่นหลัก [vZ, vN, vE]
    hor[:m] = np.sqrt(v1[:, 1] ** 2 + v1[:, 2] ** 2)
    hor[m:] = hor[m - 1]
    return hor


def _detect_s_usgs(z, n, e, fs, p_index, start_time, sta_sec, threshold, band, min_lta_sec,
                   min_sp_sec, max_sp_sec, hv_min, min_trigger_sec, aic_pre_sec, aic_post_sec,
                   origin_time):
    zf, nf, ef = (bandpass(x, fs, *band) if band else np.asarray(x, float) for x in (z, n, e))

    # 1) back-azimuth + หมุนแกน
    baz, baz_conf = estimate_backazimuth(zf, nf, ef, p_index, fs)
    rf, tf = rotate_ne_rt(nf, ef, baz)

    # 2) polarization -> CF
    hor = horizontality(zf, nf, ef, fs)
    cf = (rf ** 2 + tf ** 2) * hor ** 2

    # 3) หน้าต่างค้นหา
    nsta = max(int(sta_sec * fs), 1)
    lo = p_index + max(int(min_sp_sec * fs), int(min_lta_sec * fs) + nsta)
    hi = min(p_index + int(max_sp_sec * fs), len(cf))
    predicted = None
    if origin_time is not None:
        tp = start_time + p_index / fs
        predicted = origin_time + S_VPVS * (tp - origin_time)
        ip = int(round((predicted - start_time) * fs))
        lo = max(lo, ip - int(S_TOL_SEC * fs))
        hi = min(hi, ip + int(S_TOL_SEC * fs))
    if lo >= hi:
        return {"detected": False, "method": "usgs", "baz": baz}

    # 4) STA/LTA (LTA นับตั้งแต่ P มาถึง) + เงื่อนไข H/V
    c1 = np.concatenate(([0.0], np.cumsum(cf)))
    hh = np.concatenate(([0.0], np.cumsum(nf ** 2 + ef ** 2)))
    vv = np.concatenate(([0.0], np.cumsum(zf ** 2)))
    i = np.arange(lo, hi)
    sta = (c1[i + 1] - c1[i + 1 - nsta]) / nsta
    lta = (c1[i + 1 - nsta] - c1[p_index]) / (i + 1 - nsta - p_index)
    hv = (hh[i + 1] - hh[i + 1 - nsta]) / np.maximum(vv[i + 1] - vv[i + 1 - nsta], 1e-20)
    ratio = sta / np.maximum(lta, 1e-20)
    score = np.where(hv >= hv_min, ratio, 0.0)
    k = first_sustained(score, threshold, int(min_trigger_sec * fs))
    if k is None:
        return {"detected": False, "method": "usgs", "baz": baz, "max_ratio": float(ratio.max())}
    trig = int(i[k])

    # AIC บน R + T (รอบแรก: ข้อมูลกรอง)
    i0 = max(trig - int(aic_pre_sec * fs), p_index + int(min_sp_sec * fs))
    i1 = min(trig + int(aic_post_sec * fs), len(cf))
    _, a_r = aic_pick(rf[i0:i1])
    _, a_t = aic_pick(tf[i0:i1])
    pick = i0 + int(np.argmin(a_r + a_t))

    # AIC รอบสองบนข้อมูลดิบ (ตัวกรอง zero-phase ทำให้ pick เร็วไป) ในหน้าต่างแคบ
    rr, tt = rotate_ne_rt(np.asarray(n, float), np.asarray(e, float), baz)
    n2 = int(_cfg("P_AIC2_SEC", 1.0) * fs)
    j0, j1 = max(pick - n2, p_index + int(min_sp_sec * fs)), min(pick + n2, len(rr))
    if j1 - j0 > 10:
        seg_r = rr[j0:j1] - rr[j0:max(pick, j0 + 2)].mean()
        seg_t = tt[j0:j1] - tt[j0:max(pick, j0 + 2)].mean()
        _, b_r = aic_pick(seg_r)
        _, b_t = aic_pick(seg_t)
        pick = j0 + int(np.argmin(b_r + b_t))

    # 5) คุณภาพ pick (noise = ช่วงก่อน S ซึ่งคือ P coda) บนแกนแนวราบที่แรงกว่า
    comp = tf if np.sum(tf[pick:pick + nsta] ** 2) >= np.sum(rf[pick:pick + nsta] ** 2) else rf
    q = pick_quality(comp, pick, fs)
    weight = max(q["weight"], 1)          # S ระบุเวลายากกว่า P จึงไม่ให้ weight 0
    uncert = float(_cfg("P_WEIGHT_UNCERT", (0.05, 0.1, 0.2, 0.4, 1.0))[weight])

    return {
        "detected": True,
        "method": "usgs",
        "index": pick,
        "time": start_time + pick / fs,
        "trigger_index": trig,
        "sp_time": (pick - p_index) / fs,
        "max_ratio": float(ratio.max()),
        "baz": baz,
        "baz_conf": baz_conf,
        "predicted_time": predicted,
        "onset": q["onset"],
        "weight": weight,
        "uncertainty": uncert,
        "snr": q["snr"],
        "label": f"{q['onset']}S{weight}",
        "r": rf, "t": tf, "horizontality": hor,
    }


def predict_s_arrival(tp, origin_time, vp=C.VP_KM_S, vs=C.VS_KM_S):
    """ทำนายเวลา S จากเวลา P:  Ts = T0 + (Tp - T0) * Vp/Vs"""
    return origin_time + (tp - origin_time) * vp / vs


def sp_time_from_distance(hypo_dist_km, vp=C.VP_KM_S, vs=C.VS_KM_S):
    """ช่วงเวลา S-P จากระยะทาง:  Ts - Tp = R (1/Vs - 1/Vp)"""
    return hypo_dist_km * (1.0 / vs - 1.0 / vp)


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
