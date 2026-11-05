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
"""
import numpy as np

import config as C
from signal_utils import bandpass
from p_wave import aic_pick, first_sustained


def detect_s_wave(z, n, e, fs, p_index, start_time=0.0,
                  sta_sec=C.S_STA_SEC, threshold=C.S_TRIGGER_ON, band=C.P_BANDPASS,
                  min_lta_sec=C.S_MIN_LTA_SEC, min_sp_sec=C.S_MIN_SP_SEC,
                  max_sp_sec=C.S_MAX_SP_SEC, hv_min=C.S_HV_RATIO_MIN,
                  min_trigger_sec=0.0, aic_pre_sec=1.0, aic_post_sec=0.5,
                  method="first", jump_after_sec=5.0, jump_before_sec=15.0, hv_weight=1.0,
                  rel_score=0.5, min_score=0.5):
    """
    z, n, e : ความเร่ง 3 แกน, p_index : index ของคลื่น P (จาก p_wave.detect_p_wave)
    method  : "first" (STA/LTA จุดแรก) หรือ "jump" (จุดกระโดดของพลังงานแนวราบ + H/V)
    return dict: detected, index, time, sp_time (+ candidates เมื่อใช้ "jump")
    """
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


def predict_s_arrival(tp, origin_time, vp=C.VP_KM_S, vs=C.VS_KM_S):
    """ทำนายเวลา S จากเวลา P:  Ts = T0 + (Tp - T0) * Vp/Vs"""
    return origin_time + (tp - origin_time) * vp / vs


def sp_time_from_distance(hypo_dist_km, vp=C.VP_KM_S, vs=C.VS_KM_S):
    """ช่วงเวลา S-P จากระยะทาง:  Ts - Tp = R (1/Vs - 1/Vp)"""
    return hypo_dist_km * (1.0 / vs - 1.0 / vp)
