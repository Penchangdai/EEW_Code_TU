"""
p_wave.py
ตรวจจับคลื่น P (P-wave picking)
  1) Characteristic Function ของ Allen (1978)
  2) STA/LTA trigger
  3) ปรับเวลาให้แม่นยำด้วย AIC picker (Maeda, 1985)
  4) หาทิศการเคลื่อนที่แรก (first-motion polarity)
"""
import numpy as np

import config as C
from signal_utils import bandpass


def allen_cf(x, K=None):
    """CF(i) = x(i)^2 + K * (x(i) - x(i-1))^2"""
    x = np.asarray(x, dtype=float)
    dx = np.diff(x, prepend=x[0])
    if K is None:
        K = np.sum(np.abs(x)) / max(np.sum(np.abs(dx)), 1e-12)
    return x ** 2 + K * dx ** 2


def sta_lta(cf, fs, sta_sec, lta_sec):
    """
    Classic STA/LTA
      STA(i) = ค่าเฉลี่ย CF ใน nsta จุดล่าสุด
      LTA(i) = ค่าเฉลี่ย CF ใน nlta จุดก่อนหน้าหน้าต่าง STA
      R(i)   = STA(i) / LTA(i)
    """
    cf = np.asarray(cf, dtype=float)
    nsta, nlta = int(sta_sec * fs), int(lta_sec * fs)
    n = len(cf)
    ratio = np.zeros(n)
    if n < nsta + nlta:
        return ratio
    c = np.concatenate(([0.0], np.cumsum(cf)))
    i = np.arange(nsta + nlta - 1, n)
    sta = (c[i + 1] - c[i + 1 - nsta]) / nsta
    lta = (c[i + 1 - nsta] - c[i + 1 - nsta - nlta]) / nlta
    ratio[i] = sta / np.maximum(lta, 1e-20)
    return ratio


def aic_pick(x):
    """
    AIC(k) = k*log(var(x[0:k])) + (N-k-1)*log(var(x[k:N]))
    ตำแหน่งที่ AIC ต่ำสุด = จุดเริ่มคลื่น
    """
    x = np.asarray(x, dtype=float)
    N = len(x)
    if N < 6:
        return 0, np.zeros(N)
    c1 = np.cumsum(x)
    c2 = np.cumsum(x ** 2)
    k = np.arange(2, N - 2)
    var1 = c2[k - 1] / k - (c1[k - 1] / k) ** 2
    n2 = N - k
    var2 = (c2[-1] - c2[k - 1]) / n2 - ((c1[-1] - c1[k - 1]) / n2) ** 2
    aic = np.full(N, np.inf)
    aic[k] = k * np.log(np.maximum(var1, 1e-30)) + (N - k - 1) * np.log(np.maximum(var2, 1e-30))
    return int(np.argmin(aic)), aic


def first_motion(z, p_index, fs, window_sec=0.1):
    """ขั้วการเคลื่อนที่แรกของคลื่น P: +1 = ขึ้น (Compression), -1 = ลง (Dilatation)"""
    n = max(int(window_sec * fs), 2)
    seg = np.asarray(z[p_index:p_index + n], dtype=float) - z[p_index]
    return int(np.sign(seg[np.argmax(np.abs(seg))])) if len(seg) else 0


def first_sustained(ratio, threshold, n_min):
    """index แรกที่ ratio >= threshold ต่อเนื่องอย่างน้อย n_min จุด (กัน spike ของ noise)"""
    above = ratio >= threshold
    if n_min <= 1:
        idx = np.where(above)[0]
        return int(idx[0]) if len(idx) else None
    run = np.convolve(above.astype(int), np.ones(n_min, dtype=int), mode="valid")
    idx = np.where(run == n_min)[0]
    return int(idx[0]) if len(idx) else None


def detect_p_wave(z, fs, start_time=0.0, sta_sec=C.P_STA_SEC, lta_sec=C.P_LTA_SEC,
                  threshold=C.P_TRIGGER_ON, band=C.P_BANDPASS, refine=True,
                  min_trigger_sec=C.P_MIN_TRIGGER_SEC, aic_pre_sec=C.AIC_PRE_SEC,
                  aic_post_sec=C.AIC_POST_SEC):
    """
    ตรวจหาคลื่น P จากสัญญาณแนวดิ่ง (Z)
    z          : ความเร่งแนวดิ่ง (array)
    fs         : sampling rate (Hz)
    start_time : เวลาของ sample แรก (s, เช่น epoch time)
    return dict: detected, index, time, trigger_index, max_ratio, polarity
    """
    zf = bandpass(z, fs, *band) if band else np.asarray(z, dtype=float)
    ratio = sta_lta(allen_cf(zf), fs, sta_sec, lta_sec)
    trig = first_sustained(ratio, threshold, int(min_trigger_sec * fs))
    if trig is None:
        return {"detected": False, "max_ratio": float(ratio.max()), "ratio": ratio}

    pick = trig
    if refine:
        i0 = max(trig - int(aic_pre_sec * fs), 0)
        i1 = min(trig + int(aic_post_sec * fs), len(zf))
        k, _ = aic_pick(zf[i0:i1])
        pick = i0 + k
    return {
        "detected": True,
        "index": pick,
        "time": start_time + pick / fs,
        "trigger_index": trig,
        "max_ratio": float(ratio.max()),
        "polarity": first_motion(zf, pick, fs),
        "ratio": ratio,
    }
