"""
p_wave.py
ตรวจจับคลื่น P (P-wave picking) ตามขั้นตอนแบบเครือข่าย USGS/ANSS
(เฉพาะขั้นที่ทำได้จากสถานีเดียว)

ขั้น 1  Detection
  1.1 กรองสัญญาณ band-pass
  1.2 Characteristic Function ของ Allen (1978)
  1.3 STA/LTA trigger (ต้องเกิน threshold ต่อเนื่อง เพื่อกัน spike ของ noise)

ขั้น 2  Onset refinement + คุณภาพ pick
  2.1 ปรับเวลาให้แม่นด้วย AIC picker (Maeda, 1985) 2 รอบ:
      รอบแรกบนข้อมูลกรอง (หน้าต่างกว้าง) รอบสองบนข้อมูลดิบ (หน้าต่างแคบ)
  2.2 ลักษณะ onset : i = impulsive (ชัด), e = emergent (ค่อย ๆ โผล่)
  2.3 First motion : U = ขึ้น, D = ลง, ? = ไม่ชัด (อ่านจากข้อมูลที่ไม่ได้กรอง)
  2.4 น้ำหนัก pick 0-4 แบบ Hypoinverse (0 = ดีที่สุด, 4 = ไม่ใช้) + ความไม่แน่นอน (s)

กัน noise burst: trigger ที่หลังจากนั้น 2-10 s พลังงานกลับลงไปเท่า noise (เช่น รถวิ่งผ่าน ประตูปิด)
ไม่ใช่คลื่น P จริง (คลื่น P จริงตามด้วยคลื่นระลอกหลังยาวหลายสิบวินาที) -> ข้ามไปหา trigger ถัดไป

ขั้น 3  ระบุเฟส (Pg / Pn ...) ต้องใช้หลายสถานีและตารางเวลาเดินทาง
        จึงทำตอนหาจุดศูนย์กลาง ไม่ได้ทำในไฟล์นี้
"""
import numpy as np

import config as C
from signal_utils import bandpass


def _cfg(name, default):
    """อ่านค่าจาก config.py ถ้าไม่มีใช้ค่าเริ่มต้น"""
    return getattr(C, name, default)


# ----------------------------------------------------------------------
#  ค่าสำหรับประเมินคุณภาพ pick (ขั้น 2) - เพิ่มชื่อเดียวกันใน config.py เพื่อปรับได้
# ----------------------------------------------------------------------
P_NOISE_SEC = _cfg("P_NOISE_SEC", 5.0)            # หน้าต่าง noise ก่อน pick (s)
P_SIGNAL_SEC = _cfg("P_SIGNAL_SEC", 2.0)          # หน้าต่างสัญญาณหลัง pick (s)
P_IMPULSIVE_SEC = _cfg("P_IMPULSIVE_SEC", 0.2)    # ช่วงแรกหลัง pick ที่ใช้ตัดสิน i/e (s)
P_IMPULSIVE_SNR = _cfg("P_IMPULSIVE_SNR", 8.0)    # SNR ในช่วงนั้น >= ค่านี้ = impulsive
P_WEIGHT_SNR = _cfg("P_WEIGHT_SNR", (20.0, 10.0, 5.0, 3.0))       # SNR ขั้นต่ำของ weight 0,1,2,3
P_WEIGHT_UNCERT = _cfg("P_WEIGHT_UNCERT", (0.05, 0.1, 0.2, 0.4, 1.0))  # ความไม่แน่นอน (s) ของ weight 0-4
P_POLARITY_SEC = _cfg("P_POLARITY_SEC", 0.1)      # หน้าต่างอ่าน first motion (s)
P_POLARITY_MIN_SNR = _cfg("P_POLARITY_MIN_SNR", 3.0)  # first motion ต้องสูงกว่า noise กี่เท่า
P_AIC2_SEC = _cfg("P_AIC2_SEC", 1.0)              # หน้าต่าง ± (s) ของ AIC รอบสองบนข้อมูลดิบ
P_SUSTAIN_SEC = _cfg("P_SUSTAIN_SEC", (2.0, 10.0))  # ช่วงหลัง trigger ที่ต้องยังมีพลังงาน (s)
P_SUSTAIN_RATIO = _cfg("P_SUSTAIN_RATIO", 2.0)    # RMS ช่วงนั้น / RMS noise ก่อน trigger ขั้นต่ำ


# ======================================================================
#  ขั้น 1  Detection
# ======================================================================
def allen_cf(x, K=None):
    """CF(i) = x(i)^2 + K * (x(i) - x(i-1))^2"""
    x = np.asarray(x, dtype=float)
    dx = np.diff(x, prepend=x[0])
    if K is None:
        K = np.sum(np.abs(x)) / max(np.sum(np.abs(dx)), 1e-12)
    return x ** 2 + K * dx ** 2


LTA_MIN_SEC = getattr(C, "LTA_MIN_SEC", 10.0)


def sta_lta(cf, fs, sta_sec, lta_sec, lta_min_sec=None):
    """
    Classic STA/LTA
      STA(i) = ค่าเฉลี่ย CF ใน nsta จุดล่าสุด
      LTA(i) = ค่าเฉลี่ย CF ใน nlta จุดก่อนหน้าหน้าต่าง STA
      R(i)   = STA(i) / LTA(i)
    ช่วงต้นไฟล์ที่ข้อมูลยังไม่ถึง nlta จุด ใช้ LTA จากข้อมูลที่มีทั้งหมด (อย่างน้อย lta_min_sec วินาที)
      ไฟล์ที่เริ่มก่อนเกิดเหตุไม่นาน (เช่น 8 s) ถ้ารอให้ครบ LTA 60 s จะจับคลื่น P ที่มาก่อนนั้นไม่ได้เลย
      (ทดสอบ M4.6 21 ส.ค. 2568: pick ช้า 20-40 s ทุกสถานี ตำแหน่งคลาด 150 km)
    """
    cf = np.asarray(cf, dtype=float)
    nsta, nlta = int(sta_sec * fs), int(lta_sec * fs)
    nmin = min(int((LTA_MIN_SEC if lta_min_sec is None else lta_min_sec) * fs), nlta)
    n = len(cf)
    ratio = np.zeros(n)
    if n < nsta + max(nmin, 1):
        return ratio
    c = np.concatenate(([0.0], np.cumsum(cf)))
    i = np.arange(nsta + max(nmin, 1) - 1, n)
    sta = (c[i + 1] - c[i + 1 - nsta]) / nsta
    m = np.minimum(i + 1 - nsta, nlta)                     # จำนวนจุดที่ใช้ทำ LTA
    lta = (c[i + 1 - nsta] - c[i + 1 - nsta - m]) / m
    ratio[i] = sta / np.maximum(lta, 1e-20)
    return ratio


def first_sustained(ratio, threshold, n_min):
    """index แรกที่ ratio >= threshold ต่อเนื่องอย่างน้อย n_min จุด (กัน spike ของ noise)"""
    above = ratio >= threshold
    if n_min <= 1:
        idx = np.where(above)[0]
        return int(idx[0]) if len(idx) else None
    run = np.convolve(above.astype(int), np.ones(n_min, dtype=int), mode="valid")
    idx = np.where(run == n_min)[0]
    return int(idx[0]) if len(idx) else None


def is_sustained(zf, trig, fs, sustain_sec=P_SUSTAIN_SEC, ratio_min=P_SUSTAIN_RATIO, pre_sec=15.0):
    """trigger นี้ตามด้วยพลังงานต่อเนื่องไหม (RMS ช่วง sustain_sec หลัง trigger >= ratio_min x RMS ก่อนหน้า)"""
    a, b = trig + int(sustain_sec[0] * fs), trig + int(sustain_sec[1] * fs)
    if b > len(zf):                            # ใกล้ท้ายไฟล์ ตรวจไม่ได้ -> ถือว่าผ่าน
        return True
    p0, p1 = max(trig - int(pre_sec * fs), 0), max(trig - int(0.5 * fs), 1)
    pre = np.sqrt(np.mean(zf[p0:p1] ** 2)) if p1 - p0 > 2 else 0.0
    post = np.sqrt(np.mean(zf[a:b] ** 2))
    return post >= ratio_min * max(pre, 1e-30)


def first_real_trigger(zf, ratio, fs, threshold, n_min, i_from=0, i_to=None, max_tries=50):
    """trigger แรกในช่วง [i_from, i_to) ที่เป็นคลื่นจริง (ผ่าน is_sustained) -> (index, จำนวน burst ที่ข้าม)"""
    i_to = len(ratio) if i_to is None else i_to
    start, skipped = i_from, 0
    for _ in range(max_tries):
        k = first_sustained(ratio[start:i_to], threshold, n_min)
        if k is None:
            return None, skipped
        trig = start + k
        if is_sustained(zf, trig, fs):
            return trig, skipped
        skipped += 1
        above = np.where(ratio[trig:i_to] < threshold)[0]  # ข้ามไปหลังช่วงที่ยังเกิน threshold
        start = trig + (int(above[0]) if len(above) else (i_to - trig)) + 1
        if start >= i_to:
            return None, skipped
    return None, skipped


# ======================================================================
#  ขั้น 2  Onset refinement + คุณภาพ pick
# ======================================================================
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


def _noise_window(x, p_index, fs, noise_sec):
    """ช่วง noise ก่อน pick (ถ้าต้นไฟล์สั้นกว่าที่ขอ ใช้เท่าที่มี)"""
    i0 = max(p_index - int(noise_sec * fs), 0)
    seg = np.asarray(x[i0:p_index], dtype=float)
    return seg if len(seg) >= 2 else np.asarray(x[:max(p_index, 2)], dtype=float)


def pick_quality(xf, p_index, fs):
    """
    ประเมินคุณภาพ pick จากสัญญาณที่กรองแล้ว
    return: snr, snr_onset, onset ('i'/'e'), weight (0-4), uncertainty (s)
    """
    noise = np.std(_noise_window(xf, p_index, fs, P_NOISE_SEC))
    noise = max(noise, 1e-20)
    sig = np.abs(xf[p_index:p_index + max(int(P_SIGNAL_SEC * fs), 1)])
    early = np.abs(xf[p_index:p_index + max(int(P_IMPULSIVE_SEC * fs), 1)])
    snr = float(sig.max() / noise) if len(sig) else 0.0
    snr_onset = float(early.max() / noise) if len(early) else 0.0

    # i = แอมพลิจูดพุ่งขึ้นชัดภายในเสี้ยววินาทีแรก, e = ค่อย ๆ โผล่จาก noise
    onset = "i" if snr_onset >= P_IMPULSIVE_SNR else "e"

    weight = 4
    for w, s_min in enumerate(P_WEIGHT_SNR):
        if snr >= s_min:
            weight = w
            break
    if onset == "e":
        weight = max(weight, 2)     # onset ไม่ชัด ระบุเวลาได้ไม่แม่น จึงให้น้ำหนักไม่เกินระดับ 2
    return {
        "snr": snr,
        "snr_onset": snr_onset,
        "onset": onset,
        "weight": weight,
        "uncertainty": float(P_WEIGHT_UNCERT[weight]),
    }


def first_motion(z, p_index, fs, window_sec=P_POLARITY_SEC, min_snr=P_POLARITY_MIN_SNR,
                 noise_sec=P_NOISE_SEC):
    """
    ขั้วการเคลื่อนที่แรกของคลื่น P จากข้อมูลที่ไม่ได้กรอง
    (ตัวกรองโดยเฉพาะแบบ zero-phase ทำให้เกิดการสั่นก่อนคลื่นจริงและอาจกลับขั้วได้)
    +1 = ขึ้น (Compression), -1 = ลง (Dilatation), 0 = ไม่ชัด
    """
    z = np.asarray(z, dtype=float)
    pre = _noise_window(z, p_index, fs, noise_sec)
    base, noise = np.mean(pre), max(np.std(pre), 1e-20)
    n = max(int(window_sec * fs), 2)
    seg = z[p_index:p_index + n] - base
    if len(seg) == 0:
        return 0
    peak = seg[np.argmax(np.abs(seg))]
    if abs(peak) < min_snr * noise:
        return 0
    return int(np.sign(peak))


def detect_p_wave(z, fs, start_time=0.0, sta_sec=C.P_STA_SEC, lta_sec=C.P_LTA_SEC,
                  threshold=C.P_TRIGGER_ON, band=C.P_BANDPASS, refine=True,
                  min_trigger_sec=C.P_MIN_TRIGGER_SEC, aic_pre_sec=C.AIC_PRE_SEC,
                  aic_post_sec=C.AIC_POST_SEC):
    """
    ตรวจหาคลื่น P จากสัญญาณแนวดิ่ง (Z)
    z          : ความเร่งแนวดิ่ง (array)
    fs         : sampling rate (Hz)
    start_time : เวลาของ sample แรก (s, เช่น epoch time)
    return dict: detected, index, time, trigger_index, max_ratio, polarity,
                 onset, weight, uncertainty, snr, label (เช่น "iPU0"), ratio
    """
    z = np.asarray(z, dtype=float)

    # ---- ขั้น 1  Detection ----
    zf = bandpass(z, fs, *band) if band else z
    ratio = sta_lta(allen_cf(zf), fs, sta_sec, lta_sec)
    trig, n_burst = first_real_trigger(zf, ratio, fs, threshold, int(min_trigger_sec * fs))
    if trig is None:
        return {"detected": False, "max_ratio": float(ratio.max()), "ratio": ratio, "n_burst": n_burst}

    # ---- ขั้น 2.1  AIC refinement ----
    pick = trig
    if refine:
        i0 = max(trig - int(aic_pre_sec * fs), 0)
        i1 = min(trig + int(aic_post_sec * fs), len(zf))
        k, _ = aic_pick(zf[i0:i1])
        pick = i0 + k

        # ---- ขั้น 2.1b  AIC รอบสอง บนข้อมูลที่ไม่ได้กรอง ----
        # ตัวกรอง zero-phase (filtfilt) ทำให้สัญญาณ "รั่ว" มาก่อนคลื่นจริง
        # pick จากข้อมูลกรองจึงมักเร็วไปหลายร้อย ms จึงหาซ้ำในหน้าต่างแคบบนข้อมูลดิบ
        n2 = int(P_AIC2_SEC * fs)
        j0, j1 = max(pick - n2, 0), min(pick + n2, len(z))
        if j1 - j0 > 10:
            seg = z[j0:j1] - np.mean(z[j0:max(pick, j0 + 2)])
            k2, _ = aic_pick(seg)
            pick = j0 + k2

    # ---- ขั้น 2.2-2.4  คุณภาพ pick ----
    q = pick_quality(zf, pick, fs)
    polarity = first_motion(z, pick, fs)
    if q["onset"] == "e":
        polarity = 0                # เหมือนแนวปฏิบัติทั่วไป: อ่านขั้วเฉพาะ onset ที่ชัด
    pol_char = {1: "U", -1: "D"}.get(polarity, "?")

    return {
        "detected": True,
        "index": pick,
        "time": start_time + pick / fs,
        "trigger_index": trig,
        "max_ratio": float(ratio.max()),
        "polarity": polarity,
        "onset": q["onset"],
        "weight": q["weight"],
        "uncertainty": q["uncertainty"],
        "snr": q["snr"],
        "snr_onset": q["snr_onset"],
        "label": f"{q['onset']}P{pol_char}{q['weight']}",
        "n_burst": n_burst,             # จำนวน trigger ที่ข้ามไปเพราะเป็น noise burst
        "ratio": ratio,
    }


# ======================================================================
#  Pick คลื่น P ในช่วงเวลาที่คาด (ใช้หลังหาตำแหน่งคร่าว ๆ ได้แล้ว)
# ======================================================================
def pick_p_window(z, fs, rec_t0, t_pred, half_win=10.0, max_weight=3, sta_sec=C.P_STA_SEC,
                  lta_sec=C.P_LTA_SEC, threshold=C.P_TRIGGER_ON, band=C.P_BANDPASS,
                  min_trigger_sec=C.P_MIN_TRIGGER_SEC, aic_pre_sec=C.AIC_PRE_SEC,
                  aic_post_sec=C.AIC_POST_SEC, refine=True):
    """
    pick คลื่น P ในช่วง t_pred ± half_win (แบบที่ analyst ทำ: รู้ว่าคลื่นควรมาเมื่อไร แล้วดูเฉพาะช่วงนั้น)
    ใช้ขั้นตอนเดียวกับ detect_p_wave (STA/LTA -> AIC 2 รอบ -> คุณภาพ) แต่ยอมรับ trigger เฉพาะในหน้าต่าง
    จึงไม่ไป trigger บน noise ก่อนคลื่นจริง  ถ้าไม่มี trigger ในหน้าต่าง ใช้จุดที่ STA/LTA สูงสุดแทน
    รับพารามิเตอร์ชุดเดียวกับ detect_p_wave (เช่น **config.P_REGIONAL)
    return dict แบบเดียวกับ detect_p_wave (+ "in_window": True) หรือ None ถ้าคุณภาพต่ำเกิน max_weight
    """
    z = np.asarray(z, dtype=float)
    i_pred = int(round((t_pred - rec_t0) * fs))
    w = int(half_win * fs)
    if not (0 <= i_pred < len(z)):
        return None
    # ตัดช่วงข้อมูลให้พอสำหรับ LTA + ตัวกรอง แล้วคำนวณเฉพาะช่วงนั้น (เร็วกว่าทั้งไฟล์)
    s0 = max(i_pred - w - int((lta_sec + sta_sec + aic_pre_sec + 10) * fs), 0)
    s1 = min(i_pred + w + int((aic_post_sec + 5) * fs), len(z))
    seg = z[s0:s1]
    zf = bandpass(seg, fs, *band) if band else seg
    ratio = sta_lta(allen_cf(zf), fs, sta_sec, lta_sec)
    i0, i1 = max(i_pred - w - s0, 0), min(i_pred + w - s0, len(zf))
    if i1 - i0 < 10:
        return None
    t_real, _ = first_real_trigger(zf, ratio, fs, threshold, int(min_trigger_sec * fs), i0, i1)
    k = None if t_real is None else t_real - i0
    trig = i0 + (k if k is not None else int(np.argmax(ratio[i0:i1])))
    pick = trig
    if refine:
        a0 = max(trig - int(aic_pre_sec * fs), 0)
        a1 = min(trig + int(aic_post_sec * fs), len(zf))
        kk, _ = aic_pick(zf[a0:a1])
        pick = a0 + kk
        n2 = int(P_AIC2_SEC * fs)
        j0, j1 = max(pick - n2, 0), min(pick + n2, len(seg))
        if j1 - j0 > 10:
            raw = seg[j0:j1] - np.mean(seg[j0:max(pick, j0 + 2)])
            k2, _ = aic_pick(raw)
            pick = j0 + k2
    q = pick_quality(zf, pick, fs)
    if k is None:
        q["weight"] = max(q["weight"], 3)        # ไม่มี trigger ชัด -> คุณภาพไม่เกินระดับ 3
    if q["weight"] > max_weight:
        return None
    pol = first_motion(seg, pick, fs) if q["onset"] == "i" else 0
    pc = {1: "U", -1: "D"}.get(pol, "?")
    return {"detected": True, "index": s0 + pick, "time": rec_t0 + (s0 + pick) / fs,
            "trigger_index": s0 + trig, "weight": q["weight"], "onset": q["onset"], "polarity": pol,
            "snr": q["snr"], "uncertainty": q["uncertainty"], "in_window": True,
            "label": f"{q['onset']}P{pc}{q['weight']}"}
