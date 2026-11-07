# สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง

เอกสารนี้สรุปสูตรของโค้ดแต่ละไฟล์ แยกออกมาจากตัวโค้ด ค่าตัวเลขทุกค่าตรงกับโค้ดชุดใหม่
(p_wave.py, s_wave.py, epicenter.py, pipeline.py, magnitude.py, run_magnitude.py)

ลำดับการทำงานของระบบ:

```
ข้อมูลความเร่ง 3 แกน (Z, N, E) ของหลายสถานี
   → [1] ตรวจจับคลื่น P + ให้คะแนนคุณภาพ pick            (p_wave.py)
   → [2] ตรวจสอบเวลาไฟล์ + Association                     (pipeline.py)
   → [3] ตำแหน่งเริ่มต้นจาก P + S scan                       (pipeline.py)
   → [4] Re-pick P รอบเวลาที่คาด → Geiger จากคลื่น P          (p_wave.py, epicenter.py)
   → [5] ทำนายเวลา S / Lg → pick ในหน้าต่าง บนแกน R/T        (s_wave.py)
   → [6] Geiger จาก P + S + Lg → ตรวจทาน residual (วนซ้ำ)     (epicenter.py, pipeline.py)
   → [7] ขนาดทุกชนิดของแต่ละสถานี → น้ำหนัก → ขนาดหลัก        (magnitude.py)
   → [8] การ์ดแจ้งเตือนประชาชน                                 (run_magnitude.py)
   → [9] เวลาเดินทางถึงตำแหน่งเรา + เวลาเตือนล่วงหน้า            (travel_time.py, run_travel_time.py)
   → [10] ทำนายความรุนแรง (MMI) + แจ้งเตือนตามเกณฑ์ ShakeAlert  (intensity.py, run_intensity.py)
```

**สัญลักษณ์และหน่วยที่ใช้ร่วมกัน**

| สัญลักษณ์ | ความหมาย | หน่วย |
|---|---|---|
| V1 = Vpg, Vsg | ความเร็วในเปลือกโลก P = 6.0, S = 3.5 | km/s |
| V2 = Vpn, Vsn | ความเร็วใต้ Moho P = 8.0, S = 4.6 | km/s |
| V_Lg | ความเร็วคลื่น Lg = 3.5 | km/s |
| H | ความหนาเปลือกโลก = 30 | km |
| T0 | เวลาเกิดแผ่นดินไหว (origin time) | s |
| Tp, Ts, T_Lg | เวลาที่คลื่น P, S, Lg มาถึง | s |
| Δ | ระยะทางบนผิวโลกจากจุดเหนือศูนย์ (epicentral distance) | km หรือ องศา (1° = 111.195 km) |
| h | ความลึกของจุดกำเนิด | km |
| R | ระยะทางจากจุดกำเนิด = √(Δ² + h²) | km |
| fs | อัตราการสุ่มสัญญาณ (sampling rate) | Hz |
| ความเร่ง | gal = cm/s² = 0.01 m/s² | |

---

## 0. signal_utils.py — ฟังก์ชันพื้นฐาน

**ตัวกรอง band-pass** Butterworth อันดับ 4 แบบ **causal** (กรองทางเดียว `sosfilt`)
ตัวกรอง causal ไม่ทำให้สัญญาณรั่วมาก่อนคลื่นจริง แต่ทำให้เฟสช้าลงเล็กน้อย จึงต้องปรับเวลา pick ซ้ำบนข้อมูลดิบ (ข้อ 1.3)

**อินทิเกรตแบบสี่เหลี่ยมคางหมู** (ความเร่ง → ความเร็ว → การกระจัด)

    y[i] = y[i−1] + ( x[i] + x[i−1] ) / (2·fs)

หลังอินทิเกรตแต่ละครั้งใช้ high-pass Butterworth 0.075 Hz เพื่อกำจัด drift (Wu & Kanamori, 2005) ใช้กับ M_JMA, Pd, τc

**ระยะทางบนผิวโลก (Haversine)**

    a = sin²(Δφ/2) + cos φ1 · cos φ2 · sin²(Δλ/2)
    Δ = 2 · Rₑ · arcsin( √a )          Rₑ = 6371 km

**มุมทิศ (Azimuth) จากจุด 1 ไปจุด 2**

    θ = atan2( sin Δλ · cos φ2 ,  cos φ1 · sin φ2 − sin φ1 · cos φ2 · cos Δλ )

**จุดปลายทาง** เมื่อเดินจาก (φ1, λ1) ไปทิศ θ เป็นระยะ d  (δ = d / Rₑ)

    φ2 = arcsin( sin φ1 · cos δ + cos φ1 · sin δ · cos θ )
    λ2 = λ1 + atan2( sin θ · sin δ · cos φ1 ,  cos δ − sin φ1 · sin φ2 )

---

## 1. p_wave.py — ตรวจจับคลื่น P

### 1.1 Characteristic Function ของ Allen (1978)
ทำให้จุดที่สัญญาณเปลี่ยนทั้งแอมพลิจูดและความถี่เด่นชัดขึ้น (คำนวณบนแกน Z ที่กรองแล้ว)

    CF(i) = x(i)² + K · ( x(i) − x(i−1) )²
    K     = Σ|x(i)| / Σ|x(i) − x(i−1)|

### 1.2 STA/LTA

    STA(i) = (1/Nsta) · Σ CF(j)   ช่วง Nsta จุดล่าสุด
    LTA(i) = (1/Nlta) · Σ CF(j)   ช่วง Nlta จุดก่อนหน้าหน้าต่าง STA
    R(i)   = STA(i) / LTA(i)

นับเป็น trigger เมื่อ **R(i) ≥ threshold ต่อเนื่องอย่างน้อย N_min = min_trigger_sec · fs จุด**

### 1.3 กันการ trigger บน noise burst (รถวิ่งผ่าน ประตูปิด ฯลฯ)
คลื่น P จริงจะตามด้วยคลื่นระลอกหลังยาวหลายสิบวินาที ส่วน noise burst หายไปเร็ว จึงรับ trigger ก็ต่อเมื่อ

    RMS( z[t_trig + 2 s … t_trig + 10 s] )  ≥  2 × RMS( z[t_trig − 15 s … t_trig − 0.5 s] )

ถ้าไม่ผ่าน ข้ามไปหา trigger ถัดไป (สูงสุด 50 ครั้ง)

### 1.4 ปรับเวลาให้แม่นด้วย AIC (Maeda, 1985) — 2 รอบ

    AIC(k) = k · log( var(x[1..k]) ) + (N − k − 1) · log( var(x[k+1..N]) )

เวลาคลื่น P = ตำแหน่ง k ที่ **AIC ต่ำสุด**  →  Tp = t_start + k / fs

| รอบ | ข้อมูล | หน้าต่าง |
|---|---|---|
| 1 | ข้อมูลกรองแล้ว | t_trig − aic_pre ถึง t_trig + aic_post |
| 2 | ข้อมูลดิบ (ลบค่าเฉลี่ยก่อน pick) | pick รอบ 1 ± 1.0 s |

### 1.5 คุณภาพ pick (แบบ Hypoinverse, Klein 2002)

    SNR        = max|x[Tp … Tp+2 s]| / std(x[Tp−5 s … Tp])
    SNR_onset  = max|x[Tp … Tp+0.2 s]| / std(noise)
    onset      = i (impulsive) ถ้า SNR_onset ≥ 8,  มิฉะนั้น e (emergent)

| weight | SNR ขั้นต่ำ | ความไม่แน่นอน | ตัวคูณน้ำหนักใน Geiger |
|---|---|---|---|
| 0 | 20 | ±0.05 s | 1.00 |
| 1 | 10 | ±0.10 s | 0.75 |
| 2 | 5 | ±0.20 s | 0.50 |
| 3 | 3 | ±0.40 s | 0.25 |
| 4 | < 3 | ±1.0 s | 0 (ไม่ใช้) |

ถ้า onset = e ให้ weight ไม่ดีกว่า 2 ป้ายกำกับ pick เขียนแบบ `iPU0` = impulsive, P, ขึ้น, weight 0

### 1.6 ทิศการเคลื่อนที่แรก (First motion)
อ่านจาก **ข้อมูลดิบ** (ไม่กรอง) เครื่องหมายของแอมพลิจูดสูงสุดใน 0.1 s แรก ต้องสูงกว่า 3 × std(noise)
และอ่านเฉพาะ onset = i:  U (+1) = ขึ้น (Compression), D (−1) = ลง (Dilatation), ? = ไม่ชัด

### 1.7 Pick คลื่น P ในหน้าต่าง (re-pick หลังรู้ตำแหน่งคร่าว ๆ)

    หน้าต่าง = Tp(ทำนาย) ± 10 s

ใช้ขั้นตอนเดียวกับข้อ 1.1–1.5 แต่ยอมรับ trigger เฉพาะในหน้าต่าง ถ้าไม่มี trigger ใช้จุดที่ STA/LTA สูงสุด และให้ weight ไม่ดีกว่า 3

**ค่าพารามิเตอร์ 2 ชุด** (band-pass 1–10 Hz ทั้งสองชุด)

| ชุดค่า | ใช้กับ | STA | LTA | threshold | ต่อเนื่อง | หน้าต่าง AIC รอบ 1 |
|---|---|---|---|---|---|---|
| local | สถานีใกล้ (< ~300 km) | 0.5 s | 10 s | 4 | – | −2 ถึง +0.5 s |
| regional | สถานีไกล | 2 s | 60 s | 4 | 1 s | −8 ถึง +2 s |

### 1.8 การแปลงหน่วยข้อมูลจริง (data_loader.py)

    ความเร่ง (gal) = ความเร่ง (m/s²) × 100

---

## 2. s_wave.py — ตรวจจับคลื่น S และ Lg

### 2.1 หมุนแกน N, E → R (radial), T (transverse)
ใช้ back-azimuth จริงที่คำนวณจากตำแหน่งที่หาได้ (baz = ทิศจากสถานีไปจุดเกิด)

    R = −E · sin(baz) − N · cos(baz)
    T = −E · cos(baz) + N · sin(baz)

### 2.2 Pick ในหน้าต่างรอบเวลาที่คาด (วิธีหลักของระบบ)

    Ts(ทำนาย)  = T0 + t_S(Δ, h)         (โลก 2 ชั้น ข้อ 3.1)
    T_Lg(ทำนาย) = T0 + √(Δ² + h²) / 3.5
    หน้าต่าง    = เวลาทำนาย ± 10 s  และต้องไม่ก่อน Tp + 1 s

1. กรอง N, E, Z ช่วง 0.5–5 Hz แล้วหมุนเป็น R, T
2. AIC บน R และ T รวมกัน  `AIC(k) = AIC_R(k) + AIC_T(k)`
   ใช้ข้อมูลเลยขอบหน้าต่างออกไปอีก 10 s ด้านละข้าง แต่เลือกจุดต่ำสุดเฉพาะในหน้าต่าง
   (ผลจึงไม่ขยับมากเมื่อเวลาทำนายเปลี่ยนเล็กน้อย)
3. คุณภาพ pick เทียบกับช่วงก่อนหน้า (P coda) ไม่ใช่ noise ก่อนเหตุการณ์

        ratio   = √( mean(R² + T²)[ts … ts+3 s] / mean(R² + T²)[ts−5 s … ts] )
        hv_gain = (H/V หลัง pick) / (H/V ก่อน pick)        H = R² + T²,  V = Z²
        onset   = i ถ้า ratio ช่วง 0.5 s แรก ≥ 3,  มิฉะนั้น e

| weight | ratio ขั้นต่ำ |
|---|---|
| 1 | 4.0 |
| 2 | 2.5 |
| 3 | 1.7 |
| 4 (ไม่ใช้) | < 1.7 |

คลื่น S ไม่ให้ weight 0 เพราะระบุเวลายากกว่า P
ถ้าเป็นเฟส S และ hv_gain < 1 (พลังงานเพิ่มแต่ไม่เด่นแนวราบ อาจเป็น P ระลอกหลัง) ให้ weight แย่ลง 1 ระดับ
(Lg มีพลังงานแนวดิ่งมากด้วย จึงไม่ใช้เงื่อนไขนี้)
ถ้า S กับ Lg ที่ pick ได้ห่างกันไม่ถึง 3 s ถือเป็นคลื่นเดียวกัน เก็บตัวที่ weight ดีกว่า

### 2.3 วิธีแบบไม่รู้ตำแหน่ง (blind) — ใช้เมื่อมีสถานีเดียว
**วิธี "jump"** (ค่าเริ่มต้นของ regional) หาจุดที่พลังงานแนวราบกระโดดขึ้นและ H/V เพิ่มขึ้น

    score(k) = log10( Ha / Hb ) + w · log10[ (Ha/Va) / (Hb/Vb) ]       w = 2
    Ha, Va = ค่าเฉลี่ย H, V ช่วง 5 s หลังจุด k
    Hb, Vb = ค่าเฉลี่ย H, V ช่วง 15 s ก่อนจุด k  (ไม่ย้อนไปก่อนคลื่น P)

เลือก **ยอดแรก** ที่ score ≥ 0.7 × score สูงสุด (Sn มาก่อน Lg) แล้วปรับเวลาด้วย AIC_N + AIC_E

**วิธี "usgs"**: back-azimuth จากคลื่น P (`baz = atan2(Σ E·Z, Σ N·Z) + 180°`) → หมุน R/T →
CF = (R² + T²) × horizontality² → STA/LTA (LTA นับตั้งแต่คลื่น P) → AIC บน R + T

    horizontality = √(vN² + vE²)       v = eigenvector หลักของ covariance 3 แกน (หน้าต่าง 1 s)

**วิธี "first"** (เดิม): จุดแรกที่ STA_H / LTA_H ≥ threshold และ STA_H / STA_V ≥ 1

สถานีไกลที่แผ่นดินไหวใหญ่และรอยเลื่อนแตกยาว วิธี blind มักจับคลื่น P ระลอกหลังเป็นคลื่น S
จึงควรใช้ข้อ 2.2 ทุกครั้งที่มีตั้งแต่ 3 สถานีขึ้นไป

---

## 3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว

### 3.1 เวลาเดินทาง: แบบจำลองโลก 2 ชั้น (เปลือกโลกหนา H บนเนื้อโลก)

    คลื่นตรง  (Pg / Sg):  t = √(Δ² + h²) / V1
    คลื่นหัว  (Pn / Sn):  t = Δ / V2 + 2H · √( 1/V1² − 1/V2² )
    t_P (หรือ t_S) = min( คลื่นตรง , คลื่นหัว )        (เลือกคลื่นที่มาถึงก่อน)
    t_Lg = √(Δ² + h²) / V_Lg

| เฟส | V1 | V2 | 2H·√(1/V1² − 1/V2²) | ระยะที่คลื่นหัวเริ่มมาก่อน |
|---|---|---|---|---|
| P | 6.0 | 8.0 | 6.6 s | ~ 150 km |
| S | 3.5 | 4.6 | 11.1 s | ~ 150 km |

**ระยะจากเวลาต่าง (S − P หรือ Lg − P)** สร้างตาราง Δ = 1…4000 km แล้วเปิดตารางกลับ

    Δ = f⁻¹( Ts − Tp ),   f(Δ) = t_S(Δ, h) − t_P(Δ, h)

เทียบเท่าค่าคงที่ Omori k = Δ / (Ts − Tp) ซึ่งไม่คงที่ตามระยะ
(สูตร Omori เดิม R = k (Ts − Tp), k = 8.4 km/s (Pg/Sg) หรือ 10.82 km/s (Pn/Sn) ยังมีใน `omori_distance` ไว้เทียบ)

**เวลาเกิดจากสถานีเดียว**

    T0 = Tp − t_P(Δ, h)

### 3.2 ตรวจเวลาไฟล์ (Time coverage)

    ไม่ใช้สถานี ถ้า | t_start(ไฟล์) − median( t_start ทุกไฟล์ ) | > 600 s

(คนละเหตุการณ์ หรือเวลาในไฟล์เป็นเวลาท้องถิ่นแทน UTC)

### 3.3 Association — เวลา P ของสองสถานีต้องเข้ากันได้ทางกายภาพ

    สถานี a, b ขัดกัน ถ้า  | Tp_a − Tp_b | > d_ab / Vpn + 5 s        (d_ab = ระยะระหว่างสถานี)

ตัด P ของสถานีที่ขัดกับสถานีอื่นมากที่สุด (≥ 2 สถานี และมากที่สุดเพียงสถานีเดียว) ทำซ้ำจนไม่มีสถานีขัด

### 3.4 ตำแหน่งเริ่มต้น: scan จาก P + ทิศ + เวลา Lg + พลังงาน S (grid search)

ข้อมูล 2 อย่างที่วัดจากแต่ละสถานีได้เอง **โดยไม่ต้องรู้ตำแหน่งก่อน** (กันปัญหาวนกลับ: ตำแหน่งเริ่มต้นผิด → pick S ผิด → ตำแหน่งผิด)

    baz_obs = ทิศจากคลื่น P (PCA ข้อ 3.9) ช่วงคลื่น P 20 s แรก กรอง 0.05–0.2 Hz    (สถานีไกล)
    t_env   = เวลาที่ √(N² + E²) กรอง 0.5–5 Hz เฉลี่ย 5 s สูงสุดทั้งไฟล์            (≈ คลื่น Lg, ใช้ได้แม้ไม่มี pick P)

ทดสอบกับข้อมูลจริง: ทิศคลาด 1–11°, t_env ช้ากว่าคลื่น Lg ที่ควรเป็น 3–21 s

ลองทุกจุด (φ, λ) บนตาราง (ความลึกตรึง h = 10 km)

    rᵢ   = Tpᵢ − t_P(Δᵢ) − T0,        T0 = median_j( Tpⱼ − t_P(Δⱼ) )
    e_P  = Σ wᵢ · min(|rᵢ|, 5 s) / Σ wᵢ  ÷ 2 s
    e_baz = mean( min(|baz_obs − baz(φ,λ)|, 60°) ) ÷ 15°
    e_Lg = mean( min(| t_env − T0 − t_Lg(Δ) − 8 s |, 60 s) ) ÷ 15 s
    CFᵢ(τ) = log10( RMS_H[t, t+3 s] / RMS_H[max(t−5 s, Tp+1 s), t] ),  t = Tp + τ,  CF = 0 ถ้า τ < 6 s
    sc   = mean( clip( max_เฟส CFᵢ( t_S หรือ t_Lg − t_P ), 0, 1 ) )
    cost = e_P + 1.0·e_baz + 2.0·e_Lg − 1.0·sc

ตาราง 0.1° รอบสถานี ±12° แล้วละเอียด 0.02° รอบจุดที่ดีที่สุด ±0.6°
ใช้ได้ตั้งแต่มีคลื่น P 2 สถานี

CF ตัดช่วง τ < 6 s และช่วงเทียบไม่ย้อนไปก่อนคลื่น P เพราะถ้าไม่ตัด คลื่น P เองทำให้ CF สูงมากที่ τ = 1–4 s
และ scan จะเลือกตำแหน่งที่ติดกับสถานี (ข้อผิดพลาดเดิม: จุดศูนย์กลางไปตกที่เชียงใหม่)

**ผลทดสอบ** (ทุกชุด 3–5 สถานีจาก BKKA, CMMT, NAN, PHEA, PHIT เทียบกับ USGS)

| วิธี | ระยะคลาดกลาง ๆ | แย่ที่สุด |
|---|---|---|
| เดิม (P + พลังงาน S) | ~450 km | 1,103 km |
| ใหม่ (P + ทิศ + Lg + S, ทิศใน Geiger) | ~12 km | 74 km |

### 3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`)
ไม่ทราบค่า **m = (φ, λ, h, T0)** (ถ้าตรึงความลึก เหลือ 3 ตัว) ใช้เฟส P, S, Lg พร้อมกัน

    t_calc,i = T0 + t_เฟส( Δ(φ, λ, สถานี i), h )
    rᵢ       = t_obs,i − t_calc,i

อนุพันธ์คำนวณเชิงตัวเลขแบบ central difference (ก้าว 0.01°, 0.01°, 0.5 km, 0.01 s)

    Gᵢₖ = [ t_calc,i(m + δₖ) − t_calc,i(m − δₖ) ] / (2δₖ)

**ทิศจากคลื่น P เป็นข้อมูลเพิ่ม** (สถานีที่ใช้คลื่น P และวัดทิศได้)

    r_baz,i = ( baz_obs,i − baz(φ, λ → สถานี i) ) × 0.2 s/องศา        น้ำหนัก 0.5
    ไม่ใช้ทิศที่ต่างจากตำแหน่งปัจจุบันเกิน 45° (เช่น แก้ความกำกวม 180° ผิด)

เมื่อสถานีอยู่ด้านเดียว เวลาคลื่นอย่างเดียวบอกตำแหน่งตามแนวโค้งรอบสถานีได้ไม่ดี ทิศช่วยตรึงตำแหน่งตามแนวนั้น

แก้สมการแบบ weighted, damped least squares (Levenberg–Marquardt)

    A  = Gᵀ W G                    W = diag(wᵢ),  wᵢ = ตัวคูณน้ำหนักจาก weight ของ pick (ข้อ 1.5)
    Δm = ( A + λ · diag(A) )⁻¹ · Gᵀ W r          λ = 0.01
    m  = m + Δm                    (จำกัดก้าวไม่เกิน 2° ต่อรอบ, สูงสุด 30 รอบ)

T0 เริ่มต้น = median( t_obs − t_เฟส ) จุดเริ่มต้น = ผลจากข้อ 3.4

**ความคลาดเคลื่อน**

    wRMS  = √( Σ wᵢ rᵢ² / Σ wᵢ )
    n_eff = (Σ wᵢ)² / Σ wᵢ²                       (จำนวนข้อมูลที่มีผลจริงเมื่อถ่วงน้ำหนัก)
    σ²    = wRMS² · n_eff / (n_eff − p)            p = จำนวนตัวแปรที่หา (3 หรือ 4)
    Cov   = σ² · ( Gᵀ W̃ G )⁻¹                      W̃ = W / mean(w)
    err   = √diag(Cov)    แปลงเป็น km:  err_y = err_φ × 111.195,  err_x = err_λ × 111.195 · cos φ

### 3.6 ขั้นตอนทั้งหมดของ pipeline

1. Auto-pick P ทุกสถานี (ข้อ 1) → ตรวจเวลาไฟล์ (3.2) → Association (3.3)
2. P + S scan (3.4)
3. Re-pick P ในหน้าต่าง ±10 s (1.7) → Geiger จากคลื่น P อย่างเดียว → ตรวจทาน (3.7)
4. ทำนายเวลา S / Lg → pick ในหน้าต่าง (2.2) → Geiger จาก P + S + Lg → ตรวจทาน
5. ทำขั้น 4 ซ้ำอีกรอบด้วยตำแหน่งใหม่ (รวม 2 รอบ)

ค่าเริ่มต้น: ใช้ pick ที่ weight ≤ 3, ความลึกตรึง 10 km (`fix_depth`)

### 3.7 ตรวจทาน residual (ตัด pick ที่ไม่เข้ากับตำแหน่ง)

| profile | P | S | Lg |
|---|---|---|---|
| regional | 3 s | 6 s | 6 s |
| local | 1 s | 2 s | 2 s |

**เครือข่ายเล็ก** (เวลาคลื่นที่ใช้ได้ ≤ 15 ค่า): leave-one-out ตัดได้สูงสุด 2 pick
ตัด pick ที่ residual เทียบกับตำแหน่งที่หา **โดยไม่ใช้มัน** เกินเกณฑ์ และ RMS ลดลงเหลือไม่ถึงครึ่งเมื่อตัดออก

**เครือข่ายใหญ่**: รอบแรกเกณฑ์กว้าง (P 10 s, S/Lg 15 s) แล้วปรับเกณฑ์แต่ละรอบ

    MAD   = 1.4826 · median( | rᵢ − median(r) | )        (เทียบเท่า σ ของการแจกแจงปกติ)
    เกณฑ์ใหม่ = min( max( เกณฑ์ตามตาราง , 3 · MAD ) , เกณฑ์รอบก่อน )

ตัดได้ไม่เกิน 40% ของ pick และต้องเหลือ P อย่างน้อย 3 ค่า
หยุดเมื่อตำแหน่งขยับ < 0.5 km และเกณฑ์เปลี่ยน < 0.05 s
pick ที่เคยตัดกลับมาใช้ได้ถ้าเข้ากับตำแหน่งใหม่

### 3.8 Azimuthal gap

    gap = มุมที่กว้างที่สุดระหว่างทิศของสถานีที่อยู่ติดกัน (มองจากจุดศูนย์กลาง)

gap > 180° หมายถึงสถานีอยู่ด้านเดียว ระยะตามแนวทิศของสถานีคลาดง่าย

### 3.9 วิธีสถานีเดียว (เก็บไว้เทียบ ไม่ใช่วิธีหลัก)
ทิศจาก PCA ของคลื่น P ช่วงแรก (regional 5 s, 0.1–0.5 Hz; local 1 s, 0.5–10 Hz)

    C = | cov(N,N)  cov(N,E) |       eigenvector หลัก (pN, pE) = ทิศสั่นแนวราบ (กำกวม ±180°)
        | cov(E,N)  cov(E,E) |
    ถ้า Σ Z · (pN·N + pE·E) < 0  →  (pN, pE) = −(pN, pE)
    back_azimuth = atan2(pE, pN) + 180°

ระยะจากข้อ 3.1 แล้วหาตำแหน่งด้วยสูตรจุดปลายทาง (ข้อ 0)
สถานีไกลที่คลื่น P ค่อย ๆ โผล่ การแก้ความกำกวม 180° ไม่น่าเชื่อถือ

---

## 4. magnitude.py — ขนาดแผ่นดินไหว

### 4.0 จำลองเครื่องวัดมาตรฐานในโดเมนความถี่ (s = 2πif)
ลบค่าเฉลี่ยก่อนคลื่น P, taper 5%, FFT, คูณ H(s), inverse FFT

    oscillator อันดับ 2:  D(s; T, h) = s² + 2hω₀s + ω₀²,   ω₀ = 2π/T

| เครื่องวัด | H (ความเร่ง → เอาต์พุต) | ค่าคงที่ | ใช้กับ |
|---|---|---|---|
| Wood-Anderson กำลังขยาย 1 | 1 / D(s; 0.8, 0.7) | T₀ 0.8 s, h 0.7 | ML |
| WWSSN-SP | s / [D(s; 1.0, 1.0) · D(s; 0.75, 1.0)] ปรับให้ \|H\| = 1 ที่ 1 Hz | seismometer 1 s, galvanometer 0.75 s | mb_Lg |
| ความเร็วคาบยาว | (1/s) × Butterworth high-pass 1/60 Hz × low-pass 1/3 Hz (อันดับ 2) | คาบ 3–60 s | Ms_BB |

**SNR ของการวัดแอมพลิจูด**

    SNR = |x(ยอดที่วัด)| / max|x(30 s ก่อนคลื่น P)|           ต้อง ≥ 3

### 4.1 ML — Local magnitude (IASPEI 2013 = Hutton & Boore, 1987)

    ML = log10(A) + 1.11 · log10(R) + 0.00189 · R − 2.09           R ≤ 1000 km

- A = แอมพลิจูดสูงสุด (nm) บน Wood-Anderson จำลอง หลังคลื่น P
- วัดแยกแกน N และ E แล้วเฉลี่ย (แบบ USGS / CISN) ใช้เฉพาะแกนที่ SNR ≥ 3

### 4.2 mb_Lg — Lg magnitude (Nuttli, 1973)

    mb_Lg = 3.75 + 0.90 · log10(Δ°) + log10(A/T)        0.5° ≤ Δ < 4°
    mb_Lg = 3.30 + 1.66 · log10(Δ°) + log10(A/T)        4° ≤ Δ ≤ 30°

- หน้าต่าง Lg: T0 + Δ/3.6 ถึง T0 + Δ/3.2 (ความเร็วกลุ่ม 3.6–3.2 km/s) แกน Z บน WWSSN-SP จำลอง
- A (µm) = ยอดที่ใหญ่ **อันดับ 3** ในหน้าต่าง แปลงกลับเป็นการกระจัดพื้นด้วย A = A_เครื่อง / |H(T)|
- T = 2 × ช่วงระหว่างจุดตัดศูนย์ที่ล้อมยอดนั้น จำกัดในช่วง 0.7–1.3 s

### 4.3 Ms_BB — Broadband surface-wave magnitude (IASPEI 2013)

    Ms_BB = log10( Vmax / 2π ) + 1.66 · log10(Δ°) + 0.3          2° ≤ Δ ≤ 160°,  h < 60 km

- Vmax (nm/s) = ความเร็วแนวดิ่งสูงสุด คาบ 3–60 s ในหน้าต่างคลื่น Rayleigh (T0 + Δ/4.5 ถึง T0 + Δ/2.5)
- ถ้าทุกสถานีอยู่ใกล้กว่า 20° (นอกช่วงมาตรฐานของ Ms_20) ลดน้ำหนักเหลือ 0.7 และเพิ่มช่วงคลาด +0.2

### 4.4 Mwp — Moment magnitude จากคลื่น P (Tsuboi et al., 1995; Whitmore et al., 2002)

    uz(t)  = การกระจัดแนวดิ่ง (อินทิเกรตความเร่ง 2 ครั้ง)
    M₀     = max| ∫ uz dt | · 4πρα³r / Fp         อินทิเกรตตั้งแต่ Tp ถึงก่อนคลื่น S   (N·m)
    Mwp    = ( log10 M₀ − 9.1 ) / 1.5
    Mwp*   = ( Mwp − 1.03 ) / 0.843                (แก้ค่าแบบ Whitmore et al., 2002)

- ρ = 3400 kg/m³, α = 7900 m/s, Fp = 0.52 (ค่าเฉลี่ย radiation pattern), r = ระยะ (m), ใช้ได้ 5°–105°
- ถ้าไม่มีเวลา S ใช้ Ts ≈ Tp + Δ (1/4.6 − 1/8.0)

**ตรวจความเสถียรของ Mwp** (เครื่องวัดความเร่ง baseline ลอยง่าย)
คำนวณซ้ำ 4 แบบ: ไม่กรอง และ high-pass 0.005, 0.01, 0.02 Hz (Butterworth อันดับ 4)

    Mwp(สถานี) = median( 4 ค่า ),   spread = max − min
    drift_ratio = | u(Tp) − u(Tp − 30 s) | / max| u − u(Tp) |      (> 0.2 = baseline ลอย)

### 4.5 ขนาดแบบเดิม / สำหรับเตือนภัยเร็ว (ไม่ใช่มาตรฐาน USGS ใช้เทียบ)

    M_JMA = log10(A) + 1.73 · log10(Δ) − 0.83        A = √(A_N² + A_E²) (µm),  Δ ≤ 500 km   (Tsuboi 1954)
    M_Pd  = ( log10(Pd) + 3.463 + 1.374 · log10(R) ) / 0.729    Pd = max|uz| ใน 3 s แรก (cm)  (Wu & Zhao 2006)
    M_τc  = 3.373 · log10(τc) + 5.787,   τc = 2π / √( ∫v² dt / ∫u² dt ) ใน 3 s แรก            (Wu & Kanamori 2005)

### 4.6 น้ำหนักของแต่ละสถานี (0–1)

    w = f_SNR × f_ระยะ × f_pick   (Mwp คูณเพิ่ม f_ความเสถียร)

    f_SNR      = 0.5 + 0.5 · (SNR − 3) / (10 − 3)   จำกัด 0.5–1.0     (SNR < 3 → 0, ไม่ทราบ → 0.7)
    f_ความเสถียร (Mwp) = 1 / (1 + (spread / 0.3)²) × min(1, 0.2 / drift_ratio)

| ชนิด | f_ระยะ = 1 เมื่อ | น้อยกว่านั้น |
|---|---|---|
| ML | R ≤ 600 km | 0.5 |
| Ms_BB | Δ ≥ 20° | 0.7 |
| M_Pd, M_τc | R ≤ 200 km | 0.3 (≤ 500 km), 0.1 |
| อื่น ๆ | ทุกระยะที่ใช้ได้ | – |

| weight ของ pick P | 0–1 | 2 | 3 | 4 |
|---|---|---|---|---|
| f_pick | 1.0 | 0.9 | 0.7 | 0.3 |

### 4.7 ขนาดของเครือข่ายแต่ละชนิด — มัธยฐานถ่วงน้ำหนัก

    เรียงค่า M₁ ≤ M₂ ≤ … ≤ Mₙ พร้อมน้ำหนัก
    M_net = Mₖ ตัวแรกที่ Σ_{j≤k} wⱼ ≥ 0.5 · Σ wⱼ
    MAD   = มัธยฐานถ่วงน้ำหนักของ | Mᵢ − M_net |

(USGS ใช้มัธยฐานของ station magnitude เพื่อลดผลของสถานีที่ค่าผิดปกติ)

### 4.8 คะแนนของแต่ละชนิด และขนาดหลัก

    data    = 1 − exp( −Σw / 1.5 )                 (1 สถานีน้ำหนักเต็ม ≈ 0.49, 3 สถานี ≈ 0.86)
    scatter = 1 / (1 + (MAD / 0.3)²)
    M0      = มัธยฐานถ่วงน้ำหนักของทุกชนิด  น้ำหนัก = base × data × scatter
    sat     = clip( 1 − 0.9 · (M0 − (M_sat − 0.5)) / 1.0 , 0.1 , 1 )     (ลดคะแนนเมื่อขนาดใกล้จุดอิ่มตัว)
    range   = 1 ถ้า M0 ≥ M_min,  มิฉะนั้น 0.3
    consist = 1 / (1 + ((M − M0) / 0.7)²)
    window  = clip( T_window / T_rupture(M0) , 0.2 , 1 )   (เฉพาะ Mwp),  log10 T_rupture ≈ 0.5·M − 2.0

    score   = base × sat × range × data × scatter × consist × window
    ขนาดหลัก = ชนิดที่ score สูงสุด

| ชนิด | base | อิ่มตัวที่ M_sat | M_min |
|---|---|---|---|
| Mwp | 1.00 | 8.0 | 5.0 |
| Ms_BB | 0.90 | 8.0 | 4.5 |
| mb_Lg | 0.80 | 6.5 | 3.0 |
| ML | 0.80 | 6.5 | 0.0 |
| M_JMA | 0.55 | 8.0 | 3.0 |
| M_Pd | 0.40 | 6.5 | 3.0 |
| M_τc | 0.30 | 6.5 | 3.0 |

base เรียงตามลำดับความสำคัญของ USGS (Mw > Ms > mb_Lg, ML > ขนาดเตือนภัยเร็ว)
ขนาดอิ่มตัวโดยประมาณจาก Kanamori (1983) และ Bormann (NMSOP-2)
M0 ใช้มัธยฐานแทนค่าสูงสุด เพราะบางชนิดอาจสูงเกินจริงจากข้อมูลไม่ดี (เช่น Mwp ที่ baseline ลอย)

### 4.9 ความน่าเชื่อถือ และช่วงคลาดเคลื่อน

    ระดับเริ่มต้น: score ≥ 0.5 → สูง,  ≥ 0.2 → ปานกลาง,  < 0.2 → ต่ำ
    ลดลง 1 ระดับต่อข้อ: สถานีที่ใช้ได้ < 3,  gap > 180°,  ยังไม่รู้ตำแหน่ง

    ± = max( 0.3 , MAD , 0.5 × ผลต่างระหว่างชนิดที่ score ≥ ครึ่งหนึ่งของอันดับ 1 )
        + 0.2  ถ้า gap > 180°            (directivity ทำให้แอมพลิจูดเพี้ยน)
        + 0.2  ถ้า Ms_BB และทุกสถานี < 20°
        + 0.5  ถ้ายังไม่รู้ตำแหน่ง        (ระยะจาก S−P สถานีเดียวคลาดได้มาก)

### 4.10 คำบอกขนาด (การแบ่งกลุ่มของ USGS)

| M | คำบอก | M | คำบอก |
|---|---|---|---|
| ≥ 8.0 | ใหญ่มาก (Great) | 4.0 – 4.9 | เล็ก (Light) |
| 7.0 – 7.9 | ใหญ่ (Major) | 3.0 – 3.9 | เล็กมาก (Minor) |
| 6.0 – 6.9 | รุนแรง (Strong) | < 3.0 | จิ๋ว (Micro) |
| 5.0 – 5.9 | ปานกลาง (Moderate) | | |

---

## 5. run_magnitude.py — การ์ดแจ้งเตือน

**ที่มาของตำแหน่ง**

| สถานการณ์ | วิธี | การ์ดแสดง |
|---|---|---|
| EPICENTER = None, ≥ 3 สถานี | pipeline (ข้อ 3.6) | "คำนวณจากคลื่นที่ N สถานี (P x, S y, Lg z) คลาดประมาณ ±e km" |
| EPICENTER = None, 1–2 สถานี | Δ จาก S − P (ข้อ 3.1) ของแต่ละสถานี | ไม่มีตำแหน่ง ขนาดเพิ่มช่วงคลาด +0.5 |
| ใส่ EPICENTER เอง | ไม่ได้คำนวณ | "ค่าที่ใส่เอง (EPICENTER) ไม่ได้คำนวณ" |

    e (คลาดของตำแหน่ง) = max( err_x , err_y )         (จากข้อ 3.5)

**เงื่อนไขการแสดงความลึก**

    แสดงความลึก ⇔ DEPTH_KM = None (คำนวณความลึก)  และ  0 < err_z ≤ 15 km

ถ้าตรึงความลึก (DEPTH_KM = ตัวเลข) ความลึกไม่ได้มาจากการวัด จึงไม่แสดง
ความลึกจะวัดได้ดีเมื่อมีสถานีห่างจากจุดเกิดไม่เกินราว 1–2 เท่าของความลึก

**เงื่อนไขอื่น**

    ไม่คำนวณขนาด ถ้า RMS ของตำแหน่ง > 10 s
    ใช้สถานี ถ้ามีข้อมูลก่อนคลื่น P ≥ 5 s และหลังคลื่น P ≥ 30 s
    หัวการ์ด "แจ้งเตือนแผ่นดินไหว" ถ้า M ≥ 3.5,  มิฉะนั้น "รายงานแผ่นดินไหว"
    เวลาประเทศไทย = UTC + 7 ชม.

---

## 6. travel_time.py + run_travel_time.py — เวลาที่คลื่นเดินทางมาถึงตำแหน่งของเรา และเวลาเตือนล่วงหน้า

ใช้ตำแหน่ง เวลาเกิด T0 เวลาคลื่นของแต่ละสถานี และขนาดหลัก ชุดเดียวกับ run_magnitude.py (ข้อ 3 และ 4)
ความลึก h = ค่าที่ตรึงไว้ (10 km) ถ้าไม่ได้วัด

**เวลาเดินทางถึงตำแหน่งเรา** (แบบจำลองเดียวกับข้อ 3.1)

    Δ = Haversine(จุดศูนย์กลาง, ตำแหน่งเรา)
    Tp(site)  = T0 + min( Pg , Pn )
    Ts(site)  = T0 + min( Sg , Sn )              คลื่น S ตัวแรก
    TLg(site) = T0 + √(Δ² + h²) / 3.5            ที่ระยะไกล แรงสั่นสูงสุดมักมากับ Lg และคลื่นผิวพื้น

**เวลาที่ระบบรู้ข้อมูลแต่ละขั้น** (เวลาที่ข้อมูลมาถึงสถานีจริง + เวลาประมวลผล 1 s)

| ขั้น | เวลา | ได้อะไร |
|---|---|---|
| 1 ตรวจพบคลื่น P | Tp ของสถานีที่ n + 3 s + 1 s   (n = min(3, จำนวนสถานี)) | ตำแหน่งจากคลื่น P อย่างเดียว |
| 2 ตำแหน่งสุดท้าย | max(เวลาคลื่นทุกตัวที่ใช้หาตำแหน่ง) + 1 s | ตำแหน่งที่แสดงในการ์ด |
| 3 ขนาดหลัก | max(เวลาที่ข้อมูลของขนาดชนิดนั้นครบ ทุกสถานีที่ใช้) + 1 s | ขนาดที่แสดงในการ์ด |

เวลาที่ข้อมูลของขนาดแต่ละชนิดครบ (ตามหน้าต่างที่ magnitude.py ใช้วัด)

| ชนิด | ข้อมูลครบเมื่อ |
|---|---|
| ML | เวลายอดบน Wood-Anderson |
| mb_Lg | T0 + Δ / 3.2 |
| Ms_BB | T0 + Δ / 2.5 |
| Mwp | Tp + ช่วงคลื่น P ก่อนคลื่น S |
| M_Pd, M_τc | Tp + 3 s |

ถ้าใส่ EPICENTER และไม่เลือกสถานี ใช้ T_alert = T0 + 10 s (ALERT_DELAY_S)

**เวลาเตือนล่วงหน้า (Lead time)**

    Lead time = T_คลื่น(site) − T_ขั้น          คิดทั้งก่อนคลื่น S และก่อนคลื่น Lg

- Lead time > 0 คือมีเวลาเตรียมตัว,  ≤ 0 คือคลื่นมาถึงก่อนระบบรู้ข้อมูลขั้นนั้น

**รัศมี Blind zone** ของแต่ละขั้น (เปิดตารางเวลาเดินทางคลื่น S กลับ)

    Δ_b = t_S⁻¹( T_ขั้น − T0 )

**ผลกับข้อมูลทดสอบ** (มธ. รังสิต, Δ ≈ 1000 km)

| ขั้น | หลังเกิดเหตุ | ก่อนคลื่น S | ก่อนคลื่น Lg |
|---|---|---|---|
| 1 ตรวจพบคลื่น P | 92 s | +137 s | +195 s |
| 2 ตำแหน่งสุดท้าย | 202 s | +27 s | +84 s |
| 3 ขนาดหลัก (Ms_BB) | 418 s | −189 s | −131 s |

ขั้น 1 เตือนได้เร็วแต่ตำแหน่งคลาด ~140 km และขนาดยังไม่รู้ (Pd/τc อิ่มตัว)
ขนาดหลักของแผ่นดินไหวใหญ่ (Ms_BB ใช้คลื่นผิวพื้น) จะได้หลังคลื่นมาถึงกรุงเทพฯ แล้ว
จึงควรเตือนขั้นแรกด้วย "แผ่นดินไหวขนาดใหญ่ ≥ ค่าเบื้องต้น" แล้วปรับปรุงเมื่อได้ขั้น 2–3

---

## 7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน

### 7.0 ภาพรวม

```
ตำแหน่ง + ขนาด M(t) ณ เวลานั้น ─┬─> ก) GMPE: PGV บนหิน × อัตราขยายดิน ─> MMI (Worden 2012)
                                 └─> ข) IPE: MMI จาก M และระยะโดยตรง (Allen 2012) + ผลของดิน
        MMI ทำนาย = (ก + ข)/2  + ค่าปรับจากสถานีที่สั่นผ่านไปแล้ว
        ─> เกณฑ์ ShakeAlert ─> ส่งคำเตือนไปยังตำแหน่งที่คลื่นยังไม่ถึง
```

### 7.1 MMI ที่วัดได้ที่สถานี — GMICE ของ Worden et al. (2012) (ที่ USGS ShakeMap ใช้)

    Y = PGA (gal) หรือ PGV (cm/s) ของแกนราบที่ใหญ่กว่า
    log Y < t1b        :  MMI = c1b + c2b · log Y
    t1b ≤ log Y < t1   :  MMI = c1  + c2  · log Y
    log Y ≥ t1         :  MMI = c3  + c4  · log Y          (จำกัด 1–10)

| | c1 | c2 | c3 | c4 | t1 | c1b | c2b | t1b | σ (MMI) |
|---|---|---|---|---|---|---|---|---|---|
| PGA (gal) | 1.78 | 1.55 | −1.60 | 3.70 | 1.57 | 1.71 | 2.08 | 0.14 | 0.66 |
| PGV (cm/s) | 3.78 | 1.47 | 2.89 | 3.16 | 0.53 | 4.62 | 2.17 | −1.21 | 0.63 |

ใช้ MMI จาก PGV เป็นค่าหลัก แสดง MMI จาก PGA ไว้เทียบ
(สูตรเดิม Wald et al. 1999 `MMI = 3.47 log PGV + 2.35` ใช้ได้เฉพาะ MMI ≥ V ยังมีใน `mmi_from_pgv`)

### 7.2 JMA Instrumental Intensity (คำนวณจากคลื่นความเร่งจริงที่สถานี)

1. แปลงความเร่ง 3 แกน (gal) เป็นโดเมนความถี่ด้วย FFT
2. คูณด้วย filter  W(f) = F1(f) · F2(f) · F3(f)

        F1(f) = √(1/f)                                              (period-effect)
        F2(f) = 1 / √(1 + 0.694X² + 0.241X⁴ + 0.0557X⁶
                      + 0.009664X⁸ + 0.00134X¹⁰ + 0.000155X¹²) ,  X = f/10   (high-cut)
        F3(f) = √( 1 − exp( −(f/0.5)³ ) )                            (low-cut)

3. Inverse FFT กลับเป็นโดเมนเวลา แล้วรวมเวกเตอร์

        a(t) = √( a_Z(t)² + a_N(t)² + a_E(t)² )

4. หาค่า a₀ ที่ทำให้ **เวลารวมที่ a(t) ≥ a₀ เท่ากับ 0.3 วินาที**
5. คำนวณความรุนแรง

        I = 2 · log₁₀(a₀) + 0.94

   ปัดเศษทศนิยมตำแหน่งที่ 3 แล้วตัดให้เหลือทศนิยม 1 ตำแหน่ง

**ตารางแปลงเป็นระดับ JMA**

| ค่า I | ระดับ | ค่า I | ระดับ |
|---|---|---|---|
| I < 0.5 | 0 | 4.5 ≤ I < 5.0 | 5− (5 Lower) |
| 0.5 – 1.5 | 1 | 5.0 ≤ I < 5.5 | 5+ (5 Upper) |
| 1.5 – 2.5 | 2 | 5.5 ≤ I < 6.0 | 6− (6 Lower) |
| 2.5 – 3.5 | 3 | 6.0 ≤ I < 6.5 | 6+ (6 Upper) |
| 3.5 – 4.5 | 4 | I ≥ 6.5 | 7 |

### 7.3 ทำนาย MMI — วิธี ก) GMPE → GMICE

    log10(L) = 0.5 · M − 1.85,     X = max( R − L/2 , 3 km )            (ระยะถึงรอยเลื่อนโดยประมาณ)
    log10(PGV₆₀₀) = 0.58·M + 0.0038·D + d − 1.29 − log10( X + 0.0028 · 10^(0.50·M) ) − 0.002·X   (Si & Midorikawa 1999)
    log10(ARV)    = 1.83 − 0.66 · log10(Vs30)                             (Midorikawa et al. 1994)
    PGV_surface   = PGV₆₀₀ · ARV
    MMI_ก         = Worden 2012 (PGV_surface)
    I_JMA         = 2.68 + 1.72 · log10(PGV_surface)                      (Midorikawa et al. 1999)

d = 0 (crustal), −0.02 (interplate), +0.12 (intraplate) ใช้ M ไม่เกิน 8.3 (ช่วงที่สมการปรับเทียบไว้)

### 7.4 ทำนาย MMI — วิธี ข) IPE ของ Allen, Wald & Worden (2012) (ที่ USGS PAGER ใช้) แบบระยะ hypocentral

    rm    = m1 + m2 · exp(M − 5)
    MMI_IPE = c0 + c1·M + c2 · ln √(R² + rm²)  + c4 · ln(R/50)  [เฉพาะ R > 50 km]
    σ_IPE = s1 + s2 / (1 + (R/s3)²)

| c0 | c1 | c2 | c4 | m1 | m2 | s1 | s2 | s3 |
|---|---|---|---|---|---|---|---|---|
| 2.085 | 1.428 | −1.402 | 0.078 | −0.209 | 2.042 | 0.82 | 0.37 | 22.9 |

IPE แทนดินเฉลี่ย (ข้อมูล "Did You Feel It?") จึงเพิ่มผลของดินเทียบกับดินอ้างอิง Vs30 = 400 m/s

    ΔMMI_site = (dMMI / dlog10 PGV) × log10( ARV(Vs30) / ARV(400) )        (ความชันจาก Worden 2012 ณ ระดับ PGV นั้น)
    MMI_ข     = MMI_IPE + ΔMMI_site

### 7.5 รวม 2 วิธี และความไม่แน่นอน

    MMI_ทำนาย = ( MMI_ก + MMI_ข ) / 2
    σ = √( σ_model² + (dMMI/dM · σ_M)² + (|MMI_ก − MMI_ข| / 2)² )
    σ_model = (σ_IPE + 0.63) / 2,   σ_M = ช่วงคลาดของขนาด (ข้อ 4.9)

### 7.6 ตรวจคุณภาพข้อมูล (QC) และปรับด้วยค่าที่วัดได้จริง

**QC ก่อนวัด PGA / PGV ที่สถานี**

    ช่วงที่วัด     : Tp(คาด) − 5 s  ถึง  T0 + Δ / 2.0 km/s + 120 s      (ค่าสูงสุดก่อนคลื่น P มาถึงเป็นไปไม่ได้)
    สัญญาณกระโดด : |a| > 15 × max( ค่ากลาง |a| ±1 s , ค่ากลาง |a| ช่วง 2 s ถัดไป )  -> แทนด้วยเส้นตรง
    สัญญาณขาด    : ค่ายอด (≥ 99.9% ของค่าสูงสุด) เท่ากันติดกัน ≥ 5 จุด  -> ไม่ใช้ (ค่าสูงสุดต่ำกว่าจริง)
    SNR           : PGA ÷ ค่าสูงสุดของ noise 60 s ก่อนคลื่น P  < 3  -> ไม่ใช้

**ค่าปรับ** จากสถานีที่ผ่าน QC และสั่นแรงสุดผ่านไปแล้ว ( r = MMI วัด − MMI ทำนาย )

    ค่าผิดปกติ: | r − median(r) | > max( 3 × 1.4826 MAD , 2 MMI )  -> ไม่ใช้ปรับ

| วิธี | ค่าปรับที่ตำแหน่ง x | ต้องมีสถานี |
|---|---|---|
| median | median(r) | 2 |
| type | median(r) ของสถานีประเภทเดียวกัน | 3 ต่อประเภท |
| trend | a + b·log10 R  (Huber regression) | 6 |
| median+krig / trend+krig | ส่วนข้างบน + kriging ของ r ที่เหลือ | 8 |

Kriging แบบ simple kriging (แนวคิดเดียวกับ USGS ShakeMap) ใช้ covariance แบบ exponential

    C(d) = σ² [ (1 − ν) e^(−d/L) + ν δ ]          ν = nugget,  L = ระยะสหสัมพันธ์
    c(x) = kᵀ C⁻¹ e,    kᵢ = σ² (1 − ν) e^(−dᵢ(x)/L)

**เลือกวิธีอัตโนมัติด้วย leave-one-out**: ซ่อนค่าที่วัดได้ของสถานี i แล้วทำนายจากสถานีที่เหลือ ทำทุกสถานี

    RMSE = √( mean( (rᵢ − ĉ₋ᵢ(xᵢ))² ) )         เลือกวิธี (และ L ∈ {25, 50, 100, 200, 400} km, ν ∈ {0.2, 0.5}) ที่ RMSE ต่ำสุด
    ความไม่แน่นอนสุดท้ายของ MMI ทำนาย = RMSE นี้ (วัดจากข้อมูลจริง ไม่ใช่ค่าจากสมการ)

**ผลกับข้อมูลจริง 48 สถานี** (แผ่นดินไหวเมียนมา 28 มี.ค. 2568, ผ่าน QC 47, ใช้ปรับ 45)

| วิธี | RMSE leave-one-out (MMI) |
|---|---|
| ไม่ปรับ | 1.12 |
| ค่ากลางทุกสถานี (วิธีเดิม) | 0.69 |
| ค่ากลางแยกประเภทสถานี | 0.69 |
| แนวโน้มตามระยะ | 0.67 |
| แนวโน้มตามระยะ + kriging (L 200 km, ν 0.2) | **0.47** (เลือก) |

QC พบ: SURT สัญญาณรบกวนก่อนคลื่นมาถึง (MMI 9.1 ที่ระยะ 1,400 km) และ CRMM SNR ต่ำ ไม่ใช้ทั้งคู่,
PHRS มีสัญญาณกระโดด (PGA 665 gal) ตัดแล้ว, BKSI residual +3.2 ไม่ใช้ปรับ

### 7.7 จำลองการทำงานแบบเวลาจริง (run_intensity.py)
ทุก 5 s ตั้งแต่ขั้น 1 (ข้อ 6) จนคลื่น Lg ถึงตำแหน่งสุดท้าย:

1. ขนาด M(t) = ขนาดหลัก (ข้อ 4.8) จากเฉพาะค่าที่ **ข้อมูลมาถึงสถานีแล้ว** (เวลาพร้อมตามตารางข้อ 6)
2. ตำแหน่ง = จากคลื่น P อย่างเดียว ก่อนขั้น 2, ตำแหน่งสุดท้ายหลังจากนั้น
3. ทำนาย MMI ทุกตำแหน่ง (ข้อ 7.3–7.6)
4. ตำแหน่งที่เข้าเกณฑ์ครั้งแรก = เวลาส่งคำเตือน

    Lead time = Ts(site) − T_ส่งเตือน   และ   T_Lg(site) − T_ส่งเตือน

### 7.8 คำอธิบายระดับ MMI (ตามตาราง USGS ShakeMap)  ระดับ = ค่าที่ปัดเป็นจำนวนเต็ม

| MMI | การสั่น | ความเสียหาย | MMI | การสั่น | ความเสียหาย |
|---|---|---|---|---|---|
| I | ไม่รู้สึก | ไม่มี | VI | แรง | เล็กน้อย |
| II–III | อ่อน | ไม่มี | VII | แรงมาก | ปานกลาง |
| IV | เบา | ไม่มี | VIII | รุนแรง | ปานกลางถึงหนัก |
| V | ปานกลาง | น้อยมาก | IX | รุนแรงมาก | หนัก |
| | | | X+ | รุนแรงที่สุด | หนักมาก |

---

## 8. เกณฑ์การออกคำเตือน — ShakeAlert (USGS)

| ระดับ | เงื่อนไข | เทียบกับ ShakeAlert |
|---|---|---|
| แจ้งเตือนสาธารณะ | M ≥ 5.0 และ MMI ทำนาย ≥ 3.5 (ระดับ IV) | Wireless Emergency Alert |
| แจ้งเตือนผ่านแอป | M ≥ 4.5 และ MMI ทำนาย ≥ 2.5 (ระดับ III) | แอปมือถือ เช่น MyShake |
| ไม่แจ้ง | อื่น ๆ | |

เกณฑ์ JMA เดิม (I_JMA ทำนาย ≥ 4.5 → WARNING) ยังมีค่าใน config.py ไว้เทียบ

**ผลกับข้อมูลทดสอบ** (M 8.3 Ms_BB, หลังปรับด้วยสถานี +1.1 MMI)

| ตำแหน่ง | Δ (km) | MMI ทำนาย | MMI วัดได้ | ส่งเตือน | ก่อนคลื่น S | ก่อนคลื่น Lg |
|---|---|---|---|---|---|---|
| CMMT | 470 | VI | V | +147 s | −34 s | −13 s |
| NAN | 603 | V | IV | +147 s | −5 s | +25 s |
| PHEA | 615 | V | VII | +147 s | −2 s | +29 s |
| PHIT | 708 | V | VI | +147 s | +18 s | +56 s |
| มธ. รังสิต | 1003 | V | – | +192 s | +37 s | +95 s |
| BKKA | 1043 | IV | IV | +192 s | +46 s | +106 s |

ช่วงแรก (+92 s) มีแค่ M_Pd 4.4 จากคลื่น P 3 s แรก (อิ่มตัว) จึงยังไม่เข้าเกณฑ์
เตือนได้เมื่อมี ML / mb_Lg (+147 s) และยกระดับเมื่อได้ Ms_BB (+192 s)

---

## เอกสารอ้างอิง

**การตรวจจับคลื่นและคุณภาพ pick**
- Allen, R. V. (1978). Automatic earthquake recognition and timing from single traces. *BSSA*, 68(5), 1521–1532.
- Maeda, N. (1985). A method for reading and checking phase times in auto-processing system of seismic wave data. *Zisin*, 38, 365–379.
- Klein, F. W. (2002). User's guide to HYPOINVERSE-2000. *USGS Open-File Report* 02-171.
- Diehl, T., Kissling, E. & Bormann, P. (2012). Tutorial for consistent phase picking at local to regional distances. *NMSOP-2*, IS 11.4.

**การหาตำแหน่ง**
- Geiger, L. (1912). Probability method for the determination of earthquake epicenters from the arrival time only. *Bull. St. Louis Univ.*, 8, 60–71.
- Lienert, B. R., Berg, E. & Frazer, L. N. (1986). HYPOCENTER: An earthquake location method using centered, scaled, and adaptively damped least squares. *BSSA*, 76(3), 771–783.
- Lomax, A., Michelini, A. & Curtis, A. (2009). Earthquake location, direct, global-search methods. *Encyclopedia of Complexity and Systems Science*, Springer.
- Kennett, B. L. N. (1986). Lg waves and structural boundaries. *BSSA*, 76(4), 1133–1141.

**ขนาดแผ่นดินไหว**
- IASPEI (2013). Summary of magnitude working group recommendations on standard procedures for determining earthquake magnitudes from digital data.
- Hutton, L. K. & Boore, D. M. (1987). The ML scale in Southern California. *BSSA*, 77(6), 2074–2094.
- Nuttli, O. W. (1973). Seismic wave attenuation and magnitude relations for eastern North America. *JGR*, 78(5), 876–885.
- Tsuboi, S., Abe, K., Takano, K. & Yamanaka, Y. (1995). Rapid determination of Mw from broadband P waveforms. *BSSA*, 85(2), 606–613.
- Whitmore, P. M., Tsuboi, S., Hirshorn, B. & Sokolowski, T. J. (2002). Magnitude-dependent correction for Mwp. *Science of Tsunami Hazards*, 20(4), 187–192.
- Kanamori, H. (1983). Magnitude scale and quantification of earthquakes. *Tectonophysics*, 93, 185–199.
- Bormann, P. (ed.) (2012). *New Manual of Seismological Observatory Practice (NMSOP-2)*, บทที่ 3. GFZ Potsdam.
- Tsuboi, C. (1954). Determination of the Gutenberg-Richter's magnitude of earthquakes occurring in and near Japan. *Zisin*, 7.
- Kanamori, H. (2005). Real-time seismology and earthquake damage mitigation. *Annu. Rev. Earth Planet. Sci.*, 33.
- Wu, Y.-M. & Kanamori, H. (2005). Experiment on an onsite early warning method for the Taiwan early warning system. *BSSA*, 95(1).
- Wu, Y.-M. & Zhao, L. (2006). Magnitude estimation using the first three seconds P-wave amplitude in earthquake early warning. *GRL*, 33.

**ความรุนแรงและสภาพดิน**
- Japan Meteorological Agency (1996). Note on the JMA seismic intensity.
- Si, H. & Midorikawa, S. (1999). New attenuation relationships for peak ground acceleration and velocity considering effects of fault type and site condition. *J. Struct. Constr. Eng. AIJ*, 523.
- Midorikawa, S., Matsuoka, M. & Sakugawa, K. (1994). Site effects on strong-motion records observed during the 1987 Chiba-ken-toho-oki, Japan earthquake. *Proc. 9th Japan Earthquake Engineering Symposium*, 3.
- Midorikawa, S., Fujimoto, K. & Muramatsu, I. (1999). Correlation of new J.M.A. instrumental seismic intensity with former J.M.A. seismic intensity and ground motion parameters. *J. Inst. Social Safety Science*, 1.
- Wald, D. J., Quitoriano, V., Heaton, T. H. & Kanamori, H. (1999). Relationships between PGA, PGV, and MMI in California. *Earthquake Spectra*, 15(3).
- Worden, C. B., Gerstenberger, M. C., Rhoades, D. A. & Wald, D. J. (2012). Probabilistic relationships between ground-motion parameters and Modified Mercalli intensity in California. *BSSA*, 102(1), 204–221.
- Allen, T. I., Wald, D. J. & Worden, C. B. (2012). Intensity attenuation for active crustal regions. *J. Seismology*, 16, 409–433.
- Worden, C. B. & Wald, D. J. (2016). ShakeMap Manual Online. USGS.
- Kohler, M. D. et al. (2020). Earthquake Early Warning ShakeAlert 2.0: Public Rollout. *Seismol. Res. Lett.*, 91(3), 1763–1775.
- Hoshiba, M. & Aoki, S. (2015). Numerical shake prediction for earthquake early warning: data assimilation, real-time shake mapping, and simulation of wave propagation. *BSSA*, 105(3).

**ภาพรวม EEW (เอกสารแนบ)**
- Hoshiba, M. (2021). Real-time prediction of impending ground shaking: Review of wavefield-based (ground-motion-based) method for earthquake early warning. *Front. Earth Sci.*, 9:722784.
- Allen, R. M. & Melgar, D. (2019). Earthquake early warning: Advances, scientific challenges, and societal needs. *Annu. Rev. Earth Planet. Sci.*, 47.
- Multi-Channels Echo State Extreme Learning Machine for EEW. *J. Curr. Sci. Technol.*, 15(4), Article 139 (2025).

## ข้อจำกัดที่ควรทราบ

- **แบบจำลองความเร็ว** เป็นโลก 2 ชั้นค่าเดียวทั้งภูมิภาค (Vp 6.0/8.0, Vs 3.5/4.6, H 30 km) ถ้ามีแบบจำลองเปลือกโลกของไทย/เมียนมา (เช่น จากกรมอุตุนิยมวิทยา) ควรแก้ใน `config.py`
- **สูตรคลื่นหัว (Pn/Sn)** ในโค้ดใช้ 2H ซึ่งถูกต้องเมื่อจุดกำเนิดอยู่ที่ผิว สูตรที่รวมความลึกคือ
  `t = Δ/V2 + (2H − h)·√(1/V1² − 1/V2²)` ความต่างนี้คงที่ทุกสถานี (P ~1.1 s, S ~1.9 s ที่ h = 10 km)
  ผลต่อ P ถูกดูดซับไปใน T0 แต่ทำให้ S−P เพี้ยนประมาณ 0.8 s
- **ความลึก** ถ้าทุกสถานีอยู่ไกลหลายร้อย km หาความลึกไม่ได้ (เช่น ข้อมูลทดสอบได้ 46 ± 93 km) จึงตรึงไว้ที่ 10 km และไม่แสดง
- **สถานีอยู่ด้านเดียว** (gap > 180°) ระยะตามแนวทิศของสถานีคลาดง่าย ทั้งตำแหน่งและขนาด (directivity)
- **ขนาด** ML และ mb_Lg อิ่มตัวที่ M ~6.5, Pd / τc อิ่มตัวที่ M ~6.5 (ใช้คลื่นแค่ 3 s), Ms_BB และ Mwp อิ่มตัวที่ M ~8
  ขนาดมาตรฐานที่สุดของ USGS คือ Mww (W-phase moment tensor) ซึ่งระบบนี้ยังไม่ได้ทำ
- **ค่าคงที่ของ M_JMA, Pd, τc** ปรับเทียบจากข้อมูลญี่ปุ่นและไต้หวัน ถ้าจะใช้กับข้อมูลไทยควรทำ regression ใหม่กับ catalog ของกรมอุตุนิยมวิทยา
- **Point source** ระบบไม่ได้คำนวณรอยเลื่อนแบบ finite fault แผ่นดินไหวที่รอยเลื่อนแตกยาวหลายร้อย km (เช่น Mw 7.7 เมียนมา 2025) จุดศูนย์กลางคือจุดเริ่มแตก ไม่ใช่จุดที่สั่นแรงที่สุด
- **สมการความรุนแรง** (Si & Midorikawa, Allen 2012) ปรับเทียบไว้ที่ระยะไม่เกิน ~300 km สถานีไทยกับแผ่นดินไหวเมียนมาอยู่ 400–1000 km จึงเป็นการ extrapolate ค่าปรับจากสถานีที่วัดได้ (ข้อ 7.6) ช่วยลดความคลาดนี้ (ข้อมูลทดสอบ: ทำนายต่ำกว่าจริงเฉลี่ย 1.1 MMI ก่อนปรับ)
- **Vs30 ของสถานี** ถ้าไม่มีในไฟล์ stations.csv ใช้ 400 m/s ทุกสถานี ทำให้สถานีบนดินอ่อน (เช่น PHEA วัดได้ MMI VII แต่ทำนาย IV–V) คลาดมาก ควรเติม Vs30 จริงของแต่ละสถานี
- **แอ่งดินอ่อน** อย่างกรุงเทพฯ (และเม็กซิโกซิตี) ขยายคลื่นคาบยาว 1–5 s และทำให้สั่นนาน อาคารสูงอาจรู้สึกแรงกว่าค่า MMI ที่ทำนายจาก PGV
- **จำลองเวลาจริง** ขนาดของแต่ละช่วงเวลาใช้ระยะจากตำแหน่งสุดท้าย (ระบบจริงจะใช้ตำแหน่ง ณ ขณะนั้น) จึงอาจดีกว่าที่ระบบจริงทำได้เล็กน้อยในช่วงแรก