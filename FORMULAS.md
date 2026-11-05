# สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning)

เอกสารนี้อธิบายสูตรของโค้ดแต่ละไฟล์ แยกออกมาจากตัวโค้ด
ลำดับการทำงานของระบบ:

```
ข้อมูลความเร่ง 3 แกน (Z, N, E)
   → [1] ตรวจจับคลื่น P   (p_wave.py)
   → [2] ตรวจจับคลื่น S   (s_wave.py)
   → [3] หาจุดศูนย์กลาง   (epicenter.py)
   → [4] หาขนาด M         (magnitude.py)
   → [5] เวลาเดินทางมาถึงตำแหน่งเรา (travel_time.py)
   → [6] ทำนายความรุนแรง และตัดสินใจออกคำเตือน (intensity.py, eew_system.py)
```

**สัญลักษณ์และหน่วยที่ใช้ร่วมกัน**

| สัญลักษณ์ | ความหมาย | หน่วย |
|---|---|---|
| Vp | ความเร็วคลื่น P = 6.0 | km/s |
| Vs | ความเร็วคลื่น S = 3.5 | km/s |
| T0 | เวลาเกิดแผ่นดินไหว (origin time) | s |
| Tp, Ts | เวลาที่คลื่น P, S มาถึง | s |
| Δ | ระยะทางจากจุดเหนือศูนย์ (epicentral distance) | km |
| h | ความลึกของจุดกำเนิด | km |
| R | ระยะทางจากจุดกำเนิด (hypocentral distance) | km |
| fs | อัตราการสุ่มสัญญาณ (sampling rate) | Hz |
| ความเร่ง | gal = cm/s² | |

---

## 0. signal_utils.py — ฟังก์ชันพื้นฐาน

**อินทิเกรตแบบสี่เหลี่ยมคางหมู** (ความเร่ง → ความเร็ว → การกระจัด)

    y[i] = y[i−1] + ( x[i] + x[i−1] ) / (2·fs)

หลังอินทิเกรตแต่ละครั้งใช้ high-pass Butterworth 0.075 Hz เพื่อกำจัด drift (ตาม Wu & Kanamori)

**ระยะทางบนผิวโลก (Haversine)**

    a = sin²(Δφ/2) + cos φ1 · cos φ2 · sin²(Δλ/2)
    Δ = 2 · Rₑ · arcsin( √a )          Rₑ = 6371 km

**มุมทิศ (Azimuth) จากจุด 1 ไปจุด 2**

    θ = atan2( sin Δλ · cos φ2 ,  cos φ1 · sin φ2 − sin φ1 · cos φ2 · cos Δλ )

**จุดปลายทาง** เมื่อเดินจาก (φ1, λ1) ไปทิศ θ เป็นระยะ d  (δ = d / Rₑ)

    φ2 = arcsin( sin φ1 · cos δ + cos φ1 · sin δ · cos θ )
    λ2 = λ1 + atan2( sin θ · sin δ · cos φ1 ,  cos δ − sin φ1 · sin φ2 )

**แปลง lat/lon เป็นพิกัดระนาบ (km)** รอบจุดอ้างอิง (φ0, λ0)

    x = Rₑ · (λ − λ0) · cos φ0       (เรเดียน)
    y = Rₑ · (φ − φ0)

---

## 1. p_wave.py — ตรวจจับคลื่น P

### 1.1 Characteristic Function ของ Allen (1978)
ทำให้จุดที่สัญญาณเปลี่ยนแปลงทั้งแอมพลิจูดและความถี่เด่นชัดขึ้น

    CF(i) = x(i)² + K · ( x(i) − x(i−1) )²
    K     = Σ|x(i)| / Σ|x(i) − x(i−1)|

### 1.2 STA/LTA (Short-Term Average / Long-Term Average)

    STA(i) = (1/Nsta) · Σ CF(j)   ช่วง Nsta จุดล่าสุด       (Nsta = 0.5 s · fs)
    LTA(i) = (1/Nlta) · Σ CF(j)   ช่วง Nlta จุดก่อนหน้า STA  (Nlta = 10 s · fs)
    R(i)   = STA(i) / LTA(i)

ถ้า **R(i) ≥ 4.0** ถือว่าตรวจพบคลื่น P (จุด trigger)

### 1.3 ปรับเวลาให้แม่นยำด้วย AIC (Maeda, 1985)
ในหน้าต่าง 2 s ก่อน ถึง 0.5 s หลังจุด trigger (N จุด)

    AIC(k) = k · log( var(x[1..k]) ) + (N − k − 1) · log( var(x[k+1..N]) )

เวลาคลื่น P คือตำแหน่ง k ที่ **AIC ต่ำสุด**  →  Tp = t_start + k / fs

### 1.4 เงื่อนไข trigger ต่อเนื่อง (กัน noise)
นับว่าเป็นคลื่น P ต่อเมื่อ R(i) ≥ threshold **ต่อเนื่อง** อย่างน้อย N_min จุด

    N_min = min_trigger_sec · fs

**ค่าพารามิเตอร์ 2 ชุด**

| ชุดค่า | ใช้กับ | STA | LTA | threshold | ต่อเนื่อง | หน้าต่าง AIC |
|---|---|---|---|---|---|---|
| local | สถานีใกล้ (< ~300 km) | 0.5 s | 10 s | 4 | – | −2 ถึง +0.5 s |
| regional | สถานีไกล เช่น กรุงเทพฯ กับแผ่นดินไหวเมียนมา | 2 s | 60 s | 4 | 1 s | −8 ถึง +2 s |

สถานีไกลคลื่น P จะค่อย ๆ โผล่ (emergent) และมี noise ของเมืองที่พุ่งขึ้นชั่วขณะ จึงต้องใช้หน้าต่างยาวขึ้น

### 1.5 การแปลงหน่วยข้อมูลจริง (data_loader.py)

    ความเร่ง (gal) = ความเร่ง (m/s²) × 100

### 1.6 ทิศการเคลื่อนที่แรก (First motion)
เครื่องหมายของแอมพลิจูดสูงสุดใน 0.1 s แรก: +1 = ขึ้น (Compression), −1 = ลง (Dilatation)

---

## 2. s_wave.py — ตรวจจับคลื่น S

คลื่น S สั่นในแนวราบเป็นหลัก จึงใช้แกน N และ E

    H(i) = N(i)² + E(i)²                    (CF แนวราบ)
    V(i) = Z(i)²                            (CF แนวดิ่ง)

    STA_H(i) = ค่าเฉลี่ย H ในช่วง 0.3 s ล่าสุด
    LTA_H(i) = ค่าเฉลี่ย H ตั้งแต่คลื่น P มาถึงจนถึงก่อนหน้าต่าง STA   (post-P LTA)
    STA_V(i) = ค่าเฉลี่ย V ในช่วง 0.3 s ล่าสุด

เงื่อนไขตรวจพบคลื่น S (ต้องผ่านทั้ง 2 ข้อ)

    STA_H / LTA_H ≥ 3.0      และ     STA_H / STA_V ≥ 1.0

แล้วปรับเวลาด้วย AIC (สูตรเดียวกับข้อ 1.3) โดยใช้ผลรวม AIC_N(k) + AIC_E(k)

**ค่าพารามิเตอร์คลื่น S 2 ชุด**

| ชุดค่า | STA | LTA หลัง P ขั้นต่ำ | threshold | ต่อเนื่อง | ค้นหา S−P | หน้าต่าง AIC |
|---|---|---|---|---|---|---|
| local | 0.3 s | 0.5 s | 3 | – | 0.3 – 60 s | −1 ถึง +0.5 s |
| regional | 2 s | 20 s | 4 | 1 s | 5 – 300 s | −8 ถึง +2 s |

**ทำนายเวลาคลื่น S จากคลื่น P** (เมื่อยังไม่เห็นคลื่น S)

    Ts = T0 + (Tp − T0) · Vp / Vs

**ช่วงเวลา S−P จากระยะทาง**

    Ts − Tp = R · ( 1/Vs − 1/Vp )

---

## 3. epicenter.py — หาจุดศูนย์กลางแผ่นดินไหว

### 3.1 วิธีสถานีเดียว (แนวคิดเดียวกับ JMA single-station)

**(ก) ระยะทาง — สูตรของ Omori**

    R = k · (Ts − Tp)
    k = Vp · Vs / (Vp − Vs) = 6.0 × 3.5 / 2.5 = 8.4 km/s
    Δ = √( R² − h² )          (h สมมติ 10 km)

สถานีไกล (> ~200 km) คลื่นแรกที่มาถึงวิ่งใต้เปลือกโลก (Pn, Sn) ใช้ Vp = 8.0, Vs = 4.6 km/s

    k = 8.0 × 4.6 / (8.0 − 4.6) = 10.82 km/s

**เวลาเกิดแผ่นดินไหวโดยประมาณจากสถานีเดียว**

    T0 = Tp − R / Vp

**(ข) ทิศทาง — Principal Component Analysis ของคลื่น P ช่วง 1 วินาทีแรก**

    C = | cov(N,N)  cov(N,E) |
        | cov(E,N)  cov(E,E) |

eigenvector (pN, pE) ของ eigenvalue ที่ใหญ่สุด = ทิศการสั่นในแนวราบ (ยังกำกวม ±180°)

แก้ความกำกวม: คลื่น P ที่ผิวดิน ถ้าพื้นขยับ **ขึ้น** พื้นจะขยับ **ออกจาก** แหล่งกำเนิด

    ถ้า Σ Z(i) · ( pN·N(i) + pE·E(i) ) < 0   →   (pN, pE) = −(pN, pE)
    away_az      = atan2(pE, pN)
    back_azimuth = away_az + 180°           (ทิศจากสถานีไปหาแหล่งกำเนิด)

ความเป็นเส้นตรงของการสั่น (rectilinearity) = 1 − λmin/λmax

**ค่าพารามิเตอร์ PCA**

| ชุดค่า | ช่วงคลื่น P ที่ใช้ | ช่วงความถี่ |
|---|---|---|
| local | 1 วินาทีแรก | 0.5 – 10 Hz |
| regional | 5 วินาทีแรก | 0.1 – 0.5 Hz |

**(ค) ตำแหน่ง** = จุดปลายทางจากสถานี ไปทิศ back_azimuth เป็นระยะ Δ (สูตรข้อ 0)

ในทางปฏิบัติ การแก้ความกำกวม 180° ด้วยทิศการเคลื่อนที่แรกของแกน Z ไม่น่าเชื่อถือเมื่อคลื่น P ค่อย ๆ โผล่ (สถานีไกล)
โปรแกรม `run_epicenter.py` จึงแสดงทั้ง 2 ตำแหน่ง

    ตัวเลือก A : ทิศ θ
    ตัวเลือก B : ทิศ θ + 180°

แล้วใช้สถานีที่ 2 หรือข้อมูลแนวรอยเลื่อนช่วยตัดสิน

### 3.2 วิธีหลายสถานี — Geiger's method (Iterative least squares)
ไม่ทราบค่า 4 ตัว: **m = (x, y, z, T0)**

เวลาที่คำนวณได้ของสถานี i (ใช้ v = Vp สำหรับคลื่น P, v = Vs สำหรับคลื่น S)

    Rᵢ = √( (x − xᵢ)² + (y − yᵢ)² + z² )
    tᵢ = T0 + Rᵢ / v

อนุพันธ์ย่อย (แถวของเมทริกซ์ G)

    ∂t/∂x = (x − xᵢ) / (v·Rᵢ)
    ∂t/∂y = (y − yᵢ) / (v·Rᵢ)
    ∂t/∂z =  z       / (v·Rᵢ)
    ∂t/∂T0 = 1

วนซ้ำจนกว่า |Δm| < 10⁻⁴ km

    rᵢ = tᵢ(observed) − tᵢ(calculated)          (residual)
    Δm = (GᵀG + λI)⁻¹ · Gᵀ r                    (damped least squares, λ = 10⁻³)
    m  = m + Δm

ความคลาดเคลื่อน

    RMS  = √( Σ rᵢ² / n )
    σ²   = Σ rᵢ² / (n − 4)
    Cov  = σ² · (GᵀG)⁻¹    →  error = √(diag(Cov))

เงื่อนไข: ต้องมีเวลาคลื่น (P + S) รวมกัน **อย่างน้อย 4 ค่า** และใช้ความลึกเริ่มต้นหลายค่า (5, 15, 30, 60 km) เลือกผลที่ RMS ต่ำสุด

---

## 4. magnitude.py — ขนาดแผ่นดินไหว

### 4.1 JMA Magnitude (สูตร Tsuboi, 1954) — วิธีหลักของ JMA

    M_JMA = log₁₀(A) + 1.73 · log₁₀(Δ) − 0.83
    A     = √( A_N² + A_E² )

- A_N, A_E = การกระจัดสูงสุดแกน N และ E (หน่วย **ไมครอน µm**, 1 cm = 10⁴ µm)
- Δ = ระยะทางจุดเหนือศูนย์ (km) ใช้ได้ดีที่ Δ ≤ 500 km และความลึก ≤ 60 km

### 4.2 Pd Magnitude (Wu & Zhao, 2006) — ใช้คลื่น P 3 วินาทีแรก

    Pd = max | u_z(t) |        ,   Tp ≤ t ≤ Tp + 3 s       (u_z = การกระจัดแนวดิ่ง, cm)
    log₁₀(Pd) = −3.463 + 0.729 · M − 1.374 · log₁₀(R)
    ⇒  M_Pd = ( log₁₀(Pd) + 3.463 + 1.374 · log₁₀(R) ) / 0.729

### 4.3 τc Magnitude (Kanamori, 2005; Wu & Kanamori, 2005)

    r  = ∫ v(t)² dt  /  ∫ u(t)² dt        (อินทิเกรตช่วง Tp ถึง Tp + 3 s)
    τc = 2π / √r                          (วินาที; คาบเด่นของคลื่น P)
    M_τc = 3.373 · log₁₀(τc) + 5.787

### 4.4 รวมผลหลายสถานี

    M_network = median( M₁, M₂, …, Mₙ )

(ใช้มัธยฐานเพื่อลดผลกระทบจากสถานีที่ค่าผิดปกติ)

---

## 5. travel_time.py — เวลาที่คลื่นเดินทางมาถึงตำแหน่งของเรา

    R  = √( Δ² + h² )                     Δ จากสูตร Haversine
    Tp(site) = T0 + R / Vp                 เวลาคลื่น P มาถึง
    Ts(site) = T0 + R / Vs                 เวลาคลื่น S มาถึง (คลื่นที่สร้างความเสียหาย)

**เวลาเตือนล่วงหน้า (Lead time)**

    Lead time = Ts(site) − T_alert

- T_alert = เวลาที่สถานีที่ 3 ได้ข้อมูลคลื่น P ครบ 3 วินาที + เวลาประมวลผล (1 s)
- ถ้าใช้สถานีเดียว ต้องรอให้คลื่น S มาถึงสถานีก่อน: T_alert = max(T_alert, Ts(station) + 1 s)
- Lead time > 0 คือมีเวลาเตรียมตัว,  ≤ 0 คืออยู่ใน Blind zone

**รัศมี Blind zone** (พื้นที่ที่คลื่น S มาถึงก่อนคำเตือน)

    R_b = Vs · (T_alert − T0)
    Δ_b = √( R_b² − h² )

---

## 6. intensity.py — ความรุนแรงแผ่นดินไหว (มาตรา JMA 0–7)

### 6.1 JMA Instrumental Intensity (คำนวณจากคลื่นความเร่งจริงที่สถานี)

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

### 6.2 ทำนายความรุนแรงที่ตำแหน่งของเรา (วิธีเดียวกับ JMA EEW)

**(ก) ขนาดรอยเลื่อนและระยะถึงรอยเลื่อน** (ถือว่ารอยเลื่อนเป็นทรงกลม)

    log₁₀(L) = 0.5 · M − 1.85        (L = ความยาวรอยเลื่อน, km)
    X = max( R − L/2 , 3 km )

**(ข) PGV บนชั้นหินฐานวิศวกรรม (Vs = 600 m/s) — Si & Midorikawa (1999)**

    log₁₀(PGV₆₀₀) = 0.58·Mw + 0.0038·D + d − 1.29
                    − log₁₀( X + 0.0028 · 10^(0.50·Mw) ) − 0.002·X

- D = ความลึก (km), X = ระยะถึงรอยเลื่อน (km), PGV หน่วย cm/s
- d = 0 (แผ่นดินไหวในเปลือกโลก), −0.02 (ระหว่างแผ่น), +0.12 (ภายในแผ่น)
- ใช้ M_JMA แทน Mw

**(ค) ขยายตามสภาพดิน (Midorikawa et al., 1994)**

    log₁₀(ARV) = 1.83 − 0.66 · log₁₀(Vs30)
    PGV_surface = PGV₆₀₀ · ARV

(Vs30 = ความเร็วคลื่นเฉือนเฉลี่ยของดินชั้นบน 30 m, หน่วย m/s; ดินอ่อนกรุงเทพฯ ≈ 150–250)

**(ง) แปลง PGV เป็นความรุนแรง JMA (Midorikawa et al., 1999)**

    I_JMA = 2.68 + 1.72 · log₁₀(PGV_surface)

**(จ) เทียบมาตรา Modified Mercalli (Wald et al., 1999)**

    MMI = 3.47 · log₁₀(PGV) + 2.35

### 6.3 เกณฑ์การออกคำเตือน (ตามแนวทาง JMA)

| เงื่อนไข | ผล |
|---|---|
| I_JMA ทำนาย ≥ 4.5 (ระดับ 5− ขึ้นไป) | **WARNING** – เตือนประชาชน |
| M ≥ 3.5 หรือ I_JMA ทำนาย ≥ 2.5 (ระดับ 3) | **FORECAST** – แจ้งข้อมูล |
| อื่น ๆ | ไม่ออกคำเตือน |

---

## เอกสารอ้างอิง

- Allen, R. V. (1978). Automatic earthquake recognition and timing from single traces. *BSSA*, 68(5).
- Maeda, N. (1985). A method for reading and checking phase times in auto-processing system of seismic wave data. *Zisin*, 38.
- Tsuboi, C. (1954). Determination of the Gutenberg-Richter's magnitude of earthquakes occurring in and near Japan. *Zisin*, 7.
- Kanamori, H. (2005). Real-time seismology and earthquake damage mitigation. *Annu. Rev. Earth Planet. Sci.*, 33.
- Wu, Y.-M. & Kanamori, H. (2005). Experiment on an onsite early warning method for the Taiwan early warning system. *BSSA*, 95(1).
- Wu, Y.-M. & Zhao, L. (2006). Magnitude estimation using the first three seconds P-wave amplitude in earthquake early warning. *GRL*, 33.
- Si, H. & Midorikawa, S. (1999). New attenuation relationships for peak ground acceleration and velocity considering effects of fault type and site condition. *J. Struct. Constr. Eng. AIJ*, 523.
- Midorikawa, S., Fujimoto, K. & Muramatsu, I. (1999). Correlation of new J.M.A. instrumental seismic intensity with former J.M.A. seismic intensity and ground motion parameters. *J. Inst. Social Safety Science*, 1.
- Wald, D. J. et al. (1999). Relationships between PGA, PGV, and MMI in California. *Earthquake Spectra*, 15(3).
- Japan Meteorological Agency (1996). Note on the JMA seismic intensity.
- Hoshiba, M. (2021). Real-time prediction of impending ground shaking: Review of wavefield-based (ground-motion-based) method for earthquake early warning. *Front. Earth Sci.*, 9:722784 (เอกสารแนบ).
- Allen, R. M. & Melgar, D. (2019). Earthquake early warning: Advances, scientific challenges, and societal needs. *Annu. Rev. Earth Planet. Sci.*, 47 (เอกสารแนบ).
- Multi-Channels Echo State Extreme Learning Machine for EEW. *J. Curr. Sci. Technol.*, 15(4), Article 139 (2025) (เอกสารแนบ).

## ข้อจำกัดที่ควรทราบ

- แบบจำลองความเร็วคลื่นเป็นแบบชั้นเดียว (Vp, Vs คงที่) ถ้ามีแบบจำลองชั้นเปลือกโลกของประเทศไทยควรแก้ใน `config.py`
- ค่าคงที่ของสูตร Tsuboi, Pd และ τc ถูกปรับเทียบจากข้อมูลญี่ปุ่นและไต้หวัน ถ้านำไปใช้กับข้อมูลไทยควรปรับเทียบค่าคงที่ใหม่ด้วยการทำ regression กับ catalog ของกรมอุตุนิยมวิทยา
- วิธี Pd / τc ให้ค่าอิ่มตัว (saturate) ที่ M > ~7 เพราะใช้คลื่นแค่ 3 วินาที
- ระบบเป็นแบบ point source ไม่ได้คำนวณความยาวรอยเลื่อนจริงแบบ finite fault
