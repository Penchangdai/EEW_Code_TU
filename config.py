"""
config.py
ค่าคงที่และพารามิเตอร์ทั้งหมดของระบบ EEW (แก้ไขที่ไฟล์นี้ไฟล์เดียว)
"""

# ---------------- ข้อมูลเหตุการณ์ที่ใช้ ----------------
# ชื่อโฟลเดอร์เหตุการณ์ใน station/ ที่สคริปต์ run_*.py ใช้เมื่อตั้งเป็น "auto" หรือใส่แค่ชื่อไฟล์/รหัสสถานี
# เช่น "2025-3-28_Mw7.7_station", "2025-03-29_M5.9_AllTMD", "2026-05-18_M5.2_AllTMD"
# ใส่ None = ใช้ไฟล์ที่วางไว้ในโฟลเดอร์ station/ โดยตรง (แบบเดิม)
EVENT = "2025-3-28_Mw7.7_station"

# ---------------- แบบจำลองความเร็วคลื่น (Homogeneous half-space) ----------------
VP_KM_S = 6.0          # ความเร็วคลื่น P ในเปลือกโลก (km/s)
VS_KM_S = 3.5          # ความเร็วคลื่น S ในเปลือกโลก (km/s)  (Vp/Vs ~ 1.73)
EARTH_RADIUS_KM = 6371.0

# ---------------- การตรวจจับคลื่น P (STA/LTA) ----------------
P_STA_SEC = 0.5        # ความยาวหน้าต่าง Short-Term Average (s)
P_LTA_SEC = 10.0       # ความยาวหน้าต่าง Long-Term Average (s)
P_TRIGGER_ON = 4.0     # ค่า STA/LTA ที่ถือว่าตรวจพบคลื่น P
P_BANDPASS = (1.0, 10.0)   # ช่วงความถี่กรองสัญญาณก่อนตรวจจับ (Hz)
AIC_PRE_SEC = 2.0      # หน้าต่าง AIC ก่อนจุด trigger (s)
AIC_POST_SEC = 0.5     # หน้าต่าง AIC หลังจุด trigger (s)
P_MIN_TRIGGER_SEC = 0.0    # STA/LTA ต้องเกิน threshold ต่อเนื่องนานเท่านี้ (s) จึงนับว่าเป็นคลื่น P

# ชุดค่าสำหรับสถานีที่อยู่ไกล (regional, > ~300 km) คลื่น P ค่อย ๆ โผล่และมี noise เมือง
P_REGIONAL = dict(sta_sec=2.0, lta_sec=60.0, threshold=4.0, band=(1.0, 10.0),
                  min_trigger_sec=1.0, aic_pre_sec=8.0, aic_post_sec=2.0)

# ---------------- การตรวจจับคลื่น S ----------------
S_STA_SEC = 0.3
S_MIN_LTA_SEC = 0.5    # LTA หลังคลื่น P ต้องยาวอย่างน้อยเท่านี้ (s)
S_TRIGGER_ON = 3.0
S_MIN_SP_SEC = 0.3     # เวลา S-P ขั้นต่ำ (s)
S_MAX_SP_SEC = 60.0    # เวลา S-P สูงสุดที่ค้นหา (s)
S_HV_RATIO_MIN = 1.0   # อัตราส่วนพลังงาน แนวราบ/แนวดิ่ง ขั้นต่ำของคลื่น S

# ชุดค่าคลื่น S สำหรับสถานีที่อยู่ไกล (regional) : S-P อาจยาวถึง 2-3 นาที
#   method="jump" = หาจุดที่พลังงานแนวราบกระโดดขึ้นพร้อม H/V เพิ่มขึ้น (กันการจับคลื่น P ระลอกหลังเป็นคลื่น S)
#   ถ้าต้องการวิธีเดิม (STA/LTA จุดแรก) ให้เปลี่ยนเป็น method="first"
S_REGIONAL = dict(sta_sec=2.0, threshold=4.0, band=(0.5, 10.0), min_lta_sec=20.0,
                  min_sp_sec=15.0, max_sp_sec=300.0, hv_min=1.0, min_trigger_sec=1.0,
                  aic_pre_sec=8.0, aic_post_sec=2.0,
                  method="jump", jump_after_sec=5.0, jump_before_sec=15.0, hv_weight=2.0,
                  rel_score=0.7, min_score=0.5)

# ความเร็วคลื่นที่วิ่งใต้เปลือกโลก (Pn, Sn) ใช้กับสถานีไกล > ~200 km
VPN_KM_S = 8.0
VSN_KM_S = 4.6
CRUST_THICKNESS_KM = 35.0  # ความหนาเปลือกโลก (เวลาเดินทางแบบ 2 ชั้น: Pg/Pn, Sg/Sn)  ค่าที่ปรับเทียบกับข้อมูล ไม่ใช่ค่าวัด:
                           # คลื่นหัว (2H − h)·q ที่ h = 10 km เท่าสูตรเดิม 2·30·q  ผลที่ปรับเทียบไว้จึงไม่เปลี่ยน

# ---------------- หาทิศทางด้วย PCA (back-azimuth) ----------------
BAZ_LOCAL = dict(window_sec=1.0, band=(0.5, 10.0))      # สถานีใกล้: คลื่น P ชัด ใช้ 1 วินาทีแรก
BAZ_REGIONAL = dict(window_sec=5.0, band=(0.1, 0.5))    # สถานีไกล: ใช้ความถี่ต่ำ 5 วินาทีแรก

# ---------------- การคำนวณขนาด (Magnitude) ----------------
PD_TAUC_WINDOW_SEC = 3.0   # ใช้คลื่น P ช่วงแรก 3 วินาที (Wu & Kanamori)
DISP_HIGHPASS_HZ = 0.075   # high-pass หลัง integrate เป็น displacement

# ---------------- การเตือนภัย (ตามเกณฑ์ JMA) ----------------
WARNING_INTENSITY = 4.5    # ทำนายความรุนแรง JMA >= 5- (4.5) -> ออก "Warning"
FORECAST_MAGNITUDE = 3.5   # M >= 3.5 หรือ Intensity >= 3 -> ออก "Forecast"
FORECAST_INTENSITY = 2.5
DEFAULT_VS30 = 400.0       # ค่า Vs30 ของพื้นที่ผู้ใช้ (m/s) ถ้าไม่ทราบ
