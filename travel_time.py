"""
travel_time.py
คำนวณเวลาที่คลื่นเดินทางจากจุดกำเนิดแผ่นดินไหวมาถึง "ตำแหน่งของเรา"
และเวลาเตือนภัยล่วงหน้า (Lead time / Warning time)

ใช้แบบจำลองโลก 2 ชั้นเดียวกับที่ใช้หาตำแหน่ง (epicenter.py / pipeline.py)
  P, S : t = min( √(Δ² + h²)/V1 ,  Δ/V2 + 2H·√(1/V1² − 1/V2²) )      (Pg/Pn, Sg/Sn)
  Lg   : t = √(Δ² + h²) / 3.5      (ที่ระยะไกลหลายร้อย km แรงสั่นสูงสุดมักมากับคลื่น Lg)
  V1 = 6.0 / 3.5 km/s (เปลือกโลก), V2 = 8.0 / 4.6 km/s (ใต้ Moho), H = 30 km   -> แก้ใน config.py

ใช้แทนไฟล์เดิมได้: ชื่อฟังก์ชันและ key ของผลลัพธ์เหมือนเดิม
  ถ้าส่ง vp / vs เอง (แบบเดิม) จะคำนวณแบบโลกชั้นเดียว R / V เหมือนเดิม
  ถ้าไม่ส่ง (ค่าเริ่มต้น) จะใช้โลก 2 ชั้น
"""
import numpy as np

import config as C
from epicenter import phase_travel_time
from signal_utils import haversine_km, azimuth_deg

DEFAULT_DEPTH_KM = 10.0
P_WINDOW_SEC = getattr(C, "PD_TAUC_WINDOW_SEC", 3.0)    # ต้องรอคลื่น P กี่วินาทีก่อนประมวลผล
PROCESSING_DELAY_S = getattr(C, "PROCESSING_DELAY_S", 1.0)


# ======================================================================
#  ระยะทาง และเวลาเดินทาง
# ======================================================================
def hypocentral_distance(epi_lat, epi_lon, depth_km, site_lat, site_lon):
    """R = sqrt(Δ^2 + h^2)  ->  (R, Δ)"""
    delta = haversine_km(epi_lat, epi_lon, site_lat, site_lon)
    return float(np.sqrt(delta ** 2 + depth_km ** 2)), float(delta)


def _delta_from_R(R, depth_km):
    return float(np.sqrt(max(R ** 2 - depth_km ** 2, 0.0)))


def travel_time(delta_km, depth_km, phase="P"):
    """เวลาเดินทาง (s) ของเฟส "P", "S" หรือ "Lg" (โลก 2 ชั้น)"""
    return float(phase_travel_time(delta_km, depth_km, phase))


def phase_name(delta_km, depth_km, phase):
    """คลื่นที่มาถึงก่อนเป็นคลื่นตรง (Pg/Sg) หรือคลื่นใต้เปลือกโลก (Pn/Sn)"""
    if phase == "Lg":
        return "Lg"
    v1, v2 = (C.VP_KM_S, C.VPN_KM_S) if phase == "P" else (C.VS_KM_S, C.VSN_KM_S)
    direct = np.sqrt(delta_km ** 2 + depth_km ** 2) / v1
    head = delta_km / v2 + 2 * C.CRUST_THICKNESS_KM * np.sqrt(1 / v1 ** 2 - 1 / v2 ** 2)
    return phase + ("g" if direct <= head else "n")


def p_arrival_time(origin_time, R, vp=None, depth_km=DEFAULT_DEPTH_KM):
    """Tp = T0 + tP(Δ, h)   (vp ใส่เอง = แบบเดิม T0 + R / Vp)"""
    if vp is not None:
        return origin_time + R / vp
    return origin_time + travel_time(_delta_from_R(R, depth_km), depth_km, "P")


def s_arrival_time(origin_time, R, vs=None, depth_km=DEFAULT_DEPTH_KM):
    """Ts = T0 + tS(Δ, h)   (vs ใส่เอง = แบบเดิม T0 + R / Vs)"""
    if vs is not None:
        return origin_time + R / vs
    return origin_time + travel_time(_delta_from_R(R, depth_km), depth_km, "S")


def lg_arrival_time(origin_time, R):
    """T_Lg = T0 + R / 3.5"""
    return origin_time + R / getattr(C, "VLG_KM_S", 3.5)


# ======================================================================
#  เวลาออกคำเตือน และ Blind zone
# ======================================================================
def alert_time_from_picks(p_times, n_stations=3, p_window=P_WINDOW_SEC, processing=PROCESSING_DELAY_S):
    """
    T_alert = Tp ของสถานีที่ n (เรียงตามเวลา) + เวลารอคลื่น P + เวลาประมวลผล
    p_times : เวลาคลื่น P ของแต่ละสถานี (epoch s)   n = min(n_stations, จำนวนสถานี)
    """
    t = sorted(x for x in p_times if x is not None)
    if not t:
        return None
    return t[min(n_stations, len(t)) - 1] + p_window + processing


def blind_zone_radius(alert_time, origin_time, depth_km, vs=None):
    """
    รัศมี Δ_b (km) ที่คลื่น S มาถึงก่อนการเตือน (Blind zone)
      โลก 2 ชั้น : Δ_b = tS⁻¹(T_alert − T0)   (เปิดตารางเวลาเดินทางกลับ)
      vs ใส่เอง  : R_b = Vs (T_alert − T0),  Δ_b = √(R_b² − h²)   (แบบเดิม)
    """
    dt = alert_time - origin_time
    if vs is not None:
        Rb = vs * dt
        return float(np.sqrt(max(Rb ** 2 - depth_km ** 2, 0.0)))
    d = np.linspace(0.0, 4000.0, 8001)
    ts = phase_travel_time(d, depth_km, "S")
    return float(np.interp(dt, ts, d)) if dt > ts[0] else 0.0


# ======================================================================
#  ผลรวมของตำแหน่งหนึ่ง
# ======================================================================
def site_warning(hypo, site_lat, site_lon, alert_time, vp=None, vs=None):
    """
    hypo       : dict {lat, lon, depth, origin_time}  (เช่น ผลจาก pipeline.locate_event / run_magnitude)
    alert_time : เวลาที่ระบบส่งคำเตือนออกไป (epoch s)
    vp, vs     : None = โลก 2 ชั้น (แนะนำ),  ใส่ตัวเลข = โลกชั้นเดียวแบบเดิม
    return     : ระยะทาง, เวลาคลื่น P / S / Lg ถึง, เวลาเตือนล่วงหน้าก่อนคลื่น S และ Lg, blind zone
    """
    h = hypo.get("depth")
    h = DEFAULT_DEPTH_KM if h is None else float(h)
    R, delta = hypocentral_distance(hypo["lat"], hypo["lon"], h, site_lat, site_lon)
    t0 = hypo["origin_time"]
    tp = p_arrival_time(t0, R, vp, h)
    ts = s_arrival_time(t0, R, vs, h)
    tlg = lg_arrival_time(t0, R)
    one_layer = vp is not None or vs is not None
    return {
        "epi_dist_km": delta,
        "hypo_dist_km": R,
        "azimuth_from_epicenter": float(azimuth_deg(hypo["lat"], hypo["lon"], site_lat, site_lon)),
        "p_phase": "P" if vp is not None else phase_name(delta, h, "P"),
        "s_phase": "S" if vs is not None else phase_name(delta, h, "S"),
        "p_travel_s": tp - t0,
        "s_travel_s": ts - t0,
        "lg_travel_s": tlg - t0,
        "p_arrival": tp,
        "s_arrival": ts,
        "lg_arrival": tlg,
        "lead_time_s": ts - alert_time,          # > 0 = มีเวลาเตรียมตัว, <= 0 = อยู่ใน blind zone
        "lead_time_lg_s": tlg - alert_time,      # ก่อนคลื่น Lg (แรงสั่นสูงสุดที่ระยะไกล)
        "blind_zone_km": blind_zone_radius(alert_time, t0, h, vs),
        "model": "โลกชั้นเดียว" if one_layer else "โลก 2 ชั้น",
    }