"""
travel_time.py
คำนวณเวลาที่คลื่นเดินทางจากจุดกำเนิดแผ่นดินไหวมาถึง "ตำแหน่งของเรา"
และเวลาเตือนภัยล่วงหน้า (Lead time / Warning time)
"""
import numpy as np

import config as C
from signal_utils import haversine_km, azimuth_deg


def hypocentral_distance(epi_lat, epi_lon, depth_km, site_lat, site_lon):
    """R = sqrt(Δ^2 + h^2)"""
    delta = haversine_km(epi_lat, epi_lon, site_lat, site_lon)
    return float(np.sqrt(delta ** 2 + depth_km ** 2)), float(delta)


def p_arrival_time(origin_time, R, vp=C.VP_KM_S):
    """Tp = T0 + R / Vp"""
    return origin_time + R / vp


def s_arrival_time(origin_time, R, vs=C.VS_KM_S):
    """Ts = T0 + R / Vs"""
    return origin_time + R / vs


def blind_zone_radius(alert_time, origin_time, depth_km, vs=C.VS_KM_S):
    """
    รัศมีพื้นที่ที่คลื่น S มาถึงก่อนการเตือน (Blind zone)
      R_b = Vs * (T_alert - T0),   Δ_b = sqrt(R_b^2 - h^2)
    """
    Rb = vs * (alert_time - origin_time)
    return float(np.sqrt(max(Rb ** 2 - depth_km ** 2, 0.0)))


def site_warning(hypo, site_lat, site_lon, alert_time, vp=C.VP_KM_S, vs=C.VS_KM_S):
    """
    hypo : dict {lat, lon, depth, origin_time}
    alert_time : เวลาที่ระบบส่งคำเตือนออกไป
    return: ระยะทาง, เวลาคลื่น P/S ถึง, เวลาเตือนล่วงหน้า (lead time) ของคลื่น S
    """
    R, delta = hypocentral_distance(hypo["lat"], hypo["lon"], hypo["depth"], site_lat, site_lon)
    t0 = hypo["origin_time"]
    tp = p_arrival_time(t0, R, vp)
    ts = s_arrival_time(t0, R, vs)
    return {
        "epi_dist_km": delta,
        "hypo_dist_km": R,
        "azimuth_from_epicenter": float(azimuth_deg(hypo["lat"], hypo["lon"], site_lat, site_lon)),
        "p_travel_s": R / vp,
        "s_travel_s": R / vs,
        "p_arrival": tp,
        "s_arrival": ts,
        "lead_time_s": ts - alert_time,          # > 0 = มีเวลาเตรียมตัว, <= 0 = อยู่ใน blind zone
        "blind_zone_km": blind_zone_radius(alert_time, t0, hypo["depth"], vs),
    }
