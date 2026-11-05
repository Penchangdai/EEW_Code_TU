"""
signal_utils.py
ฟังก์ชันพื้นฐานที่ไฟล์อื่นใช้ร่วมกัน: กรองสัญญาณ, อินทิเกรต, คำนวณระยะทางบนผิวโลก
"""
import numpy as np
from scipy.signal import butter, sosfilt, sosfiltfilt

from config import EARTH_RADIUS_KM


# ============================ Signal processing ============================
def remove_offset(x, n_pre=None):
    """ลบค่า offset (ค่าเฉลี่ย) ออก ถ้ากำหนด n_pre จะใช้เฉพาะช่วง noise ก่อนเหตุการณ์"""
    x = np.asarray(x, dtype=float)
    ref = x[:n_pre] if n_pre else x
    return x - np.mean(ref)


def bandpass(x, fs, fmin=None, fmax=None, order=4, causal=True):
    """Butterworth filter (causal = ใช้ได้แบบ real-time)"""
    nyq = fs / 2.0
    if fmin and fmax:
        sos = butter(order, [fmin / nyq, fmax / nyq], btype="band", output="sos")
    elif fmin:
        sos = butter(order, fmin / nyq, btype="high", output="sos")
    elif fmax:
        sos = butter(order, fmax / nyq, btype="low", output="sos")
    else:
        return np.asarray(x, dtype=float)
    return sosfilt(sos, x) if causal else sosfiltfilt(sos, x)


def integrate(x, fs):
    """อินทิเกรตเชิงตัวเลขแบบสี่เหลี่ยมคางหมู  y[i] = y[i-1] + (x[i]+x[i-1])/(2 fs)"""
    x = np.asarray(x, dtype=float)
    y = np.zeros_like(x)
    y[1:] = np.cumsum((x[1:] + x[:-1]) / (2.0 * fs))
    return y


def differentiate(x, fs):
    return np.gradient(np.asarray(x, dtype=float), 1.0 / fs)


def acc_to_vel(acc, fs, highpass_hz=0.075, n_pre=None):
    """ความเร่ง (gal = cm/s^2) -> ความเร็ว (cm/s)"""
    v = integrate(remove_offset(acc, n_pre), fs)
    return bandpass(v, fs, fmin=highpass_hz)


def acc_to_disp(acc, fs, highpass_hz=0.075, n_pre=None):
    """ความเร่ง (gal) -> การกระจัด (cm) : อินทิเกรต 2 ครั้ง + high-pass กัน drift"""
    v = acc_to_vel(acc, fs, highpass_hz, n_pre)
    d = integrate(v, fs)
    return bandpass(d, fs, fmin=highpass_hz)


# ============================ Geodesy ============================
def haversine_km(lat1, lon1, lat2, lon2):
    """ระยะทางบนผิวโลก (epicentral distance) ด้วยสูตร Haversine"""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dlmb = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))


def azimuth_deg(lat1, lon1, lat2, lon2):
    """มุมทิศ (azimuth) จากจุด 1 ไปจุด 2 วัดตามเข็มนาฬิกาจากทิศเหนือ (องศา)"""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    x = np.sin(dl) * np.cos(p2)
    y = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return (np.degrees(np.arctan2(x, y)) + 360.0) % 360.0


def destination_point(lat, lon, az_deg, dist_km):
    """หาพิกัดปลายทางเมื่อเดินจาก (lat, lon) ไปทิศ az เป็นระยะ dist_km"""
    d = dist_km / EARTH_RADIUS_KM
    az = np.radians(az_deg)
    p1, l1 = np.radians(lat), np.radians(lon)
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(az))
    l2 = l1 + np.arctan2(np.sin(az) * np.sin(d) * np.cos(p1),
                         np.cos(d) - np.sin(p1) * np.sin(p2))
    return np.degrees(p2), (np.degrees(l2) + 540.0) % 360.0 - 180.0


def latlon_to_xy(lat, lon, lat0, lon0):
    """แปลง lat/lon เป็นพิกัดระนาบ x (ตะวันออก), y (เหนือ) หน่วย km รอบจุดอ้างอิง"""
    x = EARTH_RADIUS_KM * np.radians(np.asarray(lon) - lon0) * np.cos(np.radians(lat0))
    y = EARTH_RADIUS_KM * np.radians(np.asarray(lat) - lat0)
    return x, y


def xy_to_latlon(x, y, lat0, lon0):
    lat = lat0 + np.degrees(y / EARTH_RADIUS_KM)
    lon = lon0 + np.degrees(x / (EARTH_RADIUS_KM * np.cos(np.radians(lat0))))
    return lat, lon
