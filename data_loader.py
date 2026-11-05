"""
data_loader.py
อ่านไฟล์ข้อมูลสถานี (.txt ที่แปลงจาก miniSEED) รูปแบบ:
  # Station: BKKA
  # Sampling Rate (Hz): 200.00
  # Unit: m/s²
  DateTime  Acc_X(E)  Acc_Y(N)  Acc_Z
แล้วแปลงเป็นรูปแบบเดียวกับที่ระบบ EEW ใช้ (หน่วย gal = cm/s²)

รองรับข้อมูล 2 ชนิด แล้วแปลงเป็น "ความเร่ง" ให้เหมือนกันทั้งหมด
  - Accelerometer (สถานีวัดอัตราเร่งของพื้นดิน)  หน่วย m/s²  ช่อง HN?, EN? ...  -> ใช้ได้เลย
  - Seismometer   (สถานีตรวจวัดแผ่นดินไหว)       หน่วย m/s   ช่อง HH?, BH?, EH? -> differentiate เป็นความเร่ง
"""
import os
from datetime import datetime, timezone

import numpy as np

from signal_utils import differentiate


CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(CODE_DIR, "station")     # โฟลเดอร์เก็บไฟล์ข้อมูลสถานี


def resolve_data_path(path):
    """หาไฟล์ข้อมูล: ตามที่ระบุ -> ในโฟลเดอร์ station -> ในโฟลเดอร์โค้ด"""
    if os.path.isabs(path) or os.path.exists(path):
        return path
    for folder in (DATA_DIR, CODE_DIR):
        cand = os.path.join(folder, path)
        if os.path.exists(cand):
            return cand
    return path


def find_station_file(code, folders=None):
    """
    หาไฟล์ข้อมูลของสถานีจากรหัส (ดูบรรทัด '# Station: <code>' ใน header)
    ค้นในโฟลเดอร์ station ก่อน แล้วโฟลเดอร์โค้ด  คืนค่า path หรือ None ถ้าไม่พบ
    """
    code = code.strip().upper()
    for folder in folders or (DATA_DIR, CODE_DIR):
        if not os.path.isdir(folder):
            continue
        for fn in sorted(os.listdir(folder)):
            if not fn.lower().endswith(".txt"):
                continue
            path = os.path.join(folder, fn)
            try:
                with open(path, encoding="latin-1") as f:
                    first = f.readline()
            except OSError:
                continue
            if first.startswith("# Station:") and first.split(":", 1)[1].strip().upper() == code:
                return path
    return None


def _parse_time(s):
    return datetime.fromisoformat(s.strip().replace("Z", "")).replace(tzinfo=timezone.utc)


def load_station_txt(path):
    header = {}
    n_header = 0
    with open(path, encoding="latin-1") as f:
        for line in f:
            n_header += 1
            if line.startswith("#"):
                if ":" in line:
                    k, v = line[1:].split(":", 1)
                    header[k.strip()] = v.strip()
            else:                       # บรรทัดชื่อคอลัมน์ DateTime Acc_X Acc_Y Acc_Z
                break

    data = np.loadtxt(path, skiprows=n_header, usecols=(1, 2, 3), encoding="latin-1")
    fs = float(header["Sampling Rate (Hz)"])
    start = _parse_time(header["Start Time"])

    sensor = detect_sensor_type(header)
    e, n, z = data[:, 0] * 100.0, data[:, 1] * 100.0, data[:, 2] * 100.0   # m -> cm
    if sensor == "velocity":            # cm/s -> cm/s² (gal)
        e, n, z = differentiate(e, fs), differentiate(n, fs), differentiate(z, fs)

    return {
        "station": header.get("Station", "UNKNOWN"),
        "start_datetime": start,
        "t0": start.timestamp(),        # epoch seconds ของ sample แรก
        "fs": fs,
        "sensor": sensor,               # "acceleration" หรือ "velocity" (ข้อมูลต้นฉบับ)
        "e": e,                         # คอลัมน์ที่ 1 (ตะวันออก)  หน่วย gal
        "n": n,                         # คอลัมน์ที่ 2 (เหนือ)     หน่วย gal
        "z": z,                         # คอลัมน์ที่ 3 (แนวดิ่ง)   หน่วย gal
        "header": header,
    }


def detect_sensor_type(header):
    """
    แยกชนิดข้อมูลจาก header
      1) หน่วย: "m/s²" (มีอะไรต่อท้าย m/s) = ความเร่ง,  "m/s" เฉย ๆ = ความเร็ว
      2) ถ้าหน่วยไม่ชัด ดูตัวอักษรที่ 2 ของชื่อช่อง (SEED): N = accelerometer, H/L = seismometer
    """
    unit = header.get("Unit", "").strip().lower().replace(" ", "")
    if unit.startswith("m/s"):
        return "velocity" if unit == "m/s" else "acceleration"
    for key in ("Used for Acc_Z", "Used for Vel_Z", "Channels"):
        ch = header.get(key, "").split(",")[0].strip().upper()
        if len(ch) == 3:
            if ch[1] == "N":
                return "acceleration"
            if ch[1] in ("H", "L"):
                return "velocity"
    return "acceleration"
