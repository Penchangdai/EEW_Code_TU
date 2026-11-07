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
STATION_ROOT = os.path.join(CODE_DIR, "station")  # โฟลเดอร์รวมข้อมูล (มีโฟลเดอร์ย่อยแยกตามเหตุการณ์)


def event_dir(event=None):
    """โฟลเดอร์ข้อมูลของเหตุการณ์ (station/<event>)  ไม่ระบุ = ใช้ config.EVENT  ไม่มีโฟลเดอร์นั้น = station/"""
    if event is None:
        import config
        event = getattr(config, "EVENT", None)
    if event:
        for cand in (event, os.path.join(STATION_ROOT, event), os.path.join(CODE_DIR, event)):
            if os.path.isdir(cand):
                return os.path.abspath(cand)
        print(f"⚠ ไม่พบโฟลเดอร์เหตุการณ์ '{event}' ใน station/  ใช้ไฟล์ในโฟลเดอร์ station/ แทน")
    return STATION_ROOT


DATA_DIR = event_dir()     # โฟลเดอร์เก็บไฟล์ข้อมูลสถานีของเหตุการณ์ที่เลือก (config.EVENT)


def _search_folders():
    return tuple(dict.fromkeys((DATA_DIR, STATION_ROOT, CODE_DIR)))


def resolve_data_path(path):
    """หาไฟล์ข้อมูล: ตามที่ระบุ -> โฟลเดอร์เหตุการณ์ (config.EVENT) -> โฟลเดอร์ station -> โฟลเดอร์โค้ด"""
    if os.path.isabs(path) or os.path.exists(path):
        return path
    for folder in _search_folders():
        cand = os.path.join(folder, path)
        if os.path.exists(cand):
            return cand
    return path


def find_station_file(code, folders=None):
    """
    หาไฟล์ข้อมูลของสถานีจากรหัส (ดูบรรทัด '# Station: <code>' ใน header)
    ค้นในโฟลเดอร์เหตุการณ์ (config.EVENT) ก่อน แล้ว station/ และโฟลเดอร์โค้ด  คืนค่า path หรือ None ถ้าไม่พบ
    """
    code = code.strip().upper()
    for folder in folders or _search_folders():
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

    units_ok, units_note = check_units(header, n, fs)
    return {
        "units_ok": units_ok,           # False = ค่าในไฟล์ไม่ใช่หน่วยจริง (เช่น ยังเป็น count) ไม่ควรใช้
        "units_note": units_note,
        "station": header.get("Station", "UNKNOWN"),
        "station_type": detect_station_type(header),   # "seismic" หรือ "strong_motion" (ประเภทสถานี)
        "start_datetime": start,
        "t0": start.timestamp(),        # epoch seconds ของ sample แรก
        "fs": fs,
        "sensor": sensor,               # "acceleration" หรือ "velocity" (ชนิดข้อมูลในไฟล์นี้)
        "e": e,                         # คอลัมน์ที่ 1 (ตะวันออก)  หน่วย gal
        "n": n,                         # คอลัมน์ที่ 2 (เหนือ)     หน่วย gal
        "z": z,                         # คอลัมน์ที่ 3 (แนวดิ่ง)   หน่วย gal
        "header": header,
    }


# noise ก่อนเกิดเหตุของเครื่องวัดความเร่งจริงอยู่ที่ระดับ 0.001-1 gal  ถ้าค่ากลางของ |a| ช่วงแรกเกินนี้
# แปลว่าไฟล์ยังเป็นค่าดิบ (count) ที่ยังไม่ได้หารด้วย Scale Factor  (หน่วย gal)
MAX_NOISE_GAL = 20.0


def check_units(header, x_gal, fs, first_s=30.0):
    """
    ตรวจว่าค่าในไฟล์เป็นหน่วยจริงหรือไม่
      - Scale Factor: 1 และ noise สูงผิดปกติ = ยังไม่ได้แปลงจาก count  (ไฟล์ที่แปลงแล้วมี Scale Factor ~10^5-10^6)
      - noise ช่วง first_s วินาทีแรก (ค่ากลางของ |a − ค่าเฉลี่ย|) > MAX_NOISE_GAL = เป็นไปไม่ได้สำหรับ noise
    return: (ok, หมายเหตุ)
    """
    seg = np.asarray(x_gal[: max(int(first_s * fs), 10)], float)
    noise = float(np.median(np.abs(seg - seg.mean()))) if len(seg) else 0.0
    try:
        sf = float(header.get("Scale Factor", "nan"))
    except ValueError:
        sf = float("nan")
    if noise > MAX_NOISE_GAL:
        why = "Scale Factor = 1 (ยังไม่ได้แปลงจาก count)" if sf == 1 else f"Scale Factor = {header.get('Scale Factor', '?')}"
        return False, (f"ค่าในไฟล์ไม่ใช่หน่วยจริง: noise ก่อนเกิดเหตุ {noise:,.0f} gal (ปกติ < 1 gal)  {why}"
                       " ต้องแปลงไฟล์ใหม่ (หารด้วยค่าความไวของเครื่องวัด)")
    return True, ""


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


def detect_station_type(header):
    """
    ประเภทสถานี (ดูจากช่องที่สถานีมีทั้งหมดในบรรทัด '# Channels:' ตามรหัส SEED ตัวอักษรที่ 2)
      "seismic"       = สถานีตรวจวัดแผ่นดินไหว : มีเครื่องวัดแผ่นดินไหว (seismometer)  ช่อง ?H? / ?L?
                        เช่น HH? (broadband), EH? / SH? (คาบสั้น), BH?  (ส่วนใหญ่มีเครื่องวัดความเร่ง HN? ด้วย)
      "strong_motion" = สถานีตรวจวัดอัตราเร่งของพื้นดิน : มีแค่เครื่องวัดความเร่ง (accelerometer)  ช่อง ?N?
                        (ช่องอื่นเช่น HDF = ความดันอากาศ, HCA = สถานะเครื่อง ไม่นับ)
    ถ้าไม่มีบรรทัด Channels ใช้ชนิดข้อมูลในไฟล์แทน (ความเร็ว = seismic)
    """
    chans = [c.strip().upper() for c in header.get("Channels", "").split(",") if c.strip()]
    seed = [c for c in chans if len(c) == 3 and c.isalnum()]
    if any(c[1] in ("H", "L") for c in seed):
        return "seismic"
    if any(c[1] == "N" for c in seed):
        return "strong_motion"
    return "seismic" if detect_sensor_type(header) == "velocity" else "strong_motion"
