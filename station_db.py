"""
station_db.py
อ่านพิกัดสถานีจากไฟล์ stations.csv (คอลัมน์ code, lat, lon)
เพื่อให้โปรแกรมหาพิกัดเองจากชื่อสถานีใน header ของไฟล์ข้อมูล (# Station: BKKA)

เพิ่มสถานีใหม่: เปิด stations.csv แล้วเพิ่มบรรทัด  CODE,lat,lon,table
"""
import csv
import os

STATIONS_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stations.csv")


def load_station_db(path=STATIONS_CSV):
    """คืนค่า dict {code: (lat, lon)}"""
    db = {}
    if not os.path.exists(path):
        return db
    with open(path, encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            try:
                db[row["code"].strip().upper()] = (float(row["lat"]), float(row["lon"]))
            except (KeyError, ValueError):
                continue                    # ข้ามแถวที่ไม่มีพิกัด
    return db


def get_coords(code, path=STATIONS_CSV):
    """หาพิกัดจากรหัสสถานี ถ้าไม่พบคืนค่า None"""
    return load_station_db(path).get(str(code).strip().upper())


def get_vs30(code, path=STATIONS_CSV):
    """Vs30 (m/s) ของสถานีจากคอลัมน์ vs30 ใน stations.csv (เติมได้ด้วย vs30_tool.py) ไม่มีคืนค่า None"""
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            if row["code"].strip().upper() == str(code).strip().upper():
                try:
                    return float(row.get("vs30") or "")
                except ValueError:
                    return None
    return None
