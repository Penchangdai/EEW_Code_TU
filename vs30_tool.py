"""
vs30_tool.py
ดึงค่า Vs30 จากแผนที่ Vs30 ของ USGS ตามพิกัดสถานี แล้วเขียนลงคอลัมน์ vs30 ใน stations.csv

ขั้นตอน:
  1) ดาวน์โหลดแผนที่ Vs30 จาก USGS (ทำครั้งเดียว) เลือกแบบใดแบบหนึ่ง
       ก. ไฟล์ทั่วโลก (GMT grid .grd ในไฟล์ ZIP ~631 MB)  จาก https://earthquake.usgs.gov/data/vs30/
       ข. ตัดเฉพาะพื้นที่ (แนะนำ ไฟล์เล็ก) จาก https://apps.usgs.gov/earthquakes/vs30/
          เลือกกรอบประมาณ lat 5-21, lon 97-106 (ครอบประเทศไทย) แล้วดาวน์โหลดเป็น GMT grid หรือ ASCII grid
  2) วางไฟล์ไว้ในโฟลเดอร์โค้ด แล้วแก้ VS30_FILE ด้านล่าง
  3) กด Run  -> stations.csv จะมีคอลัมน์ vs30 ให้ run_intensity.py ใช้อัตโนมัติ

รองรับไฟล์: .grd / .nc (netCDF, ต้อง pip install netCDF4) และ .asc / .txt (ESRI ASCII grid)
ค่าจากแผนที่นี้ประมาณจากความชันภูมิประเทศ (และแผนที่ภูมิภาคถ้ามี) ไม่ใช่ค่าที่วัดจริง
"""
import csv
import os

import numpy as np

from station_db import STATIONS_CSV

# ======================================================================
VS30_FILE = "global_vs30.grd"      # ชื่อไฟล์แผนที่ Vs30 ที่ดาวน์โหลดมา
OVERWRITE = False                  # False = ไม่ทับค่า vs30 ที่กรอกไว้แล้วใน stations.csv
# ======================================================================


class Vs30Grid:
    """อ่านแผนที่ Vs30 แล้วหาค่าที่พิกัดใด ๆ (ค่าของ cell ที่ใกล้ที่สุด)"""

    def __init__(self, path):
        ext = os.path.splitext(path)[1].lower()
        if ext in (".asc", ".txt"):
            self._read_ascii(path)
        else:
            self._read_netcdf(path)

    def _read_netcdf(self, path):
        try:
            from netCDF4 import Dataset
        except ImportError:
            raise SystemExit("อ่านไฟล์ .grd ต้องติดตั้งก่อน: pip install netCDF4")
        self._ds = Dataset(path)
        names = self._ds.variables.keys()
        xname = next(n for n in ("lon", "x", "longitude") if n in names)
        yname = next(n for n in ("lat", "y", "latitude") if n in names)
        zname = next(n for n in names if n not in (xname, yname))
        self.lon = np.asarray(self._ds.variables[xname][:], float)
        self.lat = np.asarray(self._ds.variables[yname][:], float)
        self._z = self._ds.variables[zname]            # อ่านเฉพาะจุดที่ต้องการ ไม่โหลดทั้งไฟล์

    def _read_ascii(self, path):
        head = {}
        with open(path) as f:
            for _ in range(6):
                k, v = f.readline().split()[:2]
                head[k.lower()] = float(v)
        nc, nr, cs = int(head["ncols"]), int(head["nrows"]), head["cellsize"]
        x0 = head.get("xllcorner", head.get("xllcenter", 0) - cs / 2)
        y0 = head.get("yllcorner", head.get("yllcenter", 0) - cs / 2)
        z = np.loadtxt(path, skiprows=6)
        nod = head.get("nodata_value")
        if nod is not None:
            z = np.where(z == nod, np.nan, z)
        self.lon = x0 + cs / 2 + cs * np.arange(nc)
        self.lat = (y0 + cs / 2 + cs * np.arange(nr))
        self._z = z[::-1, :]                           # ไฟล์ ASCII เรียงแถวจากเหนือลงใต้

    def value(self, lat, lon):
        """Vs30 (m/s) ที่พิกัด (lat, lon) หรือ None ถ้าอยู่นอกแผนที่"""
        if not (min(self.lat) <= lat <= max(self.lat) and min(self.lon) <= lon <= max(self.lon)):
            return None
        i = int(np.argmin(np.abs(self.lat - lat)))
        j = int(np.argmin(np.abs(self.lon - lon)))
        v = float(np.ma.filled(self._z[i, j], np.nan))
        return None if not np.isfinite(v) or v <= 0 else v


def site_class(vs30):
    """ประเภทดินตาม NEHRP"""
    for lim, c in ((180, "E"), (360, "D"), (760, "C"), (1500, "B")):
        if vs30 < lim:
            return c
    return "A"


def main(vs30_file=VS30_FILE, csv_path=STATIONS_CSV, overwrite=OVERWRITE):
    if not os.path.isabs(vs30_file):
        vs30_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), vs30_file)
    if not os.path.exists(vs30_file):
        print(f"ไม่พบไฟล์ {vs30_file}\nดาวน์โหลดแผนที่ Vs30 จาก USGS ก่อน (ดูวิธีที่ด้านบนของไฟล์นี้)")
        return
    grid = Vs30Grid(vs30_file)

    with open(csv_path, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    fields = list(rows[0].keys())
    if "vs30" not in fields:
        fields.append("vs30")

    n_new = n_skip = n_out = 0
    print(f"{'Station':<8}{'Lat':>10}{'Lon':>11}{'Vs30 (m/s)':>12}  ประเภทดิน")
    for r in rows:
        if r.get("vs30") and not overwrite:
            n_skip += 1
            continue
        v = grid.value(float(r["lat"]), float(r["lon"]))
        if v is None:
            r["vs30"] = r.get("vs30", "")
            n_out += 1
            continue
        r["vs30"] = f"{v:.0f}"
        n_new += 1
        print(f"{r['code']:<8}{float(r['lat']):10.4f}{float(r['lon']):11.4f}{v:12.0f}  {site_class(v)}")

    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nเขียน vs30 ลง {os.path.basename(csv_path)} แล้ว {n_new} สถานี"
          f" (มีค่าอยู่แล้ว {n_skip}, อยู่นอกแผนที่ {n_out})")


if __name__ == "__main__":
    main()
