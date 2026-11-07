"""
run_network.py
หาจุดศูนย์กลางแผ่นดินไหวจากหลายสถานี โดยใช้ทั้งคลื่น P และคลื่น S (Sg/Sn, Lg)
ขั้นตอนทั้งหมดอยู่ใน pipeline.py (แบบเครือข่าย USGS/ANSS):
  auto-pick P -> association -> P + S scan -> re-pick P -> pick S/Lg รอบเวลาที่คาด
  -> Geiger (ทรงกลม + โลก 2 ชั้น) จาก P + S + Lg -> ตรวจทาน residual -> ทำซ้ำจนนิ่ง

รัน:  โยนไฟล์ข้อมูลเข้าโฟลเดอร์ station แล้วกด Run ใน VS Code
ผลที่ได้:
  - ตารางรายสถานีใน Terminal (pick, residual, ใช้/ไม่ใช้)
  - network_location.png  ตำแหน่งสถานีและจุดศูนย์กลาง
  - network_picks.png     กราฟตรวจ pick ของแต่ละสถานี (ดูว่าโปรแกรม pick ตรงจุดที่ควรเป็นไหม)
  - network_picks.csv     ตาราง pick ทั้งหมด
  - network_map.html      แผนที่จริง (ต้องติดตั้ง folium)
"""
import csv
import os

import numpy as np

import pipeline as PL
from data_loader import DATA_DIR
from signal_utils import destination_point
from station_db import get_coords

HERE = os.path.dirname(os.path.abspath(__file__))

# ======================================================================
#  ตั้งค่า
# ======================================================================
# "auto" = ใช้ไฟล์ .txt ทุกไฟล์ที่มี "# Station:" ในโฟลเดอร์ station
# ชื่อโฟลเดอร์ = ใช้ทุกไฟล์ในโฟลเดอร์นั้น เช่น "2026-05-18_M5.2_AllTMD" (อยู่ในโฟลเดอร์โค้ด) หรือ path เต็ม
# หรือระบุเป็น list เช่น ["BKKA", "CMMT"], [{"file": "BKKA_20250328_061959.txt"}, {"file": "X.txt", "lat": 15.0, "lon": 100.0}]
STATIONS = "2025-08-21_M4.6_AllTMD"
PROFILE = "regional"      # "regional" = สถานีไกล (> ~300 km), "local" = สถานีใกล้

USE_S = True              # True = ใช้ทั้งคลื่น P และ S / Lg,  False = คลื่น P อย่างเดียว
FIX_DEPTH_KM = 10.0       # ตรึงความลึก (km)  None = หาความลึกด้วย (ต้องมีสถานีใกล้จุดศูนย์กลาง)

# จุดศูนย์กลางอ้างอิงไว้เทียบ (lat, lon, "ชื่อ")  ไม่ใช้ใส่ None
REF_EPICENTER = (21.682, 96.121, "Reference")
REF_ORIGIN_TIME = "2025-03-28T06:20:52"      # เวลาเกิดอ้างอิง (UTC)  ไม่ใช้ใส่ None

MAKE_MAP = True           # แผนที่จริง network_map.html (ต้องติดตั้ง: pip install folium)
PLOT = True               # กราฟ PNG
SHOW_RINGS = False        # True = วาดวงกลมระยะจาก S-P / Lg-P ของแต่ละสถานีในกราฟ PNG
SAVE_CSV = True

# ค่าละเอียดของ pipeline (ไม่ต้องแก้ก็ได้) ดูความหมายใน pipeline.DEFAULTS เช่น
#   {"s_win_sec": 15.0, "scan_s_weight": 3.0, "s_phases": ("S", "Lg")}
PIPELINE_SETTINGS = {}


def find_data_files(folder=None):
    """หาไฟล์ .txt ที่เป็นข้อมูลสถานี (บรรทัดแรกขึ้นต้นด้วย '# Station:') ในโฟลเดอร์ station"""
    if folder is None:
        folder = DATA_DIR if os.path.isdir(DATA_DIR) else HERE
    print(f"ค้นหาไฟล์ข้อมูลใน : {folder}")
    files = []
    for fn in sorted(os.listdir(folder)):
        if not fn.lower().endswith(".txt"):
            continue
        try:
            with open(os.path.join(folder, fn), encoding="latin-1") as f:
                if f.readline().startswith("# Station:"):
                    files.append({"file": os.path.join(folder, fn)})
        except OSError:
            pass
    return files


def save_csv(result, path="network_picks.csv"):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["station", "lat", "lon", "dist_km", "phase", "phase_name", "time_utc", "label",
                    "weight", "residual_s", "used", "note"])
        for s in sorted(result["stations"], key=lambda q: q["dist_km"]):
            for ph in ("P", "S", "Lg"):
                if ph in s["picks"]:
                    pk = s["picks"][ph]
                    w.writerow([s["code"], s["lat"], s["lon"], round(s["dist_km"], 1), ph,
                                s["phase_name"].get(ph, ph), PL.utc(pk["time"]), pk["label"], pk["weight"],
                                round(s["residual"][ph], 2), ph in s["use"], s["notes"].get(ph, "")])
    return path


def plot_location(result, out="network_location.png"):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    hypo, sts = result["hypo"], result["stations"]
    fig, ax = plt.subplots(figsize=(9, 9))
    th = np.linspace(0, 360, 361)
    for s in sts:
        used = bool(s["use"])
        ax.plot(s["lon"], s["lat"], "^", color="k" if used else "0.7", ms=7 if len(sts) > 20 else 10)
        if len(sts) <= 30:
            ax.annotate(s["code"], (s["lon"], s["lat"]), textcoords="offset points", xytext=(5, 5), fontsize=8)
        if SHOW_RINGS and s.get("delta_sp"):
            ring = [destination_point(s["lat"], s["lon"], a, s["delta_sp"]) for a in th]
            ax.plot([q[1] for q in ring], [q[0] for q in ring], "--", lw=0.6, alpha=0.6)
    ax.plot([], [], "k^", label="Stations used")
    ax.plot([], [], "^", color="0.7", label="Stations not used")
    if result.get("hypo_p"):
        hp = result["hypo_p"]
        ax.plot(hp["lon"], hp["lat"], "o", mfc="none", mec="g", ms=12, label="P only")
    ax.plot(hypo["lon"], hypo["lat"], "r*", ms=20, label="P + S (final)")
    if REF_EPICENTER is not None:
        ax.plot(REF_EPICENTER[1], REF_EPICENTER[0], "b*", ms=14, label=REF_EPICENTER[2])
    pts = [(s["lat"], s["lon"]) for s in sts] + [(hypo["lat"], hypo["lon"])]
    if REF_EPICENTER is not None:
        pts.append(REF_EPICENTER[:2])
    la = [p[0] for p in pts]; lo = [p[1] for p in pts]
    ax.set_xlim(min(lo) - 2, max(lo) + 2); ax.set_ylim(min(la) - 2, max(la) + 2)
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_aspect("equal")
    ax.grid(alpha=0.3); ax.legend(loc="lower left")
    nb = hypo.get("n_by_phase", {})
    ax.set_title(f"Network location: {sum(1 for s in sts if s['use'])}/{len(sts)} stations used "
                 f"(P {nb.get('P', 0)}, S {nb.get('S', 0)}, Lg {nb.get('Lg', 0)})\n"
                 f"lat {hypo['lat']:.3f}, lon {hypo['lon']:.3f}, depth {hypo['depth']:.1f} km, "
                 f"RMS {hypo['rms']:.2f} s")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def main(stations=STATIONS, profile=PROFILE, plot=PLOT):
    if stations == "auto":
        stations = find_data_files()
    else:
        stations = PL.expand_specs(stations)          # ชื่อโฟลเดอร์ / ชื่อไฟล์ / รหัสสถานี -> รายการไฟล์
    print(f"Profile : {profile}   จำนวนไฟล์ : {len(stations)}   "
          f"ใช้คลื่น : {'P + S / Lg' if USE_S else 'P อย่างเดียว'}\n")
    settings = dict(PIPELINE_SETTINGS, use_s=USE_S, fix_depth=FIX_DEPTH_KM)
    result = PL.locate_event(stations, profile, get_coords, settings=settings)
    PL.print_report(result, REF_EPICENTER, REF_ORIGIN_TIME)

    if SAVE_CSV:
        print(f"\nบันทึกตาราง pick : {save_csv(result)}")
    if plot:
        out = plot_location(result)
        if out:
            print(f"บันทึกกราฟ       : {out}")
        out = PL.plot_record_section(result)
        if out:
            print(f"บันทึกกราฟ pick  : {out}  (เปิดดูว่าโปรแกรม pick ตรงจุดที่คลื่นเปลี่ยนจริงไหม)")
    if MAKE_MAP:
        try:
            from map_plot import make_multi_station_map
            hypo = result["hypo"]
            rings = [{"code": s["code"], "lat": s["lat"], "lon": s["lon"], "delta": s.get("delta_sp")}
                     for s in result["stations"] if s["use"]]
            html = make_multi_station_map(
                rings, "network_map.html", geiger=(hypo["lat"], hypo["lon"], hypo["depth"]),
                ref_epicenter=REF_EPICENTER,
                subtitle="วงกลม = ระยะจาก S-P หรือ Lg-P ของแต่ละสถานี, เป้า = จุดศูนย์กลางจาก P + S")
            if html:
                print(f"บันทึกแผนที่จริง : {html}  (เปิดด้วย browser หรือพิมพ์ start {html})")
        except Exception as ex:                          # noqa: BLE001  แผนที่ไม่ควรทำให้ผลหลักหาย
            print(f"(วาดแผนที่ไม่สำเร็จ: {ex})")
    return result


if __name__ == "__main__":
    main()