"""
run_network.py
หาจุดศูนย์กลางแผ่นดินไหวจากข้อมูลหลายสถานีรวมกัน (Geiger's method)
  1) อ่านไฟล์ของทุกสถานี แล้วหาเวลาคลื่น P และ S ของแต่ละสถานี
  2) นำเวลา P, S ของทุกสถานีมาหาตำแหน่ง (lat, lon), ความลึก และเวลาเกิด พร้อมกัน
  3) เทียบกับวิธีสถานีเดียว (Omori) ของแต่ละสถานี

ต้องมีเวลาคลื่น (P+S) รวมกันอย่างน้อย 4 ค่า  แนะนำ 3 สถานีขึ้นไป
ทุกไฟล์ต้องเป็นเหตุการณ์เดียวกัน และเวลาต้องเป็น UTC เหมือนกัน
ใส่ได้ทั้งสถานีวัดอัตราเร่ง (ACC, m/s²) และสถานีตรวจวัดแผ่นดินไหว (VEL, m/s) ปนกัน
data_loader.py แปลงเป็นความเร่งให้อัตโนมัติ

รัน:  ใส่รายชื่อไฟล์และพิกัดใน STATIONS ด้านล่าง แล้วกด Run ใน VS Code
"""
import os

import numpy as np

import config as C
from data_loader import load_station_txt, resolve_data_path, DATA_DIR
from p_wave import detect_p_wave
from s_wave import detect_s_wave
from epicenter import geiger_locate, omori_distance
from signal_utils import haversine_km, destination_point
from run_p_wave import to_utc
from station_db import get_coords
from map_plot import make_multi_station_map

HERE = os.path.dirname(os.path.abspath(__file__))


def find_data_files(folder=None):
    """
    หาไฟล์ .txt ที่เป็นข้อมูลสถานี (บรรทัดแรกขึ้นต้นด้วย '# Station:')
    ค้นในโฟลเดอร์ station ก่อน ถ้าไม่มีโฟลเดอร์นี้ จะค้นในโฟลเดอร์โค้ด
    """
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

# ======================================================================
#  เลือกไฟล์ข้อมูล
#   "auto" = ใช้ไฟล์ .txt ทุกไฟล์ที่มี "# Station:" ใน header
#            ค้นในโฟลเดอร์ "station" (ถ้ามี) ไม่อย่างนั้นค้นในโฟลเดอร์โค้ด
#            (แค่โยนไฟล์เข้าโฟลเดอร์ station แล้วกด Run)
#   หรือระบุเป็น list เช่น
#   STATIONS = [
#       {"file": "BKKA_20250328_061959.txt"},                          # พิกัดจาก stations.csv
#       {"file": "XXXX_20250328_061959.txt", "lat": 15.0, "lon": 100.0}, # ใส่พิกัดเอง
#   ]
#  พิกัดดึงจาก stations.csv ตามชื่อสถานีใน header อัตโนมัติ
# ======================================================================
STATIONS = "auto"
PROFILE = "regional"      # "regional" = สถานีไกล (> ~300 km), "local" = สถานีใกล้

# แผนที่จริง (ไฟล์ network_map.html เปิดด้วย browser) ต้องติดตั้ง: pip install folium
MAKE_MAP = True
# จุดศูนย์กลางอ้างอิงไว้เทียบบนแผนที่ (lat, lon, "ชื่อ") เช่น (22.01, 95.92, "USGS")  ไม่ใช้ใส่ None
REF_EPICENTER = None


def pick_station(st, regional):
    """อ่านไฟล์หนึ่งสถานี แล้วหาเวลา P และ S"""
    path = resolve_data_path(st["file"])
    rec = load_station_txt(path)
    fs = rec["fs"]
    if st.get("lat") is None or st.get("lon") is None:
        c = get_coords(rec["station"])
        if c is None:
            print(f"⚠ ข้ามสถานี {rec['station']}: ไม่พบพิกัดใน stations.csv")
            return None
        st = {**st, "lat": c[0], "lon": c[1]}
    p = detect_p_wave(rec["z"], fs, rec["t0"], **(C.P_REGIONAL if regional else {}))
    out = {"code": rec["station"], "lat": st["lat"], "lon": st["lon"], "rec": rec,
           "sensor": "ACC" if rec["sensor"] == "acceleration" else "VEL",
           "tp": None, "ts": None}
    if not p["detected"]:
        return out
    out["tp"] = p["time"]
    s = detect_s_wave(rec["z"], rec["n"], rec["e"], fs, p["index"], rec["t0"],
                      **(C.S_REGIONAL if regional else {}))
    if s["detected"]:
        out["ts"] = s["time"]
    return out


def main(stations=STATIONS, profile=PROFILE, plot=True):
    if stations == "auto":
        stations = find_data_files()
    regional = profile == "regional"
    vp, vs = (C.VPN_KM_S, C.VSN_KM_S) if regional else (C.VP_KM_S, C.VS_KM_S)
    print(f"Profile : {profile}   (Vp = {vp}, Vs = {vs} km/s)")
    print(f"จำนวนสถานี : {len(stations)}\n")

    # ---------- 1) หาเวลา P, S ของแต่ละสถานี ----------
    picks = []
    print(f"{'Station':<8}{'Type':<5}{'Lat':>8}{'Lon':>9}   {'P arrival (UTC)':<24}{'S arrival (UTC)':<24}{'S-P':>7}{'R_Omori':>9}")
    for st in stations:
        pk = pick_station(st, regional)
        if pk is None:
            continue
        tp = to_utc(pk["tp"]) if pk["tp"] else "-"
        ts = to_utc(pk["ts"]) if pk["ts"] else "-"
        sp = (pk["ts"] - pk["tp"]) if pk["tp"] and pk["ts"] else None
        pk["R_omori"] = omori_distance(sp, vp, vs) if sp else None
        sp_txt = f"{sp:.2f}" if sp else "-"
        r_txt = f"{pk['R_omori']:.0f}" if sp else "-"
        print(f"{pk['code']:<8}{pk['sensor']:<5}{pk['lat']:8.3f}{pk['lon']:9.3f}   {tp:<24}{ts:<24}{sp_txt:>7}{r_txt:>9}")
        if pk["tp"]:
            picks.append(pk)

    n_obs = sum(1 + (p["ts"] is not None) for p in picks)
    print(f"\nสถานีที่หาคลื่น P ได้ : {len(picks)}   จำนวนเวลาคลื่น (P+S) : {n_obs}")
    if n_obs < 4:
        print("ข้อมูลไม่พอสำหรับ Geiger (ต้องมีเวลาคลื่นอย่างน้อย 4 ค่า)")
        print("ถ้ามีสถานีเดียว ให้ใช้ run_epicenter.py แทน")
        return None
    if len(picks) < 3:
        print("⚠ มีน้อยกว่า 3 สถานี ผลอาจไม่เสถียร ควรเพิ่มสถานี")

    coords = {(round(p["lat"], 3), round(p["lon"], 3)) for p in picks}
    if len(coords) < len(picks):
        print("⚠ มีสถานีที่พิกัดซ้ำกัน ตรวจสอบค่า lat/lon ใน STATIONS (สถานีที่ตำแหน่งเดียวกันไม่ช่วยหาตำแหน่ง)")

    # ---------- 2) Geiger's method ----------
    hypo = geiger_locate(picks, vp=vp, vs=vs)
    print("\n================ ผลจากหลายสถานี (Geiger) ================")
    print(f"Epicenter   : lat {hypo['lat']:.3f}, lon {hypo['lon']:.3f}")
    print(f"Depth       : {hypo['depth']:.1f} km")
    print(f"Origin time : {to_utc(hypo['origin_time'])} UTC")
    print(f"RMS residual: {hypo['rms']:.3f} s   (ยิ่งน้อยยิ่งดี)")
    print(f"Error (1σ)  : x ±{hypo['err_x_km']:.1f} km, y ±{hypo['err_y_km']:.1f} km, "
          f"z ±{hypo['err_z_km']:.1f} km, T0 ±{hypo['err_t0_s']:.2f} s")
    if hypo["depth"] > 700 or not np.isfinite(hypo["err_x_km"]):
        print("⚠ ผลนี้ไม่น่าเชื่อถือ (ความลึกเกิน 700 km หรือคำนวณความคลาดเคลื่อนไม่ได้)")
        print("  สาเหตุที่พบบ่อย: สถานีน้อยเกินไป, พิกัดสถานีผิด/ซ้ำกัน, หรือ pick คลื่นผิด")

    print("\nระยะทางจากจุดศูนย์กลางถึงแต่ละสถานี:")
    for p in picks:
        d = haversine_km(hypo["lat"], hypo["lon"], p["lat"], p["lon"])
        r_om = f"{p['R_omori']:.0f} km" if p["R_omori"] else "-"
        print(f"  {p['code']:<8} Geiger {d:7.1f} km   |   Omori (สถานีเดียว) {r_om}")
    if regional:
        print("\n* สถานีห่างกันหลายร้อยกิโลเมตร โค้ดคำนวณบนระนาบแบน (flat-earth) จึงคลาดเคลื่อนได้บ้าง")
    if REF_EPICENTER is not None:
        d = haversine_km(REF_EPICENTER[0], REF_EPICENTER[1], hypo["lat"], hypo["lon"])
        print(f"ผล Geiger ห่างจาก {REF_EPICENTER[2]} {d:.0f} km")

    # ---------- 3) แผนที่จริง: รัศมี Omori ของแต่ละสถานี + จุดศูนย์กลางจาก Geiger ----------
    if MAKE_MAP:
        rings = []
        for p in picks:
            d = None
            if p["R_omori"]:
                d = float(np.sqrt(max(p["R_omori"] ** 2 - min(hypo["depth"], p["R_omori"]) ** 2, 0.0)))
            rings.append({"code": p["code"], "lat": p["lat"], "lon": p["lon"], "delta": d})
        html = make_multi_station_map(
            rings, "network_map.html",
            geiger=(hypo["lat"], hypo["lon"], hypo["depth"]), ref_epicenter=REF_EPICENTER,
            subtitle="วงกลม = รัศมีระยะทางจาก S-P (Omori) ของแต่ละสถานี, เป้าม่วง = จุดศูนย์กลางจาก Geiger")
        if html:
            print(f"\nบันทึกแผนที่จริง: {html}  (เปิดด้วย browser หรือพิมพ์ start {html})")

    if plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            return hypo
        fig, ax = plt.subplots(figsize=(8, 8))
        th = np.linspace(0, 360, 361)
        for p in picks:
            ax.plot(p["lon"], p["lat"], "k^", ms=10)
            ax.annotate(p["code"], (p["lon"], p["lat"]), textcoords="offset points", xytext=(6, 6))
            if p["R_omori"]:
                ring = [destination_point(p["lat"], p["lon"], a, p["R_omori"]) for a in th]
                ax.plot([q[1] for q in ring], [q[0] for q in ring], "--", lw=0.8, alpha=0.7)
        ax.plot(hypo["lon"], hypo["lat"], "r*", ms=20, label="Epicenter (Geiger)")
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_aspect("equal")
        ax.grid(alpha=0.3); ax.legend(loc="lower left")
        ax.set_title(f"Network location: {len(picks)} stations\n"
                     f"lat {hypo['lat']:.3f}, lon {hypo['lon']:.3f}, depth {hypo['depth']:.1f} km "
                     f"(dashed = Omori distance)")
        fig.tight_layout()
        fig.savefig("network_location.png", dpi=120)
        print("\nบันทึกกราฟ: network_location.png")
    return hypo


if __name__ == "__main__":
    main()
