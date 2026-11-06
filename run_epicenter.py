"""
run_epicenter.py
หาตำแหน่งจุดศูนย์กลางแผ่นดินไหว (Epicenter) จากสถานีที่เลือก 1-5 สถานี
  1) หาคลื่น P และ S ของแต่ละสถานี      -> เวลา S-P
  2) ระยะทาง (สูตร Omori)               R = k (Ts - Tp)
  3) ทิศทาง (PCA ของคลื่น P)            back-azimuth (กำกวม 180°) -> ตัวเลือก A, B
  4) หลายสถานี (LOCATE_METHOD):
     "pipeline" (ค่าเริ่มต้น, 3 สถานีขึ้นไป) : ขั้นตอนแบบ USGS ใน pipeline.py ใช้ทั้งคลื่น P และ S / Lg
                             (P + S scan -> re-pick P -> pick S/Lg รอบเวลาที่คาด -> Geiger P+S -> ตรวจทาน)
     "robust" (ค่าเริ่มต้น) : Robust grid search ใช้เวลา P, S ทุกสถานี + ทิศทาง PCA
                             ตัดสถานีที่เวลาคลื่นไม่สอดคล้องออก แล้วคำนวณใหม่
     "ab"                  : วิธีเดิม เลือกตัวเลือก A/B ที่ใกล้กันที่สุดแล้วเฉลี่ย
     SHOW_GEIGER = True     : คำนวณ Geiger (ระนาบแบน) เทียบด้วย เมื่อมี 3 สถานีขึ้นไป
     ทั้ง robust และ Geiger ถ่วงน้ำหนักเวลาคลื่นตาม weight 0-4 ของ pick (แบบ Hypoinverse)
     wP / wS ในตาราง = weight ของคลื่น P / S (0 ดีที่สุด, 4 = ไม่ใช้)
  5) วาดลงแผนที่จริง (epicenter_map_<สถานี>.html) ด้วย folium

รัน:  ใส่รหัสสถานีใน SELECTED_STATIONS ด้านล่าง แล้วกด Run ใน VS Code
      โปรแกรมหาไฟล์ข้อมูลเอง (ดู '# Station:' ใน header ของไฟล์ในโฟลเดอร์ station)
      และดึงพิกัดจาก stations.csv เอง
"""
import itertools
import os
import sys

import numpy as np

import config as C
from data_loader import load_station_txt, resolve_data_path, find_station_file
from p_wave import detect_p_wave
from s_wave import detect_s_wave
from epicenter import (omori_distance, back_azimuth_pca, geiger_locate, grid_search_locate,
                       travel_time_2layer, delta_from_sp)
from signal_utils import destination_point, haversine_km, bandpass
from station_db import get_coords
from map_plot import make_epicenter_map, make_multi_station_map
from run_p_wave import to_utc

# ======================================================================
#  ตั้งค่าตรงนี้
# ======================================================================
# รหัสสถานีที่ต้องการใช้ 1-5 สถานี (ต้องมีไฟล์ข้อมูลในโฟลเดอร์ station และพิกัดใน stations.csv)
SELECTED_STATIONS = ["BKKA", "CMMT", "PHIT", "NAN", "PHEA"]
# ตัวอย่างหลายสถานี:
# SELECTED_STATIONS = ["BKKA", "BSPA", "CMAI"]

# ถ้า SELECTED_STATIONS = [] จะใช้ไฟล์นี้ไฟล์เดียวแทน
DATA_FILE = "BKKA_20250328_061959.txt"

PROFILE = "regional"      # "regional" = สถานีไกล (> ~300 km), "local" = สถานีใกล้

# พิกัดที่ไม่มีใน stations.csv ใส่เองได้ที่นี่ เช่น {"XXXX": (15.123, 100.456)}
STATION_COORDS_OVERRIDE = {}

# วิธีแปลง S-P เป็นระยะทาง
#   "2layer" = โลก 2 ชั้น (เปลือกโลก + เนื้อโลก) ใช้ได้ทั้งใกล้และไกล  <- แนะนำ
#   "omori"  = D = k (Ts - Tp) แบบเดิม (สถานีไกลจะได้ระยะยาวเกินจริง ~60-70 km)
DISTANCE_MODEL = "2layer"
DEPTH_KM = 10.0           # ความลึกสมมติ (วิธีสถานีเดียวและ robust หาความลึกไม่ได้)

# วิธีหาจุดศูนย์กลางเมื่อมี 2 สถานีขึ้นไป
LOCATE_METHOD = "pipeline"   # "pipeline" (แนะนำ, 3 สถานีขึ้นไป), "robust" หรือ "ab" (วิธีเดิม)
USE_S = True              # pipeline: True = ใช้ทั้งคลื่น P และ S / Lg, False = คลื่น P อย่างเดียว
USE_BAZ = True            # robust: ใช้ทิศทางจาก PCA ช่วย (สำคัญเมื่อสถานีอยู่ฝั่งเดียวกันของจุดศูนย์กลาง)
SHOW_GEIGER = False       # True = คำนวณ Geiger เทียบด้วย (ระนาบแบน ไม่เหมาะกับสถานีไกลหลายร้อย km)
MAX_STATIONS = 5

# แผนที่จริง (ไฟล์ .html เปิดด้วย browser) ต้องติดตั้ง: pip install folium
MAKE_MAP = True
# จุดศูนย์กลางอ้างอิงไว้เทียบบนแผนที่ (lat, lon, "ชื่อ") เช่น (22.01, 95.92, "USGS")  ไม่ใช้ใส่ None
REF_EPICENTER = (21.682, 96.121, "อ้างอิง")
# เวลาเกิดอ้างอิง (UTC) ไว้เทียบ เช่น "2025-03-28T06:20:52"  ไม่ใช้ใส่ None
REF_ORIGIN_TIME = "2025-03-28T06:20:52"


def station_coords(code):
    if code.upper() in STATION_COORDS_OVERRIDE:
        return STATION_COORDS_OVERRIDE[code.upper()]
    return get_coords(code)


def analyze_station(path, profile, plot=True):
    """วิเคราะห์สถานีเดียว: P, S, ระยะทาง, ทิศทาง, ตัวเลือก A/B  คืนค่า dict หรือ None"""
    rec = load_station_txt(path)
    code = rec["station"]
    coords = station_coords(code)
    if coords is None:
        print(f"⚠ ข้ามสถานี {code}: ไม่พบพิกัดใน stations.csv (เพิ่มใน STATION_COORDS_OVERRIDE ได้)")
        return None
    sta_lat, sta_lon = coords
    fs = rec["fs"]
    z, n, e = rec["z"], rec["n"], rec["e"]
    regional = profile == "regional"

    print(f"\n========== สถานี {code}  ({sta_lat}, {sta_lon}) ==========")
    p = detect_p_wave(z, fs, start_time=rec["t0"], **(C.P_REGIONAL if regional else {}))
    if not p["detected"]:
        print("ไม่พบคลื่น P"); return None
    s = detect_s_wave(z, n, e, fs, p["index"], start_time=rec["t0"],
                      **(C.S_REGIONAL if regional else {}))
    if not s["detected"]:
        print("ไม่พบคลื่น S"); return None
    sp = s["sp_time"]
    print(f"P-wave arrival : {to_utc(p['time'])} UTC   {p.get('label', '')}")
    print(f"S-wave arrival : {to_utc(s['time'])} UTC   {s.get('label', '')}")
    print(f"S - P time     : {sp:.2f} s")

    if DISTANCE_MODEL == "2layer":
        # โลก 2 ชั้น: คิด intercept time ของ Pn/Sn (สูตร Omori แบบเส้นตรงให้ระยะยาวเกินจริงเมื่อสถานีไกล)
        delta = delta_from_sp(sp, DEPTH_KM)
        R = float(np.hypot(delta, DEPTH_KM))
        t0 = p["time"] - float(travel_time_2layer(delta, DEPTH_KM, "P"))
        how = "โลก 2 ชั้น"
    else:
        vp, vs = (C.VPN_KM_S, C.VSN_KM_S) if regional else (C.VP_KM_S, C.VS_KM_S)
        R = omori_distance(sp, vp, vs)
        delta = np.sqrt(max(R ** 2 - DEPTH_KM ** 2, 0.0))
        t0 = p["time"] - R / vp
        how = "Omori"
    print(f"ระยะทาง ({how}) R = {R:.0f} km, Δ = {delta:.0f} km, T0 ≈ {to_utc(t0)} UTC")

    baz_params = C.BAZ_REGIONAL if regional else C.BAZ_LOCAL
    baz, rect = back_azimuth_pca(z, n, e, fs, p["index"], **baz_params)
    baz2 = (baz + 180.0) % 360.0
    lat1, lon1 = destination_point(sta_lat, sta_lon, baz, delta)
    lat2, lon2 = destination_point(sta_lat, sta_lon, baz2, delta)
    print(f"ทิศทาง (PCA)   : {baz:.0f}° หรือ {baz2:.0f}°   (rectilinearity {rect:.2f})")
    print(f"  ตัวเลือก A (ทิศ {baz:5.1f}°) : lat {lat1:7.3f}, lon {lon1:8.3f}")
    print(f"  ตัวเลือก B (ทิศ {baz2:5.1f}°) : lat {lat2:7.3f}, lon {lon2:8.3f}")

    res = {"code": code, "lat": sta_lat, "lon": sta_lon, "tp": p["time"], "ts": s["time"],
           "sp": sp, "R": R, "delta": delta, "t0": t0, "rect": rect, "chosen": None,
           "wp": p.get("weight"), "ws": s.get("weight"),          # weight 0-4 ใช้ถ่วงน้ำหนักตอนหาตำแหน่ง
           "cands": [{"name": "ตัวเลือก A", "lat": lat1, "lon": lon1, "az": baz},
                     {"name": "ตัวเลือก B", "lat": lat2, "lon": lon2, "az": baz2}]}
    if plot:
        plot_station(rec, p, res, baz_params)
    return res


def plot_station(rec, p, res, baz_params):
    """กราฟ PNG: การสั่นของคลื่น P (ซ้าย) + ตำแหน่งที่เป็นไปได้ (ขวา)"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fs, n, e = rec["fs"], rec["n"], rec["e"]
    i0, i1 = p["index"], p["index"] + int(baz_params["window_sec"] * fs)
    nf = bandpass(n, fs, *baz_params["band"])[i0:i1]
    ef = bandpass(e, fs, *baz_params["band"])[i0:i1]
    a, b = res["cands"]

    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    ax[0].plot(ef, nf, lw=0.8)
    lim = 1.1 * max(np.max(np.abs(ef)), np.max(np.abs(nf)))
    for c, col in ((a, "r"), (b, "m")):
        ax[0].plot([0, lim * np.sin(np.radians(c["az"]))], [0, lim * np.cos(np.radians(c["az"]))],
                   col, lw=2, label=f"{c['az']:.0f}°")
    ax[0].set_xlim(-lim, lim); ax[0].set_ylim(-lim, lim); ax[0].set_aspect("equal")
    ax[0].set_xlabel("E (gal)"); ax[0].set_ylabel("N (gal)"); ax[0].legend()
    ax[0].set_title("P-wave particle motion (horizontal)")
    ax[1].plot(res["lon"], res["lat"], "k^", ms=12, label=f"Station {res['code']}")
    ax[1].plot(a["lon"], a["lat"], "r*", ms=16, label="Candidate A")
    ax[1].plot(b["lon"], b["lat"], "m*", ms=16, label="Candidate B")
    ring = [destination_point(res["lat"], res["lon"], t, res["delta"]) for t in np.linspace(0, 360, 361)]
    ax[1].plot([q[1] for q in ring], [q[0] for q in ring], "k--", lw=0.7, label=f"Δ = {res['delta']:.0f} km")
    ax[1].set_xlabel("Longitude"); ax[1].set_ylabel("Latitude")
    ax[1].legend(fontsize=8, loc="lower left"); ax[1].set_aspect("equal"); ax[1].grid(alpha=0.3)
    ax[1].set_title("Possible epicenters")
    fig.tight_layout()
    out = f"epicenter_{res['code']}.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"บันทึกกราฟ: {out}")


def choose_consistent(results):
    """
    เลือกตัวเลือก A/B ของแต่ละสถานีให้อยู่ใกล้กันที่สุด
    ลองทุกแบบ (2^n แบบ, n <= 5) แล้วเลือกแบบที่ผลรวมระยะห่างระหว่างคู่น้อยที่สุด
    """
    best, best_cost = None, np.inf
    for combo in itertools.product((0, 1), repeat=len(results)):
        pts = [r["cands"][j] for r, j in zip(results, combo)]
        cost = sum(haversine_km(p1["lat"], p1["lon"], p2["lat"], p2["lon"])
                   for p1, p2 in itertools.combinations(pts, 2))
        if cost < best_cost:
            best, best_cost = combo, cost
    for r, j in zip(results, best):
        r["chosen"] = j
    pts = [r["cands"][r["chosen"]] for r in results]
    spread = max(haversine_km(p1["lat"], p1["lon"], p2["lat"], p2["lon"])
                 for p1, p2 in itertools.combinations(pts, 2))
    return (float(np.mean([q["lat"] for q in pts])), float(np.mean([q["lon"] for q in pts]))), spread


def choose_nearest(results, epi):
    """เลือกตัวเลือก A/B ของแต่ละสถานีที่อยู่ใกล้จุดศูนย์กลางที่คำนวณได้"""
    for r in results:
        d = [haversine_km(c["lat"], c["lon"], *epi) for c in r["cands"]]
        r["chosen"] = int(np.argmin(d))


def delta_from_sp_lg(sp, depth_km=DEPTH_KM):
    """Δ จาก S−P เมื่อคลื่น S เป็น Lg:  t_Lg(Δ) − t_P(Δ) = S−P  (เปิดตารางย้อนกลับ)"""
    d = np.linspace(1, 4000, 8000)
    sp_tab = np.sqrt(d ** 2 + depth_km ** 2) / getattr(C, "VLG_KM_S", 3.5) - travel_time_2layer(d, depth_km, "P")
    return float(np.interp(sp, sp_tab, d))


def relabel_as_lg(r):
    """ตีความคลื่น S เป็น Lg: คำนวณระยะ Δ และตัวเลือก A/B ใหม่ (ใช้วาดวงกลม/เส้นบนแผนที่)"""
    d = delta_from_sp_lg(r["sp"])
    print(f"  {r['code']}: คลื่น S ที่จับได้เป็น Lg  Δ {r['delta']:.0f} → {d:.0f} km")
    r["delta"] = d
    for c in r["cands"]:
        c["lat"], c["lon"] = destination_point(r["lat"], r["lon"], c["az"], d)


def _w(code):
    return "-" if code is None else str(code)


def locate_robust(results):
    """Robust grid search: เวลา P, S ทุกสถานี (+ ทิศทาง PCA) ตัด outlier แล้วคำนวณใหม่"""
    picks = [{"code": r["code"], "lat": r["lat"], "lon": r["lon"], "tp": r["tp"], "ts": r["ts"],
              "wp": r.get("wp"), "ws": r.get("ws"),
              "baz": r["cands"][0]["az"], "rect": r["rect"]} for r in results]
    if sum(2 for _ in picks) < 4:
        print("(robust ต้องมีเวลาคลื่นอย่างน้อย 4 ค่า = 2 สถานี)"); return None
    rob = grid_search_locate(picks, use_s=True, depth_km=DEPTH_KM, use_baz=USE_BAZ)
    print(f"วิธี: {rob['method']}")
    print(f"\n{'Station':<8}{'S-P (s)':>9}{'Δ จาก S-P':>11}{'Δ ถึงผล':>10}{'res P (s)':>11}{'res S (s)':>11}{'เฟส S':>7}"
          f"{'wP':>4}{'wS':>4}  สถานะ")
    res = {(x["code"], x["phase"]): x for x in rob["residuals"]}
    for r in results:
        d = haversine_km(r["lat"], r["lon"], rob["lat"], rob["lon"])
        rp, rs = res.get((r["code"], "P")), res.get((r["code"], "S"))
        bad = [ph for ph, x in (("P", rp), ("S", rs)) if x and x["outlier"]]
        r["outlier"] = len(bad) == 2                 # สีเทาบนแผนที่ เฉพาะเมื่อถูกตัดทั้ง P และ S
        r["outlier_phases"] = bad
        r["phase_s"] = rs.get("phase_name", "S") if rs else "S"
        if r["phase_s"] == "Lg" and "S" not in bad:  # คลื่น S ที่จับได้เป็น Lg -> คำนวณ Δ และตัวเลือก A/B ใหม่
            relabel_as_lg(r)
        print(f"{r['code']:<8}{r['sp']:9.1f}{r['delta']:9.0f} km{d:8.0f} km"
              f"{rp['residual'] if rp else float('nan'):11.1f}{rs['residual'] if rs else float('nan'):11.1f}"
              f"{r['phase_s']:>7}{_w(r.get('wp')):>4}{_w(r.get('ws')):>4}  "
              + ("ตัดคลื่น " + ",".join(bad) + " ออก (เวลาไม่สอดคล้อง)" if bad else "ใช้ได้"))
    if rob["baz_residuals"]:
        print("\nทิศทางจาก PCA เทียบกับทิศไปยังผลที่ได้ (ไม่สนหน้า/หลัง 180°):")
        for x in rob["baz_residuals"]:
            print(f"  {x['code']:<6} PCA {x['baz']:5.0f}°   ทิศไปผล {x['az_to_epi']:5.0f}°   ต่าง {x['diff']:4.0f}°")
    print(f"\nจุดศูนย์กลาง (Robust) : lat {rob['lat']:.3f}, lon {rob['lon']:.3f}  (ความลึกสมมติ {DEPTH_KM:.0f} km)")
    print(f"เวลาเกิด            : {to_utc(rob['origin_time'])} UTC")
    print(f"Misfit              : {rob['misfit']:.2f} s  (truncated L1 เฉลี่ย min(|residual|, 5 s), ยิ่งน้อยยิ่งดี)")
    print(f"ความไม่แน่นอนโดยประมาณ : ±{rob['uncertainty_km']:.0f} km")
    print(f"Azimuthal gap       : {rob['gap']:.0f}°" +
          ("  ⚠ เกิน 180° (สถานีอยู่ฝั่งเดียวกัน ตำแหน่งตามแนวทิศทางคลาดได้มาก)" if rob["gap"] > 180 else ""))
    if rob["passes"] == 2:
        print("(ตัดเวลาคลื่นที่ไม่สอดคล้องออกแล้วคำนวณใหม่อีกรอบ)")
    return rob


def compare_reference(combined, rob, results):
    """เทียบกับจุดศูนย์กลางอ้างอิง (เช่น USGS)"""
    from datetime import datetime, timezone
    la, lo, label = REF_EPICENTER
    print(f"\n================ เทียบกับ {label} ({la:.3f}, {lo:.3f}) ================")
    print(f"ผลที่ได้ห่างจาก {label} : {haversine_km(la, lo, *combined):.0f} km")
    if rob is not None and REF_ORIGIN_TIME:
        t_ref = datetime.fromisoformat(REF_ORIGIN_TIME).replace(tzinfo=timezone.utc).timestamp()
        print(f"เวลาเกิดต่างจาก {label}  : {rob['origin_time'] - t_ref:+.1f} s")
    print(f"\n{'Station':<8}{'Δ จาก S-P':>11}{'Δ จริง':>9}{'คลาด':>9}{'%':>6}")
    for r in results:
        d = haversine_km(r["lat"], r["lon"], la, lo)
        print(f"{r['code']:<8}{r['delta']:9.0f} km{d:7.0f} km{r['delta'] - d:+7.0f} km{100 * (r['delta'] - d) / d:+5.0f}%")


def main(stations=SELECTED_STATIONS, profile=PROFILE, plot=True):
    # ---------- หาไฟล์ของแต่ละสถานี ----------
    if stations:
        if len(stations) > MAX_STATIONS:
            print(f"⚠ เลือกได้ไม่เกิน {MAX_STATIONS} สถานี ใช้ {MAX_STATIONS} สถานีแรก")
            stations = stations[:MAX_STATIONS]
        paths = []
        for code in stations:
            path = find_station_file(code)
            if path is None:
                print(f"⚠ ไม่พบไฟล์ข้อมูลของสถานี {code} (ต้องมี '# Station: {code}' ในโฟลเดอร์ station)")
            else:
                print(f"สถานี {code:<5} -> {os.path.basename(path)}")
                paths.append(path)
    else:
        paths = [resolve_data_path(DATA_FILE)]
    print(f"Profile : {profile}")

    results = [r for r in (analyze_station(p, profile, plot) for p in paths) if r is not None]
    if not results:
        print("\nไม่มีสถานีที่วิเคราะห์ได้"); return None

    combined, geiger, rob = None, None, None
    if len(results) == 1:
        print("\n* สถานีเดียวแยกไม่ได้ว่าเป็น A หรือ B ต้องใช้สถานีที่ 2 หรือความรู้เรื่องแนวรอยเลื่อนช่วยตัดสิน")
    else:
        # ---------- รวมผลหลายสถานี ----------
        combined, spread = choose_consistent(results)
        print(f"\n================ ผลรวมจาก {len(results)} สถานี ================")
        if LOCATE_METHOD == "pipeline" and len(paths) < 3:
            print("(pipeline ต้องมี 3 สถานีขึ้นไป ใช้วิธี robust แทน)")
        if LOCATE_METHOD == "pipeline" and len(paths) >= 3:
            import pipeline as PL
            pr = PL.locate_event([{"file": p} for p in paths], profile, station_coords,
                                 settings={"use_s": USE_S, "fix_depth": DEPTH_KM})
            PL.print_report(pr, REF_EPICENTER, REF_ORIGIN_TIME)
            if plot:
                out = PL.plot_record_section(pr, "epicenter_picks.png")
                if out:
                    print(f"บันทึกกราฟ pick: {out}")
            h = pr["hypo"]
            combined = (h["lat"], h["lon"])
            pst = {q["code"]: q for q in pr["stations"]}
            for r in results:                           # ไว้แสดงบนแผนที่: ระยะและสถานะจาก pipeline
                q = pst.get(r["code"])
                if q is not None:
                    r["delta"] = q.get("delta_sp") or q["dist_km"]
                    r["outlier"] = not q["use"]
            choose_nearest(results, combined)
            rob = {"lat": h["lat"], "lon": h["lon"], "origin_time": h["origin_time"]}
        elif LOCATE_METHOD == "robust" or (LOCATE_METHOD == "pipeline" and len(paths) < 3):
            rob = locate_robust(results)
            if rob is not None:
                combined = (rob["lat"], rob["lon"])
                choose_nearest(results, combined)        # ตัวเลือก A/B ที่ใกล้ผล robust (ไว้แสดงบนแผนที่)
        if LOCATE_METHOD == "ab" or combined is None or rob is None:
            rob = None
            print(f"{'Station':<8}{'S-P (s)':>9}{'Δ (km)':>9}   {'ตัวเลือกที่สอดคล้อง':<20}{'lat':>9}{'lon':>10}")
            for r in results:
                c = r["cands"][r["chosen"]]
                print(f"{r['code']:<8}{r['sp']:9.2f}{r['delta']:9.0f}   {c['name']:<20}{c['lat']:9.3f}{c['lon']:10.3f}")
            print(f"\nตำแหน่งรวม (เฉลี่ยตัวเลือกที่สอดคล้องกัน) : lat {combined[0]:.3f}, lon {combined[1]:.3f}")
            print(f"ตัวเลือกที่เลือกห่างกันมากสุด {spread:.0f} km  (ยิ่งน้อยยิ่งสอดคล้อง)")

        if SHOW_GEIGER and len(results) >= 3:
            vp, vs = (C.VPN_KM_S, C.VSN_KM_S) if profile == "regional" else (C.VP_KM_S, C.VS_KM_S)
            g = geiger_locate([{"lat": r["lat"], "lon": r["lon"], "tp": r["tp"], "ts": r["ts"],
                                "wp": r.get("wp"), "ws": r.get("ws")} for r in results], vp=vp, vs=vs)
            geiger = (g["lat"], g["lon"], g["depth"])
            print(f"\nGeiger ({g['n_obs']} เวลาคลื่น) : lat {g['lat']:.3f}, lon {g['lon']:.3f}, "
                  f"ลึก {g['depth']:.1f} km, T0 {to_utc(g['origin_time'])} UTC, RMS {g['rms']:.3f} s")
            d = haversine_km(g["lat"], g["lon"], combined[0], combined[1])
            print(f"Geiger ห่างจากผลหลัก {d:.0f} km")
            if g["depth"] > 700 or not np.isfinite(g["err_x_km"]):
                print("⚠ ผล Geiger ไม่น่าเชื่อถือ (ความลึกเกิน 700 km หรือคำนวณความคลาดเคลื่อนไม่ได้)")

        if REF_EPICENTER is not None and LOCATE_METHOD != "pipeline":
            compare_reference(combined, rob, results)

    # ---------- แผนที่จริง ----------
    if MAKE_MAP:
        if len(results) == 1:
            r = results[0]
            html = make_epicenter_map(
                r["code"], r["lat"], r["lon"], r["delta"],
                [dict(r["cands"][0], color="red"), dict(r["cands"][1], color="purple")],
                f"epicenter_map_{r['code']}.html", ref_epicenter=REF_EPICENTER)
        else:
            html = make_multi_station_map(
                results, "epicenter_map_" + "_".join(r["code"] for r in results) + ".html",
                combined=combined, geiger=geiger, ref_epicenter=REF_EPICENTER,
                combined_label=("ผล P + S" if LOCATE_METHOD == "pipeline" else "ผล Robust")
                if rob is not None else "ตำแหน่งรวม",
                subtitle=("วงกลม = ระยะ Δ จาก S−P, ดาวแดง = ผล Robust, สีเทา = สถานีที่ถูกตัดออก"
                          if rob is not None else None))
        if html:
            print(f"\nบันทึกแผนที่จริง: {html}  (เปิดด้วย browser หรือพิมพ์ start {html})")
    return results


if __name__ == "__main__":
    if len(sys.argv) > 1:              # python run_epicenter.py BKKA BSPA CMAI [local]
        main([a for a in sys.argv[1:] if a not in ("local", "regional")],
             "local" if "local" in sys.argv else PROFILE)
    else:
        main()