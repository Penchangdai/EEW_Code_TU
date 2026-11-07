"""
run_magnitude.py
คำนวณขนาดแผ่นดินไหวจากหลายสถานี แล้วแจ้งผลเฉพาะขนาดที่น่าเชื่อถือที่สุด

  1) หาตำแหน่ง เวลาเกิด และเวลาคลื่น P / S ของแต่ละสถานี
       - ใส่ EPICENTER ที่รู้แล้ว  หรือ
       - 3 สถานีขึ้นไป: pipeline.py (P + S scan -> Geiger P + S -> ตรวจทาน) หาให้เอง  หรือ
       - 1-2 สถานี: ระยะจาก S-P ของแต่ละสถานี (คุณภาพต่ำ)
  2) คำนวณขนาดทุกชนิดของแต่ละสถานี (magnitude.station_magnitudes)
       มาตรฐาน USGS / IASPEI: ML, mb_Lg, Ms_BB, Mwp     แบบเดิม / เตือนภัยเร็ว: M_JMA, M_Pd, M_tc
  3) ให้น้ำหนัก รวมเป็นขนาดของเครือข่าย ให้คะแนนแต่ละชนิด แล้วเลือกขนาดหลัก (ดูวิธีใน magnitude.py)
  4) แสดงผล
       OUTPUT_MODE = "public" : การ์ดแจ้งเตือนสำหรับประชาชน (รายละเอียดทั้งหมดบันทึกใน DETAIL_FILE)
       OUTPUT_MODE = "detail" : ตารางละเอียดทั้งหมดใน Terminal

รัน:  กด Run ใน VS Code   หรือ   python run_magnitude.py BKKA CMMT PHIT [local]
"""
import contextlib
import io
import os
import sys
from datetime import datetime, timedelta, timezone

import numpy as np

import config as C
import magnitude as M
from data_loader import load_station_txt, resolve_data_path, find_station_file, DATA_DIR
from p_wave import detect_p_wave, pick_p_window
from s_wave import detect_s_wave
from epicenter import delta_from_sp, predict_arrival, travel_time_2layer, azimuthal_gap
from signal_utils import haversine_km, azimuth_deg, acc_to_vel, remove_offset
from station_db import get_coords

try:
    from intensity import jma_instrumental_intensity
except ImportError:                                   # ไม่มี intensity.py ก็ยังคำนวณขนาดได้
    jma_instrumental_intensity = None

# ======================================================================
#  ตั้งค่า
# ======================================================================
# รหัสสถานี เช่น ["BKKA", "CMMT", "PHIT"]  หรือ "auto" = ทุกไฟล์ในโฟลเดอร์ station
SELECTED_STATIONS = "auto"
DATA_FILE = "BKKA_20250328_061959.txt"      # ใช้เมื่อ SELECTED_STATIONS = []
PROFILE = "regional"      # "regional" = สถานีไกล (> ~300 km), "local" = สถานีใกล้

# จุดศูนย์กลางที่รู้แล้ว (lat, lon, ความลึก km, "เวลาเกิด UTC")   None = หาเองด้วย pipeline.py
EPICENTER = None
# EPICENTER = (22.011, 95.936, 10.0, "2025-03-28T06:20:52")      # USGS
# ความลึก:  ตัวเลข = ตรึงความลึกไว้ (ไม่ได้วัด จึงไม่แสดงในการ์ด/กราฟ)  -- เหมาะกับสถานีไกลหลายร้อย km
#           None   = คำนวณความลึกด้วย (ต้องมีสถานีใกล้จุดเกิด) แสดงเมื่อความคลาดเคลื่อนยอมรับได้
DEPTH_KM = 10.0
DEPTH_MAX_ERR_KM = 15.0   # แสดงความลึกที่คำนวณได้ เมื่อคลาดเคลื่อน (1σ) ไม่เกินค่านี้

OUTPUT_MODE = "public"    # "public" = การ์ดแจ้งเตือน,  "detail" = ตารางละเอียดทั้งหมด
DETAIL_FILE = "magnitude_detail.txt"
REF_MAGNITUDE = (7.7, "Mww (USGS)")     # ขนาดอ้างอิงไว้เทียบในโหมด detail   ไม่ใช้ใส่ None
PLOT = True               # กราฟ: magnitude_alert.png (ทุกโหมด) + magnitude_summary.png (โหมด detail)
ALERT_PLOT = "magnitude_alert.png"

# เมืองที่ใช้บอกระยะ (ชื่อ, lat, lon)  การ์ดแสดงกรุงเทพฯ + 2 เมืองที่ใกล้จุดศูนย์กลางที่สุด
CITIES = [("กรุงเทพฯ", 13.7563, 100.5018), ("เชียงใหม่", 18.7883, 98.9853), ("เชียงราย", 19.9105, 99.8406),
          ("แม่ฮ่องสอน", 19.3020, 97.9654), ("ตาก", 16.8840, 99.1259), ("พิษณุโลก", 16.8211, 100.2659),
          ("ขอนแก่น", 16.4419, 102.8360), ("อุบลราชธานี", 15.2287, 104.8564), ("กาญจนบุรี", 14.0228, 99.5328),
          ("ภูเก็ต", 7.8804, 98.3923), ("สงขลา", 7.1898, 100.5951)]

RMS_MAX_SEC = 10.0        # RMS ของตำแหน่งเกินนี้ = ตำแหน่งไม่น่าเชื่อถือ ไม่คำนวณขนาดต่อ
MIN_BEFORE_P_SEC = 5.0    # ต้องมีข้อมูลก่อนคลื่น P อย่างน้อยเท่านี้ (s)
MIN_AFTER_P_SEC = 30.0    # ต้องมีข้อมูลหลังคลื่น P อย่างน้อยเท่านี้ (s)


# ======================================================================
#  เครื่องมือ
# ======================================================================
def utc(t):
    return datetime.fromtimestamp(t, tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


def _epoch(utc_str):
    return datetime.fromisoformat(utc_str).replace(tzinfo=timezone.utc).timestamp()


def _folder_of(x):
    here = os.path.dirname(os.path.abspath(__file__))
    return any(os.path.isdir(c) for c in (x, os.path.join(here, x)))


def station_files(stations):
    if stations == "auto":
        folder = DATA_DIR if os.path.isdir(DATA_DIR) else os.path.dirname(os.path.abspath(__file__))
        out, seen = [], set()
        for fn in sorted(os.listdir(folder), key=lambda x: (len(x), x)):   # ชื่อสั้นก่อน (ไฟล์ _ExampleFile ทีหลัง)
            if fn.lower().endswith(".txt"):
                p = os.path.join(folder, fn)
                with open(p, encoding="latin-1") as f:
                    first = f.readline()
                if first.startswith("# Station:"):
                    code = first.split(":", 1)[1].strip().upper()
                    if code in seen:                      # ไฟล์ซ้ำของสถานีเดียวกัน ใช้ไฟล์แรกเท่านั้น
                        continue
                    seen.add(code)
                    out.append(p)
        return sorted(out)
    if isinstance(stations, str) or (stations and any(isinstance(x, str) and _folder_of(x) for x in stations)):
        import pipeline as _PL                          # ชื่อโฟลเดอร์ เช่น "2026-05-18_M5.2_AllTMD" หรือรหัส/ไฟล์
        return [d["file"] for d in _PL.expand_specs(stations)]
    if stations:
        out = []
        for code in stations:
            p = find_station_file(code)
            if p is None:
                print(f"⚠ ไม่พบไฟล์ข้อมูลของสถานี {code}")
            else:
                out.append(p)
        return out
    return [resolve_data_path(DATA_FILE)]


# ======================================================================
#  ขั้น 1  ตำแหน่ง เวลาเกิด และเวลาคลื่นของแต่ละสถานี
# ======================================================================
def _station(code, rec, lat, lon, ip, wp, s_time, delta, depth, t0):
    fs = rec["fs"]
    return dict(code=code, rec=rec, lat=lat, lon=lon, p_index=int(ip), p_weight=wp,
                s_index=int(round((s_time - rec["t0"]) * fs)), delta=float(delta), depth=depth,
                t0_rel=t0 - rec["t0"])


def locate(paths, profile):
    """return (hypo หรือ None, list ของสถานี)"""
    regional = profile == "regional"
    p_par = C.P_REGIONAL if regional else {}
    out = []

    if EPICENTER is not None:                                    # ---- ก) จุดศูนย์กลางที่รู้แล้ว ----
        la, lo, dep, t_utc = EPICENTER
        hypo = {"lat": la, "lon": lo, "depth": dep, "origin_time": _epoch(t_utc), "method": "ใส่เอง",
                "computed": False, "depth_measured": False, "source": "ค่าที่ใส่เอง (EPICENTER) ไม่ได้คำนวณ"}
        for p in paths:
            rec = load_station_txt(p)
            if not rec.get("units_ok", True):
                print(f"⚠ ข้าม {rec['station']}: {rec['units_note']}"); continue
            c = get_coords(rec["station"])
            if c is None:
                print(f"⚠ ข้าม {rec['station']}: ไม่พบพิกัด"); continue
            k = pick_p_window(rec["z"], rec["fs"], rec["t0"], predict_arrival(hypo, c[0], c[1], "P"), 10.0,
                              max_weight=3, **p_par)
            if k is None:
                print(f"⚠ ข้าม {rec['station']}: ไม่พบคลื่น P ในช่วงเวลาที่คาด"); continue
            out.append(_station(rec["station"], rec, c[0], c[1], k["index"], k["weight"],
                                predict_arrival(hypo, c[0], c[1], "S"),
                                haversine_km(la, lo, c[0], c[1]), dep, hypo["origin_time"]))
        hypo["gap"] = azimuthal_gap(la, lo, out) if out else 360.0
        return hypo, out

    if len(paths) >= 3:                                          # ---- ข) pipeline ----
        import pipeline as PL
        res = PL.locate_event([{"file": p} for p in paths], profile, get_coords,
                              settings={"fix_depth": DEPTH_KM}, verbose=False)
        hypo = res["hypo"]
        nb = hypo.get("n_by_phase", {})
        n_st = sum(1 for q in res["stations"] if q["use"])
        ez = hypo.get("err_z_km")
        hypo["computed"] = True
        hypo["depth_measured"] = (DEPTH_KM is None and ez is not None and np.isfinite(ez)
                                  and 0 < ez <= DEPTH_MAX_ERR_KM)
        exy = [v for v in (hypo.get("err_x_km"), hypo.get("err_y_km")) if v is not None and np.isfinite(v)]
        hypo["loc_err_km"] = max(exy) if exy else None
        hypo["source"] = (f"คำนวณจากคลื่นที่ {n_st} สถานี (P {nb.get('P', 0)}, S {nb.get('S', 0)}, "
                          f"Lg {nb.get('Lg', 0)})")
        # เวลาคลื่นที่ใช้หาตำแหน่ง (run_travel_time.py ใช้หาว่าระบบรู้ตำแหน่งได้เร็วที่สุดเมื่อไร)
        hypo["picks_used"] = sorted(
            ((q["code"], ph, q["picks"][ph]["time"], q["lat"], q["lon"])
             for q in res["stations"] for ph in q["use"] if ph in q["picks"]), key=lambda x: x[2])
        hypo["hypo_p"] = res.get("hypo_p")
        for line in res["log"]:
            if line.startswith("⚠"):
                print(line)
        print(f"ตำแหน่งจาก pipeline : lat {hypo['lat']:.3f}, lon {hypo['lon']:.3f}, ลึก {hypo['depth']:.0f} km, "
              f"เวลาเกิด {utc(hypo['origin_time'])} UTC, RMS {hypo['rms']:.2f} s")
        if hypo["rms"] > RMS_MAX_SEC:
            print(f"⚠ ตำแหน่งไม่น่าเชื่อถือ (RMS {hypo['rms']:.0f} s > {RMS_MAX_SEC:g} s) จึงไม่คำนวณขนาดต่อ  "
                  "แนะนำ: รัน run_network.py ดูว่าสถานีไหน pick ผิด หรือใส่ EPICENTER")
            return hypo, []
        for s in res["stations"]:
            rec, fs = s["rec"], s["rec"]["fs"]
            if "P" in s["use"]:
                ip, wp = round((s["picks"]["P"]["time"] - rec["t0"]) * fs), s["picks"]["P"]["weight"]
            else:                                                # P ไม่ผ่าน -> ใช้เวลาที่คาด (คุณภาพต่ำ)
                ip, wp = round((predict_arrival(hypo, s["lat"], s["lon"], "P") - rec["t0"]) * fs), 4
            ts = s["picks"]["S"]["time"] if "S" in s["use"] else predict_arrival(hypo, s["lat"], s["lon"], "S")
            out.append(_station(s["code"], rec, s["lat"], s["lon"], ip, wp, ts, s["dist_km"], hypo["depth"],
                                hypo["origin_time"]))
        return hypo, out

    # ---- ค) 1-2 สถานี ----
    print("⚠ มี 1-2 สถานี: ระยะมาจาก S-P ที่ pick โดยยังไม่รู้ตำแหน่ง แผ่นดินไหวใหญ่มักจับคลื่น P ระลอกหลัง"
          "เป็น S ทำให้ระยะและขนาดผิดมาก  แนะนำให้ใส่ EPICENTER หรือใช้ 3 สถานีขึ้นไป")
    for p in paths:
        rec = load_station_txt(p)
        if not rec.get("units_ok", True):
            print(f"⚠ ข้าม {rec['station']}: {rec['units_note']}"); continue
        fs = rec["fs"]
        pk = detect_p_wave(rec["z"], fs, rec["t0"], **p_par)
        if not pk["detected"]:
            print(f"⚠ ข้าม {rec['station']}: ไม่พบคลื่น P"); continue
        s = detect_s_wave(rec["z"], rec["n"], rec["e"], fs, pk["index"], rec["t0"],
                          **(C.S_REGIONAL if regional else {}))
        if not s.get("detected"):
            print(f"⚠ ข้าม {rec['station']}: ไม่พบคลื่น S จึงหาระยะไม่ได้"); continue
        dep = DEPTH_KM if DEPTH_KM is not None else 10.0
        delta = delta_from_sp(s["sp_time"], dep)
        t0 = pk["time"] - float(travel_time_2layer(delta, dep, "P"))
        c = get_coords(rec["station"]) or (np.nan, np.nan)
        out.append(_station(rec["station"], rec, c[0], c[1], pk["index"], max(pk.get("weight", 2), 2),
                            s["time"], delta, dep, t0))
    return None, out


def covers(st):
    """ไฟล์ครอบคลุมช่วงคลื่นไหม (กัน index ติดลบ / เกินท้ายไฟล์)  -> (True/False, เหตุผล)"""
    rec, fs, ip = st["rec"], st["rec"]["fs"], st["p_index"]
    n = len(rec["z"])
    if ip < int(MIN_BEFORE_P_SEC * fs):
        return False, (f"ไฟล์เริ่ม {utc(rec['t0'])} UTC หลังเวลาคลื่น P {-ip / fs + MIN_BEFORE_P_SEC:.0f} s "
                       "(ตรวจว่าเวลาในไฟล์เป็น UTC และเป็นเหตุการณ์เดียวกัน)")
    if ip > n - int(MIN_AFTER_P_SEC * fs):
        return False, f"ไฟล์จบ {utc(rec['t0'] + n / fs)} UTC เหลือข้อมูลหลังคลื่น P {max(n - ip, 0) / fs:.0f} s"
    st["s_index"] = int(min(max(st["s_index"], ip + int(5 * fs)), n - 1))
    return True, ""


# ======================================================================
#  ขั้น 2-3  ขนาดของแต่ละสถานี -> เครือข่าย -> ขนาดหลัก
# ======================================================================
def station_values(st):
    rec, ip = st["rec"], st["p_index"]
    fs = rec["fs"]
    n_pre = max(ip - int(1.0 * fs), 10)
    st["mags"] = M.station_magnitudes(rec["z"], rec["n"], rec["e"], fs, ip, st["delta"], st["depth"],
                                      t0_rel=st["t0_rel"], s_index=st["s_index"], n_pre=n_pre,
                                      legacy=True, mwp_check=True)
    za, na, ea = (remove_offset(x, n_pre) for x in (rec["z"], rec["n"], rec["e"]))
    st["pga"] = float(np.max(np.sqrt(za ** 2 + na ** 2 + ea ** 2)))
    vn = acc_to_vel(rec["n"], fs, C.DISP_HIGHPASS_HZ, n_pre)
    ve = acc_to_vel(rec["e"], fs, C.DISP_HIGHPASS_HZ, n_pre)
    st["pgv"] = float(np.max(np.sqrt(vn ** 2 + ve ** 2)))
    st["I"] = jma_instrumental_intensity(rec["z"], rec["n"], rec["e"], fs) if jma_instrumental_intensity else None


def compute(stations=SELECTED_STATIONS, profile=PROFILE):
    paths = station_files(stations)
    print(f"Profile : {profile}   จำนวนไฟล์ : {len(paths)}")
    hypo, sts = locate(paths, profile)
    good = []
    for st in sts:
        ok, why = covers(st)
        if not ok:
            print(f"⚠ ข้ามสถานี {st['code']}: {why}"); continue
        try:
            station_values(st)
            good.append(st)
        except Exception as ex:                                  # noqa: BLE001  สถานีเดียวไม่ทำให้ทั้งระบบหยุด
            print(f"⚠ ข้ามสถานี {st['code']}: คำนวณไม่สำเร็จ ({type(ex).__name__}: {ex})")
    if not good:
        print("ไม่มีสถานีที่คำนวณได้")
        return None
    net = M.network_by_type(good)
    M0, scores = M.type_scores(net)
    pick = M.choose_preferred(net, scores)
    gap = None
    if hypo is not None:
        gap = hypo.get("gap")
        if gap is None:
            gap = azimuthal_gap(hypo["lat"], hypo["lon"], good)
    rel = M.reliability(pick, net, good, gap, hypo is not None) if pick else None
    return {"hypo": hypo, "stations": good, "network": net, "M0": M0, "scores": scores,
            "preferred": pick, "reliability": rel, "gap": gap}


# ======================================================================
#  แสดงผล: การ์ดแจ้งเตือน
# ======================================================================
_TH_MONTH = ("มกราคม", "กุมภาพันธ์", "มีนาคม", "เมษายน", "พฤษภาคม", "มิถุนายน", "กรกฎาคม", "สิงหาคม",
             "กันยายน", "ตุลาคม", "พฤศจิกายน", "ธันวาคม")
_TH_DIR = ("ทิศเหนือ", "ทิศตะวันออกเฉียงเหนือ", "ทิศตะวันออก", "ทิศตะวันออกเฉียงใต้", "ทิศใต้",
           "ทิศตะวันตกเฉียงใต้", "ทิศตะวันตก", "ทิศตะวันตกเฉียงเหนือ")


def thai_time(t):
    d = datetime.fromtimestamp(t, tz=timezone.utc) + timedelta(hours=7)
    return f"{d.day} {_TH_MONTH[d.month - 1]} {d.year + 543} เวลา {d:%H:%M:%S} น."


def location_lines(hypo):
    """ข้อความตำแหน่ง (แสดงความลึกเฉพาะเมื่อคำนวณได้จริง) และที่มาของตำแหน่ง"""
    pos = f"ละติจูด {hypo['lat']:.2f}°N  ลองจิจูด {hypo['lon']:.2f}°E"
    if hypo.get("depth_measured"):
        pos += f"   ลึก {hypo['depth']:.0f} km"
    src = hypo.get("source", "คำนวณจากสถานีตรวจวัด")
    if hypo.get("loc_err_km"):
        src += f"  คลาดประมาณ ±{hypo['loc_err_km']:.0f} km"
    return pos, src


def public_alert(res):
    line = "=" * 64
    print(line)
    pick = res["preferred"] if res else None
    if pick is None:
        print("  ระบบยังคำนวณขนาดแผ่นดินไหวไม่ได้ (ข้อมูลไม่พอ)")
        print(line)
        return
    hypo = res["hypo"]
    level, why, unc = res["reliability"]
    m = pick["mag"]
    issue = m >= getattr(C, "FORECAST_MAGNITUDE", 3.5)
    print(f"  {'แจ้งเตือนแผ่นดินไหว' if issue else 'รายงานแผ่นดินไหว'}  (ข้อมูลเบื้องต้นจากระบบอัตโนมัติ)")
    print(line)
    print(f"  ขนาด (Mw)    : {mw_text(pick)}   ({M.size_class(m)})")
    print(f"  ที่มาของขนาด  : {mw_source(pick)}")
    if hypo is not None:
        print(f"  เวลาเกิด      : {thai_time(hypo['origin_time'])} (เวลาประเทศไทย)")
        pos, src = location_lines(hypo)
        print(f"  ตำแหน่ง       : {pos}")
        print(f"  ที่มาของตำแหน่ง: {src}")
        d = [(float(haversine_km(la, lo, hypo["lat"], hypo["lon"])), name, la, lo) for name, la, lo in CITIES]
        show = sorted([x for x in d if x[1] == "กรุงเทพฯ"] + sorted(x for x in d if x[1] != "กรุงเทพฯ")[:2])
        for k, (dist, name, la, lo) in enumerate(show):
            az = float(azimuth_deg(la, lo, hypo["lat"], hypo["lon"]))
            head = "  ห่างจาก       :" if k == 0 else "                  "
            print(f"{head} {name} {dist:,.0f} km ไป{_TH_DIR[int((az + 22.5) // 45) % 8]}")
    else:
        print("  ตำแหน่ง       : ยังระบุไม่ได้ (สถานีไม่พอ)")
    print(f"  ความน่าเชื่อถือ : {level}   (ขนาดอาจคลาด ±{unc:.1f})" + (f"  เพราะ{', '.join(why)}" if why else ""))
    print(line)
    if issue:
        print("  ถ้ารู้สึกสั่น: หมอบลง ป้องศีรษะ เกาะสิ่งที่มั่นคง อยู่ห่างกระจกและสิ่งของที่อาจล้ม")
        print("  ห้ามใช้ลิฟต์ระหว่างสั่น  ติดตามประกาศทางการจากกรมอุตุนิยมวิทยา")
    print("  * ค่านี้ประมาณจากสถานีตรวจวัด อาจปรับปรุงเมื่อมีข้อมูลเพิ่ม")
    print(line)


# ======================================================================
#  แสดงผล: ตารางละเอียด
# ======================================================================
def mw_text(pick):
    """ข้อความขนาด Mw  (ขนาดที่ถึงช่วงอิ่มตัวแล้ว = ค่าต่ำสุด)"""
    return f"{pick['mag']:.1f}" + (" ขึ้นไป" if pick.get("lower_bound") else "")


def mw_source(pick):
    t, raw = pick["type"], pick.get("mag_raw", pick["mag"])
    if t in ("Mwp", "Mw_P"):
        src = {"Mwp": "Mwp จากคลื่น P (Tsuboi 1995)", "Mw_P": "คลื่น P 3 วินาทีแรก (Pd + τc)"}[t]
    else:
        src = f"{t} {raw:.1f} แปลงเป็น Mw: {pick.get('mw_note', '')}"
    if pick.get("lower_bound"):
        src += "  ·  ถึงช่วงอิ่มตัวของวิธีนี้ ขนาดจริงอาจใหญ่กว่า"
    return src


def _fmt(m, w):
    if m is None or not m.get("valid"):
        return f"{'-':>6}{'':>5}"
    return f"{m['mag']:6.2f}({w:.1f})"


def detail_report(res):
    if res is None:
        return
    sts, net, M0, scores, pick = res["stations"], res["network"], res["M0"], res["scores"], res["preferred"]
    bar = "=" * 112
    print("\n" + bar)
    print("ขนาดของแต่ละสถานี   ค่าในวงเล็บ = น้ำหนักของสถานี (0-1)   - = ใช้ไม่ได้ (นอกช่วงระยะ / SNR ต่ำ / ข้อมูลไม่พอ)")
    print(bar)
    print(f"{'Station':<7}{'Δ (km)':>7}{'Δ (°)':>6}  " + "".join(f"{t:>11}" for t in M.TYPES) + f"{'PGA':>9}{'I_JMA':>7}")
    for st in sorted(sts, key=lambda s: s["delta"]):
        row = "".join(f"{_fmt(st['mags'].get(t), st['w'].get(t, 0)):>11}" for t in M.TYPES)
        I = f"{st['I']:.1f}" if st.get("I") is not None else "-"
        print(f"{st['code']:<7}{st['delta']:7.0f}{st['delta'] / M.KM_PER_DEG:6.1f}  {row}{st['pga']:9.2f}{I:>7}")

    print("\n" + bar)
    print(f"ขนาดของเครือข่ายแต่ละชนิด (มัธยฐานถ่วงน้ำหนัก) และคะแนนความน่าเชื่อถือ   ค่ากลางของทุกชนิด M0 = {M0:.2f}")
    print(bar)
    print(f"{'ชนิด':<7}{'ค่าเดิม':>7}{'Mw':>6}{'n':>4}{'MAD':>6}  {'พื้นฐาน':>7}{'อิ่มตัว':>7}{'ช่วงM':>6}{'ข้อมูล':>7}"
          f"{'กระจาย':>7}{'สอดคล้อง':>9}{'หน้าต่าง':>8}{'SCORE':>7}   หมายเหตุ")
    for t in M.TYPES:
        v = net[t]
        if v["mag"] is None:
            print(f"{t:<7}{'-':>7}{'-':>6}{0:>4}{'':>64}{'-':>7}   ใช้ไม่ได้: {v['reason']}")
            continue
        if t not in scores:                                 # แสดงไว้เทียบ (รวมอยู่ใน Mw_P แล้ว)
            mad = f"{v['mad']:.2f}" if v["mad"] is not None else "-"
            print(f"{t:<7}{v['mag_raw']:7.2f}{v['mag']:6.2f}{v['n']:>4}{mad:>6}  {'':>57}{'-':>7}   "
                  f"{M.TYPE_INFO[t][3]} (เทียบ: รวมอยู่ใน Mw_P)")
            continue
        s = scores[t]
        mad = f"{v['mad']:.2f}" if v["mad"] is not None else "-"
        note = M.TYPE_INFO[t][3]
        if s["sat"] <= 0.11:
            note = f"เกินช่วงที่สูตรใช้ได้ดี (อิ่มตัวที่ M ~{M.TYPE_INFO[t][1]:g})"
        elif s["sat"] < 0.99:
            note = f"ใกล้ช่วงอิ่มตัว (M ~{M.TYPE_INFO[t][1]:g})"
        if t == "Mwp" and s["window"] < 0.99:
            note = f"หน้าต่าง P {v['min_window']:.0f} s สั้นกว่าเวลาแตกของรอยเลื่อน ~{M.rupture_duration(M0):.0f} s"
        if t == "Mwp" and v.get("spread") is not None and v["spread"] > M.MWP_SPREAD_REF:
            note = f"ไม่เสถียร: เปลี่ยนตัวกรองแล้วค่าต่างกัน {v['spread']:.1f} (ข้อมูลคาบยาวไม่พอ)"
        if s["consist"] < 0.5:
            note += f"; ต่างจากค่ากลาง {v['mag'] - M0:+.1f}"
        mark = "  <== เลือก" if pick and t == pick["type"] else ""
        if not v.get("mw_in_range", True) and t not in ("Mwp", "Mw_P"):
            note += "; สูตรแปลงเป็น Mw ใช้นอกช่วงที่ปรับเทียบไว้"
        print(f"{t:<7}{v['mag_raw']:7.2f}{v['mag']:6.2f}{v['n']:>4}{mad:>6}  {s['base']:7.2f}{s['sat']:7.2f}{s['range']:6.2f}"
              f"{s['data']:7.2f}{s['scatter']:7.2f}{s['consist']:9.2f}{s['window']:8.2f}{s['score']:7.2f}"
              f"   {note}{mark}")
    if pick is None:
        print("\nไม่มีขนาดที่ใช้ได้"); return
    level, why, unc = res["reliability"]
    if res["hypo"] is not None:
        h = res["hypo"]
        ez = h.get("err_z_km")
        if h.get("depth_measured"):
            dtxt = f"ลึก {h['depth']:.1f} ± {ez:.1f} km (คำนวณได้)"
        elif DEPTH_KM is None and h.get("computed"):
            dtxt = (f"ความลึกคำนวณได้ {h['depth']:.0f} km แต่คลาดเคลื่อนมาก"
                    + (f" (±{ez:.0f} km)" if ez is not None and np.isfinite(ez) else "") + " จึงไม่แสดงในการ์ด")
        else:
            dtxt = f"ความลึกตรึงไว้ที่ {h['depth']:.0f} km (ไม่ได้วัด ไม่แสดงในการ์ด)"
        print(f"\nตำแหน่ง: lat {h['lat']:.3f}, lon {h['lon']:.3f}, {dtxt}  — {location_lines(h)[1]}")
    print("\n" + bar)
    print(f"ขนาดหลัก (preferred) : Mw {mw_text(pick)}  ({mw_source(pick)})")
    print(f"                       ชนิด {pick['type']}   (score {pick['score']:.2f}, "
          f"ความน่าเชื่อถือ {level}, ±{unc:.1f})")
    print(f"ค่าเฉลี่ยถ่วงคะแนน   : {pick['combined']:.1f}   จาก {', '.join(pick['combined_types'])}")
    if REF_MAGNITUDE:
        print(f"เทียบกับ {REF_MAGNITUDE[1]} = {REF_MAGNITUDE[0]:.1f} : ต่าง {pick['mag'] - REF_MAGNITUDE[0]:+.1f}")
    first = [t for t in M.USGS_TYPES if net[t]["mag"] is not None]
    if first and first[0] != pick["type"]:
        print(f"* ลำดับปกติของ USGS จะเลือก {first[0]} ก่อน แต่ข้อมูลของ {first[0]} ครั้งนี้มีคะแนนต่ำ "
              f"({scores[first[0]]['score']:.2f}) จึงเลือก {pick['type']} แทน")
    if pick["type"] == "Ms_BB" and all(st["delta"] / M.KM_PER_DEG < 20 for st in sts):
        print("* Ms_BB จากสถานีระยะ < 20° ทั้งหมด (มาตรฐาน Ms_20 ใช้ 20°-160°) อาจคลาด ±0.3-0.5 "
              "โดยเฉพาะเมื่อรอยเลื่อนแตกพุ่งเข้าหาสถานี (directivity)")
    if pick["mag"] >= 7.0 and pick["type"] not in ("Mwp", "Mw_P"):
        print("* M >= 7: Mw นี้แปลงมาจากขนาดแอมพลิจูด (ML, mb_Lg, Ms) ด้วยสูตรเชิงประสบการณ์ ถ้ามี Mww ของ USGS / GCMT "
              "ควรใช้ค่านั้น")


def _thai_font():
    """เลือกฟอนต์ภาษาไทยที่มีในเครื่อง (Windows: Tahoma / Leelawadee, Linux: Loma / Noto Sans Thai)"""
    from matplotlib import font_manager
    have = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Leelawadee UI", "Tahoma", "Noto Sans Thai", "Loma", "Garuda", "TH Sarabun New"):
        if name in have:
            return name
    return None


MAG_CLASSES = [  # (ต่ำสุด, สูงสุด, คำบอกขนาด, ชื่อกลุ่มของ USGS, สี)
    (0.0, 3.0, "จิ๋ว", "Micro", "#d9f0d3"),
    (3.0, 4.0, "เล็กมาก", "Minor", "#a6dba0"),
    (4.0, 5.0, "เล็ก", "Light", "#fee08b"),
    (5.0, 6.0, "ปานกลาง", "Moderate", "#fdae61"),
    (6.0, 7.0, "รุนแรง", "Strong", "#f46d43"),
    (7.0, 8.0, "ใหญ่", "Major", "#d73027"),
    (8.0, 10.0, "ใหญ่มาก", "Great", "#a50026"),
]


def plot_alert(res, out=ALERT_PLOT):
    """กราฟแจ้งเตือน: มาตรวัดขนาด (แบ่งสีตามกลุ่มของ USGS) + เข็มชี้ขนาด + ช่วงคลาดเคลื่อน"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle, Polygon
    except ImportError:
        return None
    if not res or res.get("preferred") is None:
        return None
    pick, hypo = res["preferred"], res["hypo"]
    level, why, unc = res["reliability"]
    m = pick["mag"]
    font = _thai_font()
    rc = {"axes.unicode_minus": False}
    if font:
        rc["font.family"] = [font, "DejaVu Sans"]
    x0, x1 = 2.0, 9.5
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=(10, 4.6))
        ax.set_xlim(x0 - 0.2, x1 + 0.2)
        ax.set_ylim(-2.3, 3.1)
        ax.axis("off")
        # แถบสีแต่ละกลุ่มขนาด
        for lo, hi, word, eng, col in MAG_CLASSES:
            a, b = max(lo, x0), min(hi, x1)
            if b <= a:
                continue
            here = lo <= m < hi
            ax.add_patch(Rectangle((a, 0), b - a, 0.8, facecolor=col, edgecolor="white", lw=2))
            ax.text((a + b) / 2, -0.22, word, ha="center", va="top", fontsize=10,
                    color="0.1", fontweight="bold" if here else "normal")
            ax.text((a + b) / 2, -0.58, eng, ha="center", va="top", fontsize=8, color="0.45")
        for v in range(int(np.ceil(x0)), int(x1) + 1):
            ax.plot([v, v], [0.8, 0.95], color="0.3", lw=1)
            ax.text(v, 1.0, f"{v}", ha="center", va="bottom", fontsize=10, color="0.3")
        # ช่วงคลาดเคลื่อน + เข็มชี้
        lo_u, hi_u = max(m - unc, x0), min(m + unc, x1)
        ax.add_patch(Rectangle((lo_u, -0.08), hi_u - lo_u, 0.96, facecolor="none", edgecolor="k",
                               lw=1.5, ls="--"))
        mx = min(max(m, x0), x1)
        ax.add_patch(Polygon([[mx, 0.85], [mx - 0.12, 1.45], [mx + 0.12, 1.45]], closed=True, color="k"))
        ax.plot([mx, mx], [-0.08, 0.85], color="k", lw=2.5)
        ax.text(mx, 1.55, f"{m:.1f}" + ("+" if pick.get("lower_bound") else ""), ha="center", va="bottom",
                fontsize=26, fontweight="bold")
        ax.text(x0, 2.75, f"ขนาดแผ่นดินไหว Mw {mw_text(pick)}  ·  {M.size_class(m)}", fontsize=16,
                fontweight="bold", va="center")
        ax.text(x1, 2.75, f"ความน่าเชื่อถือ: {level}  (±{unc:.1f})", fontsize=11, va="center", ha="right")
        info = []
        if hypo is not None:
            info.append(f"เวลาเกิด {thai_time(hypo['origin_time'])} (เวลาประเทศไทย)")
            pos, src = location_lines(hypo)
            info.append(f"ตำแหน่ง {pos}   ({src})")
        info.append(f"ที่มาของขนาด: {mw_source(pick)}")
        info.append("กรอบเส้นประ = ช่วงที่ขนาดจริงอาจอยู่   ข้อมูลเบื้องต้นจากระบบอัตโนมัติ อาจปรับปรุงเมื่อมีข้อมูลเพิ่ม")
        ax.text(x0, -1.05, "\n".join(info), fontsize=9, va="top", color="0.25", linespacing=1.5)
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        plt.close(fig)
    return out


def plot_summary(res, out="magnitude_summary.png"):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    net, scores, pick = res["network"], res["scores"], res["preferred"]
    ts = [t for t in M.TYPES if net[t]["mag"] is not None and t in scores]
    if not ts or pick is None:
        return None
    mags = [net[t]["mag"] for t in ts]
    errs = [net[t]["mad"] or 0 for t in ts]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(range(len(ts)), mags, yerr=errs, capsize=4, alpha=0.8,
           color=["tab:red" if t == pick["type"] else "tab:blue" for t in ts])
    for i, t in enumerate(ts):
        ax.text(i, mags[i] + errs[i] + 0.08, f"{mags[i]:.1f}\nscore {scores[t]['score']:.2f}", ha="center", fontsize=8)
    if REF_MAGNITUDE:
        ax.axhline(REF_MAGNITUDE[0], color="k", ls="--", lw=1)
        ax.text(-0.45, REF_MAGNITUDE[0] - 0.18, f"{REF_MAGNITUDE[1]} = {REF_MAGNITUDE[0]:g}", fontsize=8)
    ax.set_xticks(range(len(ts)))
    ax.set_xticklabels(ts)
    ax.set_ylim(min(mags) - 1, max(mags + ([REF_MAGNITUDE[0]] if REF_MAGNITUDE else [])) + 0.8)
    ax.set_ylabel("Mw (converted)")
    ax.set_title(f"Network magnitudes as Mw (red = preferred: {pick['type']} -> Mw {pick['mag']:.1f})")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


# ======================================================================
def main(stations=SELECTED_STATIONS, profile=PROFILE, mode=None, plot=PLOT):
    mode = mode or OUTPUT_MODE
    if mode == "detail":
        res = compute(stations, profile)
        detail_report(res)
        if plot and res:
            for out in (plot_summary(res), plot_alert(res)):
                if out:
                    print(f"บันทึกกราฟ: {out}")
        return res
    buf = io.StringIO()                       # โหมด public: เก็บรายละเอียดทั้งหมดลงไฟล์ แสดงแค่การ์ด
    with contextlib.redirect_stdout(buf):
        res = compute(stations, profile)
        detail_report(res)
    with open(DETAIL_FILE, "w", encoding="utf-8") as f:
        f.write(buf.getvalue())
    public_alert(res)
    if plot:
        out = plot_alert(res)
        if out:
            print(f"กราฟขนาดแผ่นดินไหว: {out}")
    print(f"(รายละเอียดสำหรับเจ้าหน้าที่: {DETAIL_FILE})")
    return res


if __name__ == "__main__":
    if len(sys.argv) > 1:              # python run_magnitude.py BKKA CMMT PHIT [local]
        main([a for a in sys.argv[1:] if a not in ("local", "regional")],
             "local" if "local" in sys.argv else PROFILE)
    else:
        main()