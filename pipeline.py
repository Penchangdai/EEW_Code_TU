"""
pipeline.py
หาตำแหน่งแผ่นดินไหวจากหลายสถานี โดยใช้ทั้งคลื่น P และคลื่น S (Sg/Sn, Lg)
ตามขั้นตอนแบบเครือข่าย USGS/ANSS  (run_network.py และ run_epicenter.py เรียกใช้ไฟล์นี้)

  ขั้น 1  อ่านข้อมูล + auto-pick คลื่น P (STA/LTA -> AIC) พร้อมคุณภาพ (weight 0-4)
  ขั้น 2  Association: ตัดสถานีที่เวลา P ขัดกับสถานีอื่น (คนละเหตุการณ์ / pick ผิดมาก)
  ขั้น 3  ตำแหน่งเริ่มต้นจาก "P + S scan" (ไม่ต้อง pick คลื่น S ก่อน)
          ลองทุกจุดบนตาราง lat/lon:  - เวลา P ของทุกสถานีต้องเข้ากัน  (บอกทิศ)
                                      - ที่เวลา P + (S-P ที่ควรเป็น) พลังงานแนวราบต้องเพิ่มขึ้น  (บอกระยะ)
          สำคัญมากเมื่อสถานีอยู่ใกล้กันและอยู่ฝั่งเดียวกัน: คลื่น P อย่างเดียวบอกระยะไม่ได้
  ขั้น 4  Re-pick P รอบเวลาที่คาด (แก้ pick ที่ไป trigger บน noise)
  ขั้น 5  Geiger (ทรงกลม + โลก 2 ชั้น) จากคลื่น P อย่างเดียว + ตรวจทาน  (เก็บไว้เทียบ)
  ขั้น 6  ทำนายเวลาคลื่น S และ Lg จากตำแหน่งที่ได้ แล้ว pick เฉพาะในช่วงนั้น
          บนแกน R/T ที่หมุนด้วยทิศจริง (ไม่ไปจับคลื่น P ระลอกหลังเป็นคลื่น S)
  ขั้น 7  Geiger จาก P + S + Lg พร้อมกัน + ตรวจทาน residual แยกตามเฟส
  ขั้น 8  ทำขั้น 6-7 ซ้ำด้วยตำแหน่งใหม่ จนตำแหน่งนิ่ง

ค่าตั้งต้นทั้งหมดอยู่ใน DEFAULTS ด้านล่าง เปลี่ยนได้โดยส่ง settings={...} หรือเพิ่มชื่อเดียวกันใน config.py
"""
import os

import numpy as np

import config as C
from data_loader import load_station_txt, resolve_data_path, find_station_file
from p_wave import detect_p_wave, pick_p_window
from s_wave import pick_s_window
from epicenter import geiger_locate_2layer, predict_arrival, phase_name, pick_weight, azimuthal_gap
from signal_utils import haversine_km, azimuth_deg, bandpass

PHASE_KEYS = {"P": ("tp", "wp"), "S": ("ts", "ws"), "Lg": ("tlg", "wlg")}

DEFAULTS = {
    "use_s": True,             # ใช้คลื่น S / Lg ช่วยหาตำแหน่ง
    "s_phases": ("S", "Lg"),   # เฟสคลื่น S ที่ pick: "S" = Sg/Sn (มาถึงก่อน), "Lg"  (profile local ใช้แค่ "S")
    "fix_depth": 10.0,         # ตรึงความลึก (km)  None = หาความลึกด้วยเสมอ (ต้องมีสถานีใกล้ < ~2 เท่าความลึก)
    # หาความลึกแบบมีเงื่อนไข: ใช้ pick ชุดสุดท้ายหาความลึก แล้วรับค่าเมื่อผ่านเกณฑ์ทุกข้อ ไม่ผ่าน = ใช้ fix_depth ต่อ
    "solve_depth": True,       # False = ใช้ fix_depth เสมอ
    "depth_near_km": 150.0,    # ต้องมีสถานีใกล้กว่านี้ (ไกลกว่านี้คลื่นตรงไม่มาก่อน เวลาคลื่นแทบไม่บอกความลึก)
    "depth_min_obs": 8,        # จำนวนเวลาคลื่นขั้นต่ำ (หา 4 ค่า: lat, lon, ความลึก, เวลาเกิด)
    "depth_max_err_km": 15.0,  # ความคลาดของความลึก (1σ) ไม่เกินค่านี้
    "depth_min_rms_gain": 0.10,  # ถ้าความลึกต่างจาก fix_depth เกิน 1σ RMS ต้องดีขึ้นอย่างน้อยสัดส่วนนี้
    "depth_max_shift_km": 30.0,  # ตำแหน่งต้องไม่ย้ายเกินค่านี้ (ย้ายมาก = ความลึกแลกกับตำแหน่ง ไม่ใช่วัดได้จริง)
    "min_weight_used": 3,      # ใช้ pick ที่ weight <= ค่านี้
    "res_max": {"regional": {"P": 3.0, "S": 6.0, "Lg": 6.0},     # |residual| (s) ที่ยอมรับได้
                "local": {"P": 1.0, "S": 2.0, "Lg": 2.0}},
    "first_pass": {"P": 10.0, "S": 15.0, "Lg": 15.0},           # เกณฑ์รอบแรก (เครือข่ายใหญ่)
    "mad_factor": 3.0,         # เกณฑ์ = max(res_max, mad_factor x MAD ของ pick ที่ผ่าน)
    "small_net_obs": 15,       # มีเวลาคลื่นไม่เกินนี้ = เครือข่ายเล็ก ใช้ตรวจทานแบบ leave-one-out
    "max_reject_small": 2,     # เครือข่ายเล็ก: ตัดได้สูงสุดกี่ pick
    "max_reject_frac": 0.4,    # เครือข่ายใหญ่: ตัดได้ไม่เกินสัดส่วนนี้
    "repick_p": True,          # ขั้น 5
    "p_win_sec": 10.0,         # หน้าต่าง ± (s) รอบเวลา P ที่คาด
    "s_win_sec": 10.0,         # หน้าต่าง ± (s) รอบเวลา S / Lg ที่คาด
    "s_iterations": 2,         # จำนวนรอบของขั้น 6-7
    "assoc_slack_sec": 5.0,    # Association: เวลา P สองสถานีต่างกันได้ไม่เกิน d/Vp + ค่านี้
    "scan_lg_weight": 2.0,     # ขั้น 3: น้ำหนักของเวลาพลังงานแนวราบสูงสุด (≈ คลื่น Lg)    0 = ไม่ใช้
    "scan_lg_delay": 8.0,      # ขั้น 3: พลังงานสูงสุดมักมาหลังคลื่น Lg เริ่มประมาณนี้ (s)
    "env_smooth_sec": 5.0,     # ขั้น 3: ความยาวหน้าต่างเฉลี่ยพลังงานแนวราบ (s)
    "scan_band": (0.5, 5.0),   # ขั้น 3: ช่วงความถี่ของพลังงานแนวราบ (Hz)
    "scan_pad_deg": 12.0,      # ขั้น 3: ขอบเขตตารางรอบสถานี (องศา)
    "max_locate_dist_km": None,     # ใช้หาตำแหน่งเฉพาะสถานีที่ห่างจากตำแหน่งเริ่มต้นไม่เกินนี้ (km) เช่น 1000
                                    #   None = ใช้ทุกสถานี  (ต้องเหลือคลื่น P อย่างน้อย min_locate_p สถานี)
    "min_locate_p": 4,
    "max_start_offset_sec": 600.0,   # ไฟล์ที่เวลาเริ่มต่างจากค่ากลางของทุกไฟล์เกินนี้ = คนละช่วงเวลา -> ไม่ใช้
    # ตรวจแผ่นดินไหวไกลมาก (teleseism) ซึ่งโมเดล 2 ชั้นจะวางตำแหน่งผิดไว้ใกล้เครือข่าย -> ต้องไม่แจ้งเตือน
    #   นับเป็น teleseism เมื่อเข้าเกณฑ์อย่างน้อย tele_min_votes ข้อ จาก 3 ข้อ:
    "tele_vapp_kms": 15.0,     # 1) ความเร็วปรากฏของหน้าคลื่น P ข้ามเครือข่าย >= ค่านี้ (Pn ~8 km/s, teleseism > 15)
    "tele_lg_frac": 0.35,      # 2) สัดส่วนสถานีที่มีคลื่น P แล้วพบ Lg (หรือ S) ตรงเวลาที่คาด < ค่านี้
    "tele_lg_err_s": 20.0,     # 3) เวลาพลังงานแนวราบสูงสุด (≈ Lg) คลาดจากตำแหน่งเฉลี่ย > ค่านี้ (s)
    "tele_min_votes": 2,
    "tele_min_p": 5,           # ต้องมีคลื่น P อย่างน้อยเท่านี้จึงตรวจ (น้อยกว่านี้ไม่ตัดสิน)
    "tele_min_aperture_km": 300.0,   # ความกว้างเครือข่ายขั้นต่ำสำหรับเกณฑ์ข้อ 1
}


def _cfg(settings):
    st = {k: getattr(C, "PIPE_" + k.upper(), v) for k, v in DEFAULTS.items()}
    st.update(settings or {})
    return st


# ======================================================================
#  โครงสร้างข้อมูลของแต่ละสถานี
#    picks[ph] = {"time", "weight", "label", "source"}   ph = "P", "S", "Lg"
#    use       = set ของเฟสที่ใช้หาตำแหน่งตอนนี้
#    notes[ph] = เหตุผลที่ไม่ใช้
# ======================================================================
def _locator_dict(st):
    """แปลงเป็น dict ที่ geiger_locate_2layer ใช้ (เฉพาะเฟสที่ใช้)"""
    d = {"code": st["code"], "lat": st["lat"], "lon": st["lon"]}
    for ph, (tk, wk) in PHASE_KEYS.items():
        if ph in st["use"] and ph in st["picks"]:
            d[tk] = st["picks"][ph]["time"]
            d[wk] = st["picks"][ph]["weight"]
    return d


def _n_obs(stations, phases=("P", "S", "Lg")):
    return sum(1 for s in stations for ph in s["use"] if ph in phases and ph in s["picks"]
               and pick_weight(s["picks"][ph]["weight"]) > 0)


def refine_depth(stations, hypo, cfg, log=None, verbose=False):
    """
    หาความลึกแบบมีเงื่อนไข จาก pick ชุดเดียวกับผลที่ตรึงความลึก (เทียบกันตรง ๆ และไม่ทำให้ pick หาย)
    คืนค่า (hypo, info)  ผ่านเกณฑ์ = hypo ใหม่ที่มีความลึกที่คำนวณได้,  ไม่ผ่าน = hypo เดิม
    info : {"solved": bool, "reason": เหตุผล, "depth", "err_z", ...}
    """
    info = {"solved": False, "reason": "", "depth": None, "err_z": None}
    if cfg["fix_depth"] is None:
        info.update(solved=True, reason="คำนวณความลึกเสมอ (fix_depth = None)", depth=hypo.get("depth"),
                    err_z=hypo.get("err_z_km"))
        return hypo, info
    if not cfg["solve_depth"]:
        info["reason"] = "ปิดการหาความลึก (solve_depth = False)"
        return hypo, info
    used = [s for s in stations if s["use"]]
    if not used:
        info["reason"] = "ไม่มีสถานีที่ใช้"
        return hypo, info
    dmin = min(float(haversine_km(hypo["lat"], hypo["lon"], s["lat"], s["lon"])) for s in used)
    info["dmin_km"] = dmin
    if dmin > cfg["depth_near_km"]:
        info["reason"] = f"สถานีที่ใกล้ที่สุดห่าง {dmin:.0f} km (ต้องไม่เกิน {cfg['depth_near_km']:.0f} km)"
        return hypo, info
    phases = ("P",) + tuple(cfg["s_phases"]) if cfg["use_s"] else ("P",)
    n_obs = _n_obs(used, phases)
    if n_obs < cfg["depth_min_obs"]:
        info["reason"] = f"เวลาคลื่นมี {n_obs} ค่า (ต้องอย่างน้อย {cfg['depth_min_obs']})"
        return hypo, info
    try:
        free = _locate(stations, dict(cfg, fix_depth=None), init=(hypo["lat"], hypo["lon"]), phases=phases)
    except (ValueError, np.linalg.LinAlgError) as ex:
        info["reason"] = f"คำนวณไม่สำเร็จ ({ex})"
        return hypo, info
    ez, z = free.get("err_z_km"), free["depth"]
    shift = float(haversine_km(hypo["lat"], hypo["lon"], free["lat"], free["lon"]))
    info.update(depth=z, err_z=ez, rms_fixed=hypo["rms"], rms_free=free["rms"], shift_km=shift)
    if ez is None or not np.isfinite(ez) or ez <= 0 or ez > cfg["depth_max_err_km"]:
        info["reason"] = f"ความลึก {z:.0f} km คลาด ±{ez:.0f} km (ต้องไม่เกิน ±{cfg['depth_max_err_km']:.0f} km)"
    elif shift > cfg["depth_max_shift_km"]:
        info["reason"] = f"ตำแหน่งย้าย {shift:.0f} km (ต้องไม่เกิน {cfg['depth_max_shift_km']:.0f} km)"
    elif free["rms"] > hypo["rms"]:
        info["reason"] = f"RMS แย่ลง ({hypo['rms']:.2f} -> {free['rms']:.2f} s)"
    elif abs(z - cfg["fix_depth"]) > ez and free["rms"] > hypo["rms"] * (1.0 - cfg["depth_min_rms_gain"]):
        info["reason"] = (f"ความลึก {z:.0f} km ต่างจากค่าตรึง แต่ RMS ดีขึ้นไม่ถึง "
                          f"{cfg['depth_min_rms_gain']:.0%} ({hypo['rms']:.2f} -> {free['rms']:.2f} s)")
    if info["reason"]:
        return hypo, info
    info["solved"] = True
    info["reason"] = f"ความลึก {z:.0f} ±{ez:.0f} km, RMS {hypo['rms']:.2f} -> {free['rms']:.2f} s"
    return dict(hypo, **free), info


def _locate(stations, cfg, init=None, phases=("P", "S", "Lg")):
    used = [_locator_dict(s) for s in stations if s["use"]]
    return geiger_locate_2layer(used, fix_depth=cfg["fix_depth"], init=init, phases=phases)


def _residuals(stations, hypo, only_used=False):
    """[(index สถานี, เฟส, residual)] ของทุก pick (หรือเฉพาะที่ใช้)"""
    out = []
    for j, s in enumerate(stations):
        for ph, pk in s["picks"].items():
            if ph not in PHASE_KEYS or (only_used and ph not in s["use"]):
                continue
            out.append((j, ph, pk["time"] - predict_arrival(hypo, s["lat"], s["lon"], ph)))
    return out


def _say(log, verbose, msg):
    log.append(msg)
    if verbose:
        print(msg)


# ======================================================================
#  ขั้น 1  อ่านข้อมูล + auto-pick
# ======================================================================
def station_files_in(folder):
    """ไฟล์ข้อมูลสถานีทั้งหมดในโฟลเดอร์ (บรรทัดแรก '# Station:')  สถานีซ้ำใช้ไฟล์แรก"""
    out, seen = [], set()
    for fn in sorted(os.listdir(folder), key=lambda x: (len(x), x)):
        p = os.path.join(folder, fn)
        if not (fn.lower().endswith(".txt") and os.path.isfile(p)):
            continue
        try:
            with open(p, encoding="latin-1") as f:
                first = f.readline()
        except OSError:
            continue
        if first.startswith("# Station:"):
            code = first.split(":", 1)[1].strip().upper()
            if code not in seen:
                seen.add(code)
                out.append({"file": p})
    return sorted(out, key=lambda d: d["file"])


def _as_folder(s):
    """ชื่อโฟลเดอร์ (path เต็ม หรือเทียบกับโฟลเดอร์โค้ด / โฟลเดอร์ station)  ไม่ใช่โฟลเดอร์ = None"""
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (s, os.path.join(here, s), os.path.join(here, "station", s)):
        if os.path.isdir(cand):
            return cand
    return None


def expand_specs(specs):
    """
    แปลงสิ่งที่ส่งมาเป็นรายการสถานี
      "ชื่อโฟลเดอร์"        -> ทุกไฟล์สถานีในโฟลเดอร์นั้น (เช่น "2026-05-18_M5.2_AllTMD")
      "ไฟล์.txt" / "BKKA"   -> สถานีเดียว
      list ของแบบใดก็ได้ (dict / ชื่อไฟล์ / รหัส / โฟลเดอร์ / (ไฟล์, lat, lon))
    """
    if isinstance(specs, (str, os.PathLike)):
        specs = [specs]
    out = []
    for sp in specs or []:
        if isinstance(sp, (str, os.PathLike)):
            folder = _as_folder(os.fspath(sp))
            if folder:
                found = station_files_in(folder)
                if not found:
                    print(f"⚠ ไม่พบไฟล์ข้อมูลสถานีในโฟลเดอร์ {sp}")
                out.extend(found)
                continue
        d = normalize_spec(sp)
        if d is not None:
            out.append(d)
    return out


def normalize_spec(spec):
    """
    รับรายการสถานีได้หลายแบบ แล้วแปลงเป็น {"file": path, "lat"?, "lon"?}
      {"file": "BKKA_....txt"}  /  "BKKA_....txt" หรือ path เต็ม  /  "BKKA" (รหัสสถานี)
      ("BKKA_....txt", lat, lon)
    """
    if isinstance(spec, dict):
        return spec if spec.get("file") else None
    if isinstance(spec, (tuple, list)) and spec:
        d = {"file": str(spec[0])}
        if len(spec) >= 3:
            d.update(lat=spec[1], lon=spec[2])
        return normalize_spec(d)
    if isinstance(spec, (str, os.PathLike)):
        s = os.fspath(spec)
        if os.path.isdir(s):
            return None
        if not s.lower().endswith(".txt") and not os.path.exists(s):     # รหัสสถานี เช่น "BKKA"
            p = find_station_file(s)
            if p is None:
                print(f"⚠ ไม่พบไฟล์ข้อมูลของสถานี {s}")
                return None
            s = p
        return {"file": s}
    return None


def load_stations(specs, profile, get_coords, log, verbose):
    regional = profile == "regional"
    p_par = C.P_REGIONAL if regional else {}
    out = []
    for spec in specs:
        spec = normalize_spec(spec)
        if spec is None:
            continue
        path = resolve_data_path(spec["file"])
        rec = load_station_txt(path)
        code = rec["station"]
        if not rec.get("units_ok", True):
            _say(log, verbose, f"⚠ ข้ามสถานี {code}: {rec['units_note']}")
            continue
        lat, lon = spec.get("lat"), spec.get("lon")
        if lat is None or lon is None:
            c = get_coords(code) if get_coords else None
            if c is None:
                _say(log, verbose, f"⚠ ข้ามสถานี {code}: ไม่พบพิกัด")
                continue
            lat, lon = c
        st = {"code": code, "lat": float(lat), "lon": float(lon), "rec": rec, "file": path,
              "sensor": "ACC" if rec["sensor"] == "acceleration" else "VEL",
              "picks": {}, "use": set(), "notes": {}}
        p = detect_p_wave(rec["z"], rec["fs"], rec["t0"], **p_par)
        if p["detected"]:
            st["picks"]["P"] = {"time": p["time"], "weight": p.get("weight", 0),
                                "label": p.get("label", "P"), "source": "auto", "index": p["index"]}
        if "P" in st["picks"] and st["picks"]["P"]["weight"] <= 3:
            st["use"].add("P")
        elif "P" in st["picks"]:
            st["notes"]["P"] = f"weight {st['picks']['P']['weight']}"
        out.append(st)
    return out


def check_time_coverage(stations, cfg, log, verbose):
    """
    ไฟล์ของเหตุการณ์เดียวกันควรเริ่มใกล้ ๆ กัน ถ้าไฟล์ใดเริ่มห่างจากค่ากลางของทุกไฟล์เกิน max_start_offset_sec
    (เช่น เป็นคนละเหตุการณ์ / เวลาในไฟล์เป็นเวลาท้องถิ่นแทน UTC / นาฬิกาสถานีผิด) จะไม่ใช้สถานีนั้น
    """
    if len(stations) < 3:
        return stations
    med = float(np.median([s["rec"]["t0"] for s in stations]))
    keep = []
    for s in stations:
        off = s["rec"]["t0"] - med
        if abs(off) > cfg["max_start_offset_sec"]:
            _say(log, verbose, f"⚠ ไม่ใช้สถานี {s['code']}: ไฟล์เริ่ม {utc(s['rec']['t0'])} UTC ต่างจากไฟล์อื่น "
                               f"{off / 60:+.0f} นาที (คนละเหตุการณ์ หรือเวลาในไฟล์ไม่ใช่ UTC?)")
            continue
        keep.append(s)
    return keep


# ======================================================================
#  ขั้น 2  Association
# ======================================================================
def associate(stations, cfg, log, verbose):
    vp = C.VPN_KM_S
    while True:
        act = [s for s in stations if "P" in s["use"]]
        if len(act) < 3:
            return
        bad = np.zeros(len(act), int)
        for a in range(len(act)):
            for b in range(a + 1, len(act)):
                d = haversine_km(act[a]["lat"], act[a]["lon"], act[b]["lat"], act[b]["lon"])
                if abs(act[a]["picks"]["P"]["time"] - act[b]["picks"]["P"]["time"]) > d / vp + cfg["assoc_slack_sec"]:
                    bad[a] += 1
                    bad[b] += 1
        j = int(np.argmax(bad))
        if bad[j] < 2 or np.sum(bad == bad[j]) > 1:
            return
        act[j]["use"].discard("P")
        act[j]["notes"]["P"] = f"association (ขัดกับ {bad[j]} สถานี)"
        _say(log, verbose, f"⚠ Association: ไม่ใช้ P ของ {act[j]['code']} (เวลา P ขัดกับสถานีอื่น {bad[j]} สถานี)")


# ======================================================================
#  ขั้น 4 / 7  ตรวจทาน residual
# ======================================================================
def review(stations, hypo, cfg, profile, phases, log, verbose, candidates=None):
    """
    เลือก pick ที่ใช้ใหม่จาก residual (pick ที่เคยตัดไปแล้วกลับมาได้ถ้าเข้ากับตำแหน่งใหม่)
    candidates : {(index, phase)} ที่เลือกได้  (None = ทุก pick ที่ weight ผ่าน และไม่ถูก association ตัด)
    """
    res_max = cfg["res_max"][profile]
    if candidates is None:
        candidates = {(j, ph) for j, s in enumerate(stations) for ph in s["picks"]
                      if ph in phases and s["picks"][ph]["weight"] <= cfg["min_weight_used"] and not s.get("far")
                      and not str(s["notes"].get(ph, "")).startswith("association")}
    n_cand = len(candidates)
    if n_cand <= cfg["small_net_obs"]:
        return _review_loo(stations, hypo, cfg, res_max, phases, candidates, log, verbose)

    thr = {ph: max(res_max[ph], cfg["first_pass"][ph]) for ph in phases}
    for it in range(12):
        res = _residuals(stations, hypo)
        keep = {(j, ph) for j, ph, r in res if (j, ph) in candidates and abs(r) <= thr[ph]}
        # ตัดได้ไม่เกินสัดส่วนที่กำหนด (เก็บตัวที่ residual เล็กที่สุดไว้ก่อน)
        n_min = int(np.ceil((1 - cfg["max_reject_frac"]) * n_cand))
        if len(keep) < n_min:
            order = sorted(((abs(r) / res_max[ph], j, ph) for j, ph, r in res if (j, ph) in candidates))
            keep = {(j, ph) for _, j, ph in order[:n_min]}
        if sum(1 for j, ph in keep if ph == "P") < 3:
            break
        for j, s in enumerate(stations):
            s["use"] = {ph for ph in s["use"] if ph not in phases} | {ph for (jj, ph) in keep if jj == j}
        new = _locate(stations, cfg, init=(hypo["lat"], hypo["lon"]), phases=phases)
        moved = haversine_km(new["lat"], new["lon"], hypo["lat"], hypo["lon"])
        hypo = new
        # เกณฑ์รอบถัดไป: MAD ของ pick ที่ผ่าน (แยก P กับ S/Lg)
        new_thr = {}
        for ph in phases:
            r_in = np.array([r for j, p2, r in _residuals(stations, hypo) if p2 == ph and (j, p2) in keep])
            mad = 1.4826 * np.median(np.abs(r_in - np.median(r_in))) if len(r_in) >= 3 else 0.0
            new_thr[ph] = min(max(res_max[ph], cfg["mad_factor"] * mad), thr[ph])
        stable = moved < 0.5 and all(abs(new_thr[ph] - thr[ph]) < 0.05 for ph in phases)
        thr = new_thr
        if stable and it > 0:
            break
    _mark_notes(stations, hypo, thr, phases)
    return hypo, thr


def _review_loo(stations, hypo, cfg, res_max, phases, candidates, log, verbose):
    """เครือข่ายเล็ก: ตัด pick ที่ห่างจากตำแหน่งที่หา "โดยไม่ใช้มัน" มากที่สุด ทีละตัว"""
    for j, s in enumerate(stations):
        s["use"] = {ph for ph in s["use"] if ph not in phases} | {ph for (jj, ph) in candidates if jj == j}
    hypo = _locate(stations, cfg, init=(hypo["lat"], hypo["lon"]), phases=phases)
    for _ in range(cfg["max_reject_small"]):
        n_free = 3 if cfg["fix_depth"] is not None else 4
        if _n_obs(stations, phases) <= n_free + 1:
            break
        worst = None
        for j, s in enumerate(stations):
            for ph in list(s["use"]):
                if ph not in phases:
                    continue
                s["use"].discard(ph)
                ok = _n_obs(stations, phases) >= n_free and _n_obs(stations, ("P",)) >= 2
                if ok:
                    h2 = _locate(stations, cfg, init=(hypo["lat"], hypo["lon"]), phases=phases)
                    r = s["picks"][ph]["time"] - predict_arrival(h2, s["lat"], s["lon"], ph)
                    if abs(r) > res_max[ph] and h2["rms"] < 0.5 * hypo["rms"] and \
                            (worst is None or abs(r) / res_max[ph] > abs(worst[2]) / res_max[worst[1]]):
                        worst = (j, ph, r, h2)
                s["use"].add(ph)
        if worst is None:
            break
        j, ph, r, h2 = worst
        stations[j]["use"].discard(ph)
        hypo = h2
        _say(log, verbose, f"⚠ ตรวจทาน: ไม่ใช้ {ph} ของ {stations[j]['code']} "
                           f"(ห่างจากตำแหน่งที่หาโดยไม่ใช้มัน {r:+.1f} s)")
    thr = {ph: res_max[ph] for ph in phases}
    _mark_notes(stations, hypo, thr, phases)
    return hypo, thr


def _mark_notes(stations, hypo, thr, phases):
    for j, s in enumerate(stations):
        for ph in phases:
            if ph in s["picks"] and ph not in s["use"] and ph not in s["notes"]:
                r = s["picks"][ph]["time"] - predict_arrival(hypo, s["lat"], s["lon"], ph)
                s["notes"][ph] = f"residual {r:+.1f} s"
            if ph in s["use"]:
                s["notes"].pop(ph, None)


# ======================================================================
#  ขั้น 5 / 6  pick ในหน้าต่างรอบเวลาที่คาด
# ======================================================================
def repick_p(stations, hypo, cfg, profile, log, verbose):
    regional = profile == "regional"
    p_par = C.P_REGIONAL if regional else {}
    res_max = cfg["res_max"][profile]["P"]
    n_fix = n_new = 0
    for s in stations:
        if s.get("far"):
            continue
        t_pred = predict_arrival(hypo, s["lat"], s["lon"], "P")
        old = s["picks"].get("P")
        if old is not None and abs(old["time"] - t_pred) <= res_max and old["weight"] <= cfg["min_weight_used"]:
            continue
        r = s["rec"]
        k = pick_p_window(r["z"], r["fs"], r["t0"], t_pred, cfg["p_win_sec"],
                          max_weight=cfg["min_weight_used"], **p_par)
        if k is None:
            continue
        if old is None:
            n_new += 1
        else:
            n_fix += 1
            s["picks"]["P_auto"] = old
        s["picks"]["P"] = {"time": k["time"], "weight": k["weight"], "label": k["label"] + "*",
                           "source": "window", "index": k["index"]}
        s["use"].add("P")
        s["notes"].pop("P", None)
    _say(log, verbose, f"Re-pick P รอบเวลาที่คาด (± {cfg['p_win_sec']:.0f} s): แก้ {n_fix} สถานี, "
                       f"ได้ P ใหม่ {n_new} สถานี")


def pick_s_all(stations, hypo, cfg, log, verbose):
    n = {ph: 0 for ph in cfg["s_phases"]}
    for s in stations:
        if s.get("far"):
            continue
        r = s["rec"]
        baz = float(azimuth_deg(s["lat"], s["lon"], hypo["lat"], hypo["lon"]))
        s["baz_true"] = baz
        tp = s["picks"]["P"]["time"] if "P" in s["use"] else predict_arrival(hypo, s["lat"], s["lon"], "P")
        got = {}
        for ph in cfg["s_phases"]:
            t_pred = predict_arrival(hypo, s["lat"], s["lon"], ph)
            k = pick_s_window(r["z"], r["n"], r["e"], r["fs"], r["t0"], t_pred, baz,
                              cfg["s_win_sec"], phase=ph, max_weight=cfg["min_weight_used"],
                              not_before=tp + 1.0)
            s["picks"].pop(ph, None)
            s["use"].discard(ph)
            s["notes"].pop(ph, None)
            if k is not None:
                got[ph] = k
            else:
                s["notes"][ph] = "ไม่พบ (พลังงานไม่เพิ่มชัดในหน้าต่าง)"
        # S กับ Lg ใกล้กันเกินไป = คลื่นเดียวกัน -> เก็บตัวที่คุณภาพดีกว่า
        if "S" in got and "Lg" in got and abs(got["S"]["time"] - got["Lg"]["time"]) < 3.0:
            drop = "S" if got["S"]["weight"] > got["Lg"]["weight"] else "Lg"
            s["notes"][drop] = "ซ้ำกับอีกเฟส"
            got.pop(drop)
        for ph, k in got.items():
            s["picks"][ph] = {"time": k["time"], "weight": k["weight"], "label": k["label"],
                              "source": "window", "ratio": k["ratio"]}
            s["use"].add(ph)
            n[ph] += 1
    _say(log, verbose, "Pick คลื่น S รอบเวลาที่คาด (± {:.0f} s): ".format(cfg["s_win_sec"])
         + ", ".join(f"{ph} {n[ph]} สถานี" for ph in cfg["s_phases"]))


# ======================================================================
#  ขั้น 3  P + S scan
# ======================================================================


def station_constraints(stations, cfg):
    """
    t_env : เวลาที่พลังงานแนวราบ (0.5-5 Hz, เฉลี่ย 5 s) สูงสุดทั้งไฟล์ ≈ เวลาคลื่น Lg (ช้ากว่าเล็กน้อย)
            ไม่ต้องรู้ตำแหน่งก่อน จึงใช้ช่วยหาจุดเริ่มต้นได้  ใช้ได้แม้สถานีที่ไม่มี pick คลื่น P
    """
    for st in stations:
        r, fs = st["rec"], st["rec"]["fs"]
        try:
            nf, ef = (bandpass(r[k], fs, *cfg["scan_band"]) for k in ("n", "e"))
            w = max(int(cfg["env_smooth_sec"] * fs), 1)
            sm = np.convolve(np.sqrt(nf ** 2 + ef ** 2), np.ones(w) / w, "same")
            st["t_env"] = float(r["t0"] + int(np.argmax(sm)) / fs)
        except Exception:                                          # noqa: BLE001
            st.pop("t_env", None)


def scan_locate(stations, cfg, use_s=True):
    """
    ขั้น 3: grid search หาตำแหน่งเริ่มต้นของ Geiger   cost = ความคลาดของเวลา P + ความคลาดของเวลา Lg
      P  : |เวลา P − ที่คาด| / 2 s            (ตัดที่ 5 s)                        -> บอกทิศ/ตำแหน่งตามแนวขวาง
      Lg : |t_env − T0 − tLg − delay| / 15 s  (ตัดที่ 60 s) x scan_lg_weight     -> บอกระยะ
    """
    from epicenter import travel_time_2layer, phase_travel_time
    act = [s for s in stations if "P" in s["use"]]
    depth = cfg["fix_depth"] or 10.0
    slat = np.array([s["lat"] for s in act]); slon = np.array([s["lon"] for s in act])
    tp = np.array([s["picks"]["P"]["time"] for s in act])
    wp = np.array([pick_weight(s["picks"]["P"]["weight"]) for s in act])
    env_st = [s for s in stations if s.get("t_env") is not None] if use_s else []
    elat = np.array([s["lat"] for s in env_st]); elon = np.array([s["lon"] for s in env_st])
    tenv = np.array([s["t_env"] for s in env_st])
    wl = cfg["scan_lg_weight"]

    def cost(LA, LO):
        D = haversine_km(LA[..., None], LO[..., None], slat, slon)
        off = tp - travel_time_2layer(D, depth, "P")
        t0 = np.median(off, axis=-1)
        mis = np.sum(wp * np.minimum(np.abs(off - t0[..., None]), 5.0), axis=-1) / np.sum(wp)
        el = np.zeros(LA.shape)
        if len(env_st) and wl > 0:
            DE = haversine_km(LA[..., None], LO[..., None], elat, elon)
            rl = tenv - t0[..., None] - phase_travel_time(DE, depth, "Lg")
            el = np.mean(np.minimum(np.abs(rl - cfg["scan_lg_delay"]), 60.0) / 15.0, axis=-1)
        return mis / 2.0 + wl * el, mis, el

    pad = cfg["scan_pad_deg"]
    lats = np.arange(slat.min() - pad, slat.max() + pad, 0.1)
    lons = np.arange(slon.min() - pad, slon.max() + pad, 0.1)
    for step, half in ((0.1, None), (0.02, 0.6)):            # ตารางหยาบ แล้วละเอียดรอบจุดที่ดีที่สุด
        if half is not None:
            lats = np.arange(la - half, la + half, step)     # noqa: F821  (la, lo ได้จากรอบแรก)
            lons = np.arange(lo - half, lo + half, step)     # noqa: F821
        LA, LO = np.meshgrid(lats, lons, indexing="ij")
        tot, mis, el = cost(LA, LO)
        i, j = np.unravel_index(np.argmin(tot), tot.shape)
        la, lo = float(lats[i]), float(lons[j])
    D = haversine_km(la, lo, slat, slon)
    t0 = float(np.median(tp - travel_time_2layer(D, depth, "P")))
    return {"lat": la, "lon": lo, "depth": depth, "origin_time": t0, "model": "2layer", "rms": 99.0,
            "scan_p_misfit": float(mis[i, j]), "scan_lg_err": float(el[i, j] * 15.0), "n_env": len(env_st)}


def delta_from_differential(dt, phase="S", depth_km=10.0):
    """ระยะบนผิวโลก Δ จากเวลาต่าง (S-P หรือ Lg-P) ด้วยแบบจำลองเดียวกับที่ใช้หาตำแหน่ง"""
    from epicenter import phase_travel_time
    d = np.linspace(1, 4000, 8000)
    tab = phase_travel_time(d, depth_km, phase) - phase_travel_time(d, depth_km, "P")
    return float(np.interp(dt, tab, d))


# ======================================================================
#  ฟังก์ชันหลัก
# ======================================================================
def locate_event(specs, profile="regional", get_coords=None, settings=None, verbose=True):
    """
    specs : list ของ {"file": ไฟล์ข้อมูล, "lat"?, "lon"?}
    return dict:
      hypo     : ผลสุดท้าย (P + S ถ้า use_s)
      hypo_p   : ผลจากคลื่น P อย่างเดียว (ไว้เทียบ)
      start    : ผลจากขั้น 3 (P + S scan)
      stations : ข้อมูลทุกสถานี (picks, use, notes, residual ...)
      log      : ข้อความทุกขั้นตอน
    """
    cfg = _cfg(settings)
    if profile == "local" and "s_phases" not in (settings or {}):
        cfg["s_phases"] = ("S",)          # สถานีใกล้: Lg ยังไม่แยกจาก Sg จึง pick แค่ S
    log = []
    specs = expand_specs(specs)
    if verbose and specs:
        print(f"ไฟล์ข้อมูลที่ใช้ : {len(specs)} สถานี")
    stations = load_stations(specs, profile, get_coords, log, verbose)
    stations = check_time_coverage(stations, cfg, log, verbose)
    associate(stations, cfg, log, verbose)

    n_p = sum(1 for s in stations if "P" in s["use"])
    if n_p < 3:
        raise ValueError("มีสถานีที่หาคลื่น P ได้น้อยกว่า 3 สถานี หาตำแหน่งไม่ได้")
    station_constraints(stations, cfg)

    # ---------- ขั้น 3  จุดเริ่มต้น: เวลา P + เวลาคลื่น Lg ----------
    hypo = scan_locate(stations, cfg, use_s=cfg["use_s"])
    _say(log, verbose, f"ขั้น 3  จุดเริ่มต้น ({'P + Lg scan' if cfg['use_s'] else 'P scan'}): "
                       f"lat {hypo['lat']:.2f}, lon {hypo['lon']:.2f}"
                       + (f"  (เวลา Lg คลาด {hypo['scan_lg_err']:.0f} s จาก {hypo['n_env']} สถานี)" if cfg["use_s"] else ""))
    start = dict(hypo)

    # ---------- จำกัดระยะของสถานีที่ใช้หาตำแหน่ง ----------
    #   ไกลกว่า ~1000 km คลื่น P วิ่งลึกลงในเนื้อโลก (เร็วกว่าแบบจำลอง 2 ชั้น) และคลื่น Lg ที่สถานีไกลมาจากส่วนของ
    #   รอยเลื่อนที่แตกทีหลัง -> ดึงตำแหน่งไปตามทิศการแตก  สถานีเหล่านี้ยังใช้หาขนาด/ความรุนแรงได้ตามปกติ
    dmax = cfg["max_locate_dist_km"]
    if dmax:
        far = [s for s in stations if haversine_km(start["lat"], start["lon"], s["lat"], s["lon"]) > dmax]
        n_p_near = sum(1 for s in stations if "P" in s["use"] and s not in far)
        if far and n_p_near >= cfg["min_locate_p"]:
            for s in far:
                for ph in list(s["use"]):
                    s["notes"][ph] = f"ไกลเกิน {dmax:.0f} km (ไม่ใช้หาตำแหน่ง)"
                s["use"] = set()
                s["far"] = True
            _say(log, verbose, f"ไม่ใช้ {len(far)} สถานีที่ไกลเกิน {dmax:.0f} km ในการหาตำแหน่ง "
                               f"(เหลือคลื่น P {n_p_near} สถานี)")

    # ---------- ขั้น 4  Re-pick P รอบเวลาที่คาด ----------
    if cfg["repick_p"]:
        repick_p(stations, hypo, cfg, profile, log, verbose)

    # ---------- ขั้น 5  P อย่างเดียว (เก็บไว้เทียบ) ----------
    hypo_p = None
    if sum(1 for s in stations if "P" in s["use"]) >= 3:
        hypo_p, _ = review(stations, start, cfg, profile, ("P",), log, verbose)
        _say(log, verbose, f"ขั้น 5  Geiger (P อย่างเดียว): lat {hypo_p['lat']:.3f}, lon {hypo_p['lon']:.3f}, "
                           f"RMS {hypo_p['rms']:.2f} s, ใช้ P {hypo_p['n_by_phase']['P']} สถานี")
        hypo_p = dict(hypo_p)
    if not cfg["use_s"]:
        hypo = hypo_p if hypo_p is not None else hypo

    # ---------- ขั้น 6-8  P + S + Lg ----------
    if cfg["use_s"]:
        phases = ("P",) + tuple(cfg["s_phases"])
        hypo = start                                   # ทำนายหน้าต่าง S จากผล scan (ใช้ S ช่วยบอกระยะแล้ว)
        for it in range(cfg["s_iterations"]):
            pick_s_all(stations, hypo, cfg, log, verbose)
            prev = hypo
            hypo, _ = review(stations, hypo, cfg, profile, phases, log, verbose)
            nb = hypo["n_by_phase"]
            _say(log, verbose, f"ขั้น 7  Geiger (P+S) รอบ {it + 1}: lat {hypo['lat']:.3f}, lon {hypo['lon']:.3f}, "
                               f"RMS {hypo['rms']:.2f} s, ใช้ P {nb['P']} / S {nb['S']} / Lg {nb['Lg']}")
            if haversine_km(prev["lat"], prev["lon"], hypo["lat"], hypo["lon"]) < 2.0 and it > 0:
                break
        if hypo["n_by_phase"]["S"] + hypo["n_by_phase"]["Lg"] == 0:
            _say(log, verbose, "⚠ ไม่มีคลื่น S / Lg ที่ใช้ได้ ผลนี้มาจากคลื่น P อย่างเดียว")

    # ---------- หาความลึก (เมื่อข้อมูลบอกได้จริง) ----------
    hypo, depth_info = refine_depth(stations, hypo, cfg)
    hypo["depth_info"] = depth_info
    if cfg["fix_depth"] is not None:
        _say(log, verbose, ("ความลึก: คำนวณได้ " if depth_info["solved"] else
                            f"ความลึก: ใช้ค่าตรึง {cfg['fix_depth']:g} km เพราะ") + " " + depth_info["reason"])

    # ---------- สรุปผลรายสถานี ----------
    used = [s for s in stations if s["use"]]
    hypo["gap"] = azimuthal_gap(hypo["lat"], hypo["lon"], used) if used else 360.0
    for s in stations:
        d = float(haversine_km(hypo["lat"], hypo["lon"], s["lat"], s["lon"]))
        s["dist_km"] = d
        s["baz_true"] = float(azimuth_deg(s["lat"], s["lon"], hypo["lat"], hypo["lon"]))
        s["residual"] = {}
        s["phase_name"] = {}
        for ph, pk in s["picks"].items():
            if ph in PHASE_KEYS:
                s["residual"][ph] = pk["time"] - predict_arrival(hypo, s["lat"], s["lon"], ph)
                s["phase_name"][ph] = phase_name(d, ph, hypo["depth"])
        # ระยะจาก S-P หรือ Lg-P ของสถานีนี้ (สถานีเดียว ไม่ขึ้นกับเวลาเกิด) -- ไว้วาดวงกลมบนแผนที่
        s["delta_sp"], s["delta_sp_phase"] = None, None
        for ph in ("S", "Lg"):
            if "P" in s["use"] and ph in s["use"]:
                s["delta_sp"] = delta_from_differential(s["picks"][ph]["time"] - s["picks"]["P"]["time"],
                                                        ph, hypo["depth"])
                s["delta_sp_phase"] = ph
                break
    hypo["teleseism"] = teleseism_check(stations, start, cfg)
    if hypo["teleseism"]["flag"]:
        _say(log, verbose, "⚠ น่าจะเป็นแผ่นดินไหวไกลมาก (teleseism): " + "; ".join(hypo["teleseism"]["reasons"])
             + "  ตำแหน่งนี้ใช้ไม่ได้ และไม่ควรแจ้งเตือน")
    return {"hypo": hypo, "hypo_p": hypo_p, "start": start, "stations": stations, "log": log,
            "settings": cfg}


def plane_wave_fit(stations):
    """
    ฟิตหน้าคลื่นระนาบกับเวลาคลื่น P ที่ใช้ (robust, IRLS)  t = t0 + sx·x + sy·y
    คืนค่า (ความเร็วปรากฏ km/s, ความกว้างเครือข่าย km)  ข้อมูลไม่พอ = None
    """
    st = [s for s in stations if "P" in s["use"] and "P" in s["picks"]]
    if len(st) < 3:
        return None
    lat = np.array([s["lat"] for s in st]); lon = np.array([s["lon"] for s in st])
    tp = np.array([s["picks"]["P"]["time"] for s in st])
    lat0, lon0 = lat.mean(), lon.mean()
    x = (lon - lon0) * 111.19 * np.cos(np.radians(lat0)); y = (lat - lat0) * 111.19
    A = np.c_[np.ones_like(x), x, y]
    w = np.ones_like(tp)
    for _ in range(5):
        c = np.linalg.lstsq(A * w[:, None], tp * w, rcond=None)[0]
        r = tp - A @ c
        mad = 1.4826 * np.median(np.abs(r)) + 0.5
        w = 1.0 / np.maximum(1.0, np.abs(r) / (2.0 * mad))
    slow = float(np.hypot(c[1], c[2]))
    vapp = 1.0 / slow if slow > 0 else np.inf
    aperture = float(np.ptp(np.hypot(x, y)))
    return vapp, aperture


def teleseism_check(stations, start, cfg):
    """
    แผ่นดินไหวไกลมาก (เช่น M8.8 คัมชัตกา ห่าง ~7000 km) คลื่น P มาเกือบพร้อมกันทั้งเครือข่าย (ความเร็วปรากฏสูง)
    และไม่มีคลื่น Lg ของเปลือกโลกมาตามเวลาที่ตำแหน่งใกล้ ๆ ทำนาย  โมเดล 2 ชั้นจึงวางตำแหน่งผิดไว้ใกล้เครือข่าย
    ใช้ 3 เกณฑ์ (ดู DEFAULTS tele_*) เข้าอย่างน้อย tele_min_votes ข้อ = teleseism
    """
    out = {"flag": False, "votes": 0, "reasons": [], "v_app": None, "lg_frac": None, "lg_err": None}
    with_p = [s for s in stations if "P" in s["picks"]]
    n_p = sum(1 for s in stations if "P" in s["use"])
    if n_p < cfg["tele_min_p"]:
        return out
    pw = plane_wave_fit(stations)
    if pw is not None:
        out["v_app"], ap = pw
        if ap >= cfg["tele_min_aperture_km"] and out["v_app"] >= cfg["tele_vapp_kms"]:
            out["votes"] += 1
            out["reasons"].append(f"คลื่น P ข้ามเครือข่ายเร็ว {out['v_app']:.0f} km/s")
    ph = "Lg" if "Lg" in cfg["s_phases"] else "S"
    if with_p:
        out["lg_frac"] = sum(1 for s in with_p if ph in s["use"]) / len(with_p)
        if out["lg_frac"] < cfg["tele_lg_frac"]:
            out["votes"] += 1
            out["reasons"].append(f"พบคลื่น {ph} ตรงเวลาที่คาดเพียง {out['lg_frac']:.0%} ของสถานี")
    if start.get("n_env"):
        out["lg_err"] = float(start.get("scan_lg_err", 0.0))
        if out["lg_err"] > cfg["tele_lg_err_s"]:
            out["votes"] += 1
            out["reasons"].append(f"เวลาพลังงานแนวราบคลาดเฉลี่ย {out['lg_err']:.0f} s")
    out["flag"] = out["votes"] >= cfg["tele_min_votes"]
    return out


# ======================================================================
#  รายงานผล
# ======================================================================
def utc(t):
    from datetime import datetime, timezone
    return datetime.fromtimestamp(t, tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%f")[:-4]


def print_report(result, ref_epicenter=None, ref_origin_time=None):
    """พิมพ์ตารางรายสถานี + สรุปผล + เทียบกับจุดอ้างอิง"""
    hypo, stations = result["hypo"], result["stations"]
    print("\n" + "=" * 100)
    print("ผลรายสถานี   (residual = เวลาที่วัด - เวลาที่คำนวณ จากตำแหน่งสุดท้าย, x = ไม่ใช้, * = pick ใหม่รอบเวลาที่คาด)")
    print("=" * 100)
    print(f"{'Station':<7}{'Δ (km)':>8}  {'P':<8}{'res':>7}  {'S':<7}{'res':>7}  {'Lg':<8}{'res':>7}"
          f"{'Δ(S-P)':>9}  หมายเหตุ")
    for s in sorted(stations, key=lambda q: q["dist_km"]):
        cells = []
        for ph, wl in (("P", 8), ("S", 7), ("Lg", 8)):
            if ph in s["picks"]:
                lab = s["picks"][ph]["label"] + ("" if ph in s["use"] else " x")
                cells.append(f"{lab:<{wl}}{s['residual'][ph]:+7.1f}")
            else:
                cells.append(f"{'-':<{wl}}{'-':>7}")
        dsp = f"{s['delta_sp']:.0f}" if s.get("delta_sp") else "-"
        note = "; ".join(f"{ph}: {v}" for ph, v in s["notes"].items()
                         if ph in PHASE_KEYS and ph not in s["use"] and "ไม่พบ" not in v)
        print(f"{s['code']:<7}{s['dist_km']:8.0f}  " + "  ".join(cells) + f"{dsp:>9}  {note}")

    nb = hypo.get("n_by_phase", {})
    print("\n" + "=" * 100)
    print("ผลหาตำแหน่ง")
    print("=" * 100)
    print(f"วิธี         : {hypo.get('method', '')}")
    print(f"Epicenter    : lat {hypo['lat']:.3f}, lon {hypo['lon']:.3f}")
    di = hypo.get("depth_info") or {}
    print(f"Depth        : {hypo['depth']:.1f} km"
          + (f"  (คำนวณได้ ±{di['err_z']:.0f} km)" if di.get("solved") and di.get("err_z") is not None
             else "  (ค่าตรึง ไม่ได้วัด)" if di else ""))
    if di and not di.get("solved"):
        print(f"               เหตุผล: {di['reason']}")
    print(f"Origin time  : {utc(hypo['origin_time'])} UTC")
    print(f"RMS residual : {hypo['rms']:.2f} s   (ถ่วงน้ำหนักตาม weight)")
    ex, ey, et = hypo.get("err_x_km", np.nan), hypo.get("err_y_km", np.nan), hypo.get("err_t0_s", np.nan)
    print(f"Error (1σ)   : ตะวันออก ±{ex:.1f} km, เหนือ ±{ey:.1f} km, T0 ±{et:.2f} s")
    print(f"เวลาคลื่นที่ใช้ : P {nb.get('P', 0)}, S {nb.get('S', 0)}, Lg {nb.get('Lg', 0)}")
    gap = hypo.get("gap", 360.0)
    print(f"Azimuthal gap: {gap:.0f}°" + ("  ⚠ เกิน 180° สถานีอยู่ฝั่งเดียวกัน (ตำแหน่งตามแนวระยะคลาดได้มาก)"
                                          if gap > 180 else ""))

    rows = [("P + S (ผลสุดท้าย)", hypo), ("P + S scan (ขั้น 3)", result.get("start")),
            ("P อย่างเดียว", result.get("hypo_p"))]
    if ref_epicenter is not None:
        la, lo, name = ref_epicenter
        t_ref = None
        if ref_origin_time:
            from datetime import datetime, timezone
            t_ref = datetime.fromisoformat(ref_origin_time).replace(tzinfo=timezone.utc).timestamp()
        print(f"\nเทียบกับ {name} ({la:.3f}, {lo:.3f}):")
        for lab, h in rows:
            if h is None:
                continue
            d = haversine_km(la, lo, h["lat"], h["lon"])
            dt = f", เวลาเกิดต่าง {h['origin_time'] - t_ref:+.1f} s" if t_ref else ""
            print(f"  {lab:<20}: lat {h['lat']:.3f}, lon {h['lon']:.3f}  ห่าง {d:.0f} km{dt}")


def plot_record_section(result, out="network_picks.png", max_stations=24, ref_hypo=None):
    """
    กราฟตรวจ pick: envelope (log) ของแต่ละสถานี เรียงตามระยะ
      ดำ = แนวดิ่ง Z, แดง = แนวราบ, เส้นทึบ = pick ที่ใช้, เส้นประ = pick ที่ไม่ใช้,
      เส้นจุด = เวลาที่คำนวณจากตำแหน่งสุดท้าย (P เขียว, S ม่วง, Lg ส้ม)
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    hypo = result["hypo"]
    sts = sorted(result["stations"], key=lambda q: q["dist_km"])[:max_stations]
    col = {"P": "g", "S": "m", "Lg": "darkorange"}
    fig, axs = plt.subplots(len(sts), 1, figsize=(13, 1.6 * len(sts) + 1), squeeze=False)
    for ax, s in zip(axs[:, 0], sts):
        r = s["rec"]
        fs = r["fs"]
        tP = predict_arrival(hypo, s["lat"], s["lon"], "P")
        tL = predict_arrival(hypo, s["lat"], s["lon"], "Lg")
        i0 = max(int((tP - 40 - r["t0"]) * fs), 0)
        i1 = min(int((tL + 60 - r["t0"]) * fs), len(r["z"]))
        if i1 - i0 < 10:
            continue
        t = r["t0"] + np.arange(i0, i1) / fs - hypo["origin_time"]
        k = max(int(0.5 * fs), 1)
        zf = bandpass(r["z"][i0:i1], fs, 1.0, 10.0)
        nf, ef = (bandpass(r[c][i0:i1], fs, 0.5, 5.0) for c in ("n", "e"))
        ez = np.sqrt(np.convolve(zf ** 2, np.ones(k) / k, "same")) + 1e-12
        eh = np.sqrt(np.convolve(nf ** 2 + ef ** 2, np.ones(k) / k, "same")) + 1e-12
        ax.semilogy(t, ez, "k", lw=0.5)
        ax.semilogy(t, eh, "r", lw=0.5, alpha=0.8)
        for ph in ("P", "S", "Lg"):
            ax.axvline(predict_arrival(hypo, s["lat"], s["lon"], ph) - hypo["origin_time"],
                       color=col[ph], ls=":", lw=1)
            if ph in s["picks"]:
                ax.axvline(s["picks"][ph]["time"] - hypo["origin_time"], color=col[ph],
                           ls="-" if ph in s["use"] else "--", lw=1.6 if ph in s["use"] else 1)
        ax.set_xlim(t[0], t[-1])
        ax.text(0.005, 0.85, f"{s['code']}  {s['dist_km']:.0f} km", transform=ax.transAxes, fontsize=8,
                bbox=dict(facecolor="white", alpha=0.7, lw=0))
        ax.tick_params(labelsize=7)
    axs[-1, 0].set_xlabel("Time from origin (s)")
    fig.suptitle("Picks: solid = used, dashed = rejected, dotted = predicted  (green P, magenta S, orange Lg)"
                 "\nblack = Z envelope (1-10 Hz), red = horizontal envelope (0.5-5 Hz)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=90)
    plt.close(fig)
    return out