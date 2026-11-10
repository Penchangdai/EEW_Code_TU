"""
run_realtime_replay.py
จำลองสถานการณ์จริง: เล่นไฟล์เหตุการณ์ที่บันทึกไว้เหมือนข้อมูลที่ทยอยมาจากสถานี แล้วให้ระบบตัดสินใจทีละรอบ

  ทุก ๆ TICK วินาทีของเวลาจำลอง
    1) ตัดข้อมูลแต่ละสถานีให้เหลือเฉพาะ sample ที่มาถึงแล้ว (ระบบไม่เห็นข้อมูลอนาคตเลย)
    2) หาตำแหน่ง -> คำนวณขนาด -> ทำนาย MMI -> ตัดสินการแจ้งเตือน ด้วยฟังก์ชันเดิมของระบบ
    3) บันทึกว่า ณ เวลานั้นระบบรู้อะไร เตือนใคร และใช้เวลาคำนวณเท่าไร

ไฟล์นี้ไม่มีสูตรใหม่ และไม่แก้ไฟล์เดิม  สูตรทั้งหมดมาจาก
    pipeline.py / p_wave.py / s_wave.py / epicenter.py   ตำแหน่งและเวลาคลื่น
    magnitude.py / run_magnitude.py                      ขนาด
    intensity.py / run_intensity.py / travel_time.py     MMI การปรับค่า และเกณฑ์แจ้งเตือน

ค่าที่ปรับตอนรัน (ไม่ใช่สูตร และมีผลเฉพาะระหว่างรันไฟล์นี้)
    pipeline.load_station_txt      = ตัวอ่านที่ส่งเฉพาะข้อมูลที่มาถึงแล้ว
    run_magnitude.MIN_AFTER_P_SEC  = 3 s  (เดิม 30 s จะเลื่อน Mw_P ไปถึง Tp + 30 s ขัดกับเวลา Tp + 3 s ในตารางของระบบ)

ค่าอ้างอิงจาก run_validation.EVENTS ใช้เทียบผลเท่านั้น ไม่ได้ส่งเข้าระบบ

รัน (จากโฟลเดอร์ EEW_):
    python run_realtime_replay.py 2026-05-18_M5.2_AllTMD
    python run_realtime_replay.py 2025-3-28_Mw7.7_station --tick 10 --duration 600
    python run_realtime_replay.py 2026-05-18_M5.2_AllTMD --speed 1              (เล่นตามเวลาจริง)
    python run_realtime_replay.py 2026-05-18_M5.2_AllTMD --latency 3 --drop 0.2 (หน่วง 3 s และสถานีขาด 20%)

ผลลัพธ์อยู่ใน replay_output/<เหตุการณ์>_*
    _summary.txt        ลำดับเวลาที่ระบบพบคลื่น P, คลื่น S/Lg, ได้ตำแหน่ง, ได้ขนาด, ทำนายความรุนแรง/แจ้งเตือน
    _milestones.csv     ตารางเดียวกันในรูป CSV
    _phases.csv         เวลาที่คลื่นถึงและเวลาที่ระบบพบคลื่น P / S ของแต่ละสถานี
    _timeline.csv       ผลของทุกรอบ         _alerts.csv   การแจ้งเตือนของแต่ละตำแหน่ง
    แผนที่จากผลรอบสุดท้าย (ฟังก์ชันเดียวกับ run_magnitude.py / run_intensity.py)
    _magnitude_map.html          ขนาด + จุดกำเนิด
    _intensity_station_map.html  ความรุนแรงแบบ ShakeMap: พื้นสี MMI + สถานีสีตาม MMI + ป้ายจังหวัด
    (--all-maps เพิ่ม _intensity_map.html = MMI ทำนาย + วงกลม + การแจ้งเตือน และ _intensity_mmi_map.png ภาพนิ่ง)
"""
import argparse
import contextlib
import csv
import io
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True                  # ไม่สร้าง __pycache__ เพิ่มในโฟลเดอร์โค้ด
os.environ.setdefault("MPLBACKEND", "Agg")      # วาดภาพลงไฟล์ ไม่เปิดหน้าต่างกราฟ

import numpy as np

import config as C
import data_loader as DL
import intensity as INT
import magnitude as M
import pipeline as PL
import run_intensity as RI
import run_magnitude as RM
import run_travel_time as RTT
import run_validation as V
import travel_time as TT
from signal_utils import haversine_km

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CODE_DIR, "replay_output")

MIN_DATA_S = 10.0                               # สถานีต้องมีข้อมูลอย่างน้อยเท่านี้ จึงส่งเข้าระบบ
MAX_START_OFFSET_S = PL.DEFAULTS["max_start_offset_sec"]   # เกณฑ์เดียวกับ pipeline.check_time_coverage
STABLE_KM = 10.0                                # ตำแหน่ง "นิ่ง" = ทุกรอบหลังจากนั้นห่างจากผลสุดท้ายไม่เกินนี้ (ใช้รายงานผลเท่านั้น)


def utc(t):
    return datetime.fromtimestamp(t, tz=timezone.utc).strftime("%H:%M:%S")


# ======================================================================
#  ตัวป้อนข้อมูล
# ======================================================================
class ReplayFeed:
    """เก็บข้อมูลเต็มไว้ แต่ส่งให้ระบบเฉพาะ sample ที่มาถึงแล้ว ณ เวลาจำลอง now"""

    def __init__(self, records, latency_s=0.0):
        self.records = records                  # {path: rec จาก data_loader}
        self.latency_s = latency_s
        self.now = -np.inf
        self.last_stations = None               # สถานีที่ pipeline อ่านในรอบล่าสุด (มี picks ของคลื่น P/S/Lg)

    def n_available(self, path):
        rec = self.records[path]
        n = int(np.floor((self.now - self.latency_s - rec["t0"]) * rec["fs"])) + 1
        return int(np.clip(n, 0, len(rec["z"])))

    def available_paths(self):
        return [p for p, r in self.records.items() if self.n_available(p) >= MIN_DATA_S * r["fs"]]

    def view(self, path):
        """แทน data_loader.load_station_txt: คืน record ที่มีเฉพาะข้อมูลที่มาถึงแล้ว"""
        rec = self.records[path]
        n = self.n_available(path)
        out = dict(rec, z=rec["z"][:n], n=rec["n"][:n], e=rec["e"][:n])
        out["units_ok"], out["units_note"] = DL.check_units(rec["header"], out["n"], rec["fs"])
        return out


@contextlib.contextmanager
def replay_patches(feed, use_correction=True):
    """ปรับค่าเฉพาะระหว่างจำลอง แล้วคืนค่าเดิมเมื่อจบ"""
    saved = PL.load_station_txt, PL.load_stations, RM.MIN_AFTER_P_SEC, RI.USE_OBS_CORRECTION
    orig_load_stations = PL.load_stations

    def load_stations_kept(*a, **k):            # เก็บรายการสถานีไว้ดูว่าระบบพบคลื่นอะไรแล้ว (ไม่เปลี่ยนผล)
        feed.last_stations = orig_load_stations(*a, **k)
        return feed.last_stations

    PL.load_station_txt = feed.view
    PL.load_stations = load_stations_kept
    RM.MIN_AFTER_P_SEC = C.PD_TAUC_WINDOW_SEC
    RI.USE_OBS_CORRECTION = use_correction
    try:
        yield
    finally:
        PL.load_station_txt, PL.load_stations, RM.MIN_AFTER_P_SEC, RI.USE_OBS_CORRECTION = saved


def load_event(event, drop_frac=0.0, seed=1):
    """อ่านไฟล์ทุกสถานีของเหตุการณ์  ตัดไฟล์เวลาผิดด้วยเกณฑ์เดียวกับ pipeline และสุ่มตัดสถานีตาม drop_frac"""
    folder = PL._as_folder(event)
    if folder is None:
        raise SystemExit(f"ไม่พบโฟลเดอร์เหตุการณ์ {event}")
    records, notes = {}, []
    for spec in PL.station_files_in(folder):
        path = spec["file"]
        try:
            rec = DL.load_station_txt(path)
        except Exception as ex:                                 # noqa: BLE001  สถานีเดียวไม่ทำให้ทั้งระบบหยุด
            notes.append(f"ข้าม {os.path.basename(path)}: อ่านไม่ได้ ({type(ex).__name__})")
            continue
        if RM.get_coords(rec["station"]) is None:
            notes.append(f"ข้าม {rec['station']}: ไม่พบพิกัด")
            continue
        records[path] = rec
    if not records:
        raise SystemExit(f"ไม่มีข้อมูลสถานีที่ใช้ได้ใน {event}")
    med = float(np.median([r["t0"] for r in records.values()]))
    for path in [p for p, r in records.items() if abs(r["t0"] - med) > MAX_START_OFFSET_S]:
        notes.append(f"ข้าม {records[path]['station']}: เวลาเริ่มไฟล์ต่างจากค่ากลางเกิน {MAX_START_OFFSET_S:.0f} s")
        del records[path]
    if drop_frac > 0:
        rng = np.random.default_rng(seed)
        keys = sorted(records)
        n_drop = int(round(drop_frac * len(keys)))
        for path in rng.choice(keys, size=n_drop, replace=False):
            notes.append(f"จำลองสถานีขาดหาย: {records[path]['station']}")
            del records[path]
    return records, notes


def make_sites(records):
    """ตำแหน่งแจ้งเตือน: ทุกสถานี + run_intensity.TARGET_SITES  (ไม่มีข้อมูลอนาคต ใช้แค่พิกัดและชนิดสถานี)"""
    sites = {}
    for rec in records.values():
        code = rec["station"]
        lat, lon = RM.get_coords(code)
        sites[code] = dict(name=code, lat=lat, lon=lon, vs30=RI.station_vs30(code), kind="station", obs=None,
                           sensor=RI.station_sensor(code, rec.get("station_type")),
                           sensor_file=rec.get("station_type"), channels=rec["header"].get("Channels", ""),
                           alert=None, upgrade=None)
    for s in RI.TARGET_SITES:
        sites[s["name"]] = dict(name=s["name"], lat=s["lat"], lon=s["lon"], vs30=float(s.get("vs30", C.DEFAULT_VS30)),
                                kind="target", obs=None, alert=None, upgrade=None)
    return sites


# ======================================================================
#  การตัดสินใจหนึ่งรอบ  (ใช้แต่ฟังก์ชันเดิมของระบบ)
# ======================================================================
def update_phases(feed, phases, now):
    """
    เวลาที่ระบบได้คลื่นแต่ละชนิดของแต่ละสถานีเป็นครั้งแรก (pick ที่ pipeline ใช้ ใน "use")
      phases[code]["P"] / ["S"] = {seen: รอบที่พบครั้งแรก, first: เวลาคลื่นถึงตอนพบ, time: เวลาคลื่นถึงล่าสุด, phase: S/Lg}
    return (จำนวนสถานีที่มีคลื่น P รอบนี้, จำนวนที่มี S/Lg, สถานีที่พบ P ใหม่, สถานีที่พบ S/Lg ใหม่)
    """
    n_p = n_s = 0
    new_p, new_s = [], []
    for st in feed.last_stations or []:
        got = set()
        for ph in ("P", "S", "Lg"):
            pk = st["picks"].get(ph)
            if ph not in st["use"] or pk is None:
                continue
            key = "P" if ph == "P" else "S"
            if key in got:                                      # มีทั้ง S และ Lg ใช้ตัวแรก
                continue
            got.add(key)
            rec = phases.setdefault(st["code"], {})
            if key not in rec:
                rec[key] = dict(seen=now, first=float(pk["time"]), phase=ph)
                (new_p if key == "P" else new_s).append(st["code"])
            rec[key].update(time=float(pk["time"]), phase=ph)
        n_p += "P" in got
        n_s += "S" in got
    return n_p, n_s, new_p, new_s


def magnitude_now(hypo, good, gap, t):
    """
    ขนาดหลักจากข้อมูลที่ครบแล้ว ณ เวลา t   ลำดับเดียวกับ run_intensity.magnitude_at
    แต่คืน network และสถานี (พร้อมน้ำหนัก w) ด้วย เพื่อใช้วาดแผนที่ขนาดแบบ run_magnitude.make_map
    """
    T0 = hypo["origin_time"]
    sts = []
    for st in good:
        mags = {k: v for k, v in st["mags"].items() if v and v.get("valid")
                and RTT.magnitude_ready_time(st, k, T0) + RTT.PROCESSING_DELAY_S <= t}
        sts.append(dict(st, mags=mags, w={}))
    with contextlib.redirect_stdout(io.StringIO()):
        net = M.network_by_type(sts)
        pick = M.choose_preferred(net)
    if pick is None:
        return None
    used = [s for s in sts if s["w"].get(pick["type"], 0) > 0]
    rel = M.reliability(pick, net, used, gap, True)
    return dict(hypo=hypo, stations=sts, network=net, preferred=pick, reliability=rel, gap=gap, used=used)


def run_tick(feed, sites, phases):
    """return dict ของรอบนี้   status บอกว่าระบบไปถึงขั้นไหน"""
    now = feed.now
    row = dict(t=now, n_data=0, status="รอข้อมูล", lat=None, lon=None, t0=None, rms=None, n_p=0, tele=False,
               mtype=None, mw=None, raw=None, unc=None, level=None, n_mag=0, hypo=None, error=None,
               n_pdet=0, n_sdet=0, new_p=[], new_s=[], mmi_max=None, mmi_site=None, res=None)
    paths = feed.available_paths()
    row["n_data"] = len(paths)
    if len(paths) < 3:
        return row

    feed.last_stations = None
    hypo = sts = None
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            hypo, sts = RM.locate(paths, RM.PROFILE)
    except ValueError as ex:                                    # คลื่น P ไม่ถึง 3 สถานี = ยังหาตำแหน่งไม่ได้ (ปกติ)
        row["status"] = "ยังตรวจพบคลื่น P ไม่ถึง 3 สถานี" if "น้อยกว่า 3" in str(ex) else f"ValueError: {ex}"
        row["error"] = None if "น้อยกว่า 3" in str(ex) else str(ex)
    except Exception as ex:                                     # noqa: BLE001
        row["status"], row["error"] = "หาตำแหน่งไม่สำเร็จ", f"{type(ex).__name__}: {ex}"
    row["n_pdet"], row["n_sdet"], row["new_p"], row["new_s"] = update_phases(feed, phases, now)
    if hypo is None:
        return row

    nb = hypo.get("n_by_phase", {})
    row.update(status="ได้ตำแหน่ง", lat=float(hypo["lat"]), lon=float(hypo["lon"]), t0=float(hypo["origin_time"]),
               rms=float(hypo["rms"]), n_p=int(nb.get("P", 0)), tele=bool(RM.is_teleseism(hypo)), hypo=hypo)
    if not sts:
        row["status"] = f"ตำแหน่งไม่น่าเชื่อถือ (RMS เกิน {RM.RMS_MAX_SEC:g} s)"
        return row

    good = []
    for st in sts:
        ok, _ = RM.covers(st)
        if not ok:
            continue
        try:
            RM.station_values(st)
            good.append(st)
        except Exception:                                       # noqa: BLE001  สถานีเดียวไม่ทำให้ทั้งระบบหยุด
            continue
    if not good:
        row["status"] = "ได้ตำแหน่ง แต่ยังไม่มีสถานีที่คำนวณขนาดได้"
        return row

    res = magnitude_now(hypo, good, hypo.get("gap"), now)       # ใช้ค่าที่ข้อมูลครบแล้วเท่านั้น
    if res is None:
        row["status"] = "ได้ตำแหน่ง แต่ยังไม่มีขนาดที่ข้อมูลครบ"
        return row
    pick, rel = res["preferred"], res["reliability"]
    unc = rel[2]
    row.update(status="ได้ขนาด", mtype=pick["type"], mw=float(pick["mag"]), raw=float(pick["mag_raw"]), unc=float(unc),
               level=rel[0], lower_bound=bool(pick.get("lower_bound")), n_mag=len(res["used"]), res=res)

    # ---- ค่าที่สถานีวัดได้จนถึงตอนนี้ (Worden 2012 + QC ของระบบ) ----
    for st in good:
        if st["code"] in sites:
            try:
                sites[st["code"]]["obs"] = RI.observed(st, hypo)
            except Exception:                                   # noqa: BLE001
                pass

    # ---- ทำนาย MMI + ปรับค่า + เกณฑ์แจ้งเตือน  (ลำดับเดียวกับ run_intensity.simulate) ----
    preds = {}
    for s in sites.values():
        d = float(haversine_km(hypo["lat"], hypo["lon"], s["lat"], s["lon"]))
        preds[s["name"]] = INT.predict_mmi(pick["mag"], d, hypo["depth"], s["vs30"], unc, None)
    all_sites = list(sites.values())
    cm = RI.correction_model(all_sites, preds, now)
    for s in all_sites:
        w = TT.site_warning(hypo, s["lat"], s["lon"], hypo["origin_time"])
        mmi = float(np.clip(preds[s["name"]]["mmi"] + RI.site_corr(cm, s, preds[s["name"]]["R"]), 1.0, 10.0))
        lv, name = INT.alert_level(pick["mag"], mmi)
        if row["tele"]:
            lv, name = 0, "ไม่แจ้ง (teleseism)"
        info = dict(t=now, M=float(pick["mag"]), type=pick["type"], mmi=mmi, level=lv, name=name,
                    lead_s_then=w["s_arrival"] - now, lead_lg_then=w["lg_arrival"] - now)
        if lv > 0 and s["alert"] is None:
            s["alert"] = info
        if lv == 2 and s["upgrade"] is None and s["alert"] is not None and s["alert"]["level"] < 2:
            s["upgrade"] = info
        if not row["tele"] and (row["mmi_max"] is None or mmi > row["mmi_max"]):
            row["mmi_max"], row["mmi_site"] = mmi, s["name"]
    return row


# ======================================================================
#  แผนที่จากผลรอบสุดท้าย  (ฟังก์ชันแผนที่เดิมของ run_magnitude / run_intensity)
# ======================================================================
def intensity_sites(res, sites):
    """
    แปลงตำแหน่งของการจำลองเป็นรูปแบบของ run_intensity.build_sites + ค่าทำนายสุดท้าย (ส่วนท้ายของ run_intensity.simulate)
    การแจ้งเตือน (alert/upgrade) เป็นของการจำลองนี้  lead time คิดจากตำแหน่งสุดท้าย
    """
    hypo, pick, unc = res["hypo"], res["preferred"], res["reliability"][2]
    used = {st["code"] for st in res["used"]}
    out = []
    for s in sites.values():
        w = TT.site_warning(hypo, s["lat"], s["lon"], hypo["origin_time"])
        d = dict(s, used=s["name"] in used, delta=w["epi_dist_km"], t_p=w["p_arrival"], t_s=w["s_arrival"],
                 t_lg=w["lg_arrival"], rjb=RI.site_rjb(s["lat"], s["lon"]))
        for k in ("alert", "upgrade"):
            if d[k] is not None:
                d[k] = dict(d[k], corr=0.0, lead_s=d["t_s"] - d[k]["t"], lead_lg=d["t_lg"] - d[k]["t"])
        d["final_raw"] = INT.predict_mmi(pick["mag"], d["delta"], hypo["depth"], d["vs30"], unc, d["rjb"])
        out.append(d)
    RI.flag_gross_outliers(out)
    cm = RI.correction_model(out, {s["name"]: s["final_raw"] for s in out}, np.inf)
    res["correction"] = cm
    for s in out:
        f = dict(s["final_raw"])
        c = RI.site_corr(cm, s, f["R"])
        if cm and cm.get("rmse") is not None and cm["method"] != "none":
            sig = float(max(cm["rmse"], 0.3))
            f.update(mmi=float(np.clip(f["mmi"] + c, 1, 10)), sigma=sig)
            f.update(mmi_low=max(f["mmi"] - sig, 1.0), mmi_high=min(f["mmi"] + sig, 10.0))
        f["corr"], f["n_obs"] = c, (cm or {}).get("n", 0)
        s["final"] = f
    return out


def make_maps(res, sites, base, all_maps=False):
    """
    return [(ชื่อ, ไฟล์)]   res = ผลรอบสุดท้ายที่มีขนาด
    ค่าตั้งต้น: แผนที่ขนาด + แผนที่ความรุนแรงแบบ ShakeMap (run_intensity.make_station_map)
    all_maps = True เพิ่มแผนที่ MMI ทำนาย + วงกลม (run_intensity.make_map) และภาพนิ่ง PNG
    """
    if res is None:
        print("⚠ ระบบยังไม่ได้ขนาดในการจำลองนี้ จึงไม่สร้างแผนที่")
        return []
    RI.select_rupture(res["hypo"])                              # แนวรอยแตกที่รู้แล้ว ใช้กับแผนที่เท่านั้น (เหมือน run_intensity)
    msites = intensity_sites(res, sites)
    jobs = [("แผนที่ขนาดและจุดกำเนิด", lambda: RM.make_map(res, base + "_magnitude_map.html")),
            ("แผนที่ความรุนแรง (MMI)", lambda: RI.make_station_map(res, msites, base + "_intensity_station_map.html"))]
    if all_maps:
        jobs += [("แผนที่ MMI ทำนายและการแจ้งเตือน", lambda: RI.make_map(res, msites, base + "_intensity_map.html")),
                 ("แผนที่ MMI ภาพนิ่ง", lambda: RI.plot_mmi_map(res, msites, base + "_intensity_mmi_map.png"))]
    files = []
    for name, job in jobs:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                f = job()
        except Exception as ex:                                 # noqa: BLE001  แผนที่หนึ่งพังไม่ทำให้ส่วนอื่นหยุด
            print(f"⚠ สร้าง{name}ไม่สำเร็จ: {type(ex).__name__}: {ex}")
            continue
        if f:
            files.append((name, f))
    return files


# ======================================================================
#  เทียบผลกับค่าอ้างอิง  (ทำหลังจบรอบ ไม่ได้ส่งเข้าระบบ)
# ======================================================================
def find_reference(event):
    ev = next((e for e in V.EVENTS if e["folder"] == os.path.basename(event.rstrip("/\\"))), None)
    if ev is None or ev.get("lat") is None:
        return None
    return dict(name=ev["name"], lat=ev["lat"], lon=ev["lon"], t0=V._epoch(ev["time"]), mw=V.ref_mw(ev),
                label=f"{ev.get('mag_type', '')} {ev.get('mag', '')}")


def score_row(row, ref):
    if ref is None or row["lat"] is None:
        return
    row["err_km"] = float(haversine_km(row["lat"], row["lon"], ref["lat"], ref["lon"]))
    row["dt0"] = row["t0"] - ref["t0"]
    if row["mw"] is not None and ref["mw"] is not None:
        row["dm"] = row["mw"] - ref["mw"]


# ======================================================================
#  แสดงผลและบันทึก
# ======================================================================
def _fmt(v, f):
    return "-" if v is None else format(v, f)


def _names(codes, n=8):
    return ", ".join(codes[:n]) + (f" และอีก {len(codes) - n}" if len(codes) > n else "")


def print_header():
    print("เวลา UTC  จากต้นไฟล์  ข้อมูล  คลื่นP คลื่นS  ตำแหน่ง lat,lon    ขนาด         MMI สูงสุด  "
          "เตือน แอป/สาธารณะ  คำนวณ  สถานะ", flush=True)


def print_row(row, t_start, n_app, n_pub, phases):
    print(f"{utc(row['t'])} +{row['t'] - t_start:5.0f}s  {row['n_data']:3d}  {row['n_pdet']:4d} {row['n_sdet']:4d}   "
          f"{_fmt(row['lat'], '7.2f'):>7},{_fmt(row['lon'], '7.2f'):>7}  "
          f"{(row['mtype'] or '-'):>5} Mw {_fmt(row['mw'], '4.1f'):>4}  {_fmt(row['mmi_max'], '4.1f'):>4}  "
          f"{n_app:4d}/{n_pub:4d}  {row['compute_s']:5.1f}s  {row['status']}"
          + ("  [teleseism]" if row["tele"] else ""), flush=True)
    if row["new_p"]:
        first = min(phases[c]["P"]["first"] for c in row["new_p"])
        print(f"      ▶ พบคลื่น P ใหม่ {len(row['new_p'])} สถานี: {_names(row['new_p'])}  (คลื่นถึงเร็วสุด {utc(first)})")
    if row["new_s"]:
        print(f"      ▶ พบคลื่น S/Lg ใหม่ {len(row['new_s'])} สถานี: {_names(row['new_s'])}")


# ----------------------------------------------------------------------
#  ลำดับเวลาที่ระบบตรวจพบ  (ขั้น 1-5 ของ pipeline)
# ----------------------------------------------------------------------
def milestones(rows, sites, phases):
    """return list ของ dict(stage, what, t, detail)   t = เวลาที่ระบบรู้ (เวลาจำลอง)"""
    out = []

    def add(stage, what, t, detail=""):
        if t is not None:
            out.append(dict(stage=stage, what=what, t=float(t), detail=detail))

    for key, stage, label in (("P", "1 คลื่น P", "คลื่น P"), ("S", "2 คลื่น S/Lg", "คลื่น S/Lg")):
        seen = sorted((v[key]["seen"], v[key]["first"], c) for c, v in phases.items() if key in v)
        if not seen:
            out.append(dict(stage=stage, what=f"ไม่พบ{label}", t=None, detail=""))
            continue
        t, tp, c = seen[0]
        ph = phases[c][key]["phase"]
        add(stage, "สถานีแรก", t, f"{c} ({ph} ถึง {utc(tp)} UTC, ระบบพบหลังคลื่นถึง {t - tp:.0f} s)")
        if key == "P" and len(seen) >= 3:
            add(stage, "ครบ 3 สถานี (พอหาตำแหน่ง)", seen[2][0], ", ".join(x[2] for x in seen[:3]))
        add(stage, "ทั้งหมด", seen[-1][0], f"พบ {len(seen)} สถานี (รอบสุดท้ายที่พบสถานีใหม่)")

    locs = [r for r in rows if r["lat"] is not None]
    if locs:
        r = locs[0]
        add("3 ตำแหน่ง", "ครั้งแรก", r["t"], f"lat {r['lat']:.2f}, lon {r['lon']:.2f}, RMS {r['rms']:.1f} s"
            + (f", คลาด {r['err_km']:.0f} km" if r.get("err_km") is not None else "")
            + (" [teleseism]" if r["tele"] else ""))
        fin = locs[-1]
        k = len(locs)                                           # ตำแหน่งนิ่ง: ทุกรอบหลังจากนี้ห่างผลสุดท้าย <= STABLE_KM
        while k > 0 and haversine_km(locs[k - 1]["lat"], locs[k - 1]["lon"], fin["lat"], fin["lon"]) <= STABLE_KM:
            k -= 1
        if k < len(locs):
            r = locs[k]
            add("3 ตำแหน่ง", f"นิ่ง (ห่างผลสุดท้าย ≤ {STABLE_KM:g} km)", r["t"],
                f"lat {r['lat']:.2f}, lon {r['lon']:.2f}" + (f", คลาด {r['err_km']:.0f} km" if r.get("err_km") is not None else ""))
    else:
        out.append(dict(stage="3 ตำแหน่ง", what="หาตำแหน่งไม่ได้", t=None, detail=""))

    mags = [r for r in rows if r["mw"] is not None]
    if mags:
        for r, what in ((mags[0], "ครั้งแรก"), (next((x for x in mags if x["mtype"] == "Ms_BB"), None), "เปลี่ยนเป็น Ms_BB"),
                        (mags[-1], "ล่าสุด")):
            if r is not None:
                add("4 ขนาด", what, r["t"], f"{r['mtype']} -> Mw {r['mw']:.1f} ±{r['unc']:.1f} ({r['n_mag']} สถานี)"
                    + (f", คลาด {r['dm']:+.1f}" if r.get("dm") is not None else ""))
    else:
        out.append(dict(stage="4 ขนาด", what="ไม่ได้ขนาด", t=None, detail=""))

    pred = next((r for r in rows if r["mmi_max"] is not None), None)
    if pred:
        add("5 ความรุนแรง", "ทำนาย MMI ครั้งแรก", pred["t"], f"สูงสุด {pred['mmi_max']:.1f} ที่ {pred['mmi_site']}")
    elif mags:
        out.append(dict(stage="5 ความรุนแรง", what="ไม่ทำนาย MMI (teleseism)", t=None, detail=""))
    app = [(s["alert"]["t"], s["name"]) for s in sites.values() if s["alert"]]
    pub = [((s["alert"] if s["alert"]["level"] == 2 else s["upgrade"])["t"], s["name"]) for s in sites.values()
           if s["alert"] and (s["alert"]["level"] == 2 or s["upgrade"])]
    if app:
        t, c = min(app)
        add("5 ความรุนแรง", "แจ้งเตือนครั้งแรก", t, f"{c}  (รวม {len(app)} ตำแหน่ง)")
    if pub:
        t, c = min(pub)
        add("5 ความรุนแรง", "แจ้งเตือนสาธารณะครั้งแรก", t, f"{c}  (รวม {len(pub)} ตำแหน่ง)")
    if not app:
        out.append(dict(stage="5 ความรุนแรง", what="ไม่มีการแจ้งเตือน", t=None, detail=""))
    felt = sorted((s["obs"]["t_peak"], s["name"], s["obs"]["mmi"]) for s in sites.values()
                  if s["kind"] == "station" and s["obs"] and s["obs"].get("usable", True) and s["obs"]["mmi"] >= 2.5)
    if felt:
        t, c, m = felt[0]
        add("5 ความรุนแรง", "สถานีแรกที่วัดได้ MMI ≥ 2.5", t, f"{c} MMI {m:.1f} (เวลาที่สั่นแรงสุด)")
    return out


def milestone_lines(ms, t_start, t0, t0_label):
    out = [f"ลำดับเวลาที่ระบบตรวจพบ  (เวลา = เวลาที่ระบบรู้ผล, T0 = {t0_label})",
           f"  {'ขั้น':<14}{'เหตุการณ์':<34}{'เวลา UTC':>9}  {'หลัง T0':>8}  {'จากต้นไฟล์':>10}  รายละเอียด"]
    for m in ms:
        if m["t"] is None:
            out.append(f"  {m['stage']:<14}{m['what']:<34}{'-':>9}  {'-':>8}  {'-':>10}")
            continue
        dt0 = f"{m['t'] - t0:+.0f} s" if t0 else "-"
        early = "  ⚠ ก่อนเวลาเกิด = pick จาก noise / สัญญาณผิด" if t0 and m["t"] < t0 else ""
        out.append(f"  {m['stage']:<14}{m['what']:<34}{utc(m['t']):>9}  {dt0:>8}  {m['t'] - t_start:>8.0f} s  "
                   f"{m['detail']}{early}")
    return out


def write_milestones(path, ms, t_start, t0):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["stage", "event", "utc", "after_origin_s", "from_file_start_s", "detail"])
        for m in ms:
            t = m["t"]
            w.writerow([m["stage"], m["what"], utc(t) if t else "", _round(t - t0, 1) if t and t0 else "",
                        _round(t - t_start, 1) if t else "", m["detail"]])


def write_phases(path, phases, sites, last_hypo):
    """เวลาคลื่นถึง (pick ล่าสุด) และเวลาที่ระบบพบครั้งแรก ของแต่ละสถานี เรียงตามเวลาพบคลื่น P"""
    rows = []
    for code, s in sites.items():
        if s["kind"] != "station":
            continue
        v = phases.get(code, {})
        p, sw = v.get("P"), v.get("S")
        d = haversine_km(last_hypo["lat"], last_hypo["lon"], s["lat"], s["lon"]) if last_hypo else None
        rows.append((p["seen"] if p else np.inf, code, d, p, sw, s))
    rows.sort(key=lambda x: (x[0], x[1]))
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["station", "dist_km", "p_arrival_utc", "p_detected_utc", "p_delay_s", "s_phase", "s_arrival_utc",
                    "s_detected_utc", "s_delay_s", "mmi_observed", "peak_utc", "first_alert_utc", "alert_level"])
        for _, code, d, p, sw, s in rows:
            o, a = s["obs"], s["alert"]
            w.writerow([code, _round(d, 1),
                        utc(p["time"]) if p else "", utc(p["seen"]) if p else "", _round(p["seen"] - p["time"], 1) if p else "",
                        sw["phase"] if sw else "", utc(sw["time"]) if sw else "", utc(sw["seen"]) if sw else "",
                        _round(sw["seen"] - sw["time"], 1) if sw else "",
                        _round((o or {}).get("mmi"), 2), utc(o["t_peak"]) if o else "",
                        utc(a["t"]) if a else "", a["level"] if a else ""])


def write_timeline(path, rows, t_start):
    keys = ["utc", "t_from_start_s", "n_stations_with_data", "status", "n_p_detected", "n_s_detected", "new_p",
            "new_s", "n_p_used", "lat", "lon", "origin_utc", "rms_s", "teleseism", "mag_type", "mw", "mag_raw",
            "uncertainty", "reliability", "n_mag_stations", "mmi_max_pred", "mmi_max_site", "n_alert_app",
            "n_alert_public", "compute_s", "over_tick", "err_km", "origin_err_s", "mw_err", "error"]
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in rows:
            w.writerow([utc(r["t"]), round(r["t"] - t_start, 1), r["n_data"], r["status"], r["n_pdet"], r["n_sdet"],
                        " ".join(r["new_p"]), " ".join(r["new_s"]), r["n_p"],
                        _round(r["lat"], 3), _round(r["lon"], 3), utc(r["t0"]) if r["t0"] else "", _round(r["rms"], 2),
                        int(r["tele"]), r["mtype"] or "", _round(r["mw"], 2), _round(r["raw"], 2), _round(r["unc"], 2),
                        r["level"] or "", r["n_mag"], _round(r["mmi_max"], 2), r["mmi_site"] or "",
                        r["n_app"], r["n_pub"], round(r["compute_s"], 2), int(r["over"]),
                        _round(r.get("err_km"), 1), _round(r.get("dt0"), 1), _round(r.get("dm"), 2), r["error"] or ""])


def _round(v, n):
    return "" if v is None else round(float(v), n)


def write_alerts(path, sites, last_hypo):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["site", "kind", "first_alert_utc", "level", "alert", "M_then", "mag_type_then", "mmi_pred_then",
                    "upgrade_public_utc", "s_arrival_utc", "lead_before_S_s", "lead_before_Lg_s", "mmi_observed"])
        for s in sites.values():
            a, up = s["alert"], s["upgrade"]
            if a is None:
                continue
            lead_s = lead_lg = t_s = None
            if last_hypo is not None:
                x = TT.site_warning(last_hypo, s["lat"], s["lon"], last_hypo["origin_time"])
                t_s, lead_s, lead_lg = x["s_arrival"], x["s_arrival"] - a["t"], x["lg_arrival"] - a["t"]
            w.writerow([s["name"], s["kind"], utc(a["t"]), a["level"], a["name"], round(a["M"], 2), a["type"],
                        round(a["mmi"], 2), utc(up["t"]) if up else "", utc(t_s) if t_s else "", _round(lead_s, 1),
                        _round(lead_lg, 1), _round((s["obs"] or {}).get("mmi"), 2)])


def origin_for_report(rows, ref):
    """T0 ที่ใช้บอกเวลา "หลังเกิดเหตุ": ค่าอ้างอิงถ้ามี ไม่มีใช้เวลาเกิดที่ระบบคำนวณได้รอบสุดท้าย"""
    if ref and ref["t0"]:
        return ref["t0"], f"เวลาเกิดอ้างอิง {utc(ref['t0'])} UTC"
    last = next((r for r in reversed(rows) if r["t0"] is not None), None)
    if last and not last["tele"]:
        return last["t0"], f"เวลาเกิดที่ระบบคำนวณ {utc(last['t0'])} UTC"
    return None, "ไม่ทราบเวลาเกิด"


def summary_lines(event, rows, sites, notes, t_start, args, ref, ms):
    out = [f"จำลองสถานการณ์จริง: {event}",
           f"รอบละ {args.tick:g} s, หน่วง {args.latency:g} s, สถานีขาด {args.drop:.0%}, "
           f"ปรับด้วยค่าที่วัดได้ = {'ใช่' if not args.no_correction else 'ไม่'}",
           f"สถานีในการจำลอง {sum(1 for s in sites.values() if s['kind'] == 'station')}, จำนวนรอบ {len(rows)}"]
    if ref:
        out.append(f"ค่าอ้างอิง (ใช้เทียบผลเท่านั้น): {ref['name']}  {ref['label']}  lat {ref['lat']}, lon {ref['lon']}"
                   + (f"  เวลาเกิด {utc(ref['t0'])} UTC" if ref["t0"] else ""))
    t0_ref = ref["t0"] if ref else None
    t0, t0_label = origin_for_report(rows, ref)
    out.append("")
    out.extend(milestone_lines(ms, t_start, t0, t0_label))

    if t0_ref:                                                  # ก่อนเวลาเกิดจริง ทุกอย่างที่ระบบรายงานคือสัญญาณเตือนผิด
        pre = [r for r in rows if r["t"] < t0_ref]
        early = [s["name"] for s in sites.values() if s["alert"] and s["alert"]["t"] < t0_ref]
        out.append(f"  ก่อนเวลาเกิดจริง    : {len(pre)} รอบ, รายงานตำแหน่ง {sum(1 for r in pre if r['lat'] is not None)} รอบ, "
                   f"รายงานขนาด {sum(1 for r in pre if r['mw'] is not None)} รอบ, แจ้งเตือน {len(early)} ตำแหน่ง")

    last = next((r for r in reversed(rows) if r["hypo"] is not None), None)
    last_mag = next((r for r in reversed(rows) if r["mw"] is not None), None)
    out.append("")
    out.append("ผลล่าสุดเมื่อจบการจำลอง")
    if last:
        out.append(f"  ตำแหน่ง lat {last['lat']:.3f}, lon {last['lon']:.3f}, RMS {last['rms']:.2f} s"
                   + (f", คลาด {last['err_km']:.0f} km, เวลาเกิดคลาด {last['dt0']:+.1f} s" if last.get("err_km") is not None else "")
                   + ("  [ตรวจเป็น teleseism]" if last["tele"] else ""))
    if last_mag:
        out.append(f"  ขนาด {last_mag['mtype']} -> Mw {last_mag['mw']:.2f} ±{last_mag['unc']:.1f} "
                   f"ความน่าเชื่อถือ{last_mag['level']}, ใช้ {last_mag['n_mag']} สถานี"
                   + (f", คลาด {last_mag['dm']:+.2f}" if last_mag.get("dm") is not None else ""))

    st_sites = [s for s in sites.values() if s["kind"] == "station"]
    usable = [s for s in st_sites if s["obs"] and s["obs"].get("usable", True)]
    alerted = [s for s in usable if s["alert"]]
    felt = [s for s in usable if s["obs"]["mmi"] >= 2.5]
    out.append("")
    out.append("การแจ้งเตือนเทียบกับค่าที่สถานีวัดได้ (MMI วัด >= 2.5)")
    out.append(f"  สถานีที่มีค่าวัดผ่าน QC {len(usable)}, สั่นถึงเกณฑ์ {len(felt)}, ได้รับการแจ้งเตือน {len(alerted)}")
    out.append(f"  เตือนและสั่นจริง {sum(1 for s in felt if s['alert'])}, สั่นแต่ไม่ได้เตือน "
               f"{sum(1 for s in felt if not s['alert'])}, เตือนแต่สั่นไม่ถึงเกณฑ์ {sum(1 for s in alerted if s['obs']['mmi'] < 2.5)}")
    if last:
        for s in sites.values():
            if s["kind"] != "target":
                continue
            x = TT.site_warning(last["hypo"], s["lat"], s["lon"], last["hypo"]["origin_time"])
            if s["alert"]:
                out.append(f"  {s['name']}: {s['alert']['name']} เมื่อ {utc(s['alert']['t'])} UTC, MMI ทำนาย "
                           f"{s['alert']['mmi']:.1f}, ก่อนคลื่น S {x['s_arrival'] - s['alert']['t']:+.0f} s, "
                           f"ก่อนคลื่น Lg {x['lg_arrival'] - s['alert']['t']:+.0f} s")
            else:
                out.append(f"  {s['name']}: ไม่ได้แจ้งเตือน")

    cost = np.array([r["compute_s"] for r in rows])
    over = int(sum(r["over"] for r in rows))
    out.append("")
    out.append("ความเร็วของการคำนวณบนเครื่องนี้ (ข้อมูลจริงมาทุก ๆ รอบ ต้องคำนวณให้ทัน)")
    out.append(f"  เฉลี่ย {cost.mean():.1f} s, สูงสุด {cost.max():.1f} s ต่อรอบ, ช้ากว่ารอบ {args.tick:g} s จำนวน {over}/{len(rows)} รอบ")
    errs = [r for r in rows if r["error"]]
    if errs:
        out.append("")
        out.append(f"รอบที่เกิดข้อผิดพลาด {len(errs)} รอบ  ตัวอย่าง: {errs[0]['error']}")
    if notes:
        out.append("")
        out.append("หมายเหตุการโหลดข้อมูล")
        out.extend("  " + n for n in notes)
    return out


# ======================================================================
def main():
    for stream in (sys.stdout, sys.stderr):
        with contextlib.suppress(AttributeError, ValueError):
            stream.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description="เล่นไฟล์เหตุการณ์เหมือนข้อมูลจริงจากสถานี")
    ap.add_argument("event", nargs="?", default="2025-03-29_M5.9_AllTMD", help="ชื่อโฟลเดอร์ใน station/")
    ap.add_argument("--tick", type=float, default=10.0, help="ระบบตัดสินใจทุกกี่วินาที (เวลาจำลอง)")
    ap.add_argument("--duration", type=float, default=None, help="จำลองกี่วินาทีจากต้นไฟล์ (ไม่ใส่ = จนหมดไฟล์)")
    ap.add_argument("--speed", type=float, default=0.0, help="ความเร็วการเล่น 1 = ตามเวลาจริง, 10 = เร็ว 10 เท่า, 0 = เร็วที่สุด")
    ap.add_argument("--latency", type=float, default=0.0, help="ข้อมูลทุกสถานีมาช้ากี่วินาที")
    ap.add_argument("--drop", type=float, default=0.0, help="สัดส่วนสถานีที่ขาดหาย 0-1")
    ap.add_argument("--seed", type=int, default=1, help="ค่าสุ่มของ --drop")
    ap.add_argument("--no-correction", action="store_true", help="ไม่ปรับ MMI ด้วยค่าที่สถานีวัดได้")
    ap.add_argument("--no-map", action="store_true", help="ไม่สร้างแผนที่")
    ap.add_argument("--no-open", action="store_true", help="สร้างแผนที่แต่ไม่เปิดอัตโนมัติ")
    ap.add_argument("--all-maps", action="store_true",
                    help="สร้างแผนที่ MMI ทำนาย + วงกลม (_intensity_map.html) และภาพนิ่ง (_intensity_mmi_map.png) เพิ่ม")
    args = ap.parse_args()

    event = os.path.basename(args.event.rstrip("/\\"))
    print(f"กำลังอ่านข้อมูล {event} ...", flush=True)
    records, notes = load_event(args.event, args.drop, args.seed)
    for n in notes:
        print("⚠", n)
    feed = ReplayFeed(records, args.latency)
    sites = make_sites(records)
    ref = find_reference(args.event)
    t_start = min(r["t0"] for r in records.values())
    t_stop = max(r["t0"] + len(r["z"]) / r["fs"] for r in records.values()) + args.latency
    if args.duration:
        t_stop = min(t_stop, t_start + args.duration)
    print(f"สถานี {len(records)}   ช่วงเวลา {utc(t_start)} - {utc(t_stop)} UTC   รอบละ {args.tick:g} s", flush=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    base = os.path.join(OUT_DIR, event)
    rows, phases = [], {}
    final_res = None
    wall_start = time.perf_counter()
    print_header()
    with replay_patches(feed, use_correction=not args.no_correction):
        now = t_start + args.tick
        while now <= t_stop:
            feed.now = now
            t1 = time.perf_counter()
            row = run_tick(feed, sites, phases)
            if row["res"] is not None:
                final_res = row["res"]
            row["res"] = None
            row["compute_s"] = time.perf_counter() - t1
            row["over"] = row["compute_s"] > args.tick
            row["n_app"] = sum(1 for s in sites.values() if s["alert"] is not None)
            row["n_pub"] = sum(1 for s in sites.values() if s["alert"] is not None and
                               (s["alert"]["level"] == 2 or s["upgrade"] is not None))
            score_row(row, ref)
            rows.append(row)
            print_row(row, t_start, row["n_app"], row["n_pub"], phases)
            if args.speed > 0:                                  # รอให้ตรงกับเวลาจำลอง
                target = wall_start + (now - t_start) / args.speed
                time.sleep(max(target - time.perf_counter(), 0.0))
            now += args.tick
        print("\nกำลังสร้างแผนที่จากผลรอบสุดท้าย ...", flush=True)
        maps = [] if args.no_map else make_maps(final_res, sites, base, args.all_maps)

    last = next((r for r in reversed(rows) if r["hypo"] is not None), None)
    ms = milestones(rows, sites, phases)
    t0, _ = origin_for_report(rows, ref)
    write_timeline(base + "_timeline.csv", rows, t_start)
    write_alerts(base + "_alerts.csv", sites, last["hypo"] if last else None)
    write_phases(base + "_phases.csv", phases, sites, last["hypo"] if last else None)
    write_milestones(base + "_milestones.csv", ms, t_start, t0)
    lines = summary_lines(event, rows, sites, notes, t_start, args, ref, ms)
    if maps:
        lines += ["", "แผนที่ (จากผลรอบสุดท้ายที่ระบบได้ขนาด)"]
        lines += [f"  {name}: {os.path.relpath(f, CODE_DIR)}" for name, f in maps]
    with open(base + "_summary.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n" + "=" * 72)
    print("\n".join(lines))
    print("=" * 72)
    print(f"บันทึกผล: {base}_summary.txt, _milestones.csv, _phases.csv, _timeline.csv, _alerts.csv")
    if maps and not args.no_open:
        open_maps(maps)


def open_maps(maps):
    """เปิดแผนที่ HTML ในเบราว์เซอร์ และภาพ PNG ด้วยโปรแกรมดูภาพของเครื่อง"""
    import webbrowser
    for name, f in maps:
        path = os.path.abspath(f)
        try:
            if path.lower().endswith(".html") or not hasattr(os, "startfile"):
                webbrowser.open(Path(path).as_uri())
            else:
                os.startfile(path)                              # Windows: เปิดด้วยโปรแกรมตั้งต้น
            print(f"เปิด{name}: {path}")
        except Exception as ex:                                 # noqa: BLE001
            print(f"⚠ เปิด{name}ไม่ได้ ({type(ex).__name__}) เปิดเองที่ {path}")


if __name__ == "__main__":
    main()
