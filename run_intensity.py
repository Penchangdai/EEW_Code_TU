"""
run_intensity.py
ประมาณความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) ที่แต่ละสถานีและตำแหน่งที่สนใจ
แล้วแจ้งเตือนก่อนคลื่นมาถึง  ใช้ผลชุดเดียวกับ run_magnitude.py และ run_travel_time.py

  ขั้น 1  ตำแหน่ง เวลาเกิด ขนาด และเวลาที่ระบบรู้ข้อมูลแต่ละขั้น      (run_travel_time.analyse / system_stages)
  ขั้น 2  จำลองการทำงานแบบเวลาจริง ทุก UPDATE_INTERVAL_S วินาที:
            - ขนาด M(t) จากข้อมูลที่มาถึงสถานีแล้วเท่านั้น (เช่น ช่วงแรกมีแค่ Pd/τc จากคลื่น P 3 s)
            - ตำแหน่ง: จากคลื่น P อย่างเดียว ก่อนได้คลื่น S/Lg, ตำแหน่งสุดท้ายหลังจากนั้น
            - ทำนาย MMI ที่ทุกตำแหน่ง (intensity.predict_mmi)
            - แจ้งเตือนตำแหน่งที่เข้าเกณฑ์ ShakeAlert (M >= 4.5 และ MMI >= III / M >= 5 และ MMI >= IV)
              ครั้งแรกที่เข้าเกณฑ์ = เวลาส่งคำเตือน   Lead time = เวลาคลื่น S (และ Lg) มาถึง − เวลาส่งคำเตือน
  ขั้น 3  ความรุนแรงที่วัดได้จริงที่สถานี (MMI จาก PGV ตาม Worden 2012, JMA I) เทียบกับค่าทำนาย

  OUTPUT_MODE = "public" : การ์ดแจ้งเตือนของแต่ละตำแหน่ง (รายละเอียดบันทึกใน DETAIL_FILE)
  OUTPUT_MODE = "detail" : ลำดับเวลาการอัปเดตและตารางเทียบทั้งหมด
  กราฟ: intensity_alert.png (เวลาแจ้งเตือน vs เวลาคลื่นมาถึง), intensity_distance.png (MMI กับระยะ)

รัน:  กด Run ใน VS Code
"""
import contextlib
import csv
import io
import os

import numpy as np

import config as C
import magnitude as M
import run_magnitude as RM
import run_travel_time as RTT
import travel_time as TT
import intensity as INT
from signal_utils import haversine_km, acc_to_vel, remove_offset

try:
    from station_db import get_vs30
except ImportError:
    get_vs30 = None
try:
    from station_db import STATIONS_CSV
except ImportError:
    STATIONS_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stations.csv")

# ======================================================================
#  ตั้งค่าตรงนี้
# ======================================================================
SELECTED_STATIONS = "auto"    # "auto" = ทุกไฟล์ในโฟลเดอร์ station, ชื่อโฟลเดอร์เหตุการณ์ เช่น "2026-05-18_M5.2_AllTMD", หรือรหัส ["BKKA", "CMMT"]
PROFILE = "regional"
EPICENTER = None              # None = คำนวณเอง (แนะนำ)  หรือ (lat, lon, ความลึก km, "เวลาเกิด UTC")
DEPTH_KM = 10.0
EQ_TYPE = "crustal"           # "crustal" (ในเปลือกโลก), "interplate" (ระหว่างแผ่น), "intraplate" (ภายในแผ่น)

# Vs30 (m/s) ของสถานี  ลำดับที่ใช้: STATION_VS30 -> คอลัมน์ vs30 ใน stations.csv -> DEFAULT_VS30 ใน config.py
STATION_VS30 = {}             # เช่น {"BKKA": 180}
# ประเภทสถานี: "csv" = คอลัมน์ table ใน stations.csv (1 = ตรวจวัดแผ่นดินไหว, 2 = วัดอัตราเร่ง) ถ้าไม่มีใช้จากไฟล์
#              "file" = จากบรรทัด Channels ในไฟล์ข้อมูล (HH?/EH? = ตรวจวัดแผ่นดินไหว, มีแค่ HN? = วัดอัตราเร่ง)
STATION_TYPE_SOURCE = "csv"

# แนวรอยแตกของแผ่นดินไหว (finite fault) [(lat, lon), ...]  None = ถือเป็นจุดกำเนิด (ใช้ระยะถึงจุดศูนย์กลาง)
#   ถ้ารู้แนวรอยแตก (เช่น จาก USGS finite fault หรือรอยเลื่อนที่รู้ว่าแตก) ใส่ไว้ จะใช้ระยะถึงรอยแตกแทน
#   ตัวอย่าง (โดยประมาณ) รอยเลื่อนสะกาย 28 มี.ค. 2568 ยาว ~475 km เหนือจุดศูนย์กลาง ~100 km ใต้ ~375 km:
#   FAULT_TRACE = [(22.90, 95.98), (22.01, 95.94), (21.0, 96.03), (20.0, 96.15), (19.0, 96.25), (18.63, 96.27)]
FAULT_TRACE = None

# ตำแหน่งที่ต้องการแจ้งเตือน (พิกัดเป็นค่าโดยประมาณ)   Vs30 ดินอ่อนกรุงเทพฯ/ปริมณฑล ประมาณ 150-250 m/s
TARGET_SITES = [
    {"name": "มธ. รังสิต", "lat": 14.07, "lon": 100.61, "vs30": 250},
]
ALERT_STATIONS = True         # แจ้งเตือนที่ตำแหน่งสถานีด้วย

UPDATE_INTERVAL_S = 5.0       # ระบบอัปเดตขนาด/ความรุนแรงทุกกี่วินาที
# ปรับค่าทำนายด้วยความรุนแรงที่สถานีวัดได้แล้ว (แนวคิดเดียวกับ ShakeMap / PLUM): เมื่อสถานีอย่างน้อย
# OBS_MIN_STATIONS แห่งสั่นแรงสุดผ่านไปแล้ว ใช้ค่ากลางของ (MMI วัด − MMI ทำนาย) ปรับตำแหน่งที่คลื่นยังไม่ถึง
USE_OBS_CORRECTION = True
OBS_MIN_STATIONS = 2
OBS_MAX_CORRECTION = 1.5      # ปรับได้ไม่เกิน ± MMI (เฉพาะวิธีค่ากลาง)
# วิธีปรับ: "auto" = ลองทุกวิธีแล้วเลือกวิธีที่ leave-one-out ได้ความคลาดเคลื่อน (RMSE) ต่ำสุด
#           หรือกำหนดเอง "none", "median", "type", "trend", "median+krig", "trend+krig"
CORRECTION_METHOD = "auto"
KRIG_LENGTHS_KM = (25.0, 50.0, 100.0, 200.0, 400.0)   # ระยะสหสัมพันธ์ของ residual ที่ลอง (km)
KRIG_NUGGET = (0.2, 0.5)                        # สัดส่วน residual ที่ไม่สัมพันธ์กับตำแหน่ง (ความผันแปรเฉพาะที่)
OUTLIER_MAD = 3.0             # residual ต่างจากค่ากลางเกิน 3 x MAD (≈ 3σ) และ
OUTLIER_MIN_MMI = 2.0         # เกิน 2 MMI = ค่าผิดปกติ ไม่ใช้ปรับ (แค่ต่างจากสถานีอื่นเล็กน้อยไม่นับ)
GROSS_OUTLIER_MMI = 4.0       # ต่างจากค่ากลางเกิน 4 MMI = ข้อมูลเสีย ไม่ใช้เป็นค่าที่วัดได้เลย
# ตรวจคุณภาพข้อมูล (QC) ก่อนวัด PGA / PGV
QC_SPIKE_RATIO = 15.0         # จุดที่ |a| > 15 เท่าของค่ากลางรอบ ๆ (±1 s) = สัญญาณกระโดด (spike) ตัดทิ้ง
QC_MIN_SNR = 3.0              # PGA / ค่าสูงสุดของ noise ก่อนคลื่น P ต่ำกว่านี้ = แยกจาก noise ไม่ได้
QC_CLIP_SAMPLES = 5           # ค่ายอดเท่ากันติดกันตั้งแต่จำนวนนี้ = สัญญาณขาด (clip) ค่าสูงสุดต่ำกว่าจริง
EVENT_END_VMIN = 2.0          # วัดถึงเวลา T0 + Δ / 2.0 km/s + EVENT_TAIL_S (คลื่นผิวพื้นผ่านไปแล้ว)
EVENT_TAIL_S = 120.0
OUTPUT_MODE = "public"        # "public" หรือ "detail"
DETAIL_FILE = "intensity_detail.txt"
PLOT = True
MAKE_MAP = True               # intensity_map.html (รายละเอียด + การแจ้งเตือน) และ
                              # intensity_station_map.html (MMI ที่วัดได้ของแต่ละสถานีที่ input บนพื้นสี MMI)


# ======================================================================
#  ข้อมูลพื้นฐาน
# ======================================================================
# ประเภทสถานี (data_loader.detect_station_type: ดูจากช่องที่สถานีมี)
SENSOR_NAME = {"seismic": "สถานีตรวจวัดแผ่นดินไหว (Seismic station)",
               "strong_motion": "สถานีตรวจวัดอัตราเร่งของพื้นดิน (Strong-motion station)"}
SENSOR_SHORT = {"seismic": "▲ ตรวจวัดแผ่นดินไหว", "strong_motion": "■ วัดอัตราเร่ง"}
_SENSOR_ALIAS = {"velocity": "seismic", "acceleration": "strong_motion"}


def sensor_key(s):
    k = s.get("sensor") if isinstance(s, dict) else s
    k = _SENSOR_ALIAS.get(k, k)
    return k if k in SENSOR_NAME else "strong_motion"


def sensor_name(s, short=False):
    """ประเภทสถานี: seismic = มีเครื่องวัดแผ่นดินไหว (HH?/EH?), strong_motion = มีแค่เครื่องวัดความเร่ง (HN?)"""
    return (SENSOR_SHORT if short else SENSOR_NAME)[sensor_key(s)]


def obs_mmi(s):
    """MMI ที่วัดได้ ถ้าผ่าน QC  (ไม่ผ่าน = None แสดงเป็น "-" บนแผนที่)"""
    o = s.get("obs")
    return o["mmi"] if o and o.get("usable", True) else None


_CSV_CACHE = None


def station_csv(code):
    """แถวของสถานีใน stations.csv (code, lat, lon, table, vs30)  หรือ {} ถ้าไม่มี"""
    global _CSV_CACHE
    if _CSV_CACHE is None:
        _CSV_CACHE = {}
        try:
            with open(STATIONS_CSV, encoding="utf-8-sig", newline="") as f:
                for r in csv.DictReader(f):
                    if r.get("code"):
                        _CSV_CACHE[r["code"].strip().upper()] = r
        except OSError:
            pass
    return _CSV_CACHE.get(str(code).upper(), {})


def station_vs30(code):
    """Vs30: STATION_VS30 -> station_db.get_vs30 -> คอลัมน์ vs30 ใน stations.csv -> DEFAULT_VS30"""
    if code in STATION_VS30:
        return float(STATION_VS30[code])
    v = get_vs30(code) if get_vs30 else None
    if not v:
        v = station_csv(code).get("vs30")
    try:
        v = float(v)
    except (TypeError, ValueError):
        v = None
    return v if v and v > 0 else float(C.DEFAULT_VS30)


def station_sensor(code, from_file):
    """ประเภทสถานี (seismic / strong_motion) ตาม STATION_TYPE_SOURCE"""
    if STATION_TYPE_SOURCE == "csv":
        t = str(station_csv(code).get("table", "")).strip()
        if t in ("1", "2"):
            return "seismic" if t == "1" else "strong_motion"
    return from_file


def site_rjb(lat, lon):
    """ระยะในแนวราบถึงแนวรอยแตก (km) ถ้ากำหนด FAULT_TRACE  ไม่งั้น None (ใช้จุดกำเนิด)"""
    return INT.distance_to_trace_km(lat, lon, FAULT_TRACE) if FAULT_TRACE else None


def _trace_dist_grid(LA, LO):
    """ระยะถึงแนวรอยแตกของทุกจุดในตาราง (vectorized)"""
    tr = np.asarray(FAULT_TRACE, float)
    f = np.linspace(0.0, 1.0, 50)[:, None]
    pts = np.vstack([a + f * (b - a) for a, b in zip(tr[:-1], tr[1:])])
    out = np.full(LA.shape, np.inf)
    for la, lo in pts:
        out = np.minimum(out, haversine_km(la, lo, LA, LO))
    return out


def despike(x, fs, ratio=None):
    """
    ตัดสัญญาณกระโดด (spike): จุดที่ |x| > ratio เท่าของค่ากลาง |x| รอบ ๆ ±1 s
    (แรงสั่นจริงจากแผ่นดินไหว ค่าสูงสุดมักไม่เกิน ~3-5 เท่าของค่ากลางรอบ ๆ)  แทนด้วยเส้นตรงระหว่างจุดข้างเคียง
    return: (สัญญาณที่แก้แล้ว, จำนวนจุดที่ตัด)
    """
    from scipy.ndimage import median_filter
    ratio = QC_SPIKE_RATIO if ratio is None else ratio
    x = np.asarray(x, float)
    w = max(int(2 * fs) | 1, 3)
    ax = np.abs(x)
    med = median_filter(ax, size=w, mode="nearest")
    after = median_filter(ax, size=w, mode="nearest", origin=-(w // 2))   # ค่ากลางช่วง 2 s ถัดไป
    floor = np.median(ax) + 1e-12
    # spike = สูงกว่าทั้งรอบ ๆ และช่วงถัดไป (onset ของคลื่นจริงตามด้วยพลังงานต่อเนื่อง จึงไม่ถูกตัด)
    bad = ax > ratio * np.maximum(np.maximum(med, after), floor)
    if not bad.any():
        return x, 0
    bad = np.convolve(bad, np.ones(3), "same") > 0                  # รวมจุดข้างเคียง
    idx = np.arange(len(x))
    y = x.copy()
    y[bad] = np.interp(idx[bad], idx[~bad], x[~bad])
    return y, int(bad.sum())


def is_clipped(x, n_min=None):
    """
    สัญญาณขาด (clipped): เครื่องวัดถึงค่าสูงสุดที่วัดได้ ค่ายอดจึง "เท่ากันเป๊ะ" ติดกันหลายจุด
    (ยอดของคลื่นจริงค่าเปลี่ยนทุก sample)  ตรวจ: มีช่วงที่ |x| >= 99.9% ของค่าสูงสุดและค่าเท่ากัน >= n_min จุดติดกัน
    """
    n_min = QC_CLIP_SAMPLES if n_min is None else n_min
    x = np.asarray(x, float)
    if len(x) < n_min + 1:
        return False
    m = np.max(np.abs(x))
    if m <= 0:
        return False
    flat = (np.abs(x[1:]) >= 0.999 * m) & (x[1:] == x[:-1])
    run = best = 0
    for f in flat:
        run = run + 1 if f else 0
        best = max(best, run)
    return best + 1 >= n_min


def observed(st, hypo=None):
    """
    ความรุนแรงที่วัดได้ที่สถานี: PGA, PGV แนวราบ, MMI (Worden 2012), JMA I และเวลาที่สั่นแรงสุด
    QC 4 ขั้น (เพื่อไม่ให้ noise / สัญญาณผิดปกติกลายเป็น "แรงสั่น")
      0) สัญญาณขาด (clip): เครื่องวัดแผ่นดินไหว (HH/EH) ขาดง่ายเมื่อสั่นแรง  ไม่ใช้ค่าสูงสุดของสถานีนั้น
      1) วัดเฉพาะช่วงที่คลื่นแผ่นดินไหวผ่านสถานี: Tp(คาด) − 5 s ถึง T0 + Δ/2.0 + 120 s
         (ค่าสูงสุดก่อนคลื่น P มาถึงเป็นไปไม่ได้ทางกายภาพ)
      2) ตัดสัญญาณกระโดด (despike)
      3) SNR = PGA ในช่วงคลื่น / ค่าสูงสุดของ noise ก่อนคลื่น P   (< QC_MIN_SNR = ใช้ไม่ได้)
    """
    rec, fs = st["rec"], st["rec"]["fs"]
    n = len(rec["z"])
    if hypo is not None:
        tp = TT.p_arrival_time(hypo["origin_time"], float(np.hypot(st["delta"], hypo["depth"])),
                               depth_km=hypo["depth"])
        i0 = int(np.clip((tp - 5.0 - rec["t0"]) * fs, 0, n - 1))
        i1 = int(np.clip((hypo["origin_time"] + st["delta"] / EVENT_END_VMIN + EVENT_TAIL_S - rec["t0"]) * fs,
                         i0 + 1, n))
    else:
        i0, i1 = max(st["p_index"] - int(5 * fs), 0), n
    n_pre = max(i0 - int(1.0 * fs), 10)
    flags, n_spk = [], 0
    comp = {}
    raw = rec.get("raw") or {}                           # ข้อมูลก่อนแปลงหน่วย (ถ้า data_loader เก็บไว้)
    clipped = [k.upper() for k in ("n", "e") if is_clipped((raw.get(k) if k in raw else rec[k])[i0:i1])]
    if clipped:
        flags.append(f"สัญญาณขาด (clip) แกน {', '.join(clipped)} ค่าสูงสุดต่ำกว่าจริง")
    for k in ("n", "e", "z"):
        x = remove_offset(rec[k], n_pre)
        x, c = despike(x, fs)
        n_spk += c
        comp[k] = x
    if n_spk:
        flags.append(f"ตัดสัญญาณกระโดด {n_spk} จุด")
    win = slice(i0, i1)
    pga_h = float(max(np.max(np.abs(comp["n"][win])), np.max(np.abs(comp["e"][win]))))   # แกนราบที่ใหญ่กว่า
    vn, ve = (acc_to_vel(comp[k], fs, C.DISP_HIGHPASS_HZ, n_pre) for k in ("n", "e"))
    pgv_h = float(max(np.max(np.abs(vn[win])), np.max(np.abs(ve[win]))))
    pre = slice(max(i0 - int(60 * fs), 0), max(i0, 1))
    noise = float(max(np.max(np.abs(comp["n"][pre])), np.max(np.abs(comp["e"][pre])))) if i0 > fs else 0.0
    snr = pga_h / noise if noise > 0 else None
    if snr is not None and snr < QC_MIN_SNR:
        flags.append(f"SNR ต่ำ ({snr:.1f}) แยกแรงสั่นจาก noise ไม่ได้")
    zw, nw, ew = (comp[k][i0:i1] for k in ("z", "n", "e"))
    det = INT.jma_intensity_detail(zw, nw, ew, fs)
    mm = INT.observed_mmi(pga_h, pgv_h)
    return dict(pga=pga_h, pgv=pgv_h, jma=det["I"], jma_class=INT.jma_intensity_class(det["I"]),
                t_peak=rec["t0"] + (i0 + det["i_peak"]) / fs, snr=snr, n_spike=n_spk, flags=flags,
                usable=not (snr is not None and snr < QC_MIN_SNR) and not clipped, **mm)


def build_sites(res):
    """ตำแหน่งที่แจ้งเตือน: ทุกสถานีที่ input เข้าไป (+ ค่าที่วัดได้) และ TARGET_SITES"""
    sites = []
    if ALERT_STATIONS:
        for st in res["stations"]:
            sites.append(dict(name=st["code"], lat=st["lat"], lon=st["lon"], vs30=station_vs30(st["code"]),
                              kind="station", used=True, obs=observed(st, res["hypo"]),
                              sensor=station_sensor(st["code"], st["rec"].get("station_type")),
                              sensor_file=st["rec"].get("station_type"),
                              channels=st["rec"]["header"].get("Channels", "")))
        # สถานีที่ input เข้าไปแต่ไม่ได้ใช้คำนวณขนาด (ไม่พบคลื่น P / ข้อมูลไม่ครอบคลุม ฯลฯ) ยังแสดงค่าที่วัดได้
        have = {s["name"] for s in sites}
        for path in RM.station_files(SELECTED_STATIONS):
            try:
                rec = RM.load_station_txt(path)
            except Exception:                                   # noqa: BLE001
                continue
            code = rec["station"]
            c = RM.get_coords(code)
            if code in have or c is None or not rec.get("units_ok", True):
                continue
            try:
                dd = float(haversine_km(res["hypo"]["lat"], res["hypo"]["lon"], c[0], c[1]))
                obs = observed({"rec": rec, "p_index": int(10 * rec["fs"]), "delta": dd}, res["hypo"])
            except Exception:                                   # noqa: BLE001
                obs = None
            sites.append(dict(name=code, lat=c[0], lon=c[1], vs30=station_vs30(code), kind="station",
                              used=False, obs=obs, sensor=station_sensor(code, rec.get("station_type")),
                              sensor_file=rec.get("station_type"),
                              channels=rec["header"].get("Channels", "")))
            have.add(code)
    for s in TARGET_SITES:
        sites.append(dict(name=s["name"], lat=s["lat"], lon=s["lon"], vs30=float(s.get("vs30", C.DEFAULT_VS30)),
                          kind="target", used=False, obs=None))
    hypo = res["hypo"]
    for s in sites:                                          # เวลาคลื่นมาถึง (ตำแหน่งสุดท้าย)
        w = TT.site_warning(hypo, s["lat"], s["lon"], hypo["origin_time"])
        s.update(delta=w["epi_dist_km"], t_s=w["s_arrival"], t_lg=w["lg_arrival"], t_p=w["p_arrival"],
                 rjb=site_rjb(s["lat"], s["lon"]))
    return sites


# ======================================================================
#  ขั้น 2  จำลองเวลาจริง
# ======================================================================
def magnitude_at(res, t):
    """ขนาดหลักจากข้อมูลที่มาถึงสถานีแล้ว ณ เวลา t (epoch)  -> (pick, reliability) หรือ (None, None)"""
    T0 = res["hypo"]["origin_time"]
    sts = []
    for st in res["stations"]:
        mags = {k: v for k, v in st["mags"].items() if v and v.get("valid")
                and RTT.magnitude_ready_time(st, k, T0) + RTT.PROCESSING_DELAY_S <= t}
        sts.append(dict(st, mags=mags, w={}))
    with contextlib.redirect_stdout(io.StringIO()):
        net = M.network_by_type(sts)
        _, scores = M.type_scores(net)
        pick = M.choose_preferred(net, scores)
    if pick is None:
        return None, None
    used = [s for s in sts if s["w"].get(pick["type"], 0) > 0]
    rel = M.reliability(pick, net, used, res.get("gap"), True)
    return pick, rel


def simulate(res, stages, sites):
    hypo = res["hypo"]
    hp = hypo.get("hypo_p") or hypo
    t_final_loc = next((g["time"] for g in stages if g["name"].startswith("ขั้น 2")), stages[0]["time"])
    t_end = max(stages[-1]["time"], max(s["t_lg"] for s in sites)) + UPDATE_INTERVAL_S
    times = sorted(set(np.arange(stages[0]["time"], t_end, UPDATE_INTERVAL_S).tolist()
                       + [g["time"] for g in stages]))
    updates, last_key = [], None
    for s in sites:
        s.update(alert=None, upgrade=None)
    for t in times:
        pick, rel = magnitude_at(res, t)
        if pick is None:
            continue
        loc = hypo if t >= t_final_loc else hp
        unc = rel[2]
        preds = {}
        for s in sites:
            d = float(haversine_km(loc["lat"], loc["lon"], s["lat"], s["lon"]))
            preds[s["name"]] = INT.predict_mmi(pick["mag"], d, loc["depth"], s["vs30"], EQ_TYPE, unc, s.get("rjb"))
        cm = correction_model(sites, preds, t)
        n_obs = cm["n"] if cm else 0
        corrs = {s["name"]: site_corr(cm, s, preds[s["name"]]["R"]) for s in sites}
        corr = float(np.median(list(corrs.values()))) if corrs else 0.0
        key = (pick["type"], round(pick["mag"], 1), loc is hypo, (cm or {}).get("method"), round(corr, 1))
        row = dict(t=t, type=pick["type"], M=pick["mag"], unc=unc, level=rel[0], corr=corr, n_obs=n_obs,
                   method=method_text(cm),
                   loc="สุดท้าย (P+S+Lg)" if loc is hypo else "คลื่น P อย่างเดียว", pred={})
        for s in sites:
            p = preds[s["name"]]
            mmi = float(np.clip(p["mmi"] + corrs[s["name"]], 1.0, 10.0))
            lv, name = INT.alert_level(pick["mag"], mmi)
            if RM.is_teleseism(hypo):                           # แผ่นดินไหวไกลมาก: ตำแหน่ง/MMI ใช้ไม่ได้ ไม่แจ้งเตือน
                lv, name = 0, "ไม่แจ้ง (teleseism)"
            row["pred"][s["name"]] = (mmi, lv)
            info = dict(t=t, M=pick["mag"], type=pick["type"], mmi=mmi, level=lv, name=name, corr=corrs[s["name"]],
                        lead_s=s["t_s"] - t, lead_lg=s["t_lg"] - t)
            if lv > 0 and s["alert"] is None:
                s["alert"] = info
            if lv == 2 and s["upgrade"] is None and s["alert"] is not None and s["alert"]["level"] < 2:
                s["upgrade"] = info
        if key != last_key:                                     # เก็บเฉพาะการอัปเดตที่ผลเปลี่ยน
            updates.append(row)
            last_key = key
    # ค่าทำนายสุดท้าย (ขนาดหลัก + ตำแหน่งสุดท้าย + ปรับด้วยสถานีทั้งหมด)
    pick, (_, _, unc) = res["preferred"], res["reliability"]
    for s in sites:
        s["final_raw"] = INT.predict_mmi(pick["mag"], s["delta"], hypo["depth"], s["vs30"], EQ_TYPE, unc,
                                         s.get("rjb"))
    flag_gross_outliers(sites)
    cm = correction_model(sites, {s["name"]: s["final_raw"] for s in sites}, np.inf)
    res["correction"] = cm
    for s in sites:
        f = dict(s["final_raw"])
        c = site_corr(cm, s, f["R"])
        if cm and cm.get("rmse") is not None and cm["method"] != "none":
            sig = float(max(cm["rmse"], 0.3))          # ความไม่แน่นอน = RMSE แบบ leave-one-out (วัดจากข้อมูลจริง)
            f.update(mmi=float(np.clip(f["mmi"] + c, 1, 10)), sigma=sig)
            f.update(mmi_low=max(f["mmi"] - sig, 1.0), mmi_high=min(f["mmi"] + sig, 10.0))
        f["corr"], f["n_obs"] = c, (cm or {}).get("n", 0)
        s["final"] = f
    return updates


# ======================================================================
#  ปรับค่าทำนายด้วยความรุนแรงที่วัดได้ (เลือกวิธีด้วย leave-one-out)
# ======================================================================
def _robust_trend(r, logR):
    """r = a + b log10 R   (Huber IRLS กันค่าผิดปกติ)"""
    A = np.c_[np.ones(len(r)), logR]
    w = np.ones(len(r))
    c = np.array([np.median(r), 0.0])
    for _ in range(6):
        c = np.linalg.lstsq(A * w[:, None], r * w, rcond=None)[0]
        e = r - A @ c
        sc = 1.4826 * np.median(np.abs(e)) + 1e-6
        w = np.sqrt(np.minimum(1.0, 1.345 * sc / np.maximum(np.abs(e), 1e-9)))
    return c


def _fit(method, r, lat, lon, logR, typ, L=None, nug=None):
    """สร้างแบบจำลองค่าปรับจาก residual ของสถานี (r = MMI วัด − MMI ทำนาย)"""
    m = {"method": method, "L": L, "nug": nug}
    if method == "none":
        m["base"] = ("const", 0.0)
        return m
    if method.startswith("trend") and len(r) >= 4:
        m["base"] = ("trend",) + tuple(_robust_trend(r, logR))
    elif method == "type":
        g = float(np.median(r))
        by = {t: float(np.median(r[typ == t])) for t in set(typ) if np.sum(typ == t) >= 3}
        m["base"] = ("type", by, g)
    else:
        m["base"] = ("const", float(np.median(r)))
    if method.endswith("krig") and L:
        e = r - _base(m, logR, typ)
        s2 = max(float(np.var(e)), 1e-3)
        D = haversine_km(lat[:, None], lon[:, None], lat[None], lon[None])
        Cm = s2 * ((1 - nug) * np.exp(-D / L) + nug * np.eye(len(r)))
        m["krig"] = dict(lat=lat, lon=lon, alpha=np.linalg.solve(Cm, e), s2=s2)
    return m


def _base(m, logR, typ):
    b = m["base"]
    if b[0] == "trend":
        return b[1] + b[2] * np.asarray(logR)
    if b[0] == "type":
        return np.array([b[1].get(t, b[2]) for t in np.atleast_1d(typ)], float)
    return np.full(np.shape(logR), b[1], float)


def model_corr(m, lat, lon, R, typ=None):
    """ค่าปรับ (MMI) ที่ตำแหน่ง (lat, lon) ระยะ R (km)  รับได้ทั้งค่าเดียวและ array"""
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    logR = np.log10(np.maximum(np.asarray(R, float), 1.0))
    typ = np.full(np.shape(logR), "target", object) if typ is None else typ
    c = _base(m, logR, typ) if m else np.zeros(np.shape(logR))
    k = (m or {}).get("krig")
    if k is not None:
        D = haversine_km(lat[..., None], lon[..., None], k["lat"], k["lon"])
        c = c + (k["s2"] * (1 - m["nug"]) * np.exp(-D / m["L"])) @ k["alpha"]
    return np.clip(c, -3.0, 3.0)


MIN_N_TREND, MIN_N_KRIG = 6, 8   # จำนวนสถานีขั้นต่ำของวิธีที่มีพารามิเตอร์มากขึ้น (กัน overfit เมื่อสถานีน้อย)


def _candidates(n):
    out = [("none", None, None), ("median", None, None), ("type", None, None)]
    if n >= MIN_N_TREND:
        out.append(("trend", None, None))
    if n >= MIN_N_KRIG:
        for L in KRIG_LENGTHS_KM:
            for nug in KRIG_NUGGET:
                out += [("median+krig", L, nug), ("trend+krig", L, nug)]
    return out


def flag_gross_outliers(sites, limit=None):
    """
    สถานีที่ MMI วัด − MMI ทำนาย ต่างจากค่ากลางของทุกสถานีเกิน GROSS_OUTLIER_MMI = ข้อมูลเสีย (เช่น สัญญาณกระโดด
    ที่ despike ตัดไม่หมด, baseline เลื่อน) ไม่ใช้เป็นค่าที่วัดได้ ไม่แสดงบนแผนที่ และไม่ใช้ปรับค่าทำนาย
      ทดสอบ: M4.6 21 ส.ค. 2568 SURT วัดได้ PGV 147 cm/s (MMI IX) ที่ระยะ 831 km ซึ่งเป็นไปไม่ได้
      ค่าจริงที่ต่างมากสุดใน 3 เหตุการณ์ทดสอบ ≈ +3 MMI (สถานีบนดินอ่อน) จึงตั้งไว้ 4 MMI
    """
    limit = GROSS_OUTLIER_MMI if limit is None else limit
    pts = [s for s in sites if obs_mmi(s) is not None and s.get("final_raw")]
    if len(pts) < 5:
        return []
    r = np.array([obs_mmi(s) - s["final_raw"]["mmi"] for s in pts])
    med = float(np.median(r))
    bad = []
    for s, x in zip(pts, r):
        if abs(x - med) > limit:
            s["obs"]["usable"] = False
            s["obs"].setdefault("flags", []).append(
                f"ค่าที่วัดได้ต่างจากสถานีอื่นเกิน {limit:.0f} MMI ({x - med:+.1f}) น่าจะเป็นข้อมูลเสีย")
            bad.append(s["name"])
    return bad


def correction_model(sites, preds, t):
    """
    แบบจำลองค่าปรับจากสถานีที่ (1) ผ่าน QC และ (2) สั่นแรงสุดผ่านไปแล้ว ณ เวลา t
      1) r = MMI วัด − MMI ทำนาย ;  ตัดค่าผิดปกติ |r − median| > OUTLIER_MAD x 1.4826 MAD
      2) ลองทุกวิธี (CORRECTION_METHOD = "auto") แล้วเลือกวิธีที่ RMSE แบบ leave-one-out ต่ำสุด
           none        : ไม่ปรับ
           median      : + ค่ากลางของ r ทุกสถานี
           type        : + ค่ากลางของ r แยกตามประเภทสถานี (ต้องมี >= 3 สถานีต่อประเภท)
           trend       : + (a + b log10 R)  แก้การลดทอนตามระยะที่สมการไม่ตรงกับพื้นที่
           *+krig      : + kriging ของ r ที่เหลือ (ShakeMap): สถานีใกล้ตำแหน่งนั้นมีน้ำหนักมาก
                         covariance = σ² [ (1-nugget) exp(-d/L) + nugget δ ]
    return dict แบบจำลอง (method, rmse, n, excluded, table ...) หรือ None ถ้าสถานีไม่พอ
    """
    if not USE_OBS_CORRECTION:
        return None
    pts = [s for s in sites if s["obs"] and s["obs"].get("usable", True) and s["name"] in preds
           and s["obs"]["t_peak"] + RTT.PROCESSING_DELAY_S <= t]
    if len(pts) < OBS_MIN_STATIONS:
        return {"method": "none", "base": ("const", 0.0), "n": len(pts), "rmse": None, "excluded": [],
                "table": []}
    r = np.array([s["obs"]["mmi"] - preds[s["name"]]["mmi"] for s in pts])
    med = np.median(r)
    mad = max(1.4826 * np.median(np.abs(r - med)), 0.25)
    ok = np.abs(r - med) <= max(OUTLIER_MAD * mad, OUTLIER_MIN_MMI)
    excluded = [(s["name"], float(x)) for s, x, o in zip(pts, r, ok) if not o]
    pts = [s for s, o in zip(pts, ok) if o]
    r = r[ok]
    lat = np.array([s["lat"] for s in pts]); lon = np.array([s["lon"] for s in pts])
    logR = np.log10(np.array([max(preds[s["name"]]["R"], 1.0) for s in pts]))
    typ = np.array([sensor_key(s) for s in pts], object)
    n = len(pts)
    cands = _candidates(n)
    if CORRECTION_METHOD != "auto":
        cands = [c for c in cands if c[0] == CORRECTION_METHOD] or [("median", None, None)]
    table = []
    for meth, L, nug in cands:
        if n >= 3:
            err = []
            for i in range(n):
                tr = np.arange(n) != i
                mi = _fit(meth, r[tr], lat[tr], lon[tr], logR[tr], typ[tr], L, nug)
                err.append(r[i] - float(model_corr(mi, lat[i:i + 1], lon[i:i + 1], 10 ** logR[i:i + 1],
                                                   typ[i:i + 1])[0]))
            rmse = float(np.sqrt(np.mean(np.square(err))))
        else:
            rmse = float(np.sqrt(np.mean(np.square(r - np.median(r))))) if meth != "none" else \
                float(np.sqrt(np.mean(r ** 2)))
        table.append((meth, L, nug, rmse))
    best = min(table, key=lambda x: x[3])
    m = _fit(best[0], r, lat, lon, logR, typ, best[1], best[2])
    if m["base"][0] == "const" and not m.get("krig"):
        m["base"] = ("const", float(np.clip(m["base"][1], -OBS_MAX_CORRECTION, OBS_MAX_CORRECTION)))
    m.update(rmse=best[3], n=n, excluded=excluded, table=sorted(table, key=lambda x: x[3]),
             rmse_none=next((x[3] for x in table if x[0] == "none"), None))
    return m


def method_text(m):
    if not m or m.get("method") in (None, "none"):
        return "ไม่ปรับ"
    name = {"median": "ค่ากลางทุกสถานี", "type": "ค่ากลางแยกประเภทสถานี", "trend": "แนวโน้มตามระยะ",
            "median+krig": "ค่ากลาง + kriging", "trend+krig": "แนวโน้มตามระยะ + kriging"}[m["method"]]
    if m.get("L"):
        name += f" (L {m['L']:.0f} km, nugget {m['nug']:.1f})"
    return name


def site_corr(m, s, R):
    if not m:
        return 0.0
    typ = np.array([sensor_key(s) if s["kind"] == "station" else "target"], object)
    return float(model_corr(m, [s["lat"]], [s["lon"]], [R], typ)[0])


def lead_text(lead_s, lead_lg):
    if lead_s > 0:
        return f"ก่อนคลื่น S {lead_s:.0f} s, ก่อนคลื่น Lg {lead_lg:.0f} s"
    if lead_lg > 0:
        return f"หลังคลื่น S มาถึง {-lead_s:.0f} s แต่ก่อนคลื่น Lg {lead_lg:.0f} s"
    return "คลื่นมาถึงก่อน (อยู่ใน Blind zone)"


# ======================================================================
#  แสดงผล
# ======================================================================
def _mmi_txt(m):
    k, rom, shake, dmg = INT.mmi_class(m)
    return f"MMI {rom} ({shake})"


def public_report(res, sites):
    hypo, pick = res["hypo"], res["preferred"]
    line = "=" * 70
    pos, src = RM.location_lines(hypo)
    print(line)
    print(f"  ประมาณความรุนแรงแผ่นดินไหว  (ข้อมูลเบื้องต้นจากระบบอัตโนมัติ)")
    print(line)
    rel = "ต่ำ (teleseism: ขนาดต่ำกว่าจริงมาก)" if RM.is_teleseism(hypo) else res["reliability"][0]
    print(f"  ขนาด (Mw) : {RM.mw_text(pick)}  ({M.size_class(pick['mag'])})   ความน่าเชื่อถือ {rel}")
    print(f"  ที่มา      : {RM.mw_source(pick)}")
    print(f"  เวลาเกิด   : {RM.thai_time(hypo['origin_time'])} (เวลาประเทศไทย)")
    print(f"  ตำแหน่ง    : {pos}")
    print(f"  ที่มาตำแหน่ง: {src}")
    if RM.is_teleseism(hypo):
        print("  ⚠ ระบบไม่ออกคำเตือน: คลื่นมาจากแผ่นดินไหวไกลมาก ค่าความรุนแรงด้านล่างเป็นค่าประมาณที่ใช้ไม่ได้")
    print(line)
    groups =[("seismic", "▲ สถานีตรวจวัดแผ่นดินไหว (Seismic station)"),
              ("strong_motion", "■ สถานีตรวจวัดอัตราเร่งของพื้นดิน (Strong-motion station)"),
              (None, "★ ตำแหน่งที่แจ้งเตือน")]
    order = []
    for key, title in groups:
        g = [s for s in sites if (s["kind"] == "station" and sensor_key(s) == key)
             or (key is None and s["kind"] != "station")]
        if g:
            order.append((title, sorted(g, key=lambda q: q["delta"])))
    for title, g in order:
        print(f"  ---------- {title} ----------")
        for s in g:
            f = s["final"]
            a = s["alert"]
            head = f"[{s['name']}]  ห่าง {s['delta']:,.0f} km"
            print(f"  {head}")
            print(f"     คาดว่าสั่น : {_mmi_txt(f['mmi'])}  (อาจอยู่ระหว่าง MMI {f['mmi_low']:.1f}–{f['mmi_high']:.1f})")
            if a is None:
                print("     แจ้งเตือน  : " + ("ไม่แจ้ง (แผ่นดินไหวไกลมาก)" if RM.is_teleseism(hypo)
                                              else "ไม่แจ้ง (ต่ำกว่าเกณฑ์ MMI III)"))
            else:
                when = f"ส่งเมื่อ +{a['t'] - hypo['origin_time']:.0f} s หลังเกิดเหตุ  {lead_text(a['lead_s'], a['lead_lg'])}"
                print(f"     แจ้งเตือน  : {a['name']}  ({when})")
                if s["upgrade"]:
                    u = s["upgrade"]
                    print(f"                 ยกระดับเป็นแจ้งเตือนสาธารณะเมื่อ +{u['t'] - hypo['origin_time']:.0f} s "
                          f"({lead_text(u['lead_s'], u['lead_lg'])})")
                print(f"     ควรทำ     : {INT.action_text(f['mmi'])}")
            if s["obs"]:
                o = s["obs"]
                print(f"     วัดได้จริง : {_mmi_txt(o['mmi'])}  JMA {o['jma_class']}"
                      + ("" if o.get("usable", True) else "   (ข้อมูลใช้ไม่ได้)"))
                for fl in o.get("flags", []):
                    print(f"     ⚠ ข้อมูล: {fl}")
            for fl in f["flags"]:
                print(f"     ⚠ {fl}")
    print(line)
    cm = res.get("correction")
    if cm and cm.get("method") not in (None, "none"):
        print(f"  * ปรับค่าทำนายด้วยความรุนแรงที่วัดได้จริงจาก {cm['n']} สถานี  วิธี: {method_text(cm)}")
        print(f"    ความคลาดเคลื่อนเมื่อทดสอบแบบ leave-one-out ±{cm['rmse']:.1f} MMI"
              + (f"  (ถ้าไม่ปรับ ±{cm['rmse_none']:.1f})" if cm.get("rmse_none") else ""))
        if cm["excluded"]:
            print("    ไม่ใช้ปรับ (ค่าผิดปกติ): " + ", ".join(f"{n} ({r:+.1f})" for n, r in cm["excluded"]))
    print("  * ค่าประมาณจากขนาดและระยะ (point source) ความรุนแรงจริงขึ้นกับทิศการแตกของรอยเลื่อนและชั้นดิน")
    print(line)


def detail_report(res, stages, updates, sites):
    T0 = res["hypo"]["origin_time"]
    print("\n================ ลำดับการอัปเดตของระบบ (แสดงเมื่อขนาดหรือตำแหน่งเปลี่ยน) ================")
    show = [s for s in sites if s["kind"] != "station"] + \
        sorted([s for s in sites if s["kind"] == "station"], key=lambda q: q["delta"])[:max(0, 8 - len(TARGET_SITES))]
    names = [s["name"] for s in show]                      # แสดงตำแหน่งแจ้งเตือน + สถานีใกล้สุด (ไม่เกิน 8 คอลัมน์)
    print(f"{'เวลา':>7} {'จากชนิด':>7} {'Mw':>4} {'±':>4}  {'ตำแหน่งจาก':<18}{'ปรับ (สถานี)':>13}"
          + "".join(f"{n[:10]:>12}" for n in names))
    for u in updates:
        cells = "".join(f"{u['pred'][n][0]:7.1f}{'*' * u['pred'][n][1]:<5}" for n in names)
        cr = f"{u['corr']:+.1f} ({u['n_obs']})"
        print(f"+{u['t'] - T0:5.0f} s {u['type']:>7} {u['M']:4.1f} {u['unc']:4.1f}  {u['loc']:<18}{cr:>13}{cells}")
    print("  (ค่าในตาราง = MMI ทำนาย,  * = แจ้งเตือนผ่านแอป,  ** = แจ้งเตือนสาธารณะ,"
          "  ปรับ = ค่าที่ปรับจากสถานีที่สั่นแรงสุดผ่านไปแล้ว (จำนวนสถานี))")

    cm = res.get("correction")
    if cm and cm.get("table"):
        print("\n================ เลือกวิธีปรับค่าทำนายด้วยค่าที่วัดได้ (leave-one-out) ================")
        print(f"ใช้ {cm['n']} สถานี: ซ่อนค่าที่วัดได้ของแต่ละสถานีทีละสถานี ทำนายจากสถานีที่เหลือ แล้ววัดความคลาด")
        print(f"{'วิธี':<44}{'RMSE (MMI)':>11}")
        seen = set()
        for meth, L, nug, rmse in cm["table"]:
            label = method_text({"method": meth, "L": L, "nug": nug})
            if meth.endswith("krig") and meth in seen:
                continue                                   # แสดง kriging เฉพาะค่า L ที่ดีที่สุดของแต่ละแบบ
            seen.add(meth)
            print(f"{label:<44}{rmse:>11.2f}" + ("   <- เลือก" if (meth, L, nug) == (cm["method"], cm["L"], cm["nug"])
                                                 else ""))
        if cm["excluded"]:
            print(f"ไม่ใช้ปรับ (residual ต่างจากค่ากลางเกิน {OUTLIER_MAD:g} MAD และ {OUTLIER_MIN_MMI:g} MMI): "
                  + ", ".join(f"{n} ({r:+.1f})" for n, r in cm["excluded"]))
    bad = [(s["name"], s["obs"]["flags"]) for s in sites if s["obs"] and s["obs"].get("flags")]
    if bad:
        print("\nQC ข้อมูลสถานี:")
        for n, fl in bad:
            print(f"  {n:<6} " + "; ".join(fl))
    ok = [s for s in sites if obs_mmi(s) is not None]
    if ok:
        r = np.array([obs_mmi(s) - s["final_raw"]["mmi"] for s in ok])
        print(f"\nความแม่นของค่าทำนายก่อนปรับ (สิ่งที่ใช้แจ้งเตือนก่อนคลื่นถึง) {len(ok)} สถานี:"
              f"  bias (วัด − ทำนาย) {r.mean():+.2f}  SD {r.std():.2f}  RMSE {np.sqrt(np.mean(r ** 2)):.2f} MMI")
        print("  ระยะที่ใช้: " + ("ถึงแนวรอยแตก FAULT_TRACE" if FAULT_TRACE else "ถึงจุดศูนย์กลาง (จุดกำเนิด)")
              + f"   น้ำหนัก GMPE: 1 ที่ <= {INT.GMPE_TAPER_KM[0]:.0f} km -> 0 ที่ {INT.GMPE_TAPER_KM[1]:.0f} km")
    diff = [(s["name"], s.get("sensor_file"), s.get("sensor")) for s in sites
            if s.get("sensor_file") and s.get("sensor") and s["sensor_file"] != s["sensor"]]
    if diff:
        print("\nประเภทสถานีใน stations.csv ต่างจากช่องสัญญาณในไฟล์ (ใช้ตาม STATION_TYPE_SOURCE = "
              f"\"{STATION_TYPE_SOURCE}\"):")
        print("  " + ", ".join(f"{n} (ไฟล์: {SENSOR_SHORT.get(a, a)}, csv: {SENSOR_SHORT.get(b, b)})"
                           for n, a, b in diff))

    print("\n================ ค่าทำนายสุดท้าย เทียบกับค่าที่วัดได้ ================")
    print(f"{'ตำแหน่ง':<12}{'ชนิด':<20}{'Δ km':>6}{'Vs30':>6}{'PGV ทำนาย':>11}{'MMI GMPE':>9}{'MMI IPE':>8}{'MMI':>6}"
          f"{'±':>5}{'JMA':>5}{'| PGA วัด':>10}{'PGV วัด':>9}{'MMI วัด':>8}{'JMA วัด':>8}{'ปรับแล้ว':>9}")
    print("  (MMI GMPE / IPE / MMI = ก่อนปรับ,  ปรับแล้ว = หลังปรับด้วยค่าที่สถานีวัดได้)")
    for s in sorted(sites, key=lambda q: (q["kind"] != "station", q["delta"])):
        f, o = s["final_raw"], s["obs"]
        obs = (f"{o['pga']:10.2f}{o['pgv']:9.3f}{o['mmi']:8.1f}{o['jma_class']:>8}" if o
               else f"{'-':>10}{'-':>9}{'-':>8}{'-':>8}")
        obs += f"{s['final']['mmi']:9.1f}" + ("  QC" if o and o.get("flags") else "")
        kind = sensor_name(s, short=True) if s["kind"] == "station" else "★ ตำแหน่งแจ้งเตือน"
        print(f"{s['name'][:11]:<12}{kind:<20}{s['delta']:6.0f}{s['vs30']:6.0f}{f['pgv_surface']:11.3f}{f['mmi_gmpe']:9.1f}"
              f"{f['mmi_ipe']:8.1f}{f['mmi']:6.1f}{f['sigma']:5.1f}{f['jma_class']:>5}{obs}")
    res_ = [s["obs"]["mmi"] - s["final_raw"]["mmi"] for s in sites if s["obs"]]
    if res_:
        print(f"MMI วัด − ทำนาย ที่สถานี: ค่ากลาง {np.median(res_):+.1f}  (ช่วง {min(res_):+.1f} ถึง {max(res_):+.1f})")
    print("PGA = แกนราบที่ใหญ่กว่า (gal), PGV = แกนราบที่ใหญ่กว่า (cm/s), MMI วัดจาก PGV ตาม Worden et al. (2012)")
    print("▲ = สถานีตรวจวัดแผ่นดินไหว (มี seismometer ช่อง HH? / EH? ส่วนใหญ่มี accelerometer HN? ด้วย)")
    print("■ = สถานีตรวจวัดอัตราเร่งของพื้นดิน (มีแค่ accelerometer ช่อง HN?)")
    for key in ("seismic", "strong_motion"):
        r2 = [s["obs"]["mmi"] - s["final_raw"]["mmi"] for s in sites
              if s["obs"] and s["kind"] == "station" and sensor_key(s) == key]
        if r2:
            print(f"  {SENSOR_SHORT[key]}: {len(r2)} สถานี, MMI วัด − ทำนาย ค่ากลาง {np.median(r2):+.1f}")


# ======================================================================
#  กราฟ
# ======================================================================
def _plt():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    font = RM._thai_font()
    if font:
        plt.rcParams["font.family"] = font
    plt.rcParams["axes.unicode_minus"] = False
    return plt


MMI_COLORS = {1: "#ffffff", 2: "#bfccff", 3: "#a0e6ff", 4: "#80ffff", 5: "#7aff93", 6: "#ffff00",
              7: "#ffc800", 8: "#ff9100", 9: "#ff0000", 10: "#c80000"}      # สีมาตรฐานของ USGS ShakeMap


def plot_alert(res, sites, out="intensity_alert.png"):
    plt = _plt()
    if plt is None:
        return None
    T0 = res["hypo"]["origin_time"]
    sites = sorted(sites, key=lambda s: s["delta"])
    fig, ax = plt.subplots(figsize=(12, 1.0 + 0.75 * len(sites)))
    for y, s in enumerate(sites):
        k = INT.mmi_class(s["final"]["mmi"])[0]
        ax.barh(y, s["t_lg"] - T0, height=0.55, color=MMI_COLORS[k], edgecolor="0.4")
        ax.plot(s["t_p"] - T0, y, "|", color="#1f77b4", ms=18, mew=2)
        ax.plot(s["t_s"] - T0, y, "|", color="#d62728", ms=18, mew=2)
        ax.plot(s["t_lg"] - T0, y, "|", color="#ff7f0e", ms=18, mew=2)
        a = s["alert"]
        txt = f"{_mmi_txt(s['final']['mmi'])}"
        if a:
            ax.plot(a["t"] - T0, y, "v", color="k", ms=11)
            txt += f"   เตือน +{a['t'] - T0:.0f} s  (ก่อน S {a['lead_s']:+.0f} s, ก่อน Lg {a['lead_lg']:+.0f} s)"
        else:
            txt += "   ไม่แจ้งเตือน"
        x_txt = max(s["t_lg"], a["t"] if a else 0.0) - T0 + 8
        ax.text(x_txt, y, txt, va="center", fontsize=9)
    ax.set_yticks(range(len(sites)))
    ax.set_yticklabels([f"{s['name']}  ({s['delta']:.0f} km)" for s in sites])
    xmax = max(s["t_lg"] - T0 for s in sites)
    ax.set_xlim(0, xmax * 1.9)
    ax.set_xlabel("เวลาหลังเกิดแผ่นดินไหว (s)")
    ax.plot([], [], "v", color="k", label="ส่งคำเตือน")
    ax.plot([], [], "|", color="#1f77b4", mew=2, ms=12, label="คลื่น P มาถึง")
    ax.plot([], [], "|", color="#d62728", mew=2, ms=12, label="คลื่น S มาถึง")
    ax.plot([], [], "|", color="#ff7f0e", mew=2, ms=12, label="คลื่น Lg มาถึง")
    ax.legend(loc="lower right", fontsize=9)
    p = res["preferred"]
    ax.set_title(f"เวลาแจ้งเตือนเทียบกับเวลาคลื่นมาถึง   Mw {p['mag']:.1f} (จาก {p['type']})   สีแถบ = MMI ที่คาด")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_distance(res, sites, out="intensity_distance.png"):
    plt = _plt()
    if plt is None:
        return None
    hypo, p = res["hypo"], res["preferred"]
    unc = res["reliability"][2]
    dmax = max(s["delta"] for s in sites) * 1.3
    d = np.linspace(5, dmax, 200)
    fig, ax = plt.subplots(figsize=(10, 6))
    for vs, ls in ((C.DEFAULT_VS30, "-"), (200.0, "--")):
        pr = [INT.predict_mmi(p["mag"], x, hypo["depth"], vs, EQ_TYPE, unc) for x in d]
        m = np.array([q["mmi"] for q in pr])
        ax.plot(d, m, "k" + ls, label=f"MMI ทำนาย Vs30 {vs:.0f} m/s (ก่อนปรับ)")
        if vs == C.DEFAULT_VS30:
            ax.fill_between(d, [q["mmi_low"] for q in pr], [q["mmi_high"] for q in pr], color="0.85",
                            label="± 1 ส่วนเบี่ยงเบนมาตรฐาน (รวมความไม่แน่นอนของขนาด)")
            ax.plot(d, [q["mmi_gmpe"] for q in pr], color="#1f77b4", lw=1, label="วิธี GMPE -> GMICE (Worden 2012)")
            ax.plot(d, [q["mmi_ipe"] for q in pr], color="#2ca02c", lw=1, label="วิธี IPE (Allen 2012)")
    for s in sites:
        if s["obs"]:
            ax.plot(s["delta"], s["obs"]["mmi"], "o", color="#d62728", mec="k", ms=8)
            ax.annotate(s["name"], (s["delta"], s["obs"]["mmi"]), textcoords="offset points", xytext=(5, 5),
                        fontsize=8)
        else:
            ax.plot(s["delta"], s["final"]["mmi"], "*", color="#ff7f0e", mec="k", ms=13)   # หลังปรับด้วยสถานี
            ax.annotate(s["name"], (s["delta"], s["final"]["mmi"]), textcoords="offset points",
                        xytext=(5, 5), fontsize=8)
    ax.plot([], [], "o", color="#d62728", mec="k", label="MMI วัดได้ที่สถานี")
    ax.plot([], [], "*", color="#ff7f0e", mec="k", ms=12, label="ตำแหน่งที่แจ้งเตือน (ทำนาย)")
    for k in range(1, 10):
        ax.axhspan(k - 0.5, k + 0.5, color=MMI_COLORS[k], alpha=0.25, lw=0)
    ax.axvline(INT.SM99_MAX_X, color="0.5", ls=":", lw=1)
    ax.text(INT.SM99_MAX_X, 9.3, "ช่วงที่สมการปรับเทียบไว้ | ", ha="right", fontsize=8, color="0.4")
    ax.set_yticks(range(1, 11))
    ax.set_yticklabels([INT.MMI_TABLE[k][0] for k in range(1, 11)])
    ax.set_ylim(0.5, 9.8)
    ax.set_xlim(0, dmax)
    ax.set_xlabel("ระยะจากจุดเหนือศูนย์กลาง (km)")
    ax.set_ylabel("MMI")
    ax.set_title(f"ความรุนแรงกับระยะทาง   Mw {p['mag']:.1f} (จาก {p['type']})")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


MAP_FILE = "intensity_map.html"


def _mmi_rings(res, sites, vs30=None):
    """รัศมี (km) ที่ MMI ทำนาย (ค่าสุดท้าย รวมค่าปรับจากสถานี) ลดลงถึงขอบล่างของแต่ละระดับ"""
    hypo, p, (_, _, unc) = res["hypo"], res["preferred"], res["reliability"]
    vs30 = vs30 or C.DEFAULT_VS30
    cm = res.get("correction")
    d = np.linspace(1, 3000, 1500)
    m = np.array([INT.predict_mmi(p["mag"], x, hypo["depth"], vs30, EQ_TYPE, unc)["mmi"] for x in d])
    base = _base(cm, np.log10(np.hypot(d, hypo["depth"])), np.full(len(d), "target", object)) if cm else 0.0
    m = m + base
    corr = float(np.median(base)) if cm else 0.0
    rings = []
    for k in range(2, 11):
        below = np.where(m < k - 0.5)[0]
        if len(below) and below[0] > 0:
            rings.append((k, float(d[below[0]])))
    return rings, vs30, corr


def make_map(res, sites, out=MAP_FILE):
    """
    แผนที่ความรุนแรงหน่วย MMI
      1) map_plot.make_mmi_map (folium: ภาพดาวเทียม + ไฟล์ .kml สำหรับ Google Earth)   ถ้ามี folium
      2) ถ้าไม่มี folium ใช้แผนที่ Leaflet ในตัว (_make_map_leaflet)
    """
    try:
        from map_plot import make_mmi_map
    except (ImportError, AttributeError):
        make_mmi_map = None
    if make_mmi_map is not None:
        h, p = res["hypo"], res["preferred"]
        rings, vs_ring, corr = _mmi_rings(res, sites)
        pts = [dict(name=s["name"], lat=s["lat"], lon=s["lon"], kind=s["kind"], mmi=s["final"]["mmi"],
                    obs=obs_mmi(s), alert=(s["alert"]["level"] if s["alert"] else 0),
                    used=s.get("used", False), delta=s["delta"], sensor=sensor_key(s),
                    lead=(s["alert"]["lead_s"] if s["alert"] else None), popup=_popup_html(res, s)) for s in sites]
        epi = dict(lat=h["lat"], lon=h["lon"], popup=_epi_html(res))
        note = (f"วงกลม = ขอบเขตที่คาดว่าสั่นถึงแต่ละระดับ (Vs30 {vs_ring:.0f} m/s"
                + (f", ปรับจากสถานี {corr:+.1f} MMI" if corr else "") + ")")
        html = make_mmi_map(pts, out, epi, rings, title=f"ความรุนแรงแผ่นดินไหว (MMI)  Mw {p['mag']:.1f}", note=note)
        if html:
            return html
    return _make_map_leaflet(res, sites, out)


def _epi_html(res):
    h, p = res["hypo"], res["preferred"]
    return (f"<b>จุดศูนย์กลาง</b><br>Mw {p['mag']:.1f} (จาก {p['type']}) {M.size_class(p['mag'])}<br>"
            f"{RM.thai_time(h['origin_time'])}<br>{RM.location_lines(h)[0]}<br>{RM.location_lines(h)[1]}")


def _popup_html(res, s):
    T0 = res["hypo"]["origin_time"]
    f, o, a = s["final"], s["obs"], s["alert"]
    k, rom, shake, dmg = INT.mmi_class(f["mmi"])
    rows = [f"ห่างจากจุดศูนย์กลาง {s['delta']:,.0f} km", f"Vs30 {s['vs30']:.0f} m/s"]
    if s["kind"] == "station":
        rows.insert(0, f"<b>{sensor_name(s)}</b>" + (f" &nbsp;ช่อง {s['channels']}" if s.get("channels") else ""))
        rows.append("ใช้คำนวณตำแหน่ง/ขนาด" if s.get("used") else "ไม่ได้ใช้คำนวณขนาด (แสดงค่าที่วัดได้)")
    rows += [f"<b>คาดว่าสั่น: MMI {rom} ({shake})</b> &nbsp;ค่า {f['mmi']:.1f} "
            f"(ช่วง {f['mmi_low']:.1f}–{f['mmi_high']:.1f})", f"ความเสียหายที่คาด: {dmg}"]
    if o:
        _, ro, so, _ = INT.mmi_class(o["mmi"])
        rows += [f"<b>วัดได้จริง: MMI {ro} ({so})</b> &nbsp;ค่า {o['mmi']:.1f}"
                 + ("" if o.get("usable", True) else " (ข้อมูลใช้ไม่ได้)"),
                 f"PGA {o['pga']:.2f} gal, PGV {o['pgv']:.3f} cm/s"]
        if o.get("flags"):
            rows.append("⚠ QC: " + "; ".join(o["flags"]))
    if a:
        rows += [f"<b>{a['name']}</b> ส่งเมื่อ +{a['t'] - T0:.0f} s หลังเกิดเหตุ", lead_text(a["lead_s"], a["lead_lg"])]
    else:
        rows.append("ไม่แจ้งเตือน (ต่ำกว่า MMI III)")
    rows.append(f"คลื่น S มาถึง +{s['t_s'] - T0:.0f} s, คลื่น Lg +{s['t_lg'] - T0:.0f} s")
    return "<br>".join(rows)


def _make_map_leaflet(res, sites, out=MAP_FILE):
    """
    แผนที่ความรุนแรง หน่วย MMI (สีมาตรฐาน USGS ShakeMap) แบบไม่ต้องติดตั้ง folium
      วงกลม  = ขอบเขตที่คาดว่าสั่นถึงแต่ละระดับ MMI (Vs30 ค่าตั้งต้น)
      จุดกลม = สถานี (วงนอก = MMI ทำนาย, จุดใน = MMI วัดได้)   ดาว = ตำแหน่งที่แจ้งเตือน
    """
    import json
    h, p = res["hypo"], res["preferred"]
    T0 = h["origin_time"]
    rings, vs_ring, corr = _mmi_rings(res, sites)
    pts = []
    for s in sites:
        f, o, a = s["final"], s["obs"], s["alert"]
        k, rom, shake, dmg = INT.mmi_class(f["mmi"])
        rows = [f"ห่างจากจุดศูนย์กลาง {s['delta']:,.0f} km", f"Vs30 {s['vs30']:.0f} m/s",
                f"<b>คาดว่าสั่น: MMI {rom} ({shake})</b> &nbsp;ค่า {f['mmi']:.1f} "
                f"(ช่วง {f['mmi_low']:.1f}–{f['mmi_high']:.1f})", f"ความเสียหายที่คาด: {dmg}"]
        if o:
            ko, ro, so, _ = INT.mmi_class(o["mmi"])
            rows += [f"<b>วัดได้จริง: MMI {ro} ({so})</b> &nbsp;ค่า {o['mmi']:.1f}",
                     f"PGA {o['pga']:.2f} gal, PGV {o['pgv']:.3f} cm/s"]
        if a:
            rows += [f"<b>{a['name']}</b> ส่งเมื่อ +{a['t'] - T0:.0f} s หลังเกิดเหตุ", lead_text(a["lead_s"], a["lead_lg"])]
        else:
            rows.append("ไม่แจ้งเตือน (ต่ำกว่า MMI III)")
        rows.append(f"คลื่น S มาถึง +{s['t_s'] - T0:.0f} s, คลื่น Lg +{s['t_lg'] - T0:.0f} s")
        pts.append(dict(name=s["name"], lat=s["lat"], lon=s["lon"], kind=s["kind"],
                        mmi=round(f["mmi"], 1), k=k, rom=rom,
                        obs=(round(obs_mmi(s), 1) if obs_mmi(s) is not None else None),
                        ko=(INT.mmi_class(obs_mmi(s))[0] if obs_mmi(s) is not None else None),
                        alert=(a["level"] if a else 0), html=_popup_html(res, s),
                        delta=round(s["delta"]), used=s.get("used", False),
                        sensor=sensor_key(s),
                        lead=(round(a["lead_s"]) if a else None),
                        romo=(INT.mmi_class(obs_mmi(s))[1] if obs_mmi(s) is not None else None)))
    data = dict(
        epi=dict(lat=h["lat"], lon=h["lon"],
                 html=f"<b>จุดศูนย์กลาง</b><br>Mw {p['mag']:.1f} (จาก {p['type']}) {M.size_class(p['mag'])}<br>"
                      f"{RM.thai_time(T0)}<br>{RM.location_lines(h)[1]}"),
        pts=pts, rings=[dict(k=k, r=r, rom=INT.MMI_TABLE[k][0]) for k, r in rings],
        colors={str(k): c for k, c in MMI_COLORS.items()},
        legend=[dict(k=k, rom=INT.MMI_TABLE[k][0], shake=INT.MMI_TABLE[k][1], dmg=INT.MMI_TABLE[k][2])
                for k in range(1, 11) if k != 2],
        title=f"ความรุนแรงแผ่นดินไหว (MMI)  Mw {p['mag']:.1f}",
        note=(f"วงกลม = ขอบเขตที่คาดว่าสั่นถึงแต่ละระดับ (Vs30 {vs_ring:.0f} m/s"
              + (f", ปรับจากสถานี {corr:+.1f} MMI" if corr else "") + ")"))
    html = _MAP_TEMPLATE.replace("__DATA__", json.dumps(data, ensure_ascii=False))
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html)
    return out


_MAP_TEMPLATE = """<!doctype html>
<html lang="th"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Intensity Map (MMI)</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
 html,body{margin:0;height:100%;font-family:Tahoma,"Leelawadee UI","Noto Sans Thai",sans-serif}
 #map{height:100%}
 .box{background:#fff;padding:8px 10px;border-radius:6px;box-shadow:0 1px 5px rgba(0,0,0,.35);font-size:12px;line-height:1.35}
 .box h4{margin:0 0 6px;font-size:14px}
 .sw{display:inline-block;width:18px;height:12px;border:1px solid #777;margin-right:6px;vertical-align:middle}
 .lbl{background:rgba(255,255,255,.85);border:0;box-shadow:none;font-weight:bold;font-size:11px;color:#222;line-height:1.25;padding:1px 4px}
</style></head><body><div id="map"></div>
<script>
const D = __DATA__;
const map = L.map('map');
const ESRI = s => 'https://server.arcgisonline.com/ArcGIS/rest/services/' + s + '/MapServer/tile/{z}/{y}/{x}';
const bases = {
  'แผนที่สีเทาอ่อน (Esri)': L.tileLayer(ESRI('Canvas/World_Light_Gray_Base'), {maxZoom: 16, attribution: 'Esri, HERE, Garmin, &copy; OpenStreetMap contributors'}),
  'ภาพดาวเทียม': L.tileLayer(ESRI('World_Imagery'), {maxZoom: 19, attribution: 'Esri, Maxar, Earthstar Geographics'}),
  'แผนที่ถนน (Esri)': L.tileLayer(ESRI('World_Street_Map'), {maxZoom: 19, attribution: 'Esri'})};
bases['แผนที่สีเทาอ่อน (Esri)'].addTo(map);
// เส้นแบ่งเขตประเทศ/จังหวัด + ชื่อ อยู่บนสุด (ไม่ถูกพื้นสีบัง) และสลับชุดตามพื้นหลัง
map.createPane('labels'); map.getPane('labels').style.zIndex = 450; map.getPane('labels').style.pointerEvents = 'none';
map.createPane('labels_dark'); map.getPane('labels_dark').style.zIndex = 451; map.getPane('labels_dark').style.pointerEvents = 'none';
map.getPane('labels_dark').style.filter = 'contrast(2.2) brightness(0.55)';
const refSat = L.tileLayer(ESRI('Reference/World_Boundaries_and_Places'), {maxZoom: 16, pane: 'labels', attribution: 'Esri'});
const refMap = L.tileLayer(ESRI('Canvas/World_Light_Gray_Reference'), {maxZoom: 16, pane: 'labels_dark', attribution: 'Esri'}).addTo(map);
L.control.layers(bases, {'เส้นแบ่งเขตประเทศ/จังหวัด (ภาพดาวเทียม)': refSat, 'เส้นแบ่งเขตประเทศ/จังหวัด (แผนที่)': refMap},
  {collapsed: false}).addTo(map);
map.on('baselayerchange', e => { const sat = e.name === 'ภาพดาวเทียม';
  const on = sat ? refSat : refMap, off = sat ? refMap : refSat;
  if (map.hasLayer(off)) map.removeLayer(off); if (!map.hasLayer(on)) map.addLayer(on); });
const col = k => D.colors[String(k)] || '#fff';
const SHAPE = (sensor, c, z) => sensor === 'seismic'
  ? '<svg width="' + z + '" height="' + z + '" style="display:block"><polygon points="' + z/2 + ',1.5 ' + (z-1.5) + ',' + (z-1.5) + ' 1.5,' + (z-1.5) + '" fill="' + c + '" stroke="#222" stroke-width="1.6"/></svg>'
  : '<svg width="' + z + '" height="' + z + '" style="display:block"><rect x="2" y="2" width="' + (z-4) + '" height="' + (z-4) + '" fill="' + c + '" stroke="#222" stroke-width="1.6"/></svg>';
const SNAME = {seismic: 'สถานีตรวจวัดแผ่นดินไหว (Seismic station)', strong_motion: 'สถานีตรวจวัดอัตราเร่งของพื้นดิน (Strong-motion station)'};
const SSYM = {seismic: '&#9650;', strong_motion: '&#9632;'};
const SLEG = () => ['seismic', 'strong_motion'].filter(k => D.pts.some(p => (p.kind === undefined || p.kind === 'station') && p.sensor === k))
  .map(k => '<div style="display:flex;align-items:center;gap:6px;margin-top:2px">' + SHAPE(k, '#fff', 14) + SNAME[k] + '</div>').join('');

const bounds = [[D.epi.lat, D.epi.lon]];
// วงกลมขอบเขต MMI (ใหญ่ไปเล็ก ให้วงเล็กอยู่บน)
D.rings.slice().sort((a,b)=>b.r-a.r).forEach(g => {
  L.circle([D.epi.lat, D.epi.lon], {radius: g.r*1000, color:'#555', weight:1, dashArray:'4 4',
    fillColor: col(g.k), fillOpacity: 0.22}).addTo(map)
   .bindTooltip('MMI ' + g.rom + ' ถึงระยะ ~' + Math.round(g.r) + ' km', {sticky:true});
});
// จุดศูนย์กลาง
L.marker([D.epi.lat, D.epi.lon], {icon: L.divIcon({className:'', html:
  '<div style="font-size:26px;color:#c00;line-height:26px;text-shadow:0 0 3px #fff">&#10033;</div>',
  iconSize:[26,26], iconAnchor:[13,13]})}).addTo(map).bindPopup(D.epi.html);
// สถานีและตำแหน่งที่แจ้งเตือน
D.pts.forEach(p => {
  bounds.push([p.lat, p.lon]);
  let m;
  if (p.kind === 'station') {
    m = L.circleMarker([p.lat, p.lon], {radius: 13, color: p.alert ? '#000' : '#888', weight: p.alert ? 3 : 1,
      fillColor: col(p.k), fillOpacity: 0.95}).addTo(map);
    L.marker([p.lat, p.lon], {icon: L.divIcon({className: '', html: SHAPE(p.sensor, p.obs !== null ? col(p.ko) : '#bbb', 14),
      iconSize: [14, 14], iconAnchor: [7, 8]})}).addTo(map).bindPopup('<b>' + p.name + '</b><br>' + p.html);
  } else {
    m = L.marker([p.lat, p.lon], {icon: L.divIcon({className:'', html:
      '<div style="font-size:30px;line-height:30px;color:' + col(p.k) + ';-webkit-text-stroke:1.5px #000">&#9733;</div>',
      iconSize:[30,30], iconAnchor:[15,15]})}).addTo(map);
  }
  m.bindPopup('<b>' + p.name + '</b><br>' + p.html);
  const lab = p.obs!==null ? 'MMI วัดได้ ' + p.romo + ' (' + p.obs.toFixed(1) + ') · ทำนาย ' + p.rom + ' (' + p.mmi.toFixed(1) + ')'
                           : 'MMI ทำนาย ' + p.rom + ' (' + p.mmi.toFixed(1) + ')';
  L.tooltip({permanent:true, direction:'right', offset:[15,0], className:'lbl'})
   .setLatLng([p.lat, p.lon]).setContent('<span style="border-left:4px solid ' + col(p.obs!==null ? p.ko : p.k) +
     ';padding-left:4px">' + (p.kind === 'station' ? SSYM[p.sensor] + ' ' : '') + p.name + '<br>' + lab + '</span>').addTo(map);
});
// ตาราง MMI ของทุกสถานี
const tab = L.control({position:'bottomleft'});
tab.onAdd = () => { const d = L.DomUtil.create('div','box'); d.style.maxHeight='45vh'; d.style.overflow='auto';
  const cell = (k, rom, v) => '<td style="background:' + col(k) + ';text-align:center;padding:1px 6px"><b>' + rom + '</b> ' + v.toFixed(1) + '</td>';
  const rows = D.pts.slice().sort((a,b)=> (a.kind!=='station') - (b.kind!=='station') || a.delta - b.delta).map(p =>
    '<tr><td style="padding:1px 6px">' + (p.kind!=='station' ? '&#9733; ' : SSYM[p.sensor] + ' ') + p.name + (p.kind==='station' && !p.used ? ' *' : '') + '</td>' +
    '<td style="text-align:right;padding:1px 6px">' + p.delta.toLocaleString() + '</td>' +
    (p.obs!==null ? cell(p.ko, p.romo, p.obs) : '<td style="text-align:center">-</td>') + cell(p.k, p.rom, p.mmi) +
    '<td style="padding:1px 6px">' + (p.alert===2 ? 'สาธารณะ' : p.alert===1 ? 'แอป' : '-') +
    (p.alert && p.lead!==null ? ' (' + (p.lead>0?'+':'') + p.lead + ' s)' : '') + '</td></tr>').join('');
  const star = D.pts.some(p => p.kind==='station' && !p.used) ? '* input เข้าไปแต่ไม่ได้ใช้คำนวณขนาด, ' : '';
  d.innerHTML = '<h4>MMI ของแต่ละสถานี</h4><table style="border-collapse:collapse"><tr style="border-bottom:1px solid #999">' +
    '<th style="text-align:left;padding:1px 6px">สถานี</th><th style="padding:1px 6px">Δ km</th><th style="padding:1px 6px">MMI วัดได้</th>' +
    '<th style="padding:1px 6px">MMI ทำนาย</th><th style="padding:1px 6px">แจ้งเตือน</th></tr>' + rows + '</table>' +
    '<div style="font-size:11px;color:#555">' + star + 'เวลาในวงเล็บ = เตือนก่อนคลื่น S<br>' +
    '&#9650; = สถานีตรวจวัดแผ่นดินไหว, &#9632; = สถานีตรวจวัดอัตราเร่งของพื้นดิน, &#9733; = ตำแหน่งที่แจ้งเตือน</div>';
  return d; };
tab.addTo(map);
map.fitBounds(bounds, {padding:[40,40]});
// คำอธิบาย
const info = L.control({position:'topright'});
info.onAdd = () => { const d = L.DomUtil.create('div','box');
  d.innerHTML = '<h4>' + D.title + '</h4>' + D.note +
   '<br>วงนอกของสถานี = MMI ทำนาย, สัญลักษณ์ตรงกลาง = MMI วัดได้, ขอบดำ = ได้รับคำเตือน<br>ดาว = ตำแหน่งที่แจ้งเตือน (ทำนาย)' + SLEG();
  return d; };
info.addTo(map);
const leg = L.control({position:'bottomright'});
leg.onAdd = () => { const d = L.DomUtil.create('div','box');
  d.innerHTML = '<h4>มาตรา MMI (USGS)</h4>' + D.legend.map(r =>
    '<div><span class="sw" style="background:' + col(r.k) + '"></span><b>' +
    (r.k===3 ? 'II–III' : r.rom) + '</b> ' + r.shake + ' · เสียหาย: ' + r.dmg + '</div>').join('');
  return d; };
leg.addTo(map);
</script></body></html>
"""


STATION_MAP_FILE = "intensity_station_map.html"


def _mmi_rgba(mmi, alpha=0.75):
    """MMI (ทศนิยม) -> สี RGBA ไล่ระดับต่อเนื่องตามสีมาตรฐาน USGS ShakeMap"""
    ks = np.arange(1, 11)
    rgb = np.array([[int(MMI_COLORS[k][i:i + 2], 16) for i in (1, 3, 5)] for k in ks], float)
    x = np.clip(mmi, 1.0, 10.0)
    out = np.empty(x.shape + (4,), np.uint8)
    for c in range(3):
        out[..., c] = np.interp(x, ks, rgb[:, c]).astype(np.uint8)
    out[..., 3] = np.clip(np.asarray(alpha) * 255, 0, 255).astype(np.uint8)
    return out


def mmi_field(res, sites, pad_deg=10.0, n=420):
    """
    พื้นสี MMI ที่คาด (ค่าสุดท้าย + ค่าปรับจากสถานี, Vs30 ค่าตั้งต้น) เป็นภาพ PNG (data URI)
    แถวของภาพเว้นระยะเท่ากันใน Web Mercator จึงวางบนแผนที่ได้ตรงโดยไม่บิด
    """
    import base64
    from io import BytesIO
    try:
        from PIL import Image
    except ImportError:
        return None
    hypo, p, (_, _, unc) = res["hypo"], res["preferred"], res["reliability"]
    cm = res.get("correction")
    corr = float(np.median([s["final"].get("corr", 0.0) for s in sites])) if sites else 0.0
    la = [hypo["lat"]] + [s["lat"] for s in sites]
    lo = [hypo["lon"]] + [s["lon"] for s in sites]
    south, north = max(min(la) - pad_deg, -80), min(max(la) + pad_deg, 80)
    west, east = min(lo) - pad_deg, max(lo) + pad_deg
    merc = lambda lat: np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))
    ys = np.linspace(merc(north), merc(south), n)
    lats = np.degrees(2 * np.arctan(np.exp(ys)) - np.pi / 2)
    lons = np.linspace(west, east, n)
    LA, LO = np.meshgrid(lats, lons, indexing="ij")
    D = _trace_dist_grid(LA, LO) if FAULT_TRACE else haversine_km(hypo["lat"], hypo["lon"], LA, LO)
    d = np.linspace(0.0, float(D.max()) + 10, 500)
    m = np.array([INT.predict_mmi(p["mag"], x, hypo["depth"], C.DEFAULT_VS30, EQ_TYPE, unc,
                                  x if FAULT_TRACE else None)["mmi"] for x in d])
    if cm and cm.get("method") not in (None, "none"):     # ค่าปรับตามตำแหน่ง (แนวโน้มระยะ + kriging)
        Rg = np.hypot(D, hypo["depth"])
        corr_g = model_corr(cm, LA, LO, Rg, np.full(D.shape, "target", object))
    else:
        corr_g = 0.0
    F = np.clip(np.interp(D, d, m) + corr_g, 1.0, 10.0)
    # จางลงที่ขอบภาพ (15% ของความกว้าง) จะได้ไม่เห็นเป็นกรอบสี่เหลี่ยม
    u = np.minimum(np.arange(n), np.arange(n)[::-1]) / (0.15 * n)
    edge = np.clip(np.minimum.outer(u, u), 0.0, 1.0)
    edge = edge * edge * (3 - 2 * edge)
    buf = BytesIO()
    Image.fromarray(_mmi_rgba(F, 0.75 * edge)).save(buf, format="PNG")
    return {"uri": "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode(),
            "bounds": [[float(south), float(west)], [float(north), float(east)]], "corr": corr}


def make_station_map(res, sites, out=STATION_MAP_FILE):
    """MMI ที่วัดได้ของแต่ละสถานีที่ input เข้าไป (ตัวเลขใหญ่ + ชื่อ) บนพื้นสี MMI  แบบภาพ ShakeMap"""
    sts = [s for s in sites if s["kind"] == "station"]
    if not sts:
        return None
    h, p = res["hypo"], res["preferred"]
    field = mmi_field(res, sts)
    pts = []
    for s in sts:
        o = s["obs"]
        rows = [f"ห่างจากจุดศูนย์กลาง {s['delta']:,.0f} km"]
        if o:
            _, rom, shake, dmg = INT.mmi_class(o["mmi"])
            rows += [f"<b>MMI วัดได้ {o['mmi']:.1f} ({rom} {shake})</b>", f"ความเสียหาย: {dmg}",
                     f"PGA {o['pga']:.2f} gal, PGV {o['pgv']:.3f} cm/s", f"MMI จาก PGA {o['mmi_pga']:.1f}",
                     f"JMA ระดับ {o['jma_class']}"]
        else:
            rows.append("ไม่มีข้อมูลคลื่นที่ใช้ได้")
        rows.append(f"MMI ทำนาย {s['final']['mmi']:.1f}")
        rows.insert(0, f"<b>{sensor_name(s)}</b>" + (f" &nbsp;ช่อง {s['channels']}" if s.get("channels") else ""))
        if o and o.get("flags"):
            rows.append("⚠ QC: " + "; ".join(o["flags"]))
        pts.append(dict(name=s["name"], lat=s["lat"], lon=s["lon"], mmi=obs_mmi(s),
                        sensor=sensor_key(s), popup="<br>".join(rows)))
    epi = dict(lat=h["lat"], lon=h["lon"], popup=_epi_html(res))
    title = "ระดับความรุนแรงแผ่นดินไหว"
    sub = f"{RM.thai_time(h['origin_time'])}   Mw {p['mag']:.1f}   ที่สถานีตรวจวัด (หน่วยวัด MMI)"
    try:
        from map_plot import make_station_mmi_map
    except (ImportError, AttributeError):
        make_station_mmi_map = None
    if make_station_mmi_map is not None:
        html = make_station_mmi_map(pts, out, epi, field, title=title, subtitle=sub)
        if html:
            return html
    import json
    data = dict(epi=epi, pts=pts, field=field, title=title, sub=sub,
                colors={str(k): c for k, c in MMI_COLORS.items()},
                roman={str(k): INT.MMI_TABLE[k][0] for k in range(1, 11)})
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(_STATION_MAP_TEMPLATE.replace("__DATA__", json.dumps(data, ensure_ascii=False)))
    return out


_STATION_MAP_TEMPLATE = """<!doctype html>
<html lang="th"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Station MMI Map</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
 html,body{margin:0;height:100%;font-family:Tahoma,"Leelawadee UI","Noto Sans Thai",sans-serif}
 #map{height:100%}
 .head{position:absolute;top:10px;left:50%;transform:translateX(-50%);z-index:1000;text-align:center;
       pointer-events:none;text-shadow:0 0 4px #fff,0 0 4px #fff}
 .head .t{font-size:28px;font-weight:bold;color:#e02020}.head .s{font-size:16px;color:#222}
 .lab{background:transparent;border:0;box-shadow:none;text-align:center;line-height:1.05;
      text-shadow:0 0 3px #fff,0 0 3px #fff,0 0 6px #fff;white-space:nowrap}
 .lab:before{display:none}
 .lab .v{font-size:26px;font-weight:bold;color:#333}.lab .n{font-size:15px;color:#222}
 .box{background:rgba(255,255,255,.92);padding:6px 10px;border-radius:6px;box-shadow:0 1px 5px rgba(0,0,0,.35);font-size:12px}
 .box td{width:30px;text-align:center;border:1px solid #999}
</style></head><body><div id="map"></div>
<script>
const D = __DATA__;
document.body.insertAdjacentHTML('afterbegin', '<div class="head"><div class="t">' + D.title +
  '</div><div class="s">' + D.sub + '</div></div>');
const map = L.map('map');
const ESRI = s => 'https://server.arcgisonline.com/ArcGIS/rest/services/' + s + '/MapServer/tile/{z}/{y}/{x}';
const bases = {
  'แผนที่สีเทาอ่อน (Esri)': L.tileLayer(ESRI('Canvas/World_Light_Gray_Base'), {maxZoom: 16, attribution: 'Esri, HERE, Garmin, &copy; OpenStreetMap contributors'}),
  'ภาพดาวเทียม': L.tileLayer(ESRI('World_Imagery'), {maxZoom: 19, attribution: 'Esri, Maxar, Earthstar Geographics'}),
  'แผนที่ถนน (Esri)': L.tileLayer(ESRI('World_Street_Map'), {maxZoom: 19, attribution: 'Esri'})};
bases['แผนที่สีเทาอ่อน (Esri)'].addTo(map);
// เส้นแบ่งเขตประเทศ/จังหวัด + ชื่อ อยู่บนสุด (ไม่ถูกพื้นสีบัง) และสลับชุดตามพื้นหลัง
map.createPane('labels'); map.getPane('labels').style.zIndex = 450; map.getPane('labels').style.pointerEvents = 'none';
map.createPane('labels_dark'); map.getPane('labels_dark').style.zIndex = 451; map.getPane('labels_dark').style.pointerEvents = 'none';
map.getPane('labels_dark').style.filter = 'contrast(2.2) brightness(0.55)';
const refSat = L.tileLayer(ESRI('Reference/World_Boundaries_and_Places'), {maxZoom: 16, pane: 'labels', attribution: 'Esri'});
const refMap = L.tileLayer(ESRI('Canvas/World_Light_Gray_Reference'), {maxZoom: 16, pane: 'labels_dark', attribution: 'Esri'}).addTo(map);
L.control.layers(bases, {'เส้นแบ่งเขตประเทศ/จังหวัด (ภาพดาวเทียม)': refSat, 'เส้นแบ่งเขตประเทศ/จังหวัด (แผนที่)': refMap},
  {collapsed: false}).addTo(map);
map.on('baselayerchange', e => { const sat = e.name === 'ภาพดาวเทียม';
  const on = sat ? refSat : refMap, off = sat ? refMap : refSat;
  if (map.hasLayer(off)) map.removeLayer(off); if (!map.hasLayer(on)) map.addLayer(on); });
if (D.field) L.imageOverlay(D.field.uri, D.field.bounds, {opacity: 0.65}).addTo(map);
const col = v => v === null ? '#bbb' : D.colors[String(Math.min(10, Math.max(1, Math.floor(v + 0.5))))];
const SHAPE = (sensor, c, z) => sensor === 'seismic'
  ? '<svg width="' + z + '" height="' + z + '" style="display:block"><polygon points="' + z/2 + ',1.5 ' + (z-1.5) + ',' + (z-1.5) + ' 1.5,' + (z-1.5) + '" fill="' + c + '" stroke="#222" stroke-width="1.6"/></svg>'
  : '<svg width="' + z + '" height="' + z + '" style="display:block"><rect x="2" y="2" width="' + (z-4) + '" height="' + (z-4) + '" fill="' + c + '" stroke="#222" stroke-width="1.6"/></svg>';
const SNAME = {seismic: 'สถานีตรวจวัดแผ่นดินไหว (Seismic station)', strong_motion: 'สถานีตรวจวัดอัตราเร่งของพื้นดิน (Strong-motion station)'};
const SSYM = {seismic: '&#9650;', strong_motion: '&#9632;'};
const SLEG = () => ['seismic', 'strong_motion'].filter(k => D.pts.some(p => (p.kind === undefined || p.kind === 'station') && p.sensor === k))
  .map(k => '<div style="display:flex;align-items:center;gap:6px;margin-top:2px">' + SHAPE(k, '#fff', 14) + SNAME[k] + '</div>').join('');
const rom = v => D.roman[String(Math.min(10, Math.max(1, Math.floor(v + 0.5))))];
L.marker([D.epi.lat, D.epi.lon], {icon: L.divIcon({className:'', html:
  '<div style="font-size:28px;color:#f39c12;line-height:28px;-webkit-text-stroke:1px #a04000">&#9733;</div>',
  iconSize:[28,28], iconAnchor:[14,14]})}).addTo(map).bindPopup(D.epi.html || D.epi.popup);
const b = [[D.epi.lat, D.epi.lon]];
D.pts.forEach(p => {
  b.push([p.lat, p.lon]);
  const txt = p.mmi === null ? '-' : p.mmi.toFixed(1);
  L.marker([p.lat, p.lon], {icon: L.divIcon({className: '', html: SHAPE(p.sensor, col(p.mmi), 18),
      iconSize: [18, 18], iconAnchor: [9, 10]})})
   .addTo(map).bindPopup('<b>' + p.name + '</b><br>' + p.popup)
   .bindTooltip('<div class="v">' + txt + '</div><div class="n">' + p.name + '</div>',
     {permanent: true, direction: 'top', offset: [0, -8], className: 'lab'});
});
map.fitBounds(b, {padding: [60, 60]});
const leg = L.control({position: 'bottomleft'});
leg.onAdd = () => { const d = L.DomUtil.create('div', 'box');
  d.innerHTML = '<b>หน่วยวัด: MMI (Modified Mercalli Intensity)</b><table style="border-collapse:collapse;margin-top:4px"><tr>' +
    [1,2,3,4,5,6,7,8,9,10].map(k => '<td style="background:' + D.colors[String(k)] + '">' + D.roman[String(k)] + '</td>').join('') +
    '</tr></table>ตัวเลข = MMI ที่วัดได้ที่สถานี (จาก PGV, Worden et al. 2012)<br>' +
    'พื้นสี = MMI ที่คาดจากขนาดและระยะ' + (D.field && D.field.corr ? ' (ปรับด้วยค่าที่สถานีวัดได้ ' +
    (D.field.corr > 0 ? '+' : '') + D.field.corr.toFixed(1) + ')' : '') + SLEG();
  return d; };
leg.addTo(map);
</script></body></html>
"""


# ======================================================================
#  ฟังก์ชันหลัก
# ======================================================================
def compute():
    RTT.SELECTED_STATIONS, RTT.PROFILE = SELECTED_STATIONS, PROFILE
    RTT.EPICENTER, RTT.DEPTH_KM = EPICENTER, DEPTH_KM
    res = RTT.analyse()
    if res is None or res.get("hypo") is None or res.get("preferred") is None or not res["stations"]:
        print("ยังคำนวณตำแหน่งและขนาดไม่ได้ จึงประมาณความรุนแรงไม่ได้")
        return None
    stages = RTT.system_stages(res)
    sites = build_sites(res)
    updates = simulate(res, stages, sites)
    return dict(res=res, stages=stages, sites=sites, updates=updates)


def main(mode=None, plot=PLOT):
    mode = mode or OUTPUT_MODE
    out = compute()
    if out is None:
        return None
    res, stages, sites, updates = out["res"], out["stages"], out["sites"], out["updates"]
    if mode == "detail":
        public_report(res, sites)
        detail_report(res, stages, updates, sites)
    else:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            detail_report(res, stages, updates, sites)
        with open(DETAIL_FILE, "w", encoding="utf-8") as f:
            f.write(buf.getvalue())
        public_report(res, sites)
        print(f"(ลำดับการอัปเดตและตารางเทียบ: {DETAIL_FILE})")
    if plot:
        for f in (plot_alert(res, sites), plot_distance(res, sites)):
            if f:
                print(f"บันทึกกราฟ: {f}")
    if MAKE_MAP:
        html2 = make_station_map(res, sites)
        if html2:
            print(f"บันทึกแผนที่ MMI ของแต่ละสถานี: {html2}")
        html = make_map(res, sites)
        if html:
            print(f"บันทึกแผนที่จริง: {html}  (เปิดด้วย browser หรือพิมพ์ start {html})")
    return out


if __name__ == "__main__":
    main()