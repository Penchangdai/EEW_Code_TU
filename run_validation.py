"""
run_validation.py
ทดสอบระบบกับแผ่นดินไหวหลายเหตุการณ์ แล้วเทียบกับค่าจริง (USGS / TMD)

  แต่ละเหตุการณ์อยู่คนละโฟลเดอร์ (ไฟล์ .txt ของแต่ละสถานี เหมือนโฟลเดอร์ station)
  ระบบคำนวณเองทั้งหมด (ตำแหน่ง ขนาด MMI การแจ้งเตือน) เหมือน run_intensity.py แล้ววัดความคลาด:
    - ตำแหน่ง   : ระยะจากจุดศูนย์กลางจริง (km)  และเวลาเกิดคลาด (s)
    - ขนาด      : Mw ของระบบ − Mw จริง  (ขนาดแรกที่ใช้แจ้งเตือน และขนาดสุดท้าย)
    - ความรุนแรง : MMI วัด − ทำนาย ก่อนปรับ (bias, RMSE) และหลังปรับแบบ leave-one-out
    - แจ้งเตือน  : จำนวนสถานีที่ได้แจ้งเตือน เทียบกับจำนวนที่สั่นจริงถึงเกณฑ์ (เตือนผิด / ไม่ได้เตือน)

  ผลลัพธ์: ตารางสรุปบนหน้าจอ + validation_results.csv
  เพิ่มเหตุการณ์: เพิ่มบรรทัดใน EVENTS (ค่าจริงจากหน้าเหตุการณ์ของ USGS)
"""
import contextlib
import csv
import io
import os
from datetime import datetime, timezone

import numpy as np

import data_loader as DL
import magnitude as M
import run_magnitude as RM
import run_intensity as RI
from signal_utils import haversine_km

CODE_DIR = os.path.dirname(os.path.abspath(__file__))

# ======================================================================
#  เหตุการณ์ที่ใช้ทดสอบ  (folder = โฟลเดอร์ข้อมูล เทียบกับโฟลเดอร์โค้ด หรือ path เต็ม)
#  mag_type "mb" = USGS ไม่มี Mw จะแปลงเป็น Mw ด้วย Scordilis (2006) ก่อนเทียบ
#  lat = None = ยังไม่ได้ใส่ค่าจริง (ระบบยังรัน แต่ไม่คิดความคลาดของตำแหน่ง/ขนาด)
# ======================================================================
EVENTS = [
    dict(name="Mw 7.7 เมียนมา (สะกาย)", folder="2025-3-28_Mw7.7_station",
         lat=22.011, lon=95.936, depth=10.0, time="2025-03-28T06:20:52.7", mag=7.7, mag_type="Mww",
         ref="USGS us7000pn9s"),
    dict(name="mb 5.1 เมียนมา (เนปิดอว์)", folder="2025-03-29_M5.9_AllTMD",
         lat=19.6858, lon=96.0761, depth=10.0, time="2025-03-29T09:20:48.2", mag=5.1, mag_type="mb",
         ref="USGS us7000pnkj (ชื่อโฟลเดอร์ใช้ขนาด 5.9)"),
    dict(name="Mw 8.8 คัมชัตกา รัสเซีย", folder="2025-07-29_M8.8_Russia_AllTMD",
         lat=52.4948, lon=160.2395, depth=35.0, time="2025-07-29T23:24:52.5", mag=8.8, mag_type="Mww",
         ref="USGS us6000qw60"),
    dict(name="mb 4.6 นอกชายฝั่งเมียนมา", folder="2025-08-21_M4.6_AllTMD",
         lat=16.172, lon=96.535, depth=10.0, time="2025-08-21T02:58:07", mag=4.6, mag_type="mb",
         ref="USGS (TMD: 16.145, 96.686 ขนาด 5.4)"),
    dict(name="Mw 5.2 ชายฝั่งเมียนมา", folder="2026-05-18_M5.2_AllTMD",
         lat=16.4177, lon=96.3726, depth=10.0, time="2026-05-18T02:05:24", mag=5.2, mag_type="Mww",
         ref="USGS us6000syg0"),
]
OUT_CSV = "validation_results.csv"
USE_KNOWN_RUPTURE = False    # False = วัดแบบตอนเกิดเหตุจริง ซึ่งยังไม่รู้แนวรอยแตก
MAX_REGIONAL_KM = 2500.0     # ไกลกว่านี้ = แผ่นดินไหวไกล (teleseism) ระบบไม่ควรแจ้งเตือนในไทย
# ======================================================================


def _epoch(t):
    return datetime.fromisoformat(t).replace(tzinfo=timezone.utc).timestamp()


def ref_mw(ev):
    """ขนาดจริงเป็น Mw (แปลง mb / ML ถ้าจำเป็น)"""
    m, t = ev.get("mag"), str(ev.get("mag_type", "")).lower()
    if m is None:
        return None
    if t.startswith("mb"):
        return 0.85 * m + 1.03                        # Scordilis (2006)
    if t == "ml":
        return 0.953 * m + 0.422                      # Akkar et al. (2010)
    return float(m)


def _folder(ev):
    f = ev["folder"]
    if os.path.isabs(f):
        return f
    for cand in (os.path.join(CODE_DIR, "station", f), os.path.join(CODE_DIR, f)):   # station/<event> ก่อน
        if os.path.isdir(cand):
            return cand
    return os.path.join(CODE_DIR, "station", f)


def run_event(ev):
    folder = _folder(ev)
    if not os.path.isdir(folder):
        return dict(status=f"ไม่พบโฟลเดอร์ {ev['folder']}")
    DL.DATA_DIR = RM.DATA_DIR = folder
    RI.SELECTED_STATIONS = "auto"
    RI.EPICENTER = None
    RI.USE_KNOWN_RUPTURE = USE_KNOWN_RUPTURE
    log = io.StringIO()
    try:
        with contextlib.redirect_stdout(log):
            out = RI.compute()
    except Exception as e:                                   # noqa: BLE001
        return dict(status=f"error: {e}")
    if out is None:
        return dict(status="คำนวณไม่ได้")
    res, sites, updates = out["res"], out["sites"], out["updates"]
    h, p = res["hypo"], res["preferred"]
    r = dict(status="ok", n_files=len(RM.station_files("auto")), n_used=len(res["stations"]),
             lat=h["lat"], lon=h["lon"], M=p["mag"], M_type=p["type"], M_raw=p.get("mag_raw"),
             skipped_units=log.getvalue().count("ไม่ใช่หน่วยจริง"),
             teleseism=bool((h.get("teleseism") or {}).get("flag")))
    if updates:
        u = updates[0]
        r.update(M_first=u["M"], M_first_type=u["type"], t_first=u["t"] - h["origin_time"])
    mw = ref_mw(ev)
    if ev.get("lat") is not None:
        r["loc_err_km"] = float(haversine_km(h["lat"], h["lon"], ev["lat"], ev["lon"]))
        r["t0_err_s"] = h["origin_time"] - _epoch(ev["time"]) if ev.get("time") else None
        r["dist_min_km"] = float(min(haversine_km(ev["lat"], ev["lon"], s["lat"], s["lon"]) for s in sites
                                     if s["kind"] == "station"))
    if mw is not None:
        r["ref_mw"] = mw
        for t, v in res.get("network", {}).items():         # ความคลาดของแต่ละชนิด (ก่อนค่าปรับเครือข่ายไทย)
            if isinstance(v, dict) and v.get("mag_raw") is not None:
                r[f"dM_{t}"] = M.to_mw(t, v["mag_raw"], regional=False)[0] - mw
        r["dM_final"] = p["mag"] - mw
        if "M_first" in r:
            r["dM_first"] = r["M_first"] - mw
    ok = [s for s in sites if s["kind"] == "station" and RI.obs_mmi(s) is not None]
    if ok:
        e = np.array([RI.obs_mmi(s) - s["final_raw"]["mmi"] for s in ok])
        r.update(n_mmi=len(ok), mmi_bias=float(e.mean()), mmi_rmse=float(np.sqrt(np.mean(e ** 2))),
                 mmi_obs_max=float(max(RI.obs_mmi(s) for s in ok)))
    cm = res.get("correction")
    if cm and cm.get("rmse") is not None:
        r["mmi_loo"] = cm["rmse"]
    st = [s for s in sites if s["kind"] == "station"]
    pub = [s for s in st if s.get("alert") and max(s["alert"]["level"], (s.get("upgrade") or {}).get("level", 0)) == 2]
    app = [s for s in st if s.get("alert")]
    felt = [s for s in st if (RI.obs_mmi(s) or 0) >= 2.5]
    strong = [s for s in st if (RI.obs_mmi(s) or 0) >= 3.5]
    r.update(n_alert=len(app), n_public=len(pub), n_obs_ge3=len(felt), n_obs_ge4=len(strong),
             false_public=len([s for s in pub if RI.obs_mmi(s) is not None and RI.obs_mmi(s) < 2.5]),
             missed=len([s for s in strong if not s.get("alert")]))
    return r


def _f(x, fmt):
    return "-" if x is None else format(x, fmt)


def main():
    rows = []
    for ev in EVENTS:
        print(f"กำลังทดสอบ: {ev['name']}  ({ev['folder']})")
        r = run_event(ev)
        r.update(event=ev["name"], folder=ev["folder"], ref=ev.get("ref", ""))
        rows.append(r)
        if r["status"] != "ok":
            print(f"   ข้าม: {r['status']}")
            continue
        print(f"   ใช้ {r['n_used']} สถานี (ข้ามไฟล์ที่ยังเป็น count {r['skipped_units']})"
              f"   ตำแหน่ง {r['lat']:.2f}, {r['lon']:.2f}   คลาด {_f(r.get('loc_err_km'), '.0f')} km"
              f"   เวลาเกิดคลาด {_f(r.get('t0_err_s'), '+.1f')} s")
        print(f"   ขนาดแรก {r.get('M_first_type', '-')} {_f(r.get('M_first'), '.1f')} (+{_f(r.get('t_first'), '.0f')} s)"
              f"   สุดท้าย {r['M_type']} -> Mw {r['M']:.1f}   จริง Mw {_f(r.get('ref_mw'), '.1f')}"
              f"   คลาด {_f(r.get('dM_final'), '+.1f')}")
        print(f"   MMI ก่อนปรับ bias {_f(r.get('mmi_bias'), '+.2f')} RMSE {_f(r.get('mmi_rmse'), '.2f')}"
              f"   หลังปรับ (LOO) {_f(r.get('mmi_loo'), '.2f')}   วัดได้สูงสุด {_f(r.get('mmi_obs_max'), '.1f')}")
        print(f"   แจ้งเตือน {r['n_alert']} สถานี (สาธารณะ {r['n_public']})  สั่นจริง >= III {r['n_obs_ge3']}"
              f"  >= IV {r['n_obs_ge4']}  เตือนสาธารณะแต่สั่น < III {r['false_public']}  ไม่ได้เตือนแต่ >= IV {r['missed']}")
        far = r.get("dist_min_km", 0) > MAX_REGIONAL_KM
        if far and r["n_alert"]:
            print("   ⚠ แผ่นดินไหวไกลมาก (teleseism) แต่ระบบแจ้งเตือน")
        elif far and r["teleseism"]:
            print("   ✓ ระบบตรวจพบว่าเป็นแผ่นดินไหวไกลมาก (teleseism) และไม่แจ้งเตือน")
        elif r["teleseism"]:
            print("   ⚠ ระบบตัดสินว่าเป็น teleseism แต่เป็นแผ่นดินไหวในภูมิภาค (ไม่แจ้งเตือนทั้งที่ควร)")

    print("\n" + "=" * 100)
    print(f"{'เหตุการณ์':<28}{'ตำแหน่ง km':>11}{'T0 s':>7}{'Mw จริง':>8}{'Mw แรก':>8}{'Mw สุดท้าย':>11}"
          f"{'MMI bias':>9}{'RMSE':>6}{'LOO':>6}{'เตือนผิด':>9}")
    for r in rows:
        if r["status"] != "ok":
            print(f"{r['event'][:27]:<28}  {r['status']}")
            continue
        print(f"{r['event'][:27]:<28}{_f(r.get('loc_err_km'), '11.0f')}{_f(r.get('t0_err_s'), '+7.1f')}"
              f"{_f(r.get('ref_mw'), '8.1f')}{_f(r.get('M_first'), '8.1f')}{r['M']:11.1f}"
              f"{_f(r.get('mmi_bias'), '+9.2f')}{_f(r.get('mmi_rmse'), '6.2f')}{_f(r.get('mmi_loo'), '6.2f')}"
              f"{r['false_public']:9d}")
    print("\nความคลาดของขนาดแต่ละชนิด (Mw ของระบบ − Mw จริง ก่อนใส่ค่าปรับเครือข่ายไทย)")
    types = list(M.TYPES)
    print(f"{'เหตุการณ์':<28}" + "".join(f"{t:>8}" for t in types))
    for r in rows:
        if r.get("status") == "ok" and r.get("ref_mw") is not None:
            print(f"{r['event'][:27]:<28}" + "".join(_f(r.get(f"dM_{t}"), "+8.2f") for t in types))
    print(f"{'ค่าปรับที่แนะนำ (−ค่าเฉลี่ย)':<28}", end="")
    sug = {}
    for t in types:
        v = [r[f"dM_{t}"] for r in rows if r.get(f"dM_{t}") is not None and r.get("ref_mw", 0) < 7.0]
        sug[t] = -float(np.mean(v)) if len(v) >= 2 else None
        print(_f(sug[t], "+8.2f"), end="")
    print("\n  (ใช้เฉพาะเหตุการณ์ Mw < 7 เพราะขนาดใหญ่อิ่มตัว  '-' = ข้อมูลไม่พอ  ใส่ใน config.py: MW_REGIONAL_CORR = {...})")
    print(f"  ค่าที่ใช้อยู่ตอนนี้: {M.MW_REGIONAL_CORR}")

    keys = sorted({k for r in rows for k in r})
    with open(os.path.join(CODE_DIR, OUT_CSV), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"\nบันทึก {OUT_CSV}")
    return rows


if __name__ == "__main__":
    main()
