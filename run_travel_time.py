"""
run_travel_time.py
คำนวณเวลาที่คลื่น P และ S เดินทางจากจุดศูนย์กลางถึงตำแหน่งของเรา และเวลาเตือนล่วงหน้า (Lead time)

  1) ระยะทาง:     Δ = Haversine(จุดศูนย์กลาง, ตำแหน่ง),  R = √(Δ² + h²)
  2) เวลาเดินทาง: โลก 2 ชั้น เลือกคลื่นที่มาถึงก่อนเอง (Pg/Pn, Sg/Sn)
                  t = min( √(Δ² + h²)/V1 ,  Δ/V2 + 2H·√(1/V1² − 1/V2²) )
  3) เวลามาถึง:   Tp = T0 + tP,  Ts = T0 + tS
  4) เวลาเตือน:   T_alert = Tp(สถานีที่ n) + 3 s (รอคลื่น P) + เวลาประมวลผล     (n = min(จำนวนสถานี, 3))
                  หรือ T_alert = T0 + ALERT_DELAY_S เมื่อไม่ได้เลือกสถานี
  5) Lead time = Ts(ตำแหน่ง) − T_alert        (≤ 0 = อยู่ใน Blind zone เตือนไม่ทัน)
     Blind zone = ระยะ Δ ที่ Ts(Δ) = T_alert
  ผลลัพธ์: ตาราง, กราฟ travel_time.png, แผนที่ travel_time_map.html

รัน:  ใส่ EPICENTER, ORIGIN_TIME และ TARGET_SITES ด้านล่าง แล้วกด Run ใน VS Code
"""
from datetime import datetime, timezone, timedelta

import numpy as np

import config as C
from data_loader import load_station_txt, find_station_file
from p_wave import detect_p_wave
from epicenter import travel_time_2layer
from signal_utils import haversine_km
from station_db import get_coords
from map_plot import make_travel_map

# ======================================================================
#  ตั้งค่าตรงนี้
# ======================================================================
# จุดศูนย์กลาง (lat, lon, ความลึก km) เช่น ผลจาก run_network.py หรือจาก USGS / กรมอุตุฯ
# ค่าตัวอย่างคือแผ่นดินไหวเมียนมา 28 มี.ค. 2025 (ค่าโดยประมาณ ควรตรวจกับ USGS ก่อนใช้)
EPICENTER = (22.013, 95.922, 10.0)

# เวลาเกิดแผ่นดินไหว (UTC) รูปแบบ "YYYY-MM-DDTHH:MM:SS"
# None = ประมาณจากเวลาคลื่น P ของสถานีใน SELECTED_STATIONS
ORIGIN_TIME = "2025-03-28T06:20:52"

# สถานีที่ใช้ตรวจจับ (ใช้หาเวลาที่ระบบออกคำเตือนจริง)  [] = ใช้ ALERT_DELAY_S แทน
SELECTED_STATIONS = ["BKKA"]
ALERT_DELAY_S = 10.0          # ถ้าไม่เลือกสถานี: สมมติว่าระบบเตือนได้หลังเกิดเหตุกี่วินาที
PROCESSING_DELAY_S = 1.0      # เวลาประมวลผลและส่งคำเตือน

# ตำแหน่งที่ต้องการทราบเวลาคลื่นมาถึง (พิกัดเป็นค่าโดยประมาณ ควรตรวจจาก Google Maps)
TARGET_SITES = [
    {"name": "มธ. รังสิต", "lat": 14.07, "lon": 100.61},
    # {"name": "เชียงใหม่", "lat": 18.79, "lon": 98.98},
]

MAKE_PLOT = True
MAKE_MAP = True
THAI_TZ = timezone(timedelta(hours=7))


def fmt(t):
    """epoch seconds -> 'UTC  (เวลาไทย)'"""
    u = datetime.fromtimestamp(t, tz=timezone.utc)
    return f"{u.strftime('%H:%M:%S.%f')[:-4]} UTC ({u.astimezone(THAI_TZ).strftime('%H:%M:%S')} น.)"


def phase_name(delta, depth, phase):
    """บอกว่าคลื่นที่มาถึงก่อนเป็นคลื่นตรง (g) หรือคลื่นใต้เปลือกโลก (n)"""
    v1, v2 = (C.VP_KM_S, C.VPN_KM_S) if phase == "P" else (C.VS_KM_S, C.VSN_KM_S)
    direct = np.sqrt(delta ** 2 + depth ** 2) / v1
    head = delta / v2 + 2 * C.CRUST_THICKNESS_KM * np.sqrt(1 / v1 ** 2 - 1 / v2 ** 2)
    return phase + ("g" if direct <= head else "n")


def blind_zone_km(t_after_origin, depth):
    """ระยะ Δ ที่คลื่น S ใช้เวลาเดินทางเท่ากับ t (ภายในรัศมีนี้ คลื่น S มาถึงก่อนคำเตือน)"""
    d = np.linspace(0, 4000, 8001)
    ts = travel_time_2layer(d, depth, "S")
    return float(np.interp(t_after_origin, ts, d)) if t_after_origin > ts[0] else 0.0


def station_picks(codes):
    """หาเวลาคลื่น P ของสถานีที่เลือก"""
    out = []
    for code in codes:
        path = find_station_file(code)
        if path is None:
            print(f"⚠ ไม่พบไฟล์ข้อมูลของสถานี {code}"); continue
        coords = get_coords(code)
        if coords is None:
            print(f"⚠ ไม่พบพิกัดของสถานี {code} ใน stations.csv"); continue
        rec = load_station_txt(path)
        p = detect_p_wave(rec["z"], rec["fs"], rec["t0"], **C.P_REGIONAL)
        if not p["detected"]:
            print(f"⚠ สถานี {code}: ไม่พบคลื่น P"); continue
        out.append({"code": code, "lat": coords[0], "lon": coords[1], "tp": p["time"]})
    return sorted(out, key=lambda s: s["tp"])


def main():
    lat0, lon0, h = EPICENTER
    print(f"จุดศูนย์กลาง : lat {lat0}, lon {lon0}, ลึก {h} km")

    # ---------- เวลาเกิด T0 ----------
    picks = station_picks(SELECTED_STATIONS) if SELECTED_STATIONS else []
    if ORIGIN_TIME:
        T0 = datetime.fromisoformat(ORIGIN_TIME).replace(tzinfo=timezone.utc).timestamp()
        src = "กำหนดเอง"
    elif picks:
        T0 = float(np.median([s["tp"] - travel_time_2layer(haversine_km(lat0, lon0, s["lat"], s["lon"]), h, "P")
                              for s in picks]))
        src = f"ประมาณจากคลื่น P ของ {len(picks)} สถานี"
    else:
        print("ต้องใส่ ORIGIN_TIME หรือเลือกสถานีใน SELECTED_STATIONS"); return
    print(f"เวลาเกิด T0  : {fmt(T0)}   ({src})")

    # ---------- เวลาออกคำเตือน ----------
    if picks:
        print("\nสถานีที่ตรวจพบคลื่น P:")
        for s in picks:
            d = haversine_km(lat0, lon0, s["lat"], s["lon"])
            print(f"  {s['code']:<6} Δ = {d:6.0f} km   P มาถึง {fmt(s['tp'])}   (หลังเกิดเหตุ {s['tp'] - T0:6.1f} s)")
        n = min(len(picks), 3)
        T_alert = picks[n - 1]["tp"] + C.PD_TAUC_WINDOW_SEC + PROCESSING_DELAY_S
        how = (f"คลื่น P ถึงสถานีที่ {n} ({picks[n - 1]['code']}) + {C.PD_TAUC_WINDOW_SEC:.0f} s "
               f"+ ประมวลผล {PROCESSING_DELAY_S:.0f} s")
    else:
        T_alert = T0 + ALERT_DELAY_S
        how = f"สมมติ {ALERT_DELAY_S:.0f} s หลังเกิดเหตุ"
    blind = blind_zone_km(T_alert - T0, h)
    print(f"\nออกคำเตือน  : {fmt(T_alert)}   (หลังเกิดเหตุ {T_alert - T0:.1f} s; {how})")
    print(f"Blind zone   : รัศมี {blind:.0f} km จากจุดศูนย์กลาง (ภายในวงนี้คลื่น S มาถึงก่อนคำเตือน)")

    # ---------- เวลาเดินทางถึงแต่ละตำแหน่ง ----------
    print("\n================ เวลาคลื่นเดินทางถึงตำแหน่งที่สนใจ ================")
    rows = []
    for s in TARGET_SITES:
        d = haversine_km(lat0, lon0, s["lat"], s["lon"])
        R = float(np.sqrt(d ** 2 + h ** 2))
        tp = float(travel_time_2layer(d, h, "P"))
        ts = float(travel_time_2layer(d, h, "S"))
        lead = T0 + ts - T_alert
        rows.append(dict(s, delta=d, R=R, tp=tp, ts=ts, lead=lead))
        status = f"เตือนล่วงหน้าได้ {lead:.1f} s" if lead > 0 else "อยู่ใน Blind zone (เตือนไม่ทัน)"
        print(f"\n[{s['name']}]  ({s['lat']}, {s['lon']})")
        print(f"  ระยะทาง          : Δ = {d:.0f} km, R = {R:.0f} km")
        print(f"  คลื่น P ({phase_name(d, h, 'P')})     : ใช้เวลา {tp:7.1f} s   มาถึง {fmt(T0 + tp)}")
        print(f"  คลื่น S ({phase_name(d, h, 'S')})     : ใช้เวลา {ts:7.1f} s   มาถึง {fmt(T0 + ts)}")
        print(f"  S − P            : {ts - tp:.1f} s")
        print(f"  ระบบออกคำเตือน   : {fmt(T_alert)}")
        print(f"  ผล               : {status}")

    # ---------- กราฟเวลาเดินทาง ----------
    if MAKE_PLOT:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import warnings
            warnings.filterwarnings("ignore", message="Glyph")       # เครื่องที่ไม่มีฟอนต์ไทย
            # ฟอนต์ที่มีอักษรไทย (Windows: Tahoma / Leelawadee UI) ใช้ตัวแรกที่มีในเครื่อง
            plt.rcParams["font.family"] = ["Tahoma", "Leelawadee UI", "TH Sarabun New", "Noto Sans Thai",
                                           "DejaVu Sans"]
            dmax = max([r["delta"] for r in rows] + [haversine_km(lat0, lon0, s["lat"], s["lon"]) for s in picks]
                       + [blind, 100.0]) * 1.15
            d = np.linspace(0, dmax, 600)
            fig, ax = plt.subplots(figsize=(11, 6))
            ax.plot(d, travel_time_2layer(d, h, "P"), "b-", label="P wave (Pg/Pn)")
            ax.plot(d, travel_time_2layer(d, h, "S"), "r-", label="S wave (Sg/Sn)")
            ax.axhline(T_alert - T0, color="k", ls="--", lw=1, label=f"alert issued ({T_alert - T0:.0f} s)")
            ax.axvspan(0, blind, color="red", alpha=0.1, label=f"blind zone (< {blind:.0f} km)")
            for s in picks:
                ds = haversine_km(lat0, lon0, s["lat"], s["lon"])
                ax.plot(ds, s["tp"] - T0, "b^", ms=9)
                ax.annotate(s["code"], (ds, s["tp"] - T0), textcoords="offset points", xytext=(5, -12))
            for r in rows:
                ax.plot([r["delta"], r["delta"]], [T_alert - T0, r["ts"]], "g-", lw=3)
                ax.plot(r["delta"], r["ts"], "rs", ms=8)
                ax.annotate(f"{r['name']}\nlead {r['lead']:.0f} s", (r["delta"], r["ts"]),
                            textcoords="offset points", xytext=(6, 4), fontsize=9)
            ax.set_xlabel("Epicentral distance (km)"); ax.set_ylabel("Time after origin (s)")
            ax.set_title("Travel time and warning lead time (green bar = lead time)")
            ax.grid(alpha=0.3); ax.legend(loc="upper left")
            fig.tight_layout()
            fig.savefig("travel_time.png", dpi=120)
            plt.close(fig)
            print("\nบันทึกกราฟ: travel_time.png")
        except ImportError:
            pass

    # ---------- แผนที่จริง ----------
    if MAKE_MAP:
        sites = [{"name": r["name"], "lat": r["lat"], "lon": r["lon"], "lead": r["lead"],
                  "info": f"Δ {r['delta']:.0f} km<br>P มาถึง {fmt(T0 + r['tp'])}<br>S มาถึง {fmt(T0 + r['ts'])}"
                          f"<br>เวลาเตือนล่วงหน้า {r['lead']:.1f} s"} for r in rows]
        t_end = max([r["ts"] for r in rows] + [T_alert - T0]) + 30
        fronts = [(t, blind_zone_km(t, h)) for t in np.arange(30, t_end, 30)]
        html = make_travel_map(EPICENTER, sites, "travel_time_map.html", blind_km=blind,
                               stations=picks, s_fronts=fronts)
        if html:
            print(f"บันทึกแผนที่จริง: {html}  (เปิดด้วย browser หรือพิมพ์ start {html})")
    return rows


if __name__ == "__main__":
    main()
