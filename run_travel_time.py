"""
run_travel_time.py
คำนวณเวลาที่คลื่น P, S และ Lg เดินทางจากจุดศูนย์กลางถึงตำแหน่งของเรา และเวลาเตือนล่วงหน้า (Lead time)
ใช้ผลชุดเดียวกับ run_magnitude.py / pipeline.py (ตำแหน่ง เวลาเกิด เวลาคลื่นของแต่ละสถานี และขนาดหลัก)

  1) ตำแหน่ง + เวลาเกิด T0 : pipeline.py คำนวณจากคลื่น P + S + Lg ของทุกสถานี   (หรือใส่ EPICENTER เอง)
  2) ระยะทาง               : Δ = Haversine(จุดศูนย์กลาง, ตำแหน่ง)
  3) เวลาเดินทาง (โลก 2 ชั้น เหมือนที่ใช้หาตำแหน่ง)
        P, S : t = min( √(Δ² + h²)/V1 ,  Δ/V2 + 2H·√(1/V1² − 1/V2²) )     (Pg/Pn, Sg/Sn)
        Lg   : t = √(Δ² + h²) / 3.5        (ที่ระยะไกลหลายร้อย km คลื่น Lg มักสั่นแรงที่สุด)
  4) เวลาที่ระบบรู้ข้อมูลแต่ละขั้น (เวลาจริงที่ข้อมูลมาถึงสถานี + เวลาประมวลผล)
        ขั้น 1  ตรวจพบคลื่น P ครบ 3 สถานี + 3 s         -> ตำแหน่งจากคลื่น P อย่างเดียว (คลาดมากถ้าสถานีอยู่ด้านเดียว)
        ขั้น 2  ได้เวลาคลื่นทุกตัวที่ใช้หาตำแหน่งแล้ว   -> ตำแหน่งสุดท้ายที่แสดงในการ์ด
        ขั้น 3  ได้ข้อมูลที่ขนาดหลักต้องใช้ครบทุกสถานี  -> ขนาดหลักที่แสดงในการ์ด
  5) Lead time = เวลาคลื่นมาถึงตำแหน่ง − เวลาที่ระบบรู้ข้อมูลขั้นนั้น     (≤ 0 = เตือนไม่ทัน)
     Blind zone = รัศมีที่คลื่น S มาถึงก่อนระบบรู้ข้อมูลขั้นนั้น
  ผลลัพธ์: ตาราง, กราฟ travel_time.png, แผนที่ travel_time_map.html (ถ้ามี map_plot.make_travel_map)

รัน:  ใส่ TARGET_SITES ด้านล่าง แล้วกด Run ใน VS Code
"""
import contextlib
import io
from datetime import datetime, timezone, timedelta

import numpy as np

import config as C
import run_magnitude as RM
import travel_time as TT
from epicenter import travel_time_2layer, phase_travel_time
from signal_utils import haversine_km

# ======================================================================
#  ตั้งค่าตรงนี้
# ======================================================================
SELECTED_STATIONS = "auto"    # รหัสสถานี เช่น ["BKKA", "CMMT", "PHIT"]  หรือ "auto" = ทุกไฟล์ในโฟลเดอร์ station
PROFILE = "regional"

# จุดศูนย์กลาง   None = คำนวณเองจากข้อมูลสถานี (แนะนำ)
#                (lat, lon, ความลึก km, "เวลาเกิด UTC") = ใส่เอง เช่นค่าจาก USGS (ระบบจะบอกว่าไม่ได้คำนวณ)
EPICENTER = None
# EPICENTER = (22.011, 95.936, 10.0, "2025-03-28T06:20:52")      # USGS
DEPTH_KM = 10.0               # ตรึงความลึก (สถานีไกลหาความลึกไม่ได้)  None = ลองคำนวณความลึก

ALERT_MIN_STATIONS = 3        # ขั้น 1: ต้องตรวจพบคลื่น P กี่สถานี
PROCESSING_DELAY_S = 1.0      # เวลาประมวลผลและส่งคำเตือน (s)
ALERT_DELAY_S = 10.0          # ใช้เมื่อใส่ EPICENTER และ SELECTED_STATIONS = [] : สมมติว่าเตือนได้หลังเกิดเหตุกี่วินาที

# ตำแหน่งที่ต้องการทราบเวลาคลื่นมาถึง (พิกัดเป็นค่าโดยประมาณ ควรตรวจจาก Google Maps)
TARGET_SITES = [
    {"name": "มธ. รังสิต", "lat": 14.07, "lon": 100.61},
    # {"name": "เชียงใหม่", "lat": 18.79, "lon": 98.98},
]

MAKE_PLOT = True
MAKE_MAP = True
THAI_TZ = timezone(timedelta(hours=7))


# ======================================================================
#  เครื่องมือ
# ======================================================================
def fmt(t):
    """epoch seconds -> 'UTC  (เวลาไทย)'"""
    u = datetime.fromtimestamp(t, tz=timezone.utc)
    return f"{u.strftime('%H:%M:%S.%f')[:-4]} UTC ({u.astimezone(THAI_TZ).strftime('%H:%M:%S')} น.)"


def phase_name(delta, depth, phase):
    """คลื่นที่มาถึงก่อนเป็นคลื่นตรง (g) หรือคลื่นใต้เปลือกโลก (n)  (Lg คงชื่อเดิม)"""
    return TT.phase_name(delta, depth, phase)


def blind_zone_km(t_after_origin, depth):
    """ระยะ Δ ที่คลื่น S ใช้เวลาเดินทางเท่ากับ t (ภายในรัศมีนี้ คลื่น S มาถึงก่อนคำเตือน)"""
    return TT.blind_zone_radius(t_after_origin, 0.0, depth)


# ======================================================================
#  ขั้น 1  ตำแหน่ง + ขนาด (ผลชุดเดียวกับ run_magnitude.py)
# ======================================================================
def analyse():
    RM.EPICENTER, RM.DEPTH_KM = EPICENTER, DEPTH_KM
    if EPICENTER is not None and SELECTED_STATIONS == []:
        la, lo, dep, t_utc = EPICENTER
        hypo = {"lat": la, "lon": lo, "depth": dep, "origin_time": RM._epoch(t_utc), "computed": False,
                "depth_measured": False, "source": "ค่าที่ใส่เอง (EPICENTER) ไม่ได้คำนวณ"}
        return {"hypo": hypo, "stations": [], "preferred": None}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):            # รายละเอียดการคำนวณอยู่ใน run_magnitude.py
        res = RM.compute(SELECTED_STATIONS, PROFILE)
    for line in buf.getvalue().splitlines():         # แสดงเฉพาะคำเตือน
        if line.lstrip().startswith("⚠"):
            print(line)
    return res


def _station_time(st, sample):
    return st["rec"]["t0"] + sample / st["rec"]["fs"]


def magnitude_ready_time(st, mtype, T0):
    """เวลา (epoch) ที่สถานีนี้มีข้อมูลครบสำหรับขนาดชนิดนี้ (ตามหน้าต่างที่ magnitude.py ใช้วัด)"""
    m = st["mags"].get(mtype) or {}
    tp = _station_time(st, st["p_index"])
    d = st["delta"]
    if mtype == "ML":                                 # ยอดบน Wood-Anderson (เร็วที่สุดที่รู้ค่านี้)
        idx = [c["index"] for c in m.get("components", {}).values() if c.get("used")]
        return _station_time(st, max(idx)) if idx else tp
    if mtype == "mb_Lg":                              # หน้าต่าง Lg 3.6-3.2 km/s
        return T0 + d / M_LG_VMIN
    if mtype == "Ms_BB":                              # หน้าต่าง Rayleigh 4.5-2.5 km/s
        return T0 + d / M_RAYLEIGH_VMIN
    if mtype == "Mwp":                                # คลื่น P ถึงก่อนคลื่น S
        return tp + (m.get("window_sec") or 0.0)
    if mtype in ("M_Pd", "M_tc", "Mw_P"):             # คลื่น P 3 s แรก
        return tp + C.PD_TAUC_WINDOW_SEC
    return _station_time(st, len(st["rec"]["z"]))     # M_JMA: ค่าสูงสุดทั้งไฟล์


M_LG_VMIN = RM.M.LG_VMIN
M_RAYLEIGH_VMIN = RM.M.RAYLEIGH_VMIN


def system_stages(res):
    """เวลาที่ระบบรู้ข้อมูลแต่ละขั้น -> list ของ dict(name, time, how, note)"""
    hypo, T0 = res["hypo"], res["hypo"]["origin_time"]
    out = []
    if not res["stations"] and not hypo.get("picks_used"):
        out.append(dict(name="ออกคำเตือน", time=T0 + ALERT_DELAY_S,
                        how=f"สมมติ {ALERT_DELAY_S:.0f} s หลังเกิดเหตุ (ไม่ได้ใช้ข้อมูลสถานี)", note=""))
        return out

    picks = hypo.get("picks_used")
    if picks:
        p_times = sorted({(c, t) for c, ph, t, _, _ in picks if ph == "P"}, key=lambda x: x[1])
    else:                                             # EPICENTER ใส่เอง: เวลา P ที่ pick รอบเวลาที่คาด
        p_times = sorted(((st["code"], _station_time(st, st["p_index"])) for st in res["stations"]
                          if st["p_weight"] <= 3), key=lambda x: x[1])

    # ---- ขั้น 1  ตรวจพบคลื่น P ----
    n = min(ALERT_MIN_STATIONS, len(p_times))
    if n:
        note = ""
        hp = hypo.get("hypo_p")
        if hp:
            off = float(haversine_km(hp["lat"], hp["lon"], hypo["lat"], hypo["lon"]))
            note = f"ตำแหน่งจากคลื่น P อย่างเดียว ห่างจากตำแหน่งสุดท้าย {off:.0f} km"
        out.append(dict(name="ขั้น 1 ตรวจพบคลื่น P",
                        time=p_times[n - 1][1] + C.PD_TAUC_WINDOW_SEC + PROCESSING_DELAY_S,
                        how=f"คลื่น P ถึงสถานีที่ {n} ({p_times[n - 1][0]}) + {C.PD_TAUC_WINDOW_SEC:.0f} s "
                            f"+ ประมวลผล {PROCESSING_DELAY_S:.0f} s", note=note))

    # ---- ขั้น 2  ตำแหน่งสุดท้าย (ต้องได้คลื่น S / Lg ที่ใช้หาตำแหน่งครบ) ----
    if picks and any(ph != "P" for _, ph, _, _, _ in picks):
        last = max(picks, key=lambda x: x[2])
        out.append(dict(name="ขั้น 2 ตำแหน่งสุดท้าย", time=last[2] + PROCESSING_DELAY_S,
                        how=f"คลื่นตัวสุดท้ายที่ใช้หาตำแหน่ง: {last[1]} ที่ {last[0]} + ประมวลผล "
                            f"{PROCESSING_DELAY_S:.0f} s", note=hypo.get("source", "")))

    # ---- ขั้น 3  ขนาดหลัก ----
    pick = res.get("preferred")
    if pick:
        t = pick["type"]
        used = [st for st in res["stations"] if st.get("w", {}).get(t, 0) > 0]
        if used:
            ready = max(magnitude_ready_time(st, t, T0) for st in used)
            out.append(dict(name=f"ขั้น 3 ขนาดหลัก ({t})", time=ready + PROCESSING_DELAY_S,
                            how=f"ข้อมูลที่ {t} ต้องใช้ครบ {len(used)} สถานี + ประมวลผล {PROCESSING_DELAY_S:.0f} s",
                            note=f"Mw {pick['mag']:.1f}"))
    return out


# ======================================================================
#  ขั้น 2  เวลาเดินทางถึงตำแหน่งที่สนใจ
# ======================================================================
def site_rows(hypo, stages):
    """เวลาเดินทางถึงแต่ละตำแหน่ง (travel_time.site_warning) + เวลาเตือนล่วงหน้าของทุกขั้น"""
    rows = []
    for s in TARGET_SITES:
        w = TT.site_warning(hypo, s["lat"], s["lon"], stages[0]["time"] if stages else hypo["origin_time"])
        t = {"P": w["p_travel_s"], "S": w["s_travel_s"], "Lg": w["lg_travel_s"]}
        lead = {g["name"]: {"S": w["s_arrival"] - g["time"], "Lg": w["lg_arrival"] - g["time"]} for g in stages}
        rows.append(dict(s, delta=w["epi_dist_km"], R=w["hypo_dist_km"], t=t, lead=lead))
    return rows


def report(res, stages, rows):
    hypo, T0, h = res["hypo"], res["hypo"]["origin_time"], res["hypo"]["depth"]
    pos, src = RM.location_lines(hypo)
    print("=" * 72)
    print(f"ตำแหน่ง        : {pos}")
    print(f"ที่มาของตำแหน่ง : {src}")
    print(f"เวลาเกิด T0    : {fmt(T0)}")
    if res.get("preferred"):
        p, (lv, _, unc) = res["preferred"], res["reliability"]
        print(f"ขนาดหลัก       : Mw {RM.mw_text(p)} ({RM.mw_source(p)})   ความน่าเชื่อถือ {lv} ±{unc:.1f}")
    if not hypo.get("depth_measured"):
        print(f"(คำนวณเวลาเดินทางด้วยความลึก {h:.0f} km ที่ตรึงไว้ ไม่ได้วัด)")

    print("\n---------------- ระบบรู้ข้อมูลเมื่อไร (นับจากเวลาเกิด) ----------------")
    for g in stages:
        print(f"{g['name']:<26} +{g['time'] - T0:6.1f} s   {fmt(g['time'])}")
        print(f"{'':<26} {g['how']}")
        if g["note"]:
            print(f"{'':<26} {g['note']}")
        print(f"{'':<26} Blind zone รัศมี {blind_zone_km(g['time'] - T0, h):.0f} km")

    for r in rows:
        print(f"\n================ [{r['name']}]  ({r['lat']}, {r['lon']}) ================")
        print(f"  ระยะทาง  : Δ = {r['delta']:.0f} km, R = {r['R']:.0f} km")
        for ph in ("P", "S", "Lg"):
            print(f"  คลื่น {phase_name(r['delta'], h, ph):<3}  : ใช้เวลา {r['t'][ph]:7.1f} s   มาถึง {fmt(T0 + r['t'][ph])}")
        print(f"  S − P    : {r['t']['S'] - r['t']['P']:.1f} s")
        print("  เวลาเตือนล่วงหน้า (บวก = เตือนทัน, ลบ = คลื่นมาถึงก่อน)")
        print(f"  {'':<26}{'ก่อนคลื่น S':>14}{'ก่อนคลื่น Lg':>14}")
        for g in stages:
            ls, ll = r["lead"][g["name"]]["S"], r["lead"][g["name"]]["Lg"]
            print(f"  {g['name']:<26}{ls:+13.1f} s{ll:+13.1f} s")
    print("\nหมายเหตุ: ที่ระยะไกลหลายร้อย km แรงสั่นสูงสุดมักมากับคลื่น Lg และคลื่นผิวพื้นที่ตามมา "
          "คลื่น S ตัวแรก (Sn) มักยังไม่แรง")


# ======================================================================
#  กราฟ / แผนที่
# ======================================================================
def plot(res, stages, rows, out="travel_time.png"):
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
    hypo, T0, h = res["hypo"], res["hypo"]["origin_time"], res["hypo"]["depth"]
    picks = hypo.get("picks_used") or []
    dist = [float(haversine_km(hypo["lat"], hypo["lon"], la, lo)) for _, _, _, la, lo in picks]
    dmax = max([r["delta"] for r in rows] + dist + [100.0]) * 1.15
    d = np.linspace(0, dmax, 600)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.plot(d, travel_time_2layer(d, h, "P"), color="#1f77b4", label="คลื่น P (Pg/Pn)")
    ax.plot(d, travel_time_2layer(d, h, "S"), color="#d62728", label="คลื่น S (Sg/Sn)")
    ax.plot(d, phase_travel_time(d, h, "Lg"), color="#ff7f0e", ls="--", label="คลื่น Lg")
    colors = ("#2ca02c", "#9467bd", "#8c564b")
    for g, c in zip(stages, colors):
        ax.axhline(g["time"] - T0, color=c, ls=":", lw=1.5, label=f"{g['name']} (+{g['time'] - T0:.0f} s)")
    mk = {"P": ("^", "#1f77b4"), "S": ("s", "#d62728"), "Lg": ("D", "#ff7f0e")}
    seen = set()
    for (code, ph, t, _, _), ds in zip(picks, dist):
        m, c = mk[ph]
        ax.plot(ds, t - T0, m, color=c, ms=8, mec="k", label=f"pick {ph} ที่สถานี" if ph not in seen else None)
        seen.add(ph)
        if ph == "P":
            ax.annotate(code, (ds, t - T0), textcoords="offset points", xytext=(5, -12), fontsize=8)
    if stages:
        g0 = stages[0]
        for r in rows:
            ax.plot([r["delta"]] * 2, [g0["time"] - T0, r["t"]["S"]], color="#2ca02c", lw=4, alpha=0.6)
            ax.plot(r["delta"], r["t"]["S"], "s", color="#d62728", ms=9)
            ax.annotate(f"{r['name']}\nเตือนก่อนคลื่น S {r['lead'][g0['name']]['S']:+.0f} s",
                        (r["delta"], r["t"]["S"]), textcoords="offset points", xytext=(6, 4), fontsize=9)
    ax.set_xlabel("ระยะทางจากจุดเหนือศูนย์กลาง (km)")
    ax.set_ylabel("เวลาหลังเกิดแผ่นดินไหว (s)")
    ax.set_title("เวลาคลื่นเดินทาง และเวลาที่ระบบรู้ข้อมูลแต่ละขั้น  (แถบเขียว = เวลาเตือนล่วงหน้าจากขั้น 1)")
    ax.set_xlim(0, dmax)
    ax.set_ylim(0, None)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def make_map(res, stages, rows):
    try:
        from map_plot import make_travel_map
    except ImportError:
        print("(map_plot.py ไม่มี make_travel_map จึงไม่ทำแผนที่)")
        return None
    hypo, T0, h = res["hypo"], res["hypo"]["origin_time"], res["hypo"]["depth"]
    g0 = stages[0]
    blind = blind_zone_km(g0["time"] - T0, h)
    sites = [{"name": r["name"], "lat": r["lat"], "lon": r["lon"], "lead": r["lead"][g0["name"]]["S"],
              "info": f"Δ {r['delta']:.0f} km<br>P มาถึง {fmt(T0 + r['t']['P'])}<br>S มาถึง {fmt(T0 + r['t']['S'])}"
                      f"<br>Lg มาถึง {fmt(T0 + r['t']['Lg'])}"
                      + "".join(f"<br>{g['name']}: เตือนก่อนคลื่น S {r['lead'][g['name']]['S']:+.1f} s"
                                for g in stages)} for r in rows]
    stations = [{"code": c, "lat": la, "lon": lo, "tp": t}
                for c, ph, t, la, lo in (hypo.get("picks_used") or []) if ph == "P"]
    t_end = max([r["t"]["S"] for r in rows] + [g0["time"] - T0]) + 30
    fronts = [(t, blind_zone_km(t, h)) for t in np.arange(30, t_end, 30)]
    return make_travel_map((hypo["lat"], hypo["lon"], h), sites, "travel_time_map.html", blind_km=blind,
                           stations=stations, s_fronts=fronts)


def main():
    res = analyse()
    if res is None or res.get("hypo") is None:
        print("ยังหาตำแหน่งไม่ได้ (ต้องมีอย่างน้อย 3 สถานี) จึงคำนวณเวลาเดินทางไม่ได้  "
              "ใส่ EPICENTER หรือเพิ่มสถานี")
        return None
    stages = system_stages(res)
    rows = site_rows(res["hypo"], stages)
    report(res, stages, rows)
    if MAKE_PLOT:
        out = plot(res, stages, rows)
        if out:
            print(f"\nบันทึกกราฟ: {out}")
    if MAKE_MAP and stages:
        html = make_map(res, stages, rows)
        if html:
            print(f"บันทึกแผนที่จริง: {html}  (เปิดด้วย browser หรือพิมพ์ start {html})")
    return {"result": res, "stages": stages, "sites": rows}


if __name__ == "__main__":
    main()