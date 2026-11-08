"""
run_s_wave.py
หาเวลาคลื่น P แล้ว pick คลื่น S และ Lg รอบเวลาที่คาดจากตำแหน่งคร่าว ๆ (WINDOW_HYPO)
จากนั้นคำนวณเวลา S-P และระยะทางจากสถานีถึงจุดเกิดแผ่นดินไหว

รัน:  แก้ DATA_FILE ด้านล่าง แล้วกด Run ใน VS Code
  หรือ python run_s_wave.py <ไฟล์ข้อมูล.txt> [local|regional]
"""
import sys

import numpy as np

import config as C
from data_loader import load_station_txt, resolve_data_path
from p_wave import detect_p_wave
from epicenter import predict_arrival
from run_p_wave import to_utc
from s_wave import pick_s_window
from signal_utils import azimuth_deg

# ======================================================================
#  ใส่ชื่อไฟล์ข้อมูลคลื่นแผ่นดินไหวตรงนี้ (เหมือนใน run_p_wave.py)
# ======================================================================
DATA_FILE = "BKKA_20250328_061959.txt"
PROFILE = "regional"      # "regional" = สถานีไกล (> ~300 km), "local" = สถานีใกล้
# ตำแหน่งคร่าว ๆ ไว้ทำนายเวลาคลื่น S / Lg: (lat, lon, ความลึก km, "เวลาเกิด UTC")  ไม่ใช้ใส่ None
WINDOW_HYPO = (21.682, 96.121, 10.0, "2025-03-28T06:20:52")
WINDOW_SEC = 10.0         # หน้าต่าง ± (s) รอบเวลาที่คาด


def _pick_s_with_hypo(rec, p):
    """pick คลื่น S และ Lg ในช่วงรอบเวลาที่คาดจาก WINDOW_HYPO แล้วเลือกตัวที่คุณภาพดีกว่า"""
    from datetime import datetime, timezone
    from station_db import get_coords
    c = get_coords(rec["station"])
    if c is None:
        print(f"ไม่พบพิกัดของ {rec['station']} ใน stations.csv จึงใช้ WINDOW_HYPO ไม่ได้"); return None
    la, lo, dep, t_ref = WINDOW_HYPO
    hypo = {"lat": la, "lon": lo, "depth": dep,
            "origin_time": datetime.fromisoformat(t_ref).replace(tzinfo=timezone.utc).timestamp()}
    baz = float(azimuth_deg(c[0], c[1], la, lo))
    print(f"\nPick คลื่น S รอบเวลาที่คาดจาก ({la}, {lo}) ± {WINDOW_SEC:.0f} s   (ทิศจริง {baz:.0f}°)")
    best = None
    for ph in ("S", "Lg"):
        t_pred = predict_arrival(hypo, c[0], c[1], ph)
        k = pick_s_window(rec["z"], rec["n"], rec["e"], rec["fs"], rec["t0"], t_pred, baz, WINDOW_SEC,
                          phase=ph, max_weight=4, not_before=p["time"] + 1.0)
        if k is None:
            print(f"  {ph:<3}: ที่คาด {to_utc(t_pred)} UTC   -> ไม่พบ")
            continue
        print(f"  {ph:<3}: ที่คาด {to_utc(t_pred)} UTC   pick {to_utc(k['time'])} UTC "
              f"({k['time'] - t_pred:+.1f} s)  {k['label']}  พลังงานเพิ่ม x{k['ratio']:.1f}"
              + ("   <- weight 4 ไม่ใช้" if k["weight"] > 3 else ""))
        if k["weight"] <= 3 and (best is None or k["weight"] < best["weight"]):
            best = k
    if best is None:
        print("ไม่พบคลื่น S / Lg ที่ชัดพอในช่วงที่คาด"); return None
    best["sp_time"] = best["time"] - p["time"]
    print(f"\nใช้ {best['phase']} : {to_utc(best['time'])} UTC   (t = {best['index'] / rec['fs']:.2f} s)  {best['label']}")
    return best


def main(path, profile="regional", plot=True):
    rec = load_station_txt(path)
    fs = rec["fs"]
    z, n, e = rec["z"], rec["n"], rec["e"]
    print(f"Station      : {rec['station']}")
    print(f"Start time   : {to_utc(rec['t0'])} UTC")
    print(f"Profile      : {profile}")

    # ---------- 1) คลื่น P ----------
    p_params = C.P_REGIONAL if profile == "regional" else {}
    p = detect_p_wave(z, fs, start_time=rec["t0"], **p_params)
    if not p["detected"]:
        print("ไม่พบคลื่น P"); return
    print(f"\nP-wave arrival : {to_utc(p['time'])} UTC   (t = {p['index'] / fs:.2f} s)   {p.get('label', '')}")

    # ---------- 2) คลื่น S ----------
    s = _pick_s_with_hypo(rec, p)
    if s is None:
        return

    # ---------- 3) S-P และระยะทาง (โลก 2 ชั้น: Pg/Pn, Sg/Sn หรือ Lg) ----------
    import pipeline as PL
    sp = s["sp_time"]
    ph = s.get("phase", "S")
    print(f"\n{ph} - P time    : {sp:.2f} s")
    R = PL.delta_from_differential(sp, ph, 10.0)
    k = R / sp if sp > 0 else float("nan")
    print(f"ระยะทางถึงจุดเกิดแผ่นดินไหว : {R:.0f} km   (โลก 2 ชั้น, เทียบเท่า k = {k:.2f} km/s)")
    from epicenter import phase_travel_time
    t0 = p["time"] - float(phase_travel_time(R, 10.0, "P"))
    print(f"เวลาเกิดแผ่นดินไหวโดยประมาณ  : {to_utc(t0)} UTC   (T0 = Tp - tP(Δ))")

    if plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("(ไม่ได้ติดตั้ง matplotlib จึงไม่วาดกราฟ)"); return
        t = np.arange(len(z)) / fs
        tp, ts = p["index"] / fs, s["index"] / fs
        fig, ax = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
        comps, names = (z, n, e), ("Z", "N", "E")
        for a, x, name in zip(ax, comps, names):
            a.plot(t, x, lw=0.4, color="0.3")
            a.axvline(tp, color="r", label="P")
            a.axvline(ts, color="b", label="S")
            a.set_ylabel(f"{name} (gal)")
        ax[0].legend(loc="upper right")
        ax[0].set_title(f"{rec['station']}   P = {to_utc(p['time'])}   S = {to_utc(s['time'])}"
                        f"   S-P = {sp:.1f} s   R ~ {R:.0f} km   {s.get('label', '')}")
        ax[2].set_xlabel("Time from start (s)")
        fig.tight_layout()
        out = f"s_wave_{rec['station']}.png"
        fig.savefig(out, dpi=120)
        print(f"\nบันทึกกราฟ: {out}")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else DATA_FILE
    path = resolve_data_path(path)       # หาในโฟลเดอร์ station หรือโฟลเดอร์โค้ด
    profile = sys.argv[2] if len(sys.argv) > 2 else PROFILE
    main(path, profile)
