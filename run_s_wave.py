"""
run_s_wave.py
หาเวลาคลื่น P และคลื่น S จากไฟล์ข้อมูลจริงของสถานี
แล้วคำนวณเวลา S-P และระยะทางจากสถานีถึงจุดเกิดแผ่นดินไหว (แบบจำลองโลก 2 ชั้น)

สถานีเดียวที่ยังไม่รู้ตำแหน่ง: pick คลื่น S แบบ "blind" (S_METHOD) ซึ่งกับแผ่นดินไหวใหญ่ที่รอยเลื่อนแตกยาว
มักจับคลื่น P ระลอกหลังเป็นคลื่น S ถ้ารู้ตำแหน่งคร่าว ๆ (เช่นจาก run_network.py หรือจุดอ้างอิง)
ให้ใส่ WINDOW_HYPO แล้วโปรแกรมจะ pick คลื่น S และ Lg เฉพาะในช่วงรอบเวลาที่ควรมาถึง (แม่นกว่ามาก)

รัน:  แก้ DATA_FILE ด้านล่าง แล้วกด Run ใน VS Code
  หรือ python run_s_wave.py <ไฟล์ข้อมูล.txt> [local|regional]
"""
import os
import sys

import numpy as np

import config as C
from data_loader import load_station_txt, resolve_data_path
from p_wave import detect_p_wave
from s_wave import detect_s_wave
from epicenter import omori_distance, predict_arrival
from run_p_wave import to_utc
from s_wave import pick_s_window
from signal_utils import azimuth_deg

# ======================================================================
#  ใส่ชื่อไฟล์ข้อมูลคลื่นแผ่นดินไหวตรงนี้ (เหมือนใน run_p_wave.py)
# ======================================================================
DATA_FILE = "BKKA_20250328_061959.txt"
PROFILE = "regional"      # "regional" = สถานีไกล (> ~300 km), "local" = สถานีใกล้
S_METHOD = "usgs"         # "usgs" (R/T + polarization), "jump", "first" หรือ None = ใช้ค่าใน config
ORIGIN_TIME = None        # เวลาเกิดแผ่นดินไหว (epoch s) ถ้ารู้ ช่วยจำกัดหน้าต่างค้นหา S (Vp/Vs = 1.73)
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
    if WINDOW_HYPO is not None:
        s = _pick_s_with_hypo(rec, p)
        if s is None:
            return
    else:
        s_params = dict(C.S_REGIONAL) if profile == "regional" else {}
        if S_METHOD:
            s_params["method"] = S_METHOD
        s = detect_s_wave(z, n, e, fs, p["index"], start_time=rec["t0"], origin_time=ORIGIN_TIME, **s_params)
        if not s["detected"]:
            print(f"ไม่พบคลื่น S (STA/LTA สูงสุด = {s.get('max_ratio', 0):.2f})"); return
        print(f"S-wave arrival : {to_utc(s['time'])} UTC   (t = {s['index'] / fs:.2f} s)")
    if WINDOW_HYPO is None and s.get("method") == "usgs":
        print(f"  Back-azimuth : {s['baz']:.0f}°  (ความเชื่อมั่น {s['baz_conf']:.2f}, จากการสั่นของคลื่น P)")
        print(f"  Onset        : {'i = impulsive (ชัด)' if s['onset'] == 'i' else 'e = emergent (ค่อย ๆ โผล่)'}")
        print(f"  SNR          : {s['snr']:.1f}  (เทียบกับ P coda ก่อนหน้า)")
        print(f"  Weight       : {s['weight']}  (0 = ดีที่สุด ... 4 = ไม่ใช้), ±{s['uncertainty']:.2f} s")
        print(f"  Pick label   : {s['label']}")
        if s.get("predicted_time"):
            print(f"  S ที่ทำนาย    : {to_utc(s['predicted_time'])} UTC  (Vp/Vs = 1.73)")
    elif s.get("method") == "jump":
        print(f"  วิธี jump    : score = {s['score']:.2f}  (พลังงานแนวราบกระโดด + H/V เพิ่ม)")
        print("  ตัวเลือกคลื่น S ที่พบ (เลือกตัวแรกที่ score >= 70% ของสูงสุด):")
        for c in s["candidates"]:
            print(f"     S-P = {c['sp_time']:7.1f} s   score = {c['score']:5.2f}   H/V = {c['hv']:5.1f}"
                  f"{'   <-- เลือก' if c['chosen'] else ''}")
    elif WINDOW_HYPO is None:
        print(f"  S/LTA max    : {s['max_ratio']:.1f}")

    # ---------- 3) S-P และระยะทาง (โลก 2 ชั้น: Pg/Pn, Sg/Sn หรือ Lg) ----------
    import pipeline as PL
    sp = s["sp_time"]
    ph = s.get("phase", "S")
    print(f"\n{ph} - P time    : {sp:.2f} s")
    R = PL.delta_from_differential(sp, ph, 10.0)
    k = R / sp if sp > 0 else float("nan")
    print(f"ระยะทางถึงจุดเกิดแผ่นดินไหว : {R:.0f} km   (โลก 2 ชั้น, เทียบเท่า k = {k:.2f} km/s)")
    if profile == "regional" and ph == "S":
        print(f"  (สูตร Omori แบบเดิม k = 10.82 km/s จะได้ {omori_distance(sp, C.VPN_KM_S, C.VSN_KM_S):.0f} km)")
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
        if s.get("method") == "usgs":
            comps, names = (z, s["r"], s["t"]), ("Z", "R", "T")
        else:
            comps, names = (z, n, e), ("Z", "N", "E")
        for a, x, name in zip(ax, comps, names):
            a.plot(t, x, lw=0.4, color="0.3")
            a.axvline(tp, color="r", label="P")
            a.axvline(ts, color="b", label="S")
            for j, c in enumerate(s.get("candidates", [])):     # ตัวเลือกอื่นของวิธี jump (เส้นประเทา)
                if not c["chosen"]:
                    a.axvline(tp + c["sp_time"], color="0.5", ls=":", lw=1,
                              label="S candidates" if j == 0 or (j == 1 and s["candidates"][0]["chosen"]) else None)
            a.set_ylabel(f"{name} (gal)" if name in ("Z", "N", "E") else f"{name} (gal, filtered)")
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
