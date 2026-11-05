"""
run_s_wave.py
หาเวลาคลื่น P และคลื่น S จากไฟล์ข้อมูลจริงของสถานี
แล้วคำนวณเวลา S-P และระยะทางจากสถานีถึงจุดเกิดแผ่นดินไหว (สูตร Omori)

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
from epicenter import omori_distance
from run_p_wave import to_utc

# ======================================================================
#  ใส่ชื่อไฟล์ข้อมูลคลื่นแผ่นดินไหวตรงนี้ (เหมือนใน run_p_wave.py)
# ======================================================================
DATA_FILE = "BKKA_20250328_061959.txt"
PROFILE = "regional"      # "regional" = สถานีไกล (> ~300 km), "local" = สถานีใกล้


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
    print(f"\nP-wave arrival : {to_utc(p['time'])} UTC   (t = {p['index'] / fs:.2f} s)")

    # ---------- 2) คลื่น S ----------
    s_params = C.S_REGIONAL if profile == "regional" else {}
    s = detect_s_wave(z, n, e, fs, p["index"], start_time=rec["t0"], **s_params)
    if not s["detected"]:
        print(f"ไม่พบคลื่น S (STA/LTA สูงสุด = {s.get('max_ratio', 0):.2f})"); return
    print(f"S-wave arrival : {to_utc(s['time'])} UTC   (t = {s['index'] / fs:.2f} s)")
    if s.get("method") == "jump":
        print(f"  วิธี jump    : score = {s['score']:.2f}  (พลังงานแนวราบกระโดด + H/V เพิ่ม)")
        print("  ตัวเลือกคลื่น S ที่พบ (เลือกตัวแรกที่ score >= 70% ของสูงสุด):")
        for c in s["candidates"]:
            print(f"     S-P = {c['sp_time']:7.1f} s   score = {c['score']:5.2f}   H/V = {c['hv']:5.1f}"
                  f"{'   <-- เลือก' if c['chosen'] else ''}")
    else:
        print(f"  S/LTA max    : {s['max_ratio']:.1f}")

    # ---------- 3) S-P และระยะทาง ----------
    sp = s["sp_time"]
    print(f"\nS - P time     : {sp:.2f} s")
    if profile == "regional":
        vp, vs, label = C.VPN_KM_S, C.VSN_KM_S, "Pn/Sn"
    else:
        vp, vs, label = C.VP_KM_S, C.VS_KM_S, "crust"
    k = vp * vs / (vp - vs)
    R = omori_distance(sp, vp, vs)
    print(f"ระยะทางถึงจุดเกิดแผ่นดินไหว : {R:.0f} km   (Omori: k = {k:.2f} km/s, Vp={vp}, Vs={vs} [{label}])")
    t0 = p["time"] - R / vp
    print(f"เวลาเกิดแผ่นดินไหวโดยประมาณ  : {to_utc(t0)} UTC   (T0 = Tp - R/Vp)")

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
        for a, x, name in zip(ax, (z, n, e), ("Z", "N", "E")):
            a.plot(t, x, lw=0.4, color="0.3")
            a.axvline(tp, color="r", label="P")
            a.axvline(ts, color="b", label="S")
            for j, c in enumerate(s.get("candidates", [])):     # ตัวเลือกอื่นของวิธี jump (เส้นประเทา)
                if not c["chosen"]:
                    a.axvline(tp + c["sp_time"], color="0.5", ls=":", lw=1,
                              label="S candidates" if j == 0 or (j == 1 and s["candidates"][0]["chosen"]) else None)
            a.set_ylabel(f"{name} (gal)")
        ax[0].legend(loc="upper right")
        ax[0].set_title(f"{rec['station']}   P = {to_utc(p['time'])}   S = {to_utc(s['time'])}"
                        f"   S-P = {sp:.1f} s   R ~ {R:.0f} km")
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
