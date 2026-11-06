"""
run_p_wave.py
หาเวลาคลื่น P จากไฟล์ข้อมูลจริงของสถานี

รัน:  แก้ DATA_FILE ด้านล่าง แล้วกด Run ใน VS Code
  หรือ python run_p_wave.py <ไฟล์ข้อมูล.txt> [local|regional]
  local    = สถานีใกล้จุดศูนย์กลาง (ค่าใน config.py)
  regional = สถานีไกล > ~300 km (ค่า P_REGIONAL)  <- ค่าเริ่มต้น
"""
import os
import sys
from datetime import datetime, timezone

import numpy as np

import config as C
from data_loader import load_station_txt
from p_wave import detect_p_wave

# ======================================================================
#  ใส่ชื่อไฟล์ข้อมูลคลื่นแผ่นดินไหวตรงนี้
#  - ถ้าไฟล์อยู่ในโฟลเดอร์เดียวกับโค้ด ใส่แค่ชื่อไฟล์
#  - ถ้าอยู่ที่อื่น ใส่ที่อยู่เต็ม เช่น r"C:\Users\me\Downloads\BKKA.txt"
# ======================================================================
DATA_FILE = "BKKA_20250328_061959.txt"
PROFILE = "regional"      # "regional" = สถานีไกล (> ~300 km), "local" = สถานีใกล้
DATA_DIR = "station"      # โฟลเดอร์ที่เก็บไฟล์ข้อมูล (อยู่ในโฟลเดอร์ code)


def to_utc(t):
    return datetime.fromtimestamp(t, tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]


def main(path, profile="regional", plot=True):
    rec = load_station_txt(path)
    fs, z = rec["fs"], rec["z"]
    print(f"Station      : {rec['station']}")
    print(f"Start time   : {to_utc(rec['t0'])} UTC")
    print(f"Sampling rate: {fs} Hz, {len(z)} samples ({len(z) / fs:.1f} s)")

    params = C.P_REGIONAL if profile == "regional" else {}
    print(f"Profile      : {profile} {params}")
    p = detect_p_wave(z, fs, start_time=rec["t0"], **params)
    if not p["detected"]:
        print(f"ไม่พบคลื่น P (STA/LTA สูงสุด = {p['max_ratio']:.2f})")
        return

    i = p["index"]
    print(f"\nP-wave arrival : {to_utc(p['time'])} UTC")
    print(f"  sample index : {i}  (t = {i / fs:.2f} s จากต้นไฟล์)")
    print(f"  STA/LTA trig : {to_utc(rec['t0'] + p['trigger_index'] / fs)} UTC")
    print(f"  STA/LTA max  : {p['max_ratio']:.1f}")
    fm = {1: "U = ขึ้น (Compression)", -1: "D = ลง (Dilatation)"}.get(p["polarity"], "? = ไม่ชัด")
    onset = "i = impulsive (ชัด)" if p["onset"] == "i" else "e = emergent (ค่อย ๆ โผล่)"
    print(f"  Onset        : {onset}")
    print(f"  First motion : {fm}")
    print(f"  SNR          : {p['snr']:.1f}")
    print(f"  Weight       : {p['weight']}  (0 = ดีที่สุด ... 4 = ไม่ใช้), ±{p['uncertainty']:.2f} s")
    print(f"  Pick label   : {p['label']}")

    if plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("(ไม่ได้ติดตั้ง matplotlib จึงไม่วาดกราฟ)")
            return
        t = np.arange(len(z)) / fs
        fig, ax = plt.subplots(3, 1, figsize=(12, 8))
        ax[0].plot(t, z, lw=0.5)
        ax[0].axvline(i / fs, color="r", label=f"P pick ({p['label']})")
        ax[0].set_ylabel("Z (gal)"); ax[0].legend(); ax[0].set_title(f"{rec['station']} - full record")
        w = (t > i / fs - 20) & (t < i / fs + 20)
        ax[1].plot(t[w], z[w], lw=0.6)
        ax[1].axvline(i / fs, color="r"); ax[1].set_ylabel("Z (gal)"); ax[1].set_title("zoom ±20 s")
        ax[2].plot(t, p["ratio"], lw=0.6)
        ax[2].axhline(params.get("threshold", C.P_TRIGGER_ON), color="k", ls="--"); ax[2].axvline(i / fs, color="r")
        ax[2].set_ylabel("STA/LTA"); ax[2].set_xlabel("Time from start (s)")
        fig.tight_layout()
        out = f"p_wave_{rec['station']}.png"
        fig.savefig(out, dpi=120)
        print(f"\nบันทึกกราฟ: {out}")


if __name__ == "__main__":
    # ถ้าพิมพ์ชื่อไฟล์ตามหลังคำสั่งใน Terminal จะใช้ไฟล์นั้นแทน DATA_FILE
    path = sys.argv[1] if len(sys.argv) > 1 else DATA_FILE
    if not os.path.isabs(path) and not os.path.exists(path):
        code_dir = os.path.dirname(os.path.abspath(__file__))
        for folder in (os.path.join(code_dir, DATA_DIR), code_dir):   # หาใน station ก่อน แล้วค่อยหาในโฟลเดอร์โค้ด
            candidate = os.path.join(folder, path)
            if os.path.exists(candidate):
                path = candidate
                break
        else:
            sys.exit(f"หาไฟล์ '{path}' ไม่เจอ ทั้งใน {DATA_DIR}\\ และในโฟลเดอร์โค้ด")
    profile = sys.argv[2] if len(sys.argv) > 2 else PROFILE
    main(path, profile)
