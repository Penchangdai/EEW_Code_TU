"""
run_intensity.py
คำนวณความรุนแรงแผ่นดินไหว (Seismic Intensity) มาตรา JMA 0-7 และ MMI

  ส่วนที่ 1  ความรุนแรงที่วัดได้จริงที่สถานี (Instrumental Intensity ของ JMA)
             FFT -> filter F1·F2·F3 -> รวม 3 แกน -> หา a0 ที่เกินรวม 0.3 s -> I = 2 log10(a0) + 0.94
             พร้อม PGA, PGV และ MMI (Wald et al., 1999)
  ส่วนที่ 2  ความรุนแรงที่ทำนายได้ (วิธีเดียวกับ JMA EEW) เมื่อใส่ MAGNITUDE และ EPICENTER
             PGV บนชั้นหิน (Si & Midorikawa, 1999) x อัตราขยายของดิน (Vs30) -> I = 2.68 + 1.72 log10(PGV)
             ทำนายที่สถานีที่เลือก (เทียบกับค่าที่วัดได้) และที่ตำแหน่งใน TARGET_SITES
  ผลลัพธ์    ตาราง, กราฟ intensity_<สถานี>.png, แผนที่ intensity_map.html

รัน:  ใส่รหัสสถานีใน SELECTED_STATIONS (และ MAGNITUDE, EPICENTER ถ้าต้องการทำนาย) แล้วกด Run
"""
import sys

import numpy as np

import config as C
from data_loader import load_station_txt, resolve_data_path, find_station_file
from intensity import jma_intensity_detail, jma_intensity_class, mmi_from_pgv, predict_intensity
from signal_utils import haversine_km, acc_to_vel, remove_offset
from station_db import get_coords, get_vs30
from map_plot import make_intensity_map

# ======================================================================
#  ตั้งค่าตรงนี้
# ======================================================================
# รหัสสถานีที่ต้องการคำนวณความรุนแรงที่วัดได้ (1-5 สถานี, ต้องมีไฟล์ในโฟลเดอร์ station)
SELECTED_STATIONS = ["BKKA", "CHUM", "PAYA"]
# ถ้า SELECTED_STATIONS = [] จะใช้ไฟล์นี้ไฟล์เดียวแทน
DATA_FILE = "BKKA_20250328_061959.txt"
MAX_STATIONS = 5

# ---- สำหรับทำนายความรุนแรง (ใส่ None ทั้งคู่ = คำนวณเฉพาะค่าที่วัดได้) ----
MAGNITUDE = 7.7              # ขนาด เช่น 7.7 หรือผลจาก run_magnitude.py
EPICENTER = (22.01, 95.92, 10.0)              # (lat, lon, ความลึก km) เช่น ผลจาก run_network.py
EQ_TYPE = "crustal"           # "crustal" (ในเปลือกโลก), "interplate" (ระหว่างแผ่น), "intraplate" (ภายในแผ่น)

# Vs30 (m/s) ของสถานี  ลำดับที่ใช้: STATION_VS30 -> คอลัมน์ vs30 ใน stations.csv (เติมด้วย vs30_tool.py)
#                      -> DEFAULT_VS30 ใน config.py (400 m/s)
STATION_VS30 = {}             # เช่น {"BKKA": 180}

# ตำแหน่งที่ต้องการทำนายความรุนแรง (พิกัดตัวอย่างเป็นค่าโดยประมาณ ควรตรวจก่อนใช้)
# Vs30 ดินอ่อนกรุงเทพฯ/ปริมณฑล ประมาณ 150-250 m/s
TARGET_SITES = [
    {"name": "มธ. รังสิต", "lat": 14.07, "lon": 100.61, "vs30": 250},
]

MAKE_MAP = True


def analyze(path, plot=True):
    rec = load_station_txt(path)
    code, fs = rec["station"], rec["fs"]
    z, n, e = rec["z"], rec["n"], rec["e"]
    n_pre = max(int(10 * fs), 10)                          # ใช้ 10 วินาทีแรกของไฟล์ลบ offset
    det = jma_intensity_detail(z, n, e, fs)
    za, na, ea = (remove_offset(x, n_pre) for x in (z, n, e))
    pga = float(np.max(np.sqrt(za ** 2 + na ** 2 + ea ** 2)))
    vz, vn, ve = (acc_to_vel(x, fs, C.DISP_HIGHPASS_HZ, n_pre) for x in (z, n, e))
    pgv = float(np.max(np.sqrt(vn ** 2 + ve ** 2)))       # PGV แนวราบ
    res = {"code": code, "I": det["I"], "cls": jma_intensity_class(det["I"]), "a0": det["a0"],
           "pga": pga, "pgv": pgv, "mmi": mmi_from_pgv(pgv), "coords": get_coords(code)}
    print(f"\n========== สถานี {code} ==========")
    print(f"PGA (3 แกนรวม)       : {pga:.2f} gal   ({pga / 981 * 100:.2f} %g)")
    print(f"PGV (แนวราบ)         : {pgv:.3f} cm/s")
    print(f"a0 (เกินรวม 0.3 s)   : {det['a0']:.3f} gal  (หลังกรองด้วย F1·F2·F3)")
    print(f"ความรุนแรง JMA       : I = {det['I']:.1f}  -> ระดับ {res['cls']}")
    print(f"MMI (จาก PGV)        : {res['mmi']:.1f}")
    if plot:
        plot_station(rec, det, res)
    return res


def plot_station(rec, det, res):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fs = rec["fs"]
    t = np.arange(len(rec["z"])) / fs
    raw = np.sqrt(sum(remove_offset(rec[k], int(10 * fs)) ** 2 for k in "zne"))
    fig, ax = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    ax[0].plot(t, raw, lw=0.5, color="0.3")
    ax[0].set_ylabel("|a| raw (gal)")
    ax[0].set_title(f"{res['code']}   PGA {res['pga']:.2f} gal | PGV {res['pgv']:.3f} cm/s | "
                    f"JMA I = {res['I']:.1f} (class {res['cls']}) | MMI {res['mmi']:.1f}")
    ax[1].plot(t, det["a"], lw=0.5, color="tab:blue")
    ax[1].axhline(det["a0"], color="r", ls="--", label=f"a0 = {det['a0']:.3f} gal (exceeded for 0.3 s)")
    ax[1].set_ylabel("|a| filtered (gal)"); ax[1].set_xlabel("Time from start (s)")
    ax[1].legend(loc="upper right")
    fig.tight_layout()
    out = f"intensity_{res['code']}.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"บันทึกกราฟ: {out}")


def station_vs30(code):
    if code in STATION_VS30:
        return STATION_VS30[code]
    v = get_vs30(code)
    return v if v else C.DEFAULT_VS30


def predict_at(name, lat, lon, vs30):
    delta = haversine_km(EPICENTER[0], EPICENTER[1], lat, lon)
    R = float(np.sqrt(delta ** 2 + EPICENTER[2] ** 2))
    p = predict_intensity(MAGNITUDE, R, EPICENTER[2], vs30, EQ_TYPE)
    return dict(p, name=name, lat=lat, lon=lon, delta=delta, R=R, vs30=vs30)


def main(stations=SELECTED_STATIONS, plot=True):
    if stations:
        paths = []
        for code in stations[:MAX_STATIONS]:
            path = find_station_file(code)
            if path is None:
                print(f"⚠ ไม่พบไฟล์ข้อมูลของสถานี {code} (ต้องมี '# Station: {code}' ในโฟลเดอร์ station)")
            else:
                paths.append(path)
    else:
        paths = [resolve_data_path(DATA_FILE)]

    # ---------- ส่วนที่ 1: ความรุนแรงที่วัดได้ ----------
    obs = [analyze(p, plot) for p in paths]
    if obs:
        print("\n================ ความรุนแรงที่วัดได้ที่สถานี ================")
        print(f"{'Station':<8}{'PGA(gal)':>10}{'PGV(cm/s)':>11}{'I_JMA':>7}{'ระดับ':>7}{'MMI':>6}")
        for r in obs:
            print(f"{r['code']:<8}{r['pga']:10.2f}{r['pgv']:11.3f}{r['I']:7.1f}{r['cls']:>7}{r['mmi']:6.1f}")

    # ---------- ส่วนที่ 2: ความรุนแรงที่ทำนายได้ ----------
    preds = []
    if MAGNITUDE is not None and EPICENTER is not None:
        print(f"\n================ ทำนายความรุนแรง (M {MAGNITUDE}, จุดศูนย์กลาง {EPICENTER}) ================")
        print(f"{'ตำแหน่ง':<14}{'Δ(km)':>8}{'Vs30':>6}{'ARV':>6}{'PGV(cm/s)':>11}{'I ทำนาย':>9}{'ระดับ':>7}"
              f"{'I วัดได้':>9}")
        for r in obs:
            if r["coords"] is None:
                continue
            p = predict_at(r["code"], *r["coords"], station_vs30(r["code"]))
            p["I_obs"] = r["I"]
            preds.append(p)
        for s in TARGET_SITES:
            preds.append(predict_at(s["name"], s["lat"], s["lon"], s.get("vs30", C.DEFAULT_VS30)))
        for p in preds:
            iobs = f"{p['I_obs']:9.1f}" if "I_obs" in p else f"{'-':>9}"
            print(f"{p['name']:<14}{p['delta']:8.0f}{p['vs30']:6.0f}{p['arv']:6.2f}{p['pgv_surface']:11.3f}"
                  f"{p['intensity']:9.1f}{p['jma_class']:>7}{iobs}")
        if any(p["delta"] > 300 for p in preds):
            print("⚠ สมการ Si & Midorikawa (1999) ปรับเทียบจากข้อมูลที่ระยะไม่เกินราว 300 km "
                  "ตำแหน่งที่ไกลกว่านี้ค่าทำนายอาจคลาดมาก")
        warn = [p["name"] for p in preds if p["intensity"] >= C.WARNING_INTENSITY]
        if warn:
            print(f"*** WARNING: คาดว่าความรุนแรง ≥ 5 ล่าง ที่ {', '.join(warn)} ***")
    else:
        print("\n(ใส่ MAGNITUDE และ EPICENTER เพื่อทำนายความรุนแรงที่สถานีและที่ TARGET_SITES)")

    # ---------- แผนที่ ----------
    if MAKE_MAP and (obs or preds):
        pts = [{"name": r["code"], "lat": r["coords"][0], "lon": r["coords"][1], "I": r["I"], "cls": r["cls"],
                "kind": "obs", "info": f"PGA {r['pga']:.2f} gal<br>PGV {r['pgv']:.3f} cm/s<br>MMI {r['mmi']:.1f}"}
               for r in obs if r["coords"] is not None]
        pts += [{"name": p["name"], "lat": p["lat"], "lon": p["lon"], "I": p["intensity"], "cls": p["jma_class"],
                 "kind": "pred", "info": f"Δ {p['delta']:.0f} km<br>PGV {p['pgv_surface']:.3f} cm/s"
                                         f"<br>Vs30 {p['vs30']:.0f} m/s"}
                for p in preds if "I_obs" not in p]
        if pts:
            html = make_intensity_map(pts, "intensity_map.html", epicenter=EPICENTER)
            if html:
                print(f"\nบันทึกแผนที่จริง: {html}  (เปิดด้วย browser หรือพิมพ์ start {html})")
    return obs, preds


if __name__ == "__main__":
    main(sys.argv[1:] or SELECTED_STATIONS)
