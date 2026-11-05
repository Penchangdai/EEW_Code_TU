"""
main.py
ตัวอย่างการใช้งาน: จำลองแผ่นดินไหว M6.3 ที่ อ.แม่ลาว จ.เชียงราย (คล้ายเหตุการณ์ 5 พ.ค. 2557)
แล้วให้ระบบ EEW ตรวจจับ คำนวณ และเตือนภัยไปยัง มธ. ศูนย์รังสิต

รัน:  python main.py
"""
from synthetic import make_station_record
from eew_system import EEWSystem

# เหตุการณ์จำลอง (ค่า "จริง" ที่ระบบไม่รู้)
EQ = {"lat": 19.68, "lon": 99.68, "depth": 7.0, "origin_time": 0.0, "M": 6.3}

# สถานีตรวจวัด (พิกัดโดยประมาณ)
STATIONS = [
    ("CRAI", 19.91, 99.83),   # เชียงราย
    ("PHAY", 19.17, 99.90),   # พะเยา
    ("CMAI", 18.79, 98.98),   # เชียงใหม่
    ("LAMP", 18.29, 99.49),   # ลำปาง
    ("NAN", 18.78, 100.78),   # น่าน
]

SITE = {"name": "Thammasat Rangsit", "lat": 14.07, "lon": 100.61, "vs30": 250.0}


def main():
    stations = []
    for i, (code, lat, lon) in enumerate(STATIONS):
        rec = make_station_record(lat, lon, EQ, seed=i)
        stations.append({"code": code, "lat": lat, "lon": lon, "record": rec})

    eew = EEWSystem(stations, SITE)
    res = eew.run()

    print("\nตรวจสอบกับค่าจริง (synthetic):")
    for p in res["picks"]:
        r = p["record"]
        ts_err = (p["ts"] - r["ts_true"]) if p["ts"] else float("nan")
        print(f"  {p['code']:<5} P error {p['tp'] - r['tp_true']:+.2f} s   S error {ts_err:+.2f} s")
    h = res["hypocenter"]
    print(f"  True epicenter {EQ['lat']}, {EQ['lon']}, depth {EQ['depth']} km, M {EQ['M']}")
    print(f"  Est. epicenter {h['lat']:.3f}, {h['lon']:.3f}, depth {h['depth']:.1f} km, M {res['M']:.1f}")


if __name__ == "__main__":
    main()
