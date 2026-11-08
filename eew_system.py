"""
eew_system.py
รวมทุกโมดูลเป็นระบบ EEW สำหรับ demo ข้อมูลสังเคราะห์ (ลำดับเดียวกับระบบหลัก):
  คลื่น P -> ตำแหน่งจากคลื่น P -> pick คลื่น S รอบเวลาที่คาด -> ตำแหน่งจาก P + S
  -> ขนาด Mw_P -> เวลาเดินทาง -> ทำนาย MMI -> เกณฑ์แจ้งเตือน
"""
import config as C
from p_wave import detect_p_wave
from s_wave import pick_s_window
from epicenter import geiger_locate_2layer, predict_arrival
from magnitude import mw_p_magnitude, network_magnitude
from signal_utils import azimuth_deg
from travel_time import site_warning, hypocentral_distance
from intensity import predict_mmi, alert_level, mmi_class


class EEWSystem:
    def __init__(self, stations, site):
        """
        stations : list ของ dict {code, lat, lon, record:{t0, fs, z, n, e}}
        site     : dict {name, lat, lon, vs30}  ตำแหน่งของเรา (ผู้รับคำเตือน)
        """
        self.stations = stations
        self.site = site

    # ---------- ขั้นที่ 1: ตรวจจับคลื่น P ----------
    def pick_p(self):
        picks = []
        for s in self.stations:
            r = s["record"]
            p = detect_p_wave(r["z"], r["fs"], r["t0"])
            if p["detected"]:
                picks.append({**s, "p": p, "tp": p["time"], "wp": p.get("weight", 0), "ts": None})
        return sorted(picks, key=lambda x: x["tp"])

    # ---------- ขั้นที่ 2-3: ตำแหน่งจากคลื่น P -> pick คลื่น S -> ตำแหน่งจาก P + S ----------
    def locate(self, picks):
        if len(picks) < 3:
            raise RuntimeError("ข้อมูลไม่พอ: ต้องมีคลื่น P อย่างน้อย 3 สถานี")
        hypo = geiger_locate_2layer(picks, phases=("P",))
        for p in picks:
            r = p["record"]
            baz = float(azimuth_deg(p["lat"], p["lon"], hypo["lat"], hypo["lon"]))
            k = pick_s_window(r["z"], r["n"], r["e"], r["fs"], r["t0"],
                              predict_arrival(hypo, p["lat"], p["lon"], "S"), baz, phase="S",
                              max_weight=3, not_before=p["tp"] + 1.0)
            if k is not None:
                p.update(ts=k["time"], ws=k["weight"])
        return geiger_locate_2layer(picks, phases=("P", "S"))

    # ---------- ขั้นที่ 4: ขนาดแผ่นดินไหว ----------
    def estimate_magnitude(self, picks, hypo):
        vals = []
        for p in picks:
            r = p["record"]
            fs, ip = r["fs"], p["p"]["index"]
            n_pre = max(ip - int(1 * fs), 10)
            R, delta = hypocentral_distance(hypo["lat"], hypo["lon"], hypo["depth"], p["lat"], p["lon"])
            m = mw_p_magnitude(r["z"], fs, ip, delta, hypo["depth"], n_pre=n_pre)
            p.update(mw_p=m, epi_dist=delta, hypo_dist=R)
            vals.append(m["mag"])
        return network_magnitude(vals)["mag"]

    # ---------- ขั้นที่ 5-6: เวลาเดินทาง + ความรุนแรงที่ตำแหน่งเรา ----------
    def run(self, processing_delay=1.0, verbose=True):
        picks = self.pick_p()
        if not picks:
            print("ไม่พบแผ่นดินไหว"); return None
        hypo = self.locate(picks)
        M = self.estimate_magnitude(picks, hypo)

        # เวลาออกคำเตือน = สถานีที่ 3 ได้ข้อมูล P ครบ 3 วินาที + เวลาประมวลผล
        alert_time = picks[2]["tp"] + C.PD_TAUC_WINDOW_SEC + processing_delay
        tt = site_warning(hypo, self.site["lat"], self.site["lon"], alert_time)
        pred = predict_mmi(M, tt["epi_dist_km"], hypo["depth"], self.site.get("vs30", C.DEFAULT_VS30))
        _, level = alert_level(M, pred["mmi"])

        result = {"picks": picks, "hypocenter": hypo, "M": M, "alert_time": alert_time,
                  "site": tt, "prediction": pred, "level": level}
        if verbose:
            self.report(result)
        return result

    def report(self, res):
        h, s, pr = res["hypocenter"], res["site"], res["prediction"]
        line = "=" * 78
        print(line); print(" EARTHQUAKE EARLY WARNING REPORT"); print(line)
        print(f"{'Sta':<6}{'Tp(s)':>9}{'Ts(s)':>9}{'S-P':>6}{'Δ(km)':>8}{'Pd(cm)':>10}{'τc':>6}{'Mw_P':>6}")
        for p in res["picks"]:
            ts, m = p["ts"], p["mw_p"]
            print(f"{p['code']:<6}{p['tp']:9.2f}{(ts if ts else float('nan')):9.2f}"
                  f"{((ts - p['tp']) if ts else float('nan')):6.2f}{p['epi_dist']:8.1f}"
                  f"{m['pd_cm']:10.2e}{m['tau_c']:6.2f}{m['mag']:6.2f}")
        print(line)
        print(f" Method      : {h['method']}")
        print(f" Epicenter   : lat {h['lat']:.3f}, lon {h['lon']:.3f}, depth {h['depth']:.1f} km")
        print(f" Origin time : {h['origin_time']:.2f} s   (RMS {h['rms']:.3f} s)")
        print(f" Magnitude   : Mw_P {res['M']:.1f}")
        print(line)
        print(f" ตำแหน่งของเรา: {self.site['name']} ({self.site['lat']}, {self.site['lon']})")
        print(f"   ระยะทางจากจุดเหนือศูนย์ : {s['epi_dist_km']:.1f} km  (hypocentral {s['hypo_dist_km']:.1f} km)")
        print(f"   คลื่น P เดินทาง          : {s['p_travel_s']:.1f} s")
        print(f"   คลื่น S เดินทาง          : {s['s_travel_s']:.1f} s")
        print(f"   ออกคำเตือนหลังเกิดเหตุ   : {res['alert_time'] - h['origin_time']:.1f} s")
        print(f"   เวลาเตือนล่วงหน้า (S)    : {s['lead_time_s']:.1f} s")
        print(f"   รัศมี Blind zone         : {s['blind_zone_km']:.1f} km")
        _, rom, shake, _ = mmi_class(pr["mmi"])
        print(f"   ความรุนแรงทำนาย          : MMI {pr['mmi']:.1f} ({rom} {shake})")
        print(f"   ระดับการแจ้งเตือน        : {res['level']}")
        print(line)
