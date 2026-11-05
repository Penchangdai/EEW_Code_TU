"""
eew_system.py
รวมทุกโมดูลเป็นระบบ EEW:
  ข้อมูล 3 แกน -> P wave -> S wave -> Epicenter -> Magnitude -> Travel time -> Intensity -> ประกาศเตือน
"""
import numpy as np

import config as C
from p_wave import detect_p_wave
from s_wave import detect_s_wave
from epicenter import geiger_locate, single_station_location
from magnitude import (jma_magnitude, pd_value, magnitude_from_pd, tau_c,
                       magnitude_from_tauc, network_magnitude)
from travel_time import site_warning, hypocentral_distance
from intensity import predict_intensity, jma_instrumental_intensity, jma_intensity_class


class EEWSystem:
    def __init__(self, stations, site):
        """
        stations : list ของ dict {code, lat, lon, record:{t0, fs, z, n, e}}
        site     : dict {name, lat, lon, vs30}  ตำแหน่งของเรา (ผู้รับคำเตือน)
        """
        self.stations = stations
        self.site = site

    # ---------- ขั้นที่ 1-2: ตรวจจับคลื่น P และ S ----------
    def pick_phases(self):
        picks = []
        for s in self.stations:
            r = s["record"]
            p = detect_p_wave(r["z"], r["fs"], r["t0"])
            if not p["detected"]:
                continue
            sw = detect_s_wave(r["z"], r["n"], r["e"], r["fs"], p["index"], r["t0"])
            picks.append({**s, "p": p, "s": sw, "tp": p["time"],
                          "ts": sw["time"] if sw["detected"] else None})
        return sorted(picks, key=lambda x: x["tp"])

    # ---------- ขั้นที่ 3: หาจุดศูนย์กลาง ----------
    def locate(self, picks):
        n_obs = sum(1 + (p["ts"] is not None) for p in picks)
        if len(picks) >= 3 and n_obs >= 4:
            return geiger_locate(picks)
        first = picks[0]
        if first["ts"] is None:
            raise RuntimeError("ข้อมูลไม่พอ: มีสถานีเดียวและยังไม่พบคลื่น S")
        r = first["record"]
        loc = single_station_location(first["lat"], first["lon"], r["z"], r["n"], r["e"],
                                      r["fs"], first["p"]["index"], first["s"]["index"])
        loc["origin_time"] = first["tp"] - loc["hypo_dist"] / C.VP_KM_S
        return loc

    # ---------- ขั้นที่ 4: ขนาดแผ่นดินไหว ----------
    def estimate_magnitude(self, picks, hypo):
        m_jma, m_pd, m_tc = [], [], []
        for p in picks:
            r = p["record"]
            fs, ip = r["fs"], p["p"]["index"]
            n_pre = max(ip - int(1 * fs), 10)
            R, delta = hypocentral_distance(hypo["lat"], hypo["lon"], hypo["depth"], p["lat"], p["lon"])
            pd = pd_value(r["z"], fs, ip, n_pre=n_pre)
            tc = tau_c(r["z"], fs, ip, n_pre=n_pre)
            mj, A = jma_magnitude(r["n"], r["e"], fs, delta, ip, n_pre=n_pre)
            p.update(Pd_cm=pd, tau_c=tc, A_um=A, M_pd=magnitude_from_pd(pd, R),
                     M_tc=magnitude_from_tauc(tc), M_jma=mj, epi_dist=delta, hypo_dist=R,
                     I_obs=jma_instrumental_intensity(r["z"], r["n"], r["e"], fs))
            m_pd.append(p["M_pd"]); m_tc.append(p["M_tc"]); m_jma.append(mj)
        return {"M_jma": network_magnitude(m_jma), "M_pd": network_magnitude(m_pd),
                "M_tc": network_magnitude(m_tc)}

    # ---------- ขั้นที่ 5-6: เวลาเดินทาง + ความรุนแรงที่ตำแหน่งเรา ----------
    def run(self, processing_delay=1.0, verbose=True):
        picks = self.pick_phases()
        if not picks:
            print("ไม่พบแผ่นดินไหว"); return None
        hypo = self.locate(picks)
        mags = self.estimate_magnitude(picks, hypo)
        M = mags["M_jma"] if mags["M_jma"] is not None else mags["M_pd"]

        # เวลาออกคำเตือน = สถานีที่ 3 (ขั้นต่ำของ Geiger) ได้ข้อมูล P ครบ 3 วินาที + เวลาประมวลผล
        n_use = min(len(picks), 3)
        alert_time = picks[n_use - 1]["tp"] + C.PD_TAUC_WINDOW_SEC + processing_delay
        if hypo["method"].startswith("single"):        # สถานีเดียวต้องรอคลื่น S ก่อน
            alert_time = max(alert_time, picks[0]["ts"] + processing_delay)
        tt = site_warning(hypo, self.site["lat"], self.site["lon"], alert_time)
        pred = predict_intensity(M, tt["hypo_dist_km"], hypo["depth"],
                                 self.site.get("vs30", C.DEFAULT_VS30))

        if pred["intensity"] >= C.WARNING_INTENSITY:
            level = "WARNING (คาดความรุนแรง JMA >= 5-)"
        elif M >= C.FORECAST_MAGNITUDE or pred["intensity"] >= C.FORECAST_INTENSITY:
            level = "FORECAST (แจ้งข้อมูล)"
        else:
            level = "ไม่ออกคำเตือน"

        result = {"picks": picks, "hypocenter": hypo, "magnitude": mags, "M": M,
                  "alert_time": alert_time, "site": tt, "prediction": pred, "level": level}
        if verbose:
            self.report(result)
        return result

    def report(self, res):
        h, s, pr = res["hypocenter"], res["site"], res["prediction"]
        line = "=" * 78
        print(line); print(" EARTHQUAKE EARLY WARNING REPORT"); print(line)
        print(f"{'Sta':<6}{'Tp(s)':>9}{'Ts(s)':>9}{'S-P':>6}{'Δ(km)':>8}{'Pd(cm)':>10}"
              f"{'τc':>6}{'M_Pd':>6}{'M_τc':>6}{'M_JMA':>6}{'I_obs':>7}")
        for p in res["picks"]:
            ts = p["ts"]
            print(f"{p['code']:<6}{p['tp']:9.2f}{(ts if ts else float('nan')):9.2f}"
                  f"{((ts - p['tp']) if ts else float('nan')):6.2f}{p['epi_dist']:8.1f}"
                  f"{p['Pd_cm']:10.2e}{p['tau_c']:6.2f}{p['M_pd']:6.2f}{p['M_tc']:6.2f}"
                  f"{p['M_jma']:6.2f}{p['I_obs']:5.1f}({jma_intensity_class(p['I_obs'])})")
        print(line)
        print(f" Method      : {h['method']}")
        print(f" Epicenter   : lat {h['lat']:.3f}, lon {h['lon']:.3f}, depth {h['depth']:.1f} km")
        print(f" Origin time : {h['origin_time']:.2f} s" + (f"   (RMS {h['rms']:.3f} s)" if 'rms' in h else ""))
        m = res["magnitude"]
        print(f" Magnitude   : M_JMA {m['M_jma']:.1f} | M_Pd {m['M_pd']:.1f} | M_τc {m['M_tc']:.1f}")
        print(line)
        print(f" ตำแหน่งของเรา: {self.site['name']} ({self.site['lat']}, {self.site['lon']})")
        print(f"   ระยะทางจากจุดเหนือศูนย์ : {s['epi_dist_km']:.1f} km  (hypocentral {s['hypo_dist_km']:.1f} km)")
        print(f"   คลื่น P เดินทาง          : {s['p_travel_s']:.1f} s")
        print(f"   คลื่น S เดินทาง          : {s['s_travel_s']:.1f} s")
        print(f"   ออกคำเตือนหลังเกิดเหตุ   : {res['alert_time'] - h['origin_time']:.1f} s")
        print(f"   เวลาเตือนล่วงหน้า (S)    : {s['lead_time_s']:.1f} s")
        print(f"   รัศมี Blind zone         : {s['blind_zone_km']:.1f} km")
        print(f"   PGV ทำนาย                : {pr['pgv_surface']:.2f} cm/s")
        print(f"   ความรุนแรงทำนาย          : JMA {pr['intensity']:.1f} (ระดับ {pr['jma_class']}),"
              f" MMI ~ {pr['mmi']:.1f}")
        print(f"   ระดับการแจ้งเตือน        : {res['level']}")
        print(line)
