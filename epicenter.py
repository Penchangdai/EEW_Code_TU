"""
epicenter.py
หาจุดศูนย์กลางแผ่นดินไหว (Epicenter) และจุดกำเนิด (Hypocenter) จากหลายสถานี

  Geiger's method (Iterative least squares) บนผิวโลกทรงกลม + โลก 2 ชั้น
    ใช้เวลาคลื่น P, S (Sg/Sn) และ Lg จากทุกสถานี หา (lat, lon, ความลึก, T0)
    ให้น้ำหนักแต่ละ pick ตาม weight 0-4 แบบ Hypoinverse (weighted least squares)
      weight 0 -> 1.00, 1 -> 0.75, 2 -> 0.50, 3 -> 0.25, 4 -> 0 (ไม่ใช้)
"""
import numpy as np

import config as C


# ตัวคูณน้ำหนักของ weight code 0-4 (แบบ Hypoinverse: 1 - code/4)
WEIGHT_FACTOR = getattr(C, "WEIGHT_FACTOR", (1.0, 0.75, 0.5, 0.25, 0.0))


def pick_weight(code):
    """แปลง weight code (0-4) เป็นตัวคูณน้ำหนัก  None = น้ำหนักเต็ม"""
    if code is None:
        return 1.0
    return float(WEIGHT_FACTOR[int(np.clip(code, 0, len(WEIGHT_FACTOR) - 1))])


# ---------------------------------------------------------------------------
#  เฟสคลื่นที่ใช้หาตำแหน่ง (key ใน dict ของสถานี -> เฟส)
#    "tp"  / "wp"  : คลื่น P ที่มาถึงก่อน (Pg หรือ Pn ตามระยะ)
#    "ts"  / "ws"  : คลื่น S ที่มาถึงก่อน (Sg หรือ Sn ตามระยะ)
#    "tlg" / "wlg" : คลื่น Lg (คลื่น S ที่สะท้อนไปมาในเปลือกโลก ความเร็ว ~3.5 km/s)
#                    ที่สถานีไกล > ~300 km มักเป็นคลื่น S ที่เห็นชัดที่สุด
# ---------------------------------------------------------------------------
PHASES = (("tp", "wp", "P"), ("ts", "ws", "S"), ("tlg", "wlg", "Lg"))
VLG_KM_S = getattr(C, "VLG_KM_S", 3.5)


def phase_travel_time(delta_km, depth_km, phase):
    """เวลาเดินทางของแต่ละเฟส: P, S (โลก 2 ชั้น เลือกคลื่นที่มาถึงก่อน), Lg (ความเร็วคงที่)"""
    if phase == "Lg":
        return np.sqrt(np.asarray(delta_km) ** 2 + depth_km ** 2) / VLG_KM_S
    return travel_time_2layer(delta_km, depth_km, phase)


def phase_name(delta_km, phase, depth_km=10.0):
    """ชื่อเฟสตามระยะและความลึก เช่น Pg/Pn, Sg/Sn"""
    if phase == "Lg":
        return "Lg"
    v1 = C.VP_KM_S if phase == "P" else C.VS_KM_S
    direct = np.sqrt(delta_km ** 2 + depth_km ** 2) / v1
    return phase + ("g" if direct <= travel_time_2layer(delta_km, depth_km, phase) + 1e-9 else "n")


def observations(stations, use_weights=True, phases=("P", "S", "Lg")):
    """รวมเวลาคลื่นทุกเฟสของทุกสถานี -> list ของ (index สถานี, เฟส, lat, lon, เวลา, น้ำหนัก)"""
    obs = []
    for j, s in enumerate(stations):
        for key, wkey, ph in PHASES:
            if ph in phases and s.get(key) is not None:
                w = pick_weight(s.get(wkey)) if use_weights else 1.0
                if w > 0:
                    obs.append((j, ph, s["lat"], s["lon"], float(s[key]), w))
    return obs


def geiger_locate_2layer(stations, init_depths=(10, 30), max_iter=30, damping=1e-2, tol=1e-3,
                         use_weights=True, fix_depth=None, init=None, phases=("P", "S", "Lg")):
    """
    Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น ใช้ได้ทั้งคลื่น P, S (Sg/Sn) และ Lg พร้อมกัน
    เหมาะกับสถานีไกลหลายร้อย km: คิด intercept time ของคลื่น Pn/Sn ที่วิ่งใต้ Moho
    พารามิเตอร์ที่หา: lat, lon, depth, T0   (อนุพันธ์คำนวณเชิงตัวเลข)
    stations  : list ของ dict {"lat", "lon", "tp", "ts", "tlg" (ถ้ามี), "wp", "ws", "wlg" (weight 0-4)}
    fix_depth : ใส่ค่า (km) เพื่อตรึงความลึก เมื่อสถานีไกลมากจนหาความลึกไม่ได้
    init      : (lat, lon) จุดเริ่มต้น เช่น ผลจาก grid_search_locate  (None = ใกล้สถานีที่ P มาถึงก่อน)
    phases    : เฟสที่ใช้ เช่น ("P",) = ใช้แค่คลื่น P
    """
    from signal_utils import haversine_km
    obs = observations(stations, use_weights, phases)
    n_free = 3 if fix_depth is not None else 4
    if len(obs) < n_free:
        raise ValueError(f"ต้องมีเวลาคลื่นอย่างน้อย {n_free} ค่า")
    slat = np.array([o[2] for o in obs]); slon = np.array([o[3] for o in obs])
    ph = np.array([o[1] for o in obs]); t_only = np.array([o[4] for o in obs])
    w_t = np.array([o[5] for o in obs])
    tobs, wts = t_only, w_t

    def tcalc(m):
        d = haversine_km(m[0], m[1], slat, slon)
        z = max(m[2], 0.0)
        t = np.empty(len(d))
        for p in ("P", "S", "Lg"):
            k = ph == p
            if k.any():
                t[k] = phase_travel_time(d[k], z, p)
        return m[3] + t

    def wrms(r):
        return float(np.sqrt(np.sum(wts * r ** 2) / np.sum(wts)))

    free = [0, 1, 3] if fix_depth is not None else [0, 1, 2, 3]
    step = np.array([0.01, 0.01, 0.5, 0.01])
    tp_st = [s for s in stations if s.get("tp") is not None]
    first = min(tp_st, key=lambda s: s["tp"]) if tp_st else stations[0]
    best = None
    for z_init in ((fix_depth,) if fix_depth is not None else init_depths):
        if init is not None:
            m = np.array([init[0], init[1], float(z_init), 0.0])
        else:
            m = np.array([first["lat"] + 0.5, first["lon"] - 0.5, float(z_init), 0.0])
        m[3] = np.median(tobs - (tcalc(m) - m[3]))
        for _ in range(max_iter):
            r = tobs - tcalc(m)
            G = np.zeros((len(tobs), 4))
            for k in free:
                dm_k = np.zeros(4); dm_k[k] = step[k]
                G[:, k] = (tcalc(m + dm_k) - tcalc(m - dm_k)) / (2 * step[k])
            Gf = G[:, free]
            A = Gf.T @ (wts[:, None] * Gf)
            dm = np.linalg.solve(A + damping * np.diag(np.diag(A) + 1e-12), Gf.T @ (wts * r))
            dm[:2] = np.clip(dm[:2], -2.0, 2.0)                      # จำกัดก้าวไม่เกิน 2° ต่อรอบ
            m[free] += dm
            m[2] = abs(m[2])
            if np.all(np.abs(dm[:2]) < tol * 0.01) and (fix_depth is not None or abs(dm[2]) < tol * 10):
                break
        r = tobs - tcalc(m)
        if best is None or wrms(r) < best[1]:
            best = (m.copy(), wrms(r), r, G[:, free])
    m, rms, r, Gf = best
    n_eff = np.sum(wts) ** 2 / np.sum(wts ** 2)
    dof = max(n_eff - len(free), 1)
    sigma2 = np.sum(wts * r ** 2) / np.sum(wts) * n_eff / dof
    err = np.full(4, np.nan)
    try:
        cov = sigma2 * np.linalg.inv(Gf.T @ ((wts / wts.mean())[:, None] * Gf))
        err[free] = np.sqrt(np.abs(np.diag(cov)))
        km_lat, km_lon = 111.195, 111.195 * np.cos(np.radians(m[0]))
        err[0], err[1] = err[1] * km_lon, err[0] * km_lat      # -> x (ตะวันออก), y (เหนือ) เป็น km
    except np.linalg.LinAlgError:
        pass
    if fix_depth is not None:
        err[2] = 0.0
    n_by = {p: int(np.sum(ph == p)) for p in ("P", "S", "Lg")}
    return {"lat": float(m[0]), "lon": float(m[1]), "depth": float(m[2]), "origin_time": float(m[3]),
            "rms": rms, "err_x_km": err[0], "err_y_km": err[1], "err_z_km": err[2], "err_t0_s": err[3],
            "n_obs": len(tobs), "n_by_phase": n_by, "model": "2layer",
            "method": "Geiger (ทรงกลม + โลก 2 ชั้น)" + (", weighted" if use_weights else "")
                      + (f", ความลึกตรึง {fix_depth:g} km" if fix_depth is not None else "")
                      + f", P {n_by['P']} / S {n_by['S']} / Lg {n_by['Lg']}"}


def predict_arrival(hypo, lat, lon, phase):
    """เวลาที่เฟสนั้นควรมาถึงสถานี (epoch s) จากตำแหน่ง hypo"""
    from signal_utils import haversine_km
    d = float(haversine_km(hypo["lat"], hypo["lon"], lat, lon))
    return hypo["origin_time"] + float(phase_travel_time(d, hypo["depth"], phase))


# ============================ Robust grid search (ทรงกลม + โลก 2 ชั้น) ============================
def travel_time_2layer(delta_km, depth_km, phase="P"):
    """
    เวลาเดินทางในโลก 2 ชั้น (เปลือกโลก + เนื้อโลก) เลือกคลื่นที่มาถึงก่อนเอง
      คลื่นตรง  (Pg/Sg): t = R / V1,  R = sqrt(Δ² + h²)
      คลื่นหัว  (Pn/Sn): h < H:  t = Δ / V2 + (2H − h)·q,   q = sqrt(1/V1² − 1/V2²)   (H = ความหนาเปลือกโลก)
                         h ≥ H:  t = sqrt(Δ² + (h − H)²) / V2 + H·q   (จุดกำเนิดอยู่ในเนื้อโลก โดยประมาณ)
      t = min(คลื่นตรง, คลื่นหัว)
    """
    v1, v2 = (C.VP_KM_S, C.VPN_KM_S) if phase == "P" else (C.VS_KM_S, C.VSN_KM_S)
    H = C.CRUST_THICKNESS_KM
    q = np.sqrt(1 / v1 ** 2 - 1 / v2 ** 2)
    direct = np.sqrt(delta_km ** 2 + depth_km ** 2) / v1
    below = np.maximum(np.asarray(depth_km, float) - H, 0.0)      # ส่วนที่ลึกกว่า Moho
    head = np.sqrt(delta_km ** 2 + below ** 2) / v2 + (2 * H - np.minimum(depth_km, H)) * q
    return np.minimum(direct, head)


def azimuthal_gap(lat, lon, stations):
    """ช่องว่างมุมที่ใหญ่ที่สุด (องศา) เมื่อมองจากจุดศูนย์กลางไปยังสถานี  (ควร < 180°)"""
    from signal_utils import azimuth_deg
    az = np.sort([azimuth_deg(lat, lon, s["lat"], s["lon"]) % 360 for s in stations])
    if len(az) < 2:
        return 360.0
    gaps = np.diff(np.concatenate((az, [az[0] + 360])))
    return float(gaps.max())


