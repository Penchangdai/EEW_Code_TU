"""
epicenter.py
หาจุดศูนย์กลางแผ่นดินไหว (Epicenter) และจุดกำเนิด (Hypocenter)

  วิธีที่ 1  สถานีเดียว (Single-station, แนวคิดเดียวกับ JMA)
            - ระยะทางจากสูตร Omori:  D = k (Ts - Tp)
            - ทิศของแหล่งกำเนิด (back-azimuth) จาก Principal Component Analysis
              ของการเคลื่อนที่ของอนุภาคคลื่น P ช่วงแรก
  วิธีที่ 2  หลายสถานี (Network) ด้วย Geiger's method (Iterative least squares)
            ใช้เวลาคลื่น P และ S จากทุกสถานี หา (x, y, z, T0)
            ให้น้ำหนักแต่ละ pick ตาม weight 0-4 แบบ Hypoinverse (weighted least squares)
              weight 0 -> 1.00, 1 -> 0.75, 2 -> 0.50, 3 -> 0.25, 4 -> 0 (ไม่ใช้)
            weight อ่านจาก key "wp" (คลื่น P) และ "ws" (คลื่น S) ของแต่ละสถานี
            ถ้าไม่มี key นี้ ถือว่า weight 0 (น้ำหนักเต็ม) เหมือนเดิม
            geiger_locate_2layer = Geiger บนผิวโลกทรงกลม + โลก 2 ชั้น ใช้คลื่น P, S และ Lg พร้อมกัน
                                   (แนะนำสำหรับสถานีไกล)
"""
import numpy as np

import config as C
from signal_utils import (bandpass, acc_to_disp, destination_point,
                          latlon_to_xy, xy_to_latlon)


# ============================ Single station ============================
def omori_distance(sp_time, vp=C.VP_KM_S, vs=C.VS_KM_S):
    """D = k * (Ts - Tp),  k = Vp*Vs/(Vp - Vs)  (~8.4 km/s เมื่อ Vp=6, Vs=3.5)"""
    k = vp * vs / (vp - vs)
    return k * sp_time


def back_azimuth_pca(z, n, e, fs, p_index, window_sec=1.0, band=(0.5, 10.0)):
    """
    หาทิศจากสถานีไปยังแหล่งกำเนิด (back-azimuth, องศาจากทิศเหนือ)
    1) ใช้ข้อมูล N, E ช่วงแรกของคลื่น P ทำ covariance matrix
    2) eigenvector ที่ค่า eigenvalue สูงสุด = ทิศการสั่นในแนวราบ (ยังกำกวม 180°)
    3) แก้ความกำกวมด้วยเครื่องหมายของ cov(Z, H):
       คลื่น P ที่ผิวดิน ถ้า Z ขึ้น -> พื้นเคลื่อนที่ "ออกจาก" แหล่งกำเนิด
    """
    i0, i1 = p_index, p_index + int(window_sec * fs)
    zf = bandpass(z, fs, *band)[i0:i1]
    nf = bandpass(n, fs, *band)[i0:i1]
    ef = bandpass(e, fs, *band)[i0:i1]
    cov = np.cov(np.vstack([nf, ef]))
    w, vec = np.linalg.eigh(cov)
    pn, pe = vec[:, np.argmax(w)]
    h = pn * nf + pe * ef
    if np.sum(zf * h) < 0:            # ให้ (pn, pe) ชี้ทิศ "ออกจาก" แหล่งกำเนิดเมื่อ Z ขึ้น
        pn, pe = -pn, -pe
    away_az = np.degrees(np.arctan2(pe, pn))
    rectilinearity = 1.0 - w.min() / max(w.max(), 1e-30)
    return (away_az + 180.0) % 360.0, rectilinearity


def single_station_location(sta_lat, sta_lon, z, n, e, fs, p_index, s_index=None,
                            sp_time=None, default_depth=10.0):
    """ประมาณตำแหน่งจากสถานีเดียว (ต้องมี S-P time)"""
    if sp_time is None:
        sp_time = (s_index - p_index) / fs
    R = omori_distance(sp_time)                       # ระยะทางจากจุดกำเนิด (hypocentral)
    depth = min(default_depth, 0.9 * R)
    delta = np.sqrt(max(R ** 2 - depth ** 2, 0.0))    # ระยะทางบนผิว (epicentral)
    baz, rect = back_azimuth_pca(z, n, e, fs, p_index)
    lat, lon = destination_point(sta_lat, sta_lon, baz, delta)
    return {"lat": lat, "lon": lon, "depth": depth, "hypo_dist": R,
            "epi_dist": delta, "back_azimuth": baz, "rectilinearity": rect,
            "method": "single-station (Omori + PCA)"}


# ============================ Geiger's method ============================
def _travel_time_and_derivs(m, xi, yi, v):
    """t = T0 + R/v,  R = sqrt((x-xi)^2 + (y-yi)^2 + z^2)  และอนุพันธ์ย่อย"""
    x, y, z, t0 = m
    dx, dy = x - xi, y - yi
    R = np.sqrt(dx ** 2 + dy ** 2 + z ** 2)
    R = np.maximum(R, 1e-6)
    t = t0 + R / v
    G = np.column_stack([dx / (v * R), dy / (v * R), z / (v * R), np.ones_like(R)])
    return t, G


# ตัวคูณน้ำหนักของ weight code 0-4 (แบบ Hypoinverse: 1 - code/4)
WEIGHT_FACTOR = getattr(C, "WEIGHT_FACTOR", (1.0, 0.75, 0.5, 0.25, 0.0))


def pick_weight(code):
    """แปลง weight code (0-4) เป็นตัวคูณน้ำหนัก  None = น้ำหนักเต็ม"""
    if code is None:
        return 1.0
    return float(WEIGHT_FACTOR[int(np.clip(code, 0, len(WEIGHT_FACTOR) - 1))])


def geiger_locate(stations, vp=C.VP_KM_S, vs=C.VS_KM_S, init_depths=(5, 15, 30, 60),
                  max_iter=50, damping=1e-3, tol=1e-4, use_weights=True):
    """
    stations : list ของ dict {"lat", "lon", "tp", "ts"(ถ้ามี), "wp", "ws"(weight 0-4, ถ้ามี)}
    ต้องมีจำนวนเวลาที่อ่านได้ (P+S) อย่างน้อย 4 ค่า (ไม่นับ pick ที่ weight 4)
    use_weights : False = ทุก pick น้ำหนักเท่ากัน (แบบเดิม)
    return dict: lat, lon, depth, origin_time, rms (ถ่วงน้ำหนัก), errors (km, s)
                 + lat0, lon0, x, y (ไว้ให้ geiger_residuals คำนวณด้วยแบบจำลองเดียวกัน)
    """
    lat0 = np.mean([s["lat"] for s in stations])
    lon0 = np.mean([s["lon"] for s in stations])

    xs, ys, tobs, vel, wts = [], [], [], [], []
    for s in stations:
        x, y = latlon_to_xy(s["lat"], s["lon"], lat0, lon0)
        for key, wkey, v in (("tp", "wp", vp), ("ts", "ws", vs)):
            if s.get(key) is not None:
                w = pick_weight(s.get(wkey)) if use_weights else 1.0
                if w <= 0:
                    continue                      # weight 4 = ไม่ใช้
                xs.append(x); ys.append(y); tobs.append(s[key]); vel.append(v); wts.append(w)
    xs, ys, tobs, vel, wts = map(np.asarray, (xs, ys, tobs, vel, wts))
    if len(tobs) < 4:
        raise ValueError("ต้องมีเวลาคลื่น (P+S) อย่างน้อย 4 ค่า สำหรับหา x, y, z, T0")
    W = np.diag(wts)

    # จุดเริ่มต้น: ใกล้สถานีที่คลื่น P มาถึงก่อน
    first = min((s for s in stations if s.get("tp") is not None), key=lambda s: s["tp"])
    x_init, y_init = latlon_to_xy(first["lat"], first["lon"], lat0, lon0)

    def wrms(r):
        return float(np.sqrt(np.sum(wts * r ** 2) / np.sum(wts)))

    best = None
    for z_init in init_depths:
        m = np.array([x_init + 0.1, y_init + 0.1, float(z_init), first["tp"] - z_init / vp])
        for _ in range(max_iter):
            tcal, G = _travel_time_and_derivs(m, xs, ys, vel)
            r = tobs - tcal                                            # residual
            dm = np.linalg.solve(G.T @ W @ G + damping * np.eye(4), G.T @ W @ r)
            m = m + dm
            m[2] = abs(m[2]) if m[2] != 0 else 0.1                     # ความลึกต้องเป็นบวก
            if np.linalg.norm(dm[:3]) < tol:
                break
        tcal, G = _travel_time_and_derivs(m, xs, ys, vel)
        r = tobs - tcal
        rms = wrms(r)
        if best is None or rms < best[1]:
            best = (m.copy(), rms, G, r)

    m, rms, G, r = best
    # ความคลาดเคลื่อน: น้ำหนักแบบสัมพัทธ์ ปรับสเกลด้วยจำนวนข้อมูลที่มีผลจริง
    n_eff = np.sum(wts) ** 2 / np.sum(wts ** 2)
    dof = max(n_eff - 4, 1)
    sigma2 = np.sum(wts * r ** 2) / np.sum(wts) * n_eff / dof
    try:
        cov = sigma2 * np.linalg.inv(G.T @ (W / wts.mean()) @ G)
        err = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        err = np.full(4, np.nan)
    lat, lon = xy_to_latlon(m[0], m[1], lat0, lon0)
    return {"lat": float(lat), "lon": float(lon), "depth": float(m[2]),
            "origin_time": float(m[3]), "rms": rms,
            "err_x_km": err[0], "err_y_km": err[1], "err_z_km": err[2], "err_t0_s": err[3],
            "n_obs": len(tobs), "lat0": float(lat0), "lon0": float(lon0),
            "x": float(m[0]), "y": float(m[1]),
            "method": "Geiger (P+S)" + (", weighted" if use_weights else "")}


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
                         use_weights=True, fix_depth=None, init=None, phases=("P", "S", "Lg"),
                         baz_obs=None, baz_s_per_deg=0.2, baz_max_dev=45.0):
    """
    Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น ใช้ได้ทั้งคลื่น P, S (Sg/Sn) และ Lg พร้อมกัน
    เหมาะกับสถานีไกลหลายร้อย km: คิด intercept time ของคลื่น Pn/Sn ที่วิ่งใต้ Moho
    พารามิเตอร์ที่หา: lat, lon, depth, T0   (อนุพันธ์คำนวณเชิงตัวเลข)
    stations  : list ของ dict {"lat", "lon", "tp", "ts", "tlg" (ถ้ามี), "wp", "ws", "wlg" (weight 0-4)}
    fix_depth : ใส่ค่า (km) เพื่อตรึงความลึก เมื่อสถานีไกลมากจนหาความลึกไม่ได้
    init      : (lat, lon) จุดเริ่มต้น เช่น ผลจาก grid_search_locate  (None = ใกล้สถานีที่ P มาถึงก่อน)
    phases    : เฟสที่ใช้ เช่น ("P",) = ใช้แค่คลื่น P
    baz_obs   : list ของ (lat, lon, back-azimuth จากคลื่น P องศา, น้ำหนัก) ใช้เป็นข้อมูลเพิ่ม
                residual ของทิศแปลงเป็นวินาทีด้วย baz_s_per_deg (0.2 s/องศา = ทิศคลาด 10° เท่ากับเวลาคลาด 2 s)
                ทิศที่ต่างจากตำแหน่งปัจจุบันเกิน baz_max_dev ไม่ใช้ในรอบนั้น
                ช่วยมากเมื่อสถานีอยู่ด้านเดียว: เวลาคลื่นอย่างเดียวบอกตำแหน่งตามแนวโค้งรอบสถานีได้ไม่ดี
    """
    from signal_utils import haversine_km
    obs = observations(stations, use_weights, phases)
    n_free = 3 if fix_depth is not None else 4
    if len(obs) < n_free:
        raise ValueError(f"ต้องมีเวลาคลื่นอย่างน้อย {n_free} ค่า")
    from signal_utils import azimuth_deg
    slat = np.array([o[2] for o in obs]); slon = np.array([o[3] for o in obs])
    ph = np.array([o[1] for o in obs]); t_only = np.array([o[4] for o in obs])
    w_t = np.array([o[5] for o in obs])
    nt = len(t_only)
    bz = [b for b in (baz_obs or []) if b[2] is not None and np.isfinite(b[2])]
    blat = np.array([b[0] for b in bz]); blon = np.array([b[1] for b in bz])
    bdeg = np.array([b[2] for b in bz]); bw0 = np.array([b[3] for b in bz], float)
    tobs = np.concatenate((t_only, bdeg * baz_s_per_deg))
    wts = np.concatenate((w_t, bw0))

    def tcalc(m):
        d = haversine_km(m[0], m[1], slat, slon)
        z = max(m[2], 0.0)
        t = np.empty(len(d))
        for p in ("P", "S", "Lg"):
            k = ph == p
            if k.any():
                t[k] = phase_travel_time(d[k], z, p)
        if not len(bz):
            return m[3] + t
        bp = azimuth_deg(blat, blon, m[0], m[1])
        bp = bdeg + ((bp - bdeg + 180.0) % 360.0) - 180.0          # ให้ต่อเนื่องรอบค่าที่วัด (ไม่กระโดด 0/360)
        return np.concatenate((m[3] + t, bp * baz_s_per_deg))

    def baz_weights(m):
        """ทิศที่ต่างจากตำแหน่งปัจจุบันเกิน baz_max_dev ไม่ใช้"""
        if not len(bz):
            return
        dev = np.abs(tobs[nt:] - tcalc(m)[nt:]) / baz_s_per_deg
        wts[nt:] = np.where(dev <= baz_max_dev, bw0, 0.0)

    def wrms(r):
        return float(np.sqrt(np.sum(w_t * r[:nt] ** 2) / np.sum(w_t)))

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
        m[3] = np.median(t_only - (tcalc(m)[:nt] - m[3]))
        for _ in range(max_iter):
            baz_weights(m)
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
    baz_weights(m)
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
            "n_obs": nt, "n_by_phase": n_by, "model": "2layer", "n_baz": int(np.sum(wts[nt:] > 0)),
            "baz_res": (tobs[nt:] - tcalc(m)[nt:]) / baz_s_per_deg if len(bz) else np.array([]),
            "method": "Geiger (ทรงกลม + โลก 2 ชั้น)" + (", weighted" if use_weights else "")
                      + (f", ความลึกตรึง {fix_depth:g} km" if fix_depth is not None else "")
                      + f", P {n_by['P']} / S {n_by['S']} / Lg {n_by['Lg']}"}


def predict_arrival(hypo, lat, lon, phase):
    """เวลาที่เฟสนั้นควรมาถึงสถานี (epoch s) จากตำแหน่ง hypo"""
    from signal_utils import haversine_km
    d = float(haversine_km(hypo["lat"], hypo["lon"], lat, lon))
    return hypo["origin_time"] + float(phase_travel_time(d, hypo["depth"], phase))


def geiger_residuals(hypo, stations, vp=C.VP_KM_S, vs=C.VS_KM_S):
    """
    residual = เวลาที่วัด - เวลาที่คำนวณ ด้วยแบบจำลองเดียวกับที่ใช้หาตำแหน่ง
    hypo จาก geiger_locate_2layer -> ทรงกลม + โลก 2 ชั้น (P, S, Lg)
    hypo จาก geiger_locate        -> ระนาบแบน + ความเร็วคงที่ (P, S)
    ใช้ได้กับสถานีที่ไม่ได้อยู่ในชุดที่ใช้หาตำแหน่งด้วย (สำหรับตรวจทานแบบ leave-one-out)
    return list ของ (index สถานี, เฟส, residual)
    """
    out = []
    if hypo.get("model") == "2layer":
        for j, s in enumerate(stations):
            for key, _, ph in PHASES:
                if s.get(key) is not None:
                    out.append((j, ph, float(s[key] - predict_arrival(hypo, s["lat"], s["lon"], ph))))
        return out
    for j, s in enumerate(stations):
        x, y = latlon_to_xy(s["lat"], s["lon"], hypo["lat0"], hypo["lon0"])
        R = np.sqrt((hypo["x"] - x) ** 2 + (hypo["y"] - y) ** 2 + hypo["depth"] ** 2)
        for key, ph, v in (("tp", "P", vp), ("ts", "S", vs)):
            if s.get(key) is not None:
                out.append((j, ph, float(s[key] - (hypo["origin_time"] + R / v))))
    return out


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


def _robust_misfit(tobs, tcal):
    """(เดิม) T0 = มัธยฐาน, misfit = มัธยฐาน |residual|  -- ทิ้งข้อมูลได้เกือบครึ่ง ไม่ใช้แล้ว"""
    d = tobs[..., :] - tcal
    t0 = np.median(d, axis=-1)
    res = d - t0[..., None]
    return np.median(np.abs(res), axis=-1), t0


def _fit_t0_truncated(tobs, tP, tSn, tLg, is_p, c=5.0, allow_lg=True, n_iter=2, w=None):
    """
    หา T0 และเลือกเฟสคลื่น S (Sn/Sg หรือ Lg) ทุกจุดบนตาราง แล้วคิด misfit แบบ truncated L1
      1) T0 เริ่มจากมัธยฐานของ (tobs − tP) ของคลื่น P   (ถ้าไม่มี P ใช้ทุกค่า)
      2) คลื่น S แต่ละค่า ตีความเป็นเฟสที่ residual เล็กกว่า (Sn หรือ Lg)
      3) ปรับ T0 = ค่าเฉลี่ยของ residual ที่คลาดไม่เกิน c วินาที (n_iter รอบ)
      4) misfit = mean( min(|r − T0|, c) )   ทุกเวลาคลื่นมีผล แต่ค่าผิดมากมีผลได้ไม่เกิน c
    w : น้ำหนักต่อเวลาคลื่น (จาก weight 0-4) ถ้ามี ใช้ค่าเฉลี่ยแบบถ่วงน้ำหนักในข้อ 3 และ 4
    คืนค่า misfit, T0, residual (หลังลบ T0), is_lg (bool ต่อเวลาคลื่น)
    """
    if w is None:
        w = np.ones(tobs.shape[-1])
    rP = tobs - tP
    if is_p.any():
        t0 = np.nanmedian(np.where(is_p, rP, np.nan), axis=-1)
        pass
    else:
        t0 = np.median(tobs - tSn, axis=-1)
    for _ in range(n_iter + 1):
        rSn = tobs - tSn - t0[..., None]
        rLg = tobs - tLg - t0[..., None]
        is_lg = (~is_p) & allow_lg & (np.abs(rLg) < np.abs(rSn))
        tcal = np.where(is_p, tP, np.where(is_lg, tLg, tSn))
        r = tobs - tcal                                        # ยังไม่ลบ T0
        inl = np.abs(r - t0[..., None]) <= c
        n_in = np.sum(np.where(inl, w, 0.0), axis=-1)
        t0 = np.where(n_in > 0, np.sum(np.where(inl, w * r, 0.0), axis=-1) / np.maximum(n_in, 1e-12), t0)
    res = r - t0[..., None]
    mis = np.sum(w * np.minimum(np.abs(res), c), axis=-1) / np.sum(w)
    return mis, t0, res, is_lg


def _axial_diff_deg(a, b):
    """มุมต่างแบบไม่สนทิศหน้า/หลัง (0-90°) ใช้กับ back-azimuth จาก PCA ที่กำกวม 180°"""
    d = np.abs((np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0)
    return np.minimum(d, 180.0 - d)


def azimuthal_gap(lat, lon, stations):
    """ช่องว่างมุมที่ใหญ่ที่สุด (องศา) เมื่อมองจากจุดศูนย์กลางไปยังสถานี  (ควร < 180°)"""
    from signal_utils import azimuth_deg
    az = np.sort([azimuth_deg(lat, lon, s["lat"], s["lon"]) % 360 for s in stations])
    if len(az) < 2:
        return 360.0
    gaps = np.diff(np.concatenate((az, [az[0] + 360])))
    return float(gaps.max())


def grid_search_locate(stations, use_s=True, depth_km=10.0, step_deg=0.25, pad_deg=12.0,
                       outlier_sec=None, use_baz=False, baz_sec_per_deg=0.03, baz_min_rect=0.5,
                       relocate=True, trunc_sec=5.0, allow_lg=True):
    """
    หาตำแหน่งด้วยการลองทุกจุดบนตาราง lat/lon (ระยะทางจริงบนผิวโลก + เวลาเดินทางโลก 2 ชั้น)
    stations : list ของ dict {"code", "lat", "lon", "tp", "ts"(ถ้ามี), "baz"(ถ้ามี), "rect"(ถ้ามี),
                              "wp", "ws" (weight 0-4 ถ้ามี, ไม่มี = น้ำหนักเต็ม)}
    use_s    : ใช้เวลาคลื่น S ด้วยไหม
    use_baz  : ใช้ทิศทางจาก PCA ช่วยด้วย (แก้ปัญหาสถานีอยู่ฝั่งเดียวกัน)
               misfit = มัธยฐาน|residual เวลา| + baz_sec_per_deg × มัธยฐาน(มุมต่าง)
               ใช้เฉพาะสถานีที่ rectilinearity >= baz_min_rect
    relocate : ตัดเวลาคลื่นที่ไม่สอดคล้อง (outlier) แล้วคำนวณใหม่อีกรอบ
    คืนค่า dict: lat, lon, origin_time, misfit, residuals, outlier, ความไม่แน่นอน, azimuthal gap
    """
    from signal_utils import haversine_km, azimuth_deg
    obs = []                                  # (index สถานี, phase, เวลา, น้ำหนัก)
    for k, s in enumerate(stations):
        if s.get("tp") is not None and pick_weight(s.get("wp")) > 0:
            obs.append((k, "P", s["tp"], pick_weight(s.get("wp"))))
        if use_s and s.get("ts") is not None and pick_weight(s.get("ws")) > 0:
            obs.append((k, "S", s["ts"], pick_weight(s.get("ws"))))
    if len(obs) < 4:
        raise ValueError("ต้องมีเวลาคลื่นอย่างน้อย 4 ค่า")
    tobs_all = np.array([o[2] for o in obs])
    slat_all = np.array([stations[o[0]]["lat"] for o in obs])
    slon_all = np.array([stations[o[0]]["lon"] for o in obs])
    isp_all = np.array([o[1] == "P" for o in obs])
    w_all = np.array([o[3] for o in obs])

    bz = [s for s in stations if use_baz and s.get("baz") is not None and s.get("rect", 1.0) >= baz_min_rect]
    blat = np.array([s["lat"] for s in bz]); blon = np.array([s["lon"] for s in bz])
    bbaz = np.array([s["baz"] for s in bz])

    def predict(LA, LO, keep):
        D = haversine_km(LA[..., None], LO[..., None], slat_all[keep], slon_all[keep])
        tP = travel_time_2layer(D, depth_km, "P")
        tSn = travel_time_2layer(D, depth_km, "S")
        tLg = np.sqrt(D ** 2 + depth_km ** 2) / getattr(C, "VLG_KM_S", 3.5)   # คลื่น Lg
        return tP, tSn, tLg

    def cost(LA, LO, keep):
        tP, tSn, tLg = predict(LA, LO, keep)
        mis, t0, _, _ = _fit_t0_truncated(tobs_all[keep], tP, tSn, tLg, isp_all[keep],
                                          c=trunc_sec, allow_lg=allow_lg, w=w_all[keep])
        tot = mis.copy()
        if len(bz):
            az = azimuth_deg(blat, blon, LA[..., None], LO[..., None])   # ทิศจากสถานีไปจุดที่ลอง
            tot = tot + baz_sec_per_deg * np.median(_axial_diff_deg(az, bbaz), axis=-1)
        return tot, mis, t0

    def search(keep):
        lat_min, lat_max = slat_all.min() - pad_deg, slat_all.max() + pad_deg
        lon_min, lon_max = slon_all.min() - pad_deg, slon_all.max() + pad_deg
        lats, lons = np.arange(lat_min, lat_max, step_deg), np.arange(lon_min, lon_max, step_deg)
        LA, LO = np.meshgrid(lats, lons, indexing="ij")
        tot, _, _ = cost(LA, LO, keep)
        i, j = np.unravel_index(np.argmin(tot), tot.shape)
        lats = np.arange(lats[i] - 1, lats[i] + 1, 0.02)
        lons = np.arange(lons[j] - 1, lons[j] + 1, 0.02)
        LA, LO = np.meshgrid(lats, lons, indexing="ij")
        tot, mis, t0 = cost(LA, LO, keep)
        i, j = np.unravel_index(np.argmin(tot), tot.shape)
        # ความไม่แน่นอน: ขอบเขตของจุดที่ cost ไม่เกินค่าต่ำสุด + 1 s
        ok = tot <= tot[i, j] + 1.0
        unc = float(np.max(haversine_km(LA[ok], LO[ok], lats[i], lons[j]))) if ok.any() else 0.0
        return float(lats[i]), float(lons[j]), float(mis[i, j]), float(t0[i, j]), float(tot[i, j]), unc

    keep = np.ones(len(obs), bool)
    la, lo, mis, t0, tot, unc = search(keep)

    def residuals(la, lo, t0):
        """residual ของทุกเวลาคลื่น (รวมค่าที่ถูกตัด) และเฟสคลื่น S ที่เลือก ที่จุดคำตอบ"""
        allk = np.ones(len(obs), bool)
        tP, tSn, tLg = predict(np.array(la), np.array(lo), allk)
        _, _, _, is_lg = _fit_t0_truncated(tobs_all, tP, tSn, tLg, isp_all, c=trunc_sec, allow_lg=allow_lg,
                                           w=w_all)
        tcal = np.where(isp_all, tP, np.where(is_lg, tLg, tSn))
        global_lg[:] = is_lg
        return tobs_all - t0 - tcal

    global_lg = np.zeros(len(obs), bool)
    res = residuals(la, lo, t0)
    mad = np.median(np.abs(res)) * 1.4826
    thr = outlier_sec if outlier_sec else max(3 * mad, 3.0)
    thr = min(thr, 2 * trunc_sec) if not outlier_sec else thr     # MAD พองเมื่อเวลาผิดเกือบครึ่ง จึงจำกัดเพดาน
    out = np.abs(res) > thr
    n_pass = 1
    if relocate and out.any() and (~out).sum() >= 4:          # รอบที่ 2: ใช้เฉพาะเวลาที่สอดคล้อง
        la, lo, mis, t0, tot, unc = search(~out)
        res = residuals(la, lo, t0)
        out = np.abs(res) > thr
        n_pass = 2

    rows = [{"code": stations[o[0]]["code"], "phase": o[1], "residual": float(r), "outlier": bool(b),
             "weight_factor": float(o[3]),
             "phase_name": "P" if o[1] == "P" else ("Lg" if lg else "Sn/Sg")}
            for o, r, b, lg in zip(obs, res, out, global_lg)]
    baz_rows = []
    for s in bz:
        az = azimuth_deg(s["lat"], s["lon"], la, lo)
        baz_rows.append({"code": s["code"], "baz": float(s["baz"]), "az_to_epi": float(az),
                         "diff": float(_axial_diff_deg(az, s["baz"]))})
    return {"lat": la, "lon": lo, "depth": depth_km, "origin_time": t0,
            "misfit": mis, "cost": tot, "threshold": float(thr), "residuals": rows,
            "n_obs": len(obs), "uncertainty_km": unc, "passes": n_pass,
            "gap": azimuthal_gap(la, lo, stations), "baz_residuals": baz_rows,
            "method": "Robust grid search, truncated L1" + (" (P+S, แยก Sn/Lg)" if use_s else " (P only)")
                      + (" + ทิศทาง PCA" if len(bz) else "")}


def delta_from_sp(sp_time, depth_km=10.0):
    """ระยะบนผิวโลก Δ จากเวลา S−P ด้วยโลก 2 ชั้น (ใช้ได้ทั้งสถานีใกล้และไกล ไม่ต้องเลือก profile)"""
    d = np.linspace(0, 4000, 8001)
    sp = travel_time_2layer(d, depth_km, "S") - travel_time_2layer(d, depth_km, "P")
    return float(np.interp(sp_time, sp, d))