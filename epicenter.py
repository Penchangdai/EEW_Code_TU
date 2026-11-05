"""
epicenter.py
หาจุดศูนย์กลางแผ่นดินไหว (Epicenter) และจุดกำเนิด (Hypocenter)

  วิธีที่ 1  สถานีเดียว (Single-station, แนวคิดเดียวกับ JMA)
            - ระยะทางจากสูตร Omori:  D = k (Ts - Tp)
            - ทิศของแหล่งกำเนิด (back-azimuth) จาก Principal Component Analysis
              ของการเคลื่อนที่ของอนุภาคคลื่น P ช่วงแรก
  วิธีที่ 2  หลายสถานี (Network) ด้วย Geiger's method (Iterative least squares)
            ใช้เวลาคลื่น P และ S จากทุกสถานี หา (x, y, z, T0)
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


def geiger_locate(stations, vp=C.VP_KM_S, vs=C.VS_KM_S, init_depths=(5, 15, 30, 60),
                  max_iter=50, damping=1e-3, tol=1e-4):
    """
    stations : list ของ dict {"lat", "lon", "tp", "ts"(ถ้ามี)}
    ต้องมีจำนวนเวลาที่อ่านได้ (P+S) อย่างน้อย 4 ค่า
    return dict: lat, lon, depth, origin_time, rms, errors (km, s)
    """
    lat0 = np.mean([s["lat"] for s in stations])
    lon0 = np.mean([s["lon"] for s in stations])

    xs, ys, tobs, vel = [], [], [], []
    for s in stations:
        x, y = latlon_to_xy(s["lat"], s["lon"], lat0, lon0)
        for key, v in (("tp", vp), ("ts", vs)):
            if s.get(key) is not None:
                xs.append(x); ys.append(y); tobs.append(s[key]); vel.append(v)
    xs, ys, tobs, vel = map(np.asarray, (xs, ys, tobs, vel))
    if len(tobs) < 4:
        raise ValueError("ต้องมีเวลาคลื่น (P+S) อย่างน้อย 4 ค่า สำหรับหา x, y, z, T0")

    # จุดเริ่มต้น: ใกล้สถานีที่คลื่น P มาถึงก่อน
    first = min(stations, key=lambda s: s["tp"])
    x_init, y_init = latlon_to_xy(first["lat"], first["lon"], lat0, lon0)

    best = None
    for z_init in init_depths:
        m = np.array([x_init + 0.1, y_init + 0.1, float(z_init), first["tp"] - z_init / vp])
        for _ in range(max_iter):
            tcal, G = _travel_time_and_derivs(m, xs, ys, vel)
            r = tobs - tcal                                            # residual
            dm = np.linalg.solve(G.T @ G + damping * np.eye(4), G.T @ r)
            m = m + dm
            m[2] = abs(m[2]) if m[2] != 0 else 0.1                     # ความลึกต้องเป็นบวก
            if np.linalg.norm(dm[:3]) < tol:
                break
        tcal, G = _travel_time_and_derivs(m, xs, ys, vel)
        r = tobs - tcal
        rms = float(np.sqrt(np.mean(r ** 2)))
        if best is None or rms < best[1]:
            best = (m.copy(), rms, G, r)

    m, rms, G, r = best
    dof = max(len(r) - 4, 1)
    sigma2 = np.sum(r ** 2) / dof
    try:
        cov = sigma2 * np.linalg.inv(G.T @ G)
        err = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        err = np.full(4, np.nan)
    lat, lon = xy_to_latlon(m[0], m[1], lat0, lon0)
    return {"lat": float(lat), "lon": float(lon), "depth": float(m[2]),
            "origin_time": float(m[3]), "rms": rms,
            "err_x_km": err[0], "err_y_km": err[1], "err_z_km": err[2], "err_t0_s": err[3],
            "n_obs": len(tobs), "method": "Geiger (P+S)"}


# ============================ Robust grid search (ทรงกลม + โลก 2 ชั้น) ============================
def travel_time_2layer(delta_km, depth_km, phase="P"):
    """
    เวลาเดินทางในโลก 2 ชั้น (เปลือกโลก + เนื้อโลก) เลือกคลื่นที่มาถึงก่อนเอง
      คลื่นตรง  (Pg/Sg): t = R / V1,  R = sqrt(Δ² + h²)
      คลื่นหัว  (Pn/Sn): t = Δ / V2 + 2H·sqrt(1/V1² − 1/V2²)   (H = ความหนาเปลือกโลก)
      t = min(คลื่นตรง, คลื่นหัว)
    """
    v1, v2 = (C.VP_KM_S, C.VPN_KM_S) if phase == "P" else (C.VS_KM_S, C.VSN_KM_S)
    H = C.CRUST_THICKNESS_KM
    direct = np.sqrt(delta_km ** 2 + depth_km ** 2) / v1
    head = delta_km / v2 + 2 * H * np.sqrt(1 / v1 ** 2 - 1 / v2 ** 2)
    return np.minimum(direct, head)


def _robust_misfit(tobs, tcal):
    """T0 = มัธยฐานของ (tobs − tcal), misfit = มัธยฐานของ |residual|  (ทนต่อค่าผิดปกติ)"""
    d = tobs[..., :] - tcal
    t0 = np.median(d, axis=-1)
    res = d - t0[..., None]
    return np.median(np.abs(res), axis=-1), t0


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
                       outlier_sec=None, use_baz=False, baz_sec_per_deg=0.1, baz_min_rect=0.5,
                       relocate=True):
    """
    หาตำแหน่งด้วยการลองทุกจุดบนตาราง lat/lon (ระยะทางจริงบนผิวโลก + เวลาเดินทางโลก 2 ชั้น)
    stations : list ของ dict {"code", "lat", "lon", "tp", "ts"(ถ้ามี), "baz"(ถ้ามี), "rect"(ถ้ามี)}
    use_s    : ใช้เวลาคลื่น S ด้วยไหม
    use_baz  : ใช้ทิศทางจาก PCA ช่วยด้วย (แก้ปัญหาสถานีอยู่ฝั่งเดียวกัน)
               misfit = มัธยฐาน|residual เวลา| + baz_sec_per_deg × มัธยฐาน(มุมต่าง)
               ใช้เฉพาะสถานีที่ rectilinearity >= baz_min_rect
    relocate : ตัดเวลาคลื่นที่ไม่สอดคล้อง (outlier) แล้วคำนวณใหม่อีกรอบ
    คืนค่า dict: lat, lon, origin_time, misfit, residuals, outlier, ความไม่แน่นอน, azimuthal gap
    """
    from signal_utils import haversine_km, azimuth_deg
    obs = []                                  # (index สถานี, phase, เวลา)
    for k, s in enumerate(stations):
        if s.get("tp") is not None:
            obs.append((k, "P", s["tp"]))
        if use_s and s.get("ts") is not None:
            obs.append((k, "S", s["ts"]))
    if len(obs) < 4:
        raise ValueError("ต้องมีเวลาคลื่นอย่างน้อย 4 ค่า")
    tobs_all = np.array([o[2] for o in obs])
    slat_all = np.array([stations[o[0]]["lat"] for o in obs])
    slon_all = np.array([stations[o[0]]["lon"] for o in obs])
    isp_all = np.array([o[1] == "P" for o in obs])

    bz = [s for s in stations if use_baz and s.get("baz") is not None and s.get("rect", 1.0) >= baz_min_rect]
    blat = np.array([s["lat"] for s in bz]); blon = np.array([s["lon"] for s in bz])
    bbaz = np.array([s["baz"] for s in bz])

    def cost(LA, LO, keep):
        D = haversine_km(LA[..., None], LO[..., None], slat_all[keep], slon_all[keep])
        T = np.where(isp_all[keep], travel_time_2layer(D, depth_km, "P"), travel_time_2layer(D, depth_km, "S"))
        mis, t0 = _robust_misfit(tobs_all[keep], T)
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
        D = haversine_km(la, lo, slat_all, slon_all)
        tcal = np.where(isp_all, travel_time_2layer(D, depth_km, "P"), travel_time_2layer(D, depth_km, "S"))
        return tobs_all - t0 - tcal

    res = residuals(la, lo, t0)
    mad = np.median(np.abs(res)) * 1.4826
    thr = outlier_sec if outlier_sec else max(3 * mad, 3.0)
    out = np.abs(res) > thr
    n_pass = 1
    if relocate and out.any() and (~out).sum() >= 4:          # รอบที่ 2: ใช้เฉพาะเวลาที่สอดคล้อง
        la, lo, mis, t0, tot, unc = search(~out)
        res = residuals(la, lo, t0)
        out = np.abs(res) > thr
        n_pass = 2

    rows = [{"code": stations[o[0]]["code"], "phase": o[1], "residual": float(r), "outlier": bool(b)}
            for o, r, b in zip(obs, res, out)]
    baz_rows = []
    for s in bz:
        az = azimuth_deg(s["lat"], s["lon"], la, lo)
        baz_rows.append({"code": s["code"], "baz": float(s["baz"]), "az_to_epi": float(az),
                         "diff": float(_axial_diff_deg(az, s["baz"]))})
    return {"lat": la, "lon": lo, "depth": depth_km, "origin_time": t0,
            "misfit": mis, "cost": tot, "threshold": float(thr), "residuals": rows,
            "n_obs": len(obs), "uncertainty_km": unc, "passes": n_pass,
            "gap": azimuthal_gap(la, lo, stations), "baz_residuals": baz_rows,
            "method": "Robust grid search" + (" (P+S)" if use_s else " (P only)")
                      + (" + ทิศทาง PCA" if len(bz) else "")}


def delta_from_sp(sp_time, depth_km=10.0):
    """ระยะบนผิวโลก Δ จากเวลา S−P ด้วยโลก 2 ชั้น (ใช้ได้ทั้งสถานีใกล้และไกล ไม่ต้องเลือก profile)"""
    d = np.linspace(0, 4000, 8001)
    sp = travel_time_2layer(d, depth_km, "S") - travel_time_2layer(d, depth_km, "P")
    return float(np.interp(sp_time, sp, d))
