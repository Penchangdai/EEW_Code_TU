"""
map_plot.py
วาดผลลงบนแผนที่จริงด้วย folium: ค่าเริ่มต้นเป็นภาพดาวเทียม + เส้นแบ่งเขตประเทศ (คล้าย Google Earth)
และบันทึกไฟล์ .kml คู่กัน สำหรับเปิดใน Google Earth (ไฟล์ -> นำเข้าไฟล์ KML)
ผลลัพธ์เป็นไฟล์ .html เปิดด้วย browser แล้วซูม/คลิกดูรายละเอียดได้ (ต้องต่ออินเทอร์เน็ตตอนเปิดดู)

ติดตั้ง:  pip install folium

แผนที่ที่มี
  make_epicenter_map / make_multi_station_map : ตำแหน่งจุดศูนย์กลาง
  make_mmi_map        : ความรุนแรงมาตรา MMI แบบ USGS (สี ShakeMap) + การแจ้งเตือน   <- run_intensity.py ใช้
  make_station_mmi_map: MMI ที่วัดได้ของแต่ละสถานี บนพื้นสี MMI แบบ ShakeMap      <- run_intensity.py ใช้
  make_intensity_map  : ความรุนแรงมาตรา JMA (แบบเดิม)
  make_travel_map     : เวลาคลื่นเดินทางและเวลาเตือนล่วงหน้า
"""
import numpy as np

import html as _html
import os
import re

from signal_utils import destination_point
from station_db import load_station_db

# ======================================================================
#  แผนที่พื้นหลัง
#  ค่าเริ่มต้น = ภาพดาวเทียม + เส้นแบ่งเขตประเทศ/ชื่อประเทศ/ชื่อเมือง (คล้าย Google Earth)
#  เปลี่ยนได้ที่ปุ่มชั้นข้อมูลมุมขวาบนของแผนที่
# ======================================================================
ESRI = "https://server.arcgisonline.com/ArcGIS/rest/services/{}/MapServer/tile/{{z}}/{{y}}/{{x}}"
BASEMAP = "satellite"        # "satellite" (ดาวเทียม), "light" (สีเทาอ่อน แบบ ShakeMap), "street" (แผนที่ถนน)
#   หมายเหตุ: แผนที่ของ CARTO ตอนนี้ขึ้น "API KEY REQUIRED" ถ้าไม่มี key จึงเปลี่ยนมาใช้ของ Esri ทั้งหมด
#             (ใส่ "carto" ไว้ในโค้ดเดิม จะใช้ "light" แทนให้อัตโนมัติ)
SAVE_KML = True              # True = บันทึกไฟล์ .kml คู่กับ .html สำหรับเปิดใน Google Earth

BASE_LAYERS = {
    "satellite": ("World_Imagery", "Esri, Maxar, Earthstar Geographics", "ภาพดาวเทียม", 19),
    "light": ("Canvas/World_Light_Gray_Base", "Esri, HERE, Garmin, &copy; OpenStreetMap contributors",
              "แผนที่สีเทาอ่อน (Esri)", 16),
    "street": ("World_Street_Map", "Esri", "แผนที่ถนน (Esri)", 19),
    "topo": ("World_Topo_Map", "Esri", "แผนที่ภูมิประเทศ (Esri)", 19),
}
# ชั้นเส้นแบ่งเขตประเทศ/จังหวัด + ชื่อสถานที่ (พื้นโปร่งใส วางทับพื้นหลังและพื้นสี) ที่เข้ากับพื้นหลังแต่ละแบบ
#   ภาพดาวเทียม -> เส้นสีสว่าง,  แผนที่อื่น -> เส้นสีเทา (ทำให้เข้มขึ้นด้วย CSS ให้เห็นชัดบนพื้นสี MMI)
REFERENCE_LAYER = {"satellite": "Reference/World_Boundaries_and_Places",
                   "map": "Canvas/World_Light_Gray_Reference"}
REF_NAME = {"satellite": "เส้นแบ่งเขตประเทศ/จังหวัด + ชื่อ (สำหรับภาพดาวเทียม)",
            "map": "เส้นแบ่งเขตประเทศ/จังหวัด + ชื่อ (สำหรับแผนที่)"}


def _new_map(location, zoom_start=5, basemap=None):
    """
    สร้างแผนที่ folium พร้อมพื้นหลังหลายแบบ และชั้นเส้นแบ่งเขตประเทศ/จังหวัด + ชื่อสถานที่
    basemap : None = ใช้ค่า BASEMAP ด้านบน, หรือ "satellite" / "light" / "street" / "topo"
    ชั้นเส้นแบ่งเขตอยู่บนสุด (ไม่ถูกพื้นสีบัง) และสลับชุดให้อัตโนมัติเมื่อเปลี่ยนพื้นหลัง
    """
    import folium
    base = basemap or BASEMAP
    base = "light" if base == "carto" else base
    base = base if base in BASE_LAYERS else "satellite"
    m = folium.Map(location=location, zoom_start=zoom_start, tiles=None, control_scale=True)
    for k in [base] + [k for k in BASE_LAYERS if k != base]:     # ชั้นแรก = ค่าเริ่มต้น
        svc, attr, name, zmax = BASE_LAYERS[k]
        folium.TileLayer(tiles=ESRI.format(svc), attr=attr, name=name, max_zoom=zmax,
                         max_native_zoom=zmax, show=(k == base)).add_to(m)
    folium.map.CustomPane("labels", z_index=450, pointer_events=False).add_to(m)
    folium.map.CustomPane("labels_dark", z_index=451, pointer_events=False).add_to(m)
    ref = {}
    for kind, svc in REFERENCE_LAYER.items():
        ref[kind] = folium.TileLayer(tiles=ESRI.format(svc), attr="Esri", name=REF_NAME[kind], overlay=True,
                                     control=True, show=((kind == "satellite") == (base == "satellite")),
                                     max_zoom=16, max_native_zoom=16,
                                     pane="labels" if kind == "satellite" else "labels_dark")
        ref[kind].add_to(m)
    # เส้นสีเทาบนแผนที่: เพิ่มความเข้ม  และสลับชุดเส้นตามพื้นหลังที่เลือก
    sat_name = BASE_LAYERS["satellite"][2]
    mp = m.get_name()
    m.get_root().header.add_child(folium.Element(
        "<style>.leaflet-labels_dark-pane{filter:contrast(2.2) brightness(0.55);}</style>"))
    # (รอให้หน้าโหลดเสร็จ เพราะสคริปต์นี้ถูกเขียนก่อนตัวแปรแผนที่ของ folium)
    m.get_root().script.add_child(folium.Element(f"""
    window.addEventListener('load', function () {{ {mp}.on('baselayerchange', function (e) {{
        var sat = (e.name === '{sat_name}');
        var on = sat ? {ref['satellite'].get_name()} : {ref['map'].get_name()};
        var off = sat ? {ref['map'].get_name()} : {ref['satellite'].get_name()};
        if ({mp}.hasLayer(off)) {{ {mp}.removeLayer(off); }}
        if (!{mp}.hasLayer(on)) {{ {mp}.addLayer(on); }}
    }}); }});"""))
    return m


def _finish(m, out_html):
    """ใส่ปุ่มเลือกชั้นข้อมูล บันทึก html และ kml (Google Earth)"""
    import folium
    folium.LayerControl(collapsed=False).add_to(m)
    m.save(out_html)
    if SAVE_KML:
        try:
            kml = save_kml(m, os.path.splitext(out_html)[0] + ".kml")
            print(f"บันทึกไฟล์ Google Earth: {kml}")
        except Exception as exc:                       # KML เป็นส่วนเสริม ไม่ให้โปรแกรมหยุด
            print(f"(สร้างไฟล์ .kml ไม่สำเร็จ: {exc})")
    return out_html


# ======================================================================
#  ส่งออก KML สำหรับ Google Earth (แปลงจากวัตถุบนแผนที่ folium)
# ======================================================================
_NAMED = {"red": "#d63e2a", "darkred": "#a23336", "lightred": "#ff8e7f", "orange": "#f69730",
          "beige": "#ffcb92", "green": "#72b026", "darkgreen": "#728224", "lightgreen": "#bbf970",
          "blue": "#38aadd", "darkblue": "#0067a3", "lightblue": "#8adaff", "purple": "#d252b9",
          "darkpurple": "#5b396b", "pink": "#ff91ea", "cadetblue": "#436978", "white": "#ffffff",
          "gray": "#575757", "lightgray": "#a3a3a3", "black": "#000000"}


def _kml_color(color, alpha=1.0):
    """'#rrggbb' หรือชื่อสี -> 'aabbggrr' ของ KML"""
    c = _NAMED.get(str(color).lower(), color) if color else "#ff0000"
    c = str(c).lstrip("#")
    if len(c) == 3:
        c = "".join(ch * 2 for ch in c)
    if not re.fullmatch(r"[0-9a-fA-F]{6}", c):
        c = "ff0000"
    a = max(0, min(255, int(round(float(alpha) * 255))))
    return f"{a:02x}{c[4:6]}{c[2:4]}{c[0:2]}".lower()


def _text_of(obj, kind):
    """ดึงข้อความ tooltip / popup ของวัตถุ folium"""
    for ch in obj._children.values():
        name = type(ch).__name__
        if kind == "tooltip" and name == "Tooltip":
            return str(ch.text)
        if kind == "popup" and name == "Popup":
            for h in ch.html._children.values():
                data = getattr(h, "data", None)
                if data:
                    return str(data)
    return ""


def _plain(s):
    s = re.sub(r"<br\s*/?>", "\n", s or "")
    return _html.unescape(re.sub(r"<[^>]+>", "", s)).replace("&nbsp;", " ").strip()


def _circle_coords(lat, lon, radius_km, n=90):
    pts = [destination_point(lat, lon, az, radius_km) for az in np.linspace(0, 360, n + 1)]
    return " ".join(f"{p[1]:.5f},{p[0]:.5f},0" for p in pts)


def save_kml(m, out_kml):
    """แปลงหมุด วงกลม เส้น บนแผนที่ folium เป็นไฟล์ KML สำหรับเปิดใน Google Earth"""
    pm, styles = [], {}

    def style(key, body):
        styles.setdefault(key, f"<Style id='{key}'>{body}</Style>")
        return key

    def place(name, desc, geom, sid):
        pm.append(f"<Placemark><name>{_html.escape(name)}</name>"
                  f"<description>{_html.escape(desc)}</description>"
                  f"<styleUrl>#{sid}</styleUrl>{geom}</Placemark>")

    def walk(obj):
        kind = type(obj).__name__
        if kind == "FeatureGroup":
            start = len(pm)
            for ch in list(obj._children.values()):
                walk(ch)
            if len(pm) > start:
                pm.insert(start, f"<Folder><name>{_html.escape(obj.layer_name)}</name>")
                pm.append("</Folder>")
            return
        tip = _plain(_text_of(obj, "tooltip"))
        desc = _plain(_text_of(obj, "popup")) or tip
        o = getattr(obj, "options", {}) or {}
        if kind == "Marker":
            icon = next((c for c in obj._children.values() if type(c).__name__ in ("Icon", "DivIcon")), None)
            if icon is None or type(icon).__name__ == "DivIcon":
                return                                   # ป้ายข้อความบนแผนที่ ไม่ต้องส่งออก
            io = icon.options
            col = io.get("marker_color") or io.get("markerColor") or "blue"
            sid = style(f"mk_{col}", f"<IconStyle><color>{_kml_color(col)}</color><scale>1.2</scale>"
                        "<Icon><href>http://maps.google.com/mapfiles/kml/pushpin/wht-pushpin.png</href></Icon>"
                        "</IconStyle>")
            lat, lon = obj.location
            place(tip.split(":")[0] or "จุด", desc, f"<Point><coordinates>{lon},{lat},0</coordinates></Point>", sid)
        elif kind == "CircleMarker":
            col = o.get("fillColor") or o.get("color")
            big = o.get("radius", 5) >= 6                 # จุดใหญ่ = แสดงชื่อ, จุดเล็ก (สถานีทั้งหมด) = ซ่อนชื่อ
            sid = style(f"dot_{_kml_color(col)}_{int(big)}",
                        f"<IconStyle><color>{_kml_color(col)}</color><scale>{1.0 if big else 0.6}</scale>"
                        "<Icon><href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href></Icon>"
                        f"</IconStyle><LabelStyle><scale>{0.9 if big else 0}</scale></LabelStyle>")
            lat, lon = obj.location
            place(tip.split(" (")[0] or "จุด", desc, f"<Point><coordinates>{lon},{lat},0</coordinates></Point>", sid)
        elif kind in ("Circle", "Rectangle", "PolyLine"):
            col = o.get("color", "red")
            fill = bool(o.get("fill", False))
            fcol = o.get("fillColor") or col
            fa = o.get("fillOpacity", 0.2) if fill else 0.0
            lc, pc = _kml_color(col, o.get("opacity", 1.0)), _kml_color(fcol, fa)
            sid = style(f"ln_{lc}_{pc}_{o.get('weight', 2)}",
                        f"<LineStyle><color>{lc}</color><width>{o.get('weight', 2)}</width></LineStyle>"
                        f"<PolyStyle><color>{pc}</color><fill>{int(fill)}</fill></PolyStyle>")
            if kind == "PolyLine":
                coords = " ".join(f"{p[1]:.5f},{p[0]:.5f},0" for p in obj.locations)
                place(tip or "เส้น", desc, "<LineString><tessellate>1</tessellate>"
                                          f"<coordinates>{coords}</coordinates></LineString>", sid)
                return
            if kind == "Circle":
                coords = _circle_coords(obj.location[0], obj.location[1], o["radius"] / 1000.0)
            else:
                (a, b), (c, d) = obj.locations
                coords = " ".join(f"{x},{y},0" for y, x in ((a, b), (a, d), (c, d), (c, b), (a, b)))
            place(tip or "พื้นที่", desc, "<Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing>"
                                       f"<coordinates>{coords}</coordinates></LinearRing></outerBoundaryIs></Polygon>",
                  sid)

    for ch in list(m._children.values()):
        walk(ch)
    name = os.path.splitext(os.path.basename(out_kml))[0]
    doc = ("<?xml version='1.0' encoding='UTF-8'?>\n<kml xmlns='http://www.opengis.net/kml/2.2'><Document>"
           f"<name>{_html.escape(name)}</name>" + "".join(styles.values()) + "".join(pm) + "</Document></kml>")
    with open(out_kml, "w", encoding="utf-8") as f:
        f.write(doc)
    return out_kml


def _great_circle(lat, lon, az_deg, dist_km, n=60):
    """จุดตามเส้นทางจาก (lat, lon) ไปทิศ az เป็นระยะ dist_km (ใช้วาดเส้นบนแผนที่)"""
    return [destination_point(lat, lon, az_deg, d) for d in np.linspace(0, dist_km, n)]


def make_epicenter_map(station, sta_lat, sta_lon, delta_km, candidates, out_html,
                       ref_epicenter=None, show_all_stations=True):
    """
    station       : รหัสสถานีที่บันทึกข้อมูล
    sta_lat/lon   : พิกัดสถานี
    delta_km      : ระยะทางบนผิวโลกจากสถานีถึงจุดศูนย์กลาง (km)
    candidates    : list ของ dict {"name", "lat", "lon", "az", "color"}
    ref_epicenter : (lat, lon, "ชื่อ") จุดศูนย์กลางอ้างอิง เช่น จาก USGS (ไม่ใส่ก็ได้)
    คืนค่า path ของไฟล์ html หรือ None ถ้าไม่ได้ติดตั้ง folium
    """
    try:
        import folium
    except ImportError:
        print("(ไม่ได้ติดตั้ง folium จึงไม่วาดแผนที่จริง: pip install folium)")
        return None

    m = _new_map([sta_lat, sta_lon], zoom_start=5)

    # ---- สถานีทั้งหมดใน stations.csv ----
    if show_all_stations:
        layer = folium.FeatureGroup(name="สถานีใน stations.csv")
        for code, (la, lo) in load_station_db().items():
            folium.CircleMarker([la, lo], radius=3, color="#555555", fill=True, fill_opacity=0.8,
                                weight=1, tooltip=f"{code} ({la:.4f}, {lo:.4f})").add_to(layer)
        layer.add_to(m)

    # ---- วงกลมระยะทาง Δ (ระยะจริงเป็นเมตร) ----
    folium.Circle([sta_lat, sta_lon], radius=delta_km * 1000, color="black", weight=1.5,
                  dash_array="6 6", fill=False,
                  tooltip=f"ระยะทางจากสถานี Δ = {delta_km:.0f} km").add_to(m)

    # ---- สถานีที่บันทึกข้อมูล ----
    folium.Marker([sta_lat, sta_lon], tooltip=f"สถานี {station}",
                  popup=f"<b>สถานี {station}</b><br>lat {sta_lat:.4f}<br>lon {sta_lon:.4f}",
                  icon=folium.Icon(color="blue", icon="signal")).add_to(m)

    # ---- ตำแหน่งที่เป็นไปได้ + เส้นทิศทางจากสถานี ----
    points = [(sta_lat, sta_lon)]
    for c in candidates:
        folium.PolyLine(_great_circle(sta_lat, sta_lon, c["az"], delta_km), color=c["color"],
                        weight=2, opacity=0.8, tooltip=f"ทิศ {c['az']:.0f}°").add_to(m)
        folium.Marker([c["lat"], c["lon"]], tooltip=c["name"],
                      popup=(f"<b>{c['name']}</b><br>lat {c['lat']:.3f}<br>lon {c['lon']:.3f}"
                             f"<br>ทิศจากสถานี {c['az']:.1f}°<br>ระยะ {delta_km:.0f} km"),
                      icon=folium.Icon(color=c["color"], icon="star")).add_to(m)
        points.append((c["lat"], c["lon"]))

    # ---- จุดศูนย์กลางอ้างอิง (ถ้ามี) ----
    if ref_epicenter is not None:
        la, lo, label = ref_epicenter
        folium.Marker([la, lo], tooltip=label, popup=f"<b>{label}</b><br>lat {la:.3f}<br>lon {lo:.3f}",
                      icon=folium.Icon(color="green", icon="flag")).add_to(m)
        points.append((la, lo))

    # ---- หัวเรื่อง + ซูมให้เห็นทุกจุด ----
    title = (f"<div style='position:fixed;top:10px;left:50px;z-index:9999;background:white;"
             f"padding:6px 10px;border-radius:6px;box-shadow:0 1px 4px #888;font-size:14px'>"
             f"<b>ตำแหน่งจุดศูนย์กลางที่เป็นไปได้ จากสถานี {station}</b><br>"
             f"Δ = {delta_km:.0f} km &nbsp;|&nbsp; น้ำเงิน = สถานี, ดาว = ตำแหน่งที่เป็นไปได้"
             f"{', ธงเขียว = จุดอ้างอิง' if ref_epicenter else ''}</div>")
    m.get_root().html.add_child(folium.Element(title))
    lats = [p[0] for p in points]
    lons = [p[1] for p in points]
    m.fit_bounds([[min(lats) - 1, min(lons) - 1], [max(lats) + 1, max(lons) + 1]])
    return _finish(m, out_html)


COLORS = ["blue", "darkred", "darkgreen", "orange", "cadetblue"]   # สีของแต่ละสถานี (สูงสุด 5)


def make_multi_station_map(results, out_html, combined=None, geiger=None, ref_epicenter=None,
                           show_all_stations=True, subtitle=None, combined_label="ตำแหน่งรวม"):
    """
    วาดผลจากหลายสถานีลงแผนที่จริง
    results  : list ของ dict ต่อสถานี {"code", "lat", "lon", "delta", "cands":[{"name","lat","lon","az"}], "chosen"}
               ("cands" ไม่ใส่ก็ได้ = วาดแค่สถานีและวงกลม, "delta" = None = ไม่มีวงกลม)
               chosen = index ของตัวเลือก (0=A, 1=B) ที่สอดคล้องกับสถานีอื่น หรือ None
    combined : (lat, lon) ตำแหน่งเฉลี่ยจากตัวเลือกที่สอดคล้องกัน
    geiger   : (lat, lon, depth) ผลจาก Geiger (ถ้ามี)
    """
    try:
        import folium
    except ImportError:
        print("(ไม่ได้ติดตั้ง folium จึงไม่วาดแผนที่จริง: pip install folium)")
        return None

    m = _new_map([results[0]["lat"], results[0]["lon"]], zoom_start=5)

    if show_all_stations:
        layer = folium.FeatureGroup(name="สถานีใน stations.csv")
        for code, (la, lo) in load_station_db().items():
            folium.CircleMarker([la, lo], radius=3, color="#555555", fill=True, fill_opacity=0.8,
                                weight=1, tooltip=f"{code} ({la:.4f}, {lo:.4f})").add_to(layer)
        layer.add_to(m)

    points = []
    for k, r in enumerate(results):
        col = COLORS[k % len(COLORS)]
        layer = folium.FeatureGroup(name=f"สถานี {r['code']}")
        bad = r.get("outlier", False)          # สถานีที่ถูกตัดออก -> วงกลมสีเทาจาง
        if r.get("delta"):          # สถานีที่ไม่มีคลื่น S จะไม่มีวงกลม
            folium.Circle([r["lat"], r["lon"]], radius=r["delta"] * 1000,
                          color="#999999" if bad else col, weight=1 if bad else 1.5,
                          opacity=0.4 if bad else 1.0, dash_array="2 6" if bad else "6 6", fill=False,
                          tooltip=f"{r['code']}: Δ = {r['delta']:.0f} km"
                                  + ("  (ตัดออก: เวลาคลื่นไม่สอดคล้อง)" if bad else "")).add_to(layer)
        delta_txt = f"<br>Δ = {r['delta']:.0f} km" if r.get("delta") else "<br>ไม่มีคลื่น S"
        folium.Marker([r["lat"], r["lon"]], tooltip=f"สถานี {r['code']}" + (" (ตัดออก)" if bad else ""),
                      popup=f"<b>สถานี {r['code']}</b><br>lat {r['lat']:.4f}<br>lon {r['lon']:.4f}"
                            f"{delta_txt}",
                      icon=folium.Icon(color="lightgray" if bad else col, icon="signal")).add_to(layer)
        points.append((r["lat"], r["lon"]))
        for j, c in enumerate(r.get("cands", [])):
            chosen = r.get("chosen") == j
            folium.PolyLine(_great_circle(r["lat"], r["lon"], c["az"], r["delta"]), color=col,
                            weight=2.5 if chosen else 1, opacity=0.9 if chosen else 0.4,
                            dash_array=None if chosen else "4 6",
                            tooltip=f"{r['code']} {c['name']}: ทิศ {c['az']:.0f}°").add_to(layer)
            folium.CircleMarker([c["lat"], c["lon"]], radius=7 if chosen else 4, color=col,
                                fill=True, fill_opacity=0.9 if chosen else 0.3,
                                tooltip=f"{r['code']} {c['name']} ({c['lat']:.3f}, {c['lon']:.3f})"
                                        + ("  ← สอดคล้อง" if chosen else "")).add_to(layer)
            points.append((c["lat"], c["lon"]))
        layer.add_to(m)

    if combined is not None:
        folium.Marker(list(combined[:2]), tooltip=combined_label,
                      popup=f"<b>{combined_label}</b><br>lat {combined[0]:.3f}<br>lon {combined[1]:.3f}",
                      icon=folium.Icon(color="red", icon="star")).add_to(m)
        points.append(tuple(combined[:2]))
    if geiger is not None:
        folium.Marker([geiger[0], geiger[1]], tooltip="ผลจาก Geiger",
                      popup=f"<b>Geiger</b><br>lat {geiger[0]:.3f}<br>lon {geiger[1]:.3f}"
                            f"<br>ลึก {geiger[2]:.1f} km",
                      icon=folium.Icon(color="purple", icon="screenshot")).add_to(m)
        points.append((geiger[0], geiger[1]))
    if ref_epicenter is not None:
        la, lo, label = ref_epicenter
        folium.Marker([la, lo], tooltip=label, popup=f"<b>{label}</b><br>lat {la:.3f}<br>lon {lo:.3f}",
                      icon=folium.Icon(color="green", icon="flag")).add_to(m)
        points.append((la, lo))

    codes = ", ".join(r["code"] for r in results) if len(results) <= 10 else "ดูรายชื่อที่ปุ่มมุมขวาบน"
    title = (f"<div style='position:fixed;top:10px;left:50px;z-index:9999;background:white;"
             f"padding:6px 10px;border-radius:6px;box-shadow:0 1px 4px #888;font-size:14px'>"
             f"<b>ตำแหน่งจุดศูนย์กลางจาก {len(results)} สถานี: {codes}</b><br>"
             + (subtitle or (f"วงกลม = ระยะ Δ ของแต่ละสถานี, จุดทึบ = ตัวเลือกที่สอดคล้องกัน, "
                             f"ดาวแดง = ตำแหน่งรวม{', เป้าม่วง = Geiger' if geiger else ''}")) +
             f"{', ธงเขียว = จุดอ้างอิง' if ref_epicenter else ''}</div>")
    m.get_root().html.add_child(folium.Element(title))
    lats = [p[0] for p in points]
    lons = [p[1] for p in points]
    m.fit_bounds([[min(lats) - 1, min(lons) - 1], [max(lats) + 1, max(lons) + 1]])
    return _finish(m, out_html)


# สีตามระดับความรุนแรง JMA (ใกล้เคียงสีที่ JMA ใช้)
JMA_CLASS_COLORS = {"0": "#f2f2ff", "1": "#00aaff", "2": "#0041ff", "3": "#fae696", "4": "#ffe600",
                    "5-": "#ff9900", "5+": "#ff2800", "6-": "#a50021", "6+": "#a50021", "7": "#b40068"}


def make_intensity_map(points, out_html, epicenter=None, title="ความรุนแรงแผ่นดินไหว (มาตรา JMA)"):
    """
    points    : list ของ dict {"name", "lat", "lon", "I", "cls", "kind", "info"}
                kind = "obs" (วัดได้ที่สถานี, วงกลม) หรือ "pred" (ทำนาย, สี่เหลี่ยม)
    epicenter : (lat, lon) ถ้ามี
    """
    try:
        import folium
    except ImportError:
        print("(ไม่ได้ติดตั้ง folium จึงไม่วาดแผนที่จริง: pip install folium)")
        return None
    m = _new_map([points[0]["lat"], points[0]["lon"]], zoom_start=6)
    lats, lons = [], []
    for p in points:
        col = JMA_CLASS_COLORS.get(p["cls"], "#888888")
        label = f"{p['name']}: I = {p['I']:.1f} (ระดับ {p['cls']})"
        popup = f"<b>{p['name']}</b><br>{'วัดได้' if p['kind'] == 'obs' else 'ทำนาย'}: I = {p['I']:.1f} " \
                f"(ระดับ {p['cls']})<br>{p.get('info', '')}"
        if p["kind"] == "obs":
            folium.CircleMarker([p["lat"], p["lon"]], radius=9, color="black", weight=1, fill=True,
                                fill_color=col, fill_opacity=0.9, tooltip=label, popup=popup).add_to(m)
        else:
            d = 0.12
            folium.Rectangle([[p["lat"] - d, p["lon"] - d], [p["lat"] + d, p["lon"] + d]], color="black",
                             weight=1, fill=True, fill_color=col, fill_opacity=0.9,
                             tooltip=label + " (ทำนาย)", popup=popup).add_to(m)
        folium.Marker([p["lat"], p["lon"]], icon=folium.DivIcon(
            html=f"<div style='font-size:11px;font-weight:bold;transform:translate(10px,-8px);"
                 f"white-space:nowrap;text-shadow:0 0 3px white'>{p['name']} {p['cls']}</div>")).add_to(m)
        lats.append(p["lat"]); lons.append(p["lon"])
    if epicenter is not None:
        folium.Marker(list(epicenter[:2]), tooltip="จุดศูนย์กลาง",
                      icon=folium.Icon(color="red", icon="star")).add_to(m)
        lats.append(epicenter[0]); lons.append(epicenter[1])
    legend = "".join(f"<span style='display:inline-block;width:14px;height:14px;background:{c};"
                     f"border:1px solid #333;margin:0 3px 0 8px;vertical-align:middle'></span>{k}"
                     for k, c in JMA_CLASS_COLORS.items())
    m.get_root().html.add_child(folium.Element(
        f"<div style='position:fixed;top:10px;left:50px;z-index:9999;background:white;padding:6px 10px;"
        f"border-radius:6px;box-shadow:0 1px 4px #888;font-size:13px'><b>{title}</b><br>"
        f"วงกลม = วัดได้ที่สถานี, สี่เหลี่ยม = ทำนาย, ดาวแดง = จุดศูนย์กลาง<br>{legend}</div>"))
    m.fit_bounds([[min(lats) - 1, min(lons) - 1], [max(lats) + 1, max(lons) + 1]])
    return _finish(m, out_html)


# สีระดับ MMI มาตรฐานของ USGS ShakeMap
MMI_COLORS = {1: "#ffffff", 2: "#bfccff", 3: "#a0e6ff", 4: "#80ffff", 5: "#7aff93", 6: "#ffff00",
              7: "#ffc800", 8: "#ff9100", 9: "#ff0000", 10: "#c80000"}
MMI_ROMAN = {1: "I", 2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI", 7: "VII", 8: "VIII", 9: "IX", 10: "X+"}
MMI_TEXT = {1: ("ไม่รู้สึก", "ไม่มี"), 2: ("อ่อน", "ไม่มี"), 3: ("อ่อน", "ไม่มี"), 4: ("เบา", "ไม่มี"),
            5: ("ปานกลาง", "น้อยมาก"), 6: ("แรง", "เล็กน้อย"), 7: ("แรงมาก", "ปานกลาง"),
            8: ("รุนแรง", "ปานกลางถึงหนัก"), 9: ("รุนแรงมาก", "หนัก"), 10: ("รุนแรงที่สุด", "หนักมาก")}


def mmi_level(mmi):
    """MMI (ทศนิยม) -> ระดับจำนวนเต็ม 1-10 (ค่าที่ปัดแล้ว แบบ USGS)"""
    return int(min(max(np.floor(float(mmi) + 0.5), 1), 10))


def _mmi_label(p):
    """ข้อความป้าย: MMI วัดได้ และ MMI ทำนาย (เลขโรมัน + ค่า)"""
    k, obs = mmi_level(p["mmi"]), p.get("obs")
    pred = f"ทำนาย {MMI_ROMAN[k]} ({p['mmi']:.1f})"
    if obs is None:
        return f"MMI {pred}"
    return f"MMI วัดได้ {MMI_ROMAN[mmi_level(obs)]} ({obs:.1f}) · {pred}"


def _station_table(points):
    """ตาราง MMI ของทุกสถานีที่ input เข้าไป (+ ตำแหน่งที่แจ้งเตือน)"""
    rows = []
    for p in sorted(points, key=lambda q: (q.get("kind") != "station", q.get("delta", 0))):
        obs = p.get("obs")
        ko, kp = (mmi_level(obs) if obs is not None else None), mmi_level(p["mmi"])
        cell = lambda k, v: (f"<td style='background:{MMI_COLORS[k]};text-align:center;padding:1px 6px'>"
                             f"<b>{MMI_ROMAN[k]}</b> {v:.1f}</td>")
        lead = p.get("lead")
        alert = {2: "สาธารณะ", 1: "แอป"}.get(p.get("alert", 0), "-")
        if lead is not None and p.get("alert", 0):
            alert += f" ({lead:+.0f} s)"
        name = p["name"] + ("" if p.get("kind") != "station" or p.get("used", True) else " *")
        rows.append(f"<tr><td style='padding:1px 6px'>{'&#9733; ' if p.get('kind') != 'station' else ''}{name}</td>"
                    f"<td style='text-align:right;padding:1px 6px'>{p.get('delta', 0):,.0f}</td>"
                    + (cell(ko, obs) if obs is not None else "<td style='text-align:center'>-</td>")
                    + cell(kp, p["mmi"]) + f"<td style='padding:1px 6px'>{alert}</td></tr>")
    foot = ("<div style='font-size:11px;color:#555'>* input เข้าไปแต่ไม่ได้ใช้คำนวณขนาด, "
            "เวลาในวงเล็บ = เตือนก่อนคลื่น S</div>") if any(
        p.get("kind") == "station" and not p.get("used", True) for p in points) else (
        "<div style='font-size:11px;color:#555'>เวลาในวงเล็บ = เตือนก่อนคลื่น S</div>")
    return ("<table style='border-collapse:collapse;font-size:12px;margin-top:4px'>"
            "<tr style='border-bottom:1px solid #999'><th style='text-align:left;padding:1px 6px'>สถานี</th>"
            "<th style='padding:1px 6px'>Δ km</th><th style='padding:1px 6px'>MMI วัดได้</th>"
            "<th style='padding:1px 6px'>MMI ทำนาย</th><th style='padding:1px 6px'>แจ้งเตือน</th></tr>"
            + "".join(rows) + "</table>" + foot)


def make_mmi_map(points, out_html, epicenter, rings=None, title="ความรุนแรงแผ่นดินไหว (MMI)", note=""):
    """
    แผนที่ความรุนแรงหน่วย MMI (สีมาตรฐาน USGS ShakeMap)
    points    : list ของ dict
                {"name", "lat", "lon", "kind": "station" | "target",
                 "mmi": MMI ทำนาย, "obs": MMI วัดได้ หรือ None, "alert": 0/1/2 (ไม่แจ้ง/แอป/สาธารณะ),
                 "popup": ข้อความ html,  (ไม่ใส่ก็ได้:) "delta" km, "lead" s, "used" True/False}
    ป้ายข้างสถานีบอก MMI วัดได้ + ทำนาย และมีตาราง MMI ของทุกสถานีที่มุมซ้ายล่าง
    epicenter : dict {"lat", "lon", "popup"}
    rings     : list ของ (ระดับ MMI, รัศมี km) ขอบเขตที่คาดว่าสั่นถึงแต่ละระดับ
    คืนค่า path ของไฟล์ html หรือ None ถ้าไม่ได้ติดตั้ง folium
    """
    try:
        import folium
    except ImportError:
        print("(ไม่ได้ติดตั้ง folium จึงไม่วาดแผนที่จริง: pip install folium)")
        return None
    m = _new_map([epicenter["lat"], epicenter["lon"]], zoom_start=5)
    lats, lons = [epicenter["lat"]], [epicenter["lon"]]

    # ---- วงขอบเขตระดับ MMI (วงใหญ่ก่อน วงเล็กทับด้านบน) ----
    if rings:
        layer = folium.FeatureGroup(name="ขอบเขตระดับ MMI ที่คาด")
        for k, r in sorted(rings, key=lambda x: -x[1]):
            folium.Circle([epicenter["lat"], epicenter["lon"]], radius=r * 1000, color="#555555", weight=1,
                          dash_array="4 4", fill=True, fill_color=MMI_COLORS[k], fill_opacity=0.22,
                          tooltip=f"MMI {MMI_ROMAN[k]} ({MMI_TEXT[k][0]}) ถึงระยะประมาณ {r:,.0f} km").add_to(layer)
        layer.add_to(m)

    # ---- จุดศูนย์กลาง ----
    folium.Marker([epicenter["lat"], epicenter["lon"]], tooltip="จุดศูนย์กลาง",
                  popup=folium.Popup(epicenter.get("popup", "จุดศูนย์กลาง"), max_width=320),
                  icon=folium.Icon(color="red", icon="star")).add_to(m)

    # ---- สถานี / ตำแหน่งที่แจ้งเตือน ----
    st_layer = folium.FeatureGroup(name="สถานี (วงนอก = ทำนาย, จุดใน = วัดได้)")
    tg_layer = folium.FeatureGroup(name="ตำแหน่งที่แจ้งเตือน")
    for p in points:
        k = mmi_level(p["mmi"])
        alerted = p.get("alert", 0) > 0
        obs = p.get("obs")
        label = f"{p['name']}: MMI {MMI_ROMAN[k]} ({MMI_TEXT[k][0]})" + (
            f", วัดได้ {MMI_ROMAN[mmi_level(obs)]}" if obs is not None else "")
        popup = folium.Popup(f"<b>{p['name']}</b><br>{p.get('popup', '')}", max_width=360)
        layer = st_layer if p.get("kind") == "station" else tg_layer
        folium.CircleMarker([p["lat"], p["lon"]], radius=13 if p.get("kind") == "station" else 15,
                            color="black" if alerted else "#888888", weight=3 if alerted else 1,
                            fill=True, fill_color=MMI_COLORS[k], fill_opacity=0.95,
                            tooltip=label, popup=popup).add_to(layer)
        if obs is not None:
            folium.CircleMarker([p["lat"], p["lon"]], radius=6, color="#333333", weight=1, fill=True,
                                fill_color=MMI_COLORS[mmi_level(obs)], fill_opacity=1.0,
                                tooltip=f"{p['name']}: วัดได้ MMI {MMI_ROMAN[mmi_level(obs)]} ({obs:.1f})",
                                popup=popup).add_to(layer)
        if p.get("kind") != "station":
            folium.Marker([p["lat"], p["lon"]], icon=folium.DivIcon(
                html="<div style='font-size:18px;line-height:18px;transform:translate(-6px,-9px)'>&#9733;</div>"
            )).add_to(layer)
        folium.Marker([p["lat"], p["lon"]], icon=folium.DivIcon(
            html=f"<div style='font-size:11px;font-weight:bold;transform:translate(17px,-14px);white-space:nowrap;"
                 f"line-height:1.25;background:rgba(255,255,255,.85);padding:1px 4px;border-radius:3px;"
                 f"border-left:4px solid {MMI_COLORS[mmi_level(obs if obs is not None else p['mmi'])]}'>"
                 f"{p['name']}<br>" + _mmi_label(p) + "</div>")).add_to(layer)
        lats.append(p["lat"]); lons.append(p["lon"])
    st_layer.add_to(m)
    tg_layer.add_to(m)

    # ---- หัวเรื่อง + คำอธิบายสี ----
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;top:10px;left:50px;z-index:9999;background:white;padding:6px 10px;"
        "border-radius:6px;box-shadow:0 1px 4px #888;font-size:13px;max-width:460px'>"
        f"<b>{title}</b><br>{note}<br>"
        "สถานี: วงนอก = MMI ทำนาย, จุดใน = MMI วัดได้, ขอบดำ = ได้รับคำเตือน<br>"
        "ดาว = ตำแหน่งที่แจ้งเตือน, ดาวแดง = จุดศูนย์กลาง</div>"))
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;bottom:25px;left:10px;z-index:9999;background:white;padding:6px 10px;"
        "border-radius:6px;box-shadow:0 1px 4px #888;font-size:12px;max-height:45%;overflow:auto'>"
        "<b>MMI ของแต่ละสถานี</b>" + _station_table(points) + "</div>"))
    legend = "".join(
        f"<div><span style='display:inline-block;width:18px;height:12px;background:{MMI_COLORS[k]};"
        f"border:1px solid #777;margin-right:6px;vertical-align:middle'></span>"
        f"<b>{'II–III' if k == 3 else MMI_ROMAN[k]}</b> {MMI_TEXT[k][0]} · เสียหาย: {MMI_TEXT[k][1]}</div>"
        for k in range(1, 11) if k != 2)
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;bottom:25px;right:10px;z-index:9999;background:white;padding:6px 10px;"
        f"border-radius:6px;box-shadow:0 1px 4px #888;font-size:12px'><b>มาตรา MMI (USGS)</b>{legend}</div>"))
    m.fit_bounds([[min(lats) - 1, min(lons) - 1], [max(lats) + 1, max(lons) + 1]])
    return _finish(m, out_html)


def make_station_mmi_map(stations, out_html, epicenter, field=None, title="ระดับความรุนแรงแผ่นดินไหว",
                         subtitle="ที่สถานีตรวจวัด (หน่วยวัด MMI)", basemap="light"):
    """
    แผนที่ MMI ที่วัดได้ของแต่ละสถานี (ตัวเลขใหญ่ + ชื่อสถานี) บนพื้นสี MMI แบบ ShakeMap
    stations : list ของ dict {"name", "lat", "lon", "mmi": MMI วัดได้ (หรือ None), "popup"}
    epicenter: dict {"lat", "lon", "popup"}
    field    : dict {"uri": ภาพ PNG (data URI), "bounds": [[lat_s, lon_w], [lat_n, lon_e]]}  พื้นสี MMI ที่คาด
    """
    try:
        import folium
    except ImportError:
        print("(ไม่ได้ติดตั้ง folium จึงไม่วาดแผนที่จริง: pip install folium)")
        return None
    m = _new_map([epicenter["lat"], epicenter["lon"]], zoom_start=6, basemap=basemap)
    if field:
        folium.raster_layers.ImageOverlay(image=field["uri"], bounds=field["bounds"], opacity=0.65,
                                          name="พื้นสี MMI ที่คาด (ShakeMap)", interactive=False,
                                          zindex=1).add_to(m)
    folium.Marker([epicenter["lat"], epicenter["lon"]], tooltip="จุดศูนย์กลาง",
                  popup=folium.Popup(epicenter.get("popup", "จุดศูนย์กลาง"), max_width=320),
                  icon=folium.Icon(color="red", icon="star")).add_to(m)
    layer = folium.FeatureGroup(name="สถานี (MMI วัดได้)")
    lats, lons = [epicenter["lat"]], [epicenter["lon"]]
    for p in stations:
        v = p.get("mmi")
        col = MMI_COLORS[mmi_level(v)] if v is not None else "#bbbbbb"
        txt = f"{v:.1f}" if v is not None else "-"
        tip = f"{p['name']}: MMI {txt}" + (f" ({MMI_ROMAN[mmi_level(v)]} {MMI_TEXT[mmi_level(v)][0]})"
                                           if v is not None else " (ไม่มีข้อมูล)")
        folium.CircleMarker([p["lat"], p["lon"]], radius=7, color="#333333", weight=1.5, fill=True,
                            fill_color=col, fill_opacity=1.0, tooltip=tip,
                            popup=folium.Popup(f"<b>{p['name']}</b><br>{p.get('popup', '')}", max_width=360)
                            ).add_to(layer)
        folium.Marker([p["lat"], p["lon"]], icon=folium.DivIcon(
            icon_size=(160, 60), icon_anchor=(80, 62),
            html=("<div style='width:160px;text-align:center;font-family:Tahoma,sans-serif;"
                  "text-shadow:0 0 3px #fff,0 0 3px #fff,0 0 6px #fff;line-height:1.05'>"
                  f"<div style='font-size:26px;font-weight:bold;color:#333'>{txt}</div>"
                  f"<div style='font-size:15px;color:#222'>{p['name']}</div></div>"))).add_to(layer)
        lats.append(p["lat"]); lons.append(p["lon"])
    layer.add_to(m)
    legend = "".join(f"<td style='background:{MMI_COLORS[k]};width:34px;text-align:center;border:1px solid #999'>"
                     f"{MMI_ROMAN[k]}</td>" for k in range(1, 11))
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;top:10px;left:50%;transform:translateX(-50%);z-index:9999;text-align:center;"
        "font-family:Tahoma,sans-serif;text-shadow:0 0 4px #fff,0 0 4px #fff;pointer-events:none'>"
        f"<div style='font-size:28px;font-weight:bold;color:#e02020'>{title}</div>"
        f"<div style='font-size:17px;color:#222'>{subtitle}</div></div>"))
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;bottom:25px;left:10px;z-index:9999;background:rgba(255,255,255,.9);"
        "padding:6px 10px;border-radius:6px;box-shadow:0 1px 4px #888;font-size:12px'>"
        "<b>หน่วยวัด: MMI (Modified Mercalli Intensity)</b>"
        f"<table style='border-collapse:collapse;margin-top:4px'><tr>{legend}</tr></table>"
        "ตัวเลข = MMI ที่วัดได้ที่สถานี (จาก PGV, Worden et al. 2012)<br>"
        "พื้นสี = MMI ที่คาดจากขนาดและระยะ (ปรับด้วยค่าที่สถานีวัดได้)</div>"))
    m.fit_bounds([[min(lats) - 1, min(lons) - 1], [max(lats) + 1, max(lons) + 1]])
    return _finish(m, out_html)


def make_travel_map(epicenter, sites, out_html, blind_km=None, stations=None, s_fronts=None):
    """
    epicenter : (lat, lon, depth)
    sites     : list ของ dict {"name", "lat", "lon", "lead", "info"}  lead = เวลาเตือนล่วงหน้า (s) หรือ None
    blind_km  : รัศมี blind zone (km)
    stations  : list ของ dict {"code", "lat", "lon"} สถานีที่ใช้ตรวจจับ
    s_fronts  : list ของ (วินาทีหลังเกิดเหตุ, รัศมี km) แนวหน้าคลื่น S
    """
    try:
        import folium
    except ImportError:
        print("(ไม่ได้ติดตั้ง folium จึงไม่วาดแผนที่จริง: pip install folium)")
        return None
    m = _new_map(list(epicenter[:2]), zoom_start=5)
    lats, lons = [epicenter[0]], [epicenter[1]]
    if s_fronts:
        layer = folium.FeatureGroup(name="แนวหน้าคลื่น S")
        for t, r in s_fronts:
            folium.Circle(list(epicenter[:2]), radius=r * 1000, color="#666666", weight=1, dash_array="4 6",
                          fill=False, tooltip=f"คลื่น S ที่ {t:.0f} s หลังเกิดเหตุ (รัศมี {r:.0f} km)").add_to(layer)
        layer.add_to(m)
    if blind_km:
        folium.Circle(list(epicenter[:2]), radius=blind_km * 1000, color="red", weight=2, fill=True,
                      fill_opacity=0.12, tooltip=f"Blind zone: รัศมี {blind_km:.0f} km (เตือนไม่ทัน)").add_to(m)
    folium.Marker(list(epicenter[:2]), tooltip="จุดศูนย์กลาง",
                  popup=f"<b>จุดศูนย์กลาง</b><br>lat {epicenter[0]:.3f}<br>lon {epicenter[1]:.3f}"
                        f"<br>ลึก {epicenter[2]:.0f} km",
                  icon=folium.Icon(color="red", icon="star")).add_to(m)
    for s in stations or []:
        folium.CircleMarker([s["lat"], s["lon"]], radius=5, color="black", fill=True, fill_color="#555555",
                            fill_opacity=0.9, tooltip=f"สถานีตรวจจับ {s['code']}").add_to(m)
        lats.append(s["lat"]); lons.append(s["lon"])
    for s in sites:
        lead = s.get("lead")
        col = "gray" if lead is None else ("red" if lead <= 0 else ("orange" if lead < 10 else "green"))
        txt = "-" if lead is None else (f"{lead:.0f} s" if lead > 0 else "เตือนไม่ทัน")
        folium.Marker([s["lat"], s["lon"]], tooltip=f"{s['name']}: เวลาเตือนล่วงหน้า {txt}",
                      popup=f"<b>{s['name']}</b><br>{s.get('info', '')}",
                      icon=folium.Icon(color=col, icon="home")).add_to(m)
        lats.append(s["lat"]); lons.append(s["lon"])
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;top:10px;left:50px;z-index:9999;background:white;padding:6px 10px;"
        "border-radius:6px;box-shadow:0 1px 4px #888;font-size:13px'><b>เวลาคลื่นเดินทางและเวลาเตือนล่วงหน้า</b><br>"
        "บ้าน: เขียว = เตือนได้ ≥ 10 s, ส้ม = < 10 s, แดง = เตือนไม่ทัน; วงแดง = Blind zone; "
        "วงประ = แนวหน้าคลื่น S</div>"))
    m.fit_bounds([[min(lats) - 1, min(lons) - 1], [max(lats) + 1, max(lons) + 1]])
    return _finish(m, out_html)