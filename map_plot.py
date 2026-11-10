"""
map_plot.py
วาดผลลงบนแผนที่จริงด้วย folium: ค่าเริ่มต้นเป็นภาพดาวเทียม + เส้นแบ่งเขตประเทศ (คล้าย Google Earth)
และบันทึกไฟล์ .kml คู่กัน สำหรับเปิดใน Google Earth (ไฟล์ -> นำเข้าไฟล์ KML)
ผลลัพธ์เป็นไฟล์ .html เปิดด้วย browser แล้วซูม/คลิกดูรายละเอียดได้ (ต้องต่ออินเทอร์เน็ตตอนเปิดดู)

ติดตั้ง:  pip install folium

แผนที่ที่มี
  make_mmi_map        : ความรุนแรงมาตรา MMI แบบ USGS (สี ShakeMap) + การแจ้งเตือน   <- run_intensity.py ใช้
  make_station_mmi_map: MMI ที่วัดได้ของแต่ละสถานี บนพื้นสี MMI แบบ ShakeMap      <- run_intensity.py ใช้
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
    "natgeo": ("NatGeo_World_Map", "Esri, National Geographic, Garmin, HERE", "แผนที่ภูมิศาสตร์ (Esri)", 12),
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
                                     control=True, show=((kind == "satellite") == (base == "satellite")) and base != "natgeo",
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


def _finish(m, out_html, collapsed=False):
    """ใส่ปุ่มเลือกชั้นข้อมูล บันทึก html และ kml (Google Earth)"""
    import folium
    folium.LayerControl(collapsed=collapsed).add_to(m)
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


def make_mmi_map(points, out_html, epicenter, rings=None, rupture=None, title="ความรุนแรงแผ่นดินไหว (MMI)", note=""):
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
    for la, lo in _add_rupture(m, rupture):
        lats.append(la); lons.append(lo)

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


def _add_rupture(m, rupture):
    """วาดแนวรอยแตก หรือขอบเขตรอยแตกโดยประมาณ  คืนค่าจุดที่ต้องรวมในขอบเขตแผนที่"""
    import folium
    if not rupture:
        return []
    if rupture.get("trace"):
        tr = [list(p) for p in rupture["trace"]]
        folium.PolyLine(tr, color="#7a0000", weight=5, opacity=0.95, tooltip=rupture["label"], popup=rupture["label"],
                        name=rupture["label"]).add_to(m)
        return tr
    if rupture.get("radius_km"):
        folium.Circle(list(rupture["center"]), radius=rupture["radius_km"] * 1000, color="#7a0000", weight=2,
                      dash_array="6 6", fill=False, tooltip=rupture["label"]).add_to(m)
    return []


def _mmi_tag(folium, lat, lon, text, mmi, lift_px):
    """ป้าย "ข้อความ + ค่า MMI" เหนือจุดที่กำหนด lift_px พิกเซล  เมาส์ทะลุผ่านได้ จึงไม่บังจุดข้างใต้"""
    col = MMI_COLORS[mmi_level(mmi)]
    return folium.Marker([lat, lon], interactive=False, icon=folium.DivIcon(
        icon_size=(190, 22), icon_anchor=(95, lift_px),
        html=("<div style='text-align:center;pointer-events:none'><span style='display:inline-block;"
              "background:rgba(255,255,255,.92);border:1px solid #555;border-radius:4px;padding:0 3px 0 5px;"
              f"font:bold 12px/18px Tahoma,sans-serif;color:#222;white-space:nowrap'>{text} "
              f"<span style='background:{col};color:{_text_on(col)};border-radius:3px;padding:0 4px'>{mmi:.1f}</span>"
              "</span></div>")))


def make_station_mmi_map(stations, out_html, epicenter, field=None, title="ระดับความรุนแรงแผ่นดินไหว",
                         subtitle="", basemap="natgeo", rupture=None, reference=None):
    """
    แผนที่ MMI แบบ USGS ShakeMap: แผนที่ภูมิศาสตร์ + พื้นสี MMI โปร่งใส + แนวรอยแตก + สถานีสีตาม MMI ที่วัดได้
    stations : list ของ dict {"name", "lat", "lon", "mmi": MMI วัดได้ (หรือ None), "popup", "province": ชื่อจังหวัดบนป้าย (หรือ None)}
               สถานีที่ไม่มีป้ายแสดงเฉพาะสี  รายละเอียดดูได้เมื่อชี้เมาส์หรือคลิก
    epicenter: dict {"lat", "lon", "popup", "mmi": MMI คาดที่จุดศูนย์กลาง (ไม่มี = ไม่แสดงป้าย)}
    field    : dict {"uri": ภาพ PNG (data URI), "bounds": [[lat_s, lon_w], [lat_n, lon_e]]}  พื้นสี MMI ที่คาด
    rupture  : dict จาก run_intensity.rupture_info   reference : ค่าอ้างอิงของ USGS (lat, lon)
    """
    try:
        import folium
        from folium.plugins import MousePosition
    except ImportError:
        print("(ไม่ได้ติดตั้ง folium จึงไม่วาดแผนที่จริง: pip install folium)")
        return None
    m = _new_map([epicenter["lat"], epicenter["lon"]], zoom_start=6, basemap=basemap)
    if field:
        folium.raster_layers.ImageOverlay(image=field["uri"], bounds=field["bounds"], opacity=0.72,
                                          name="พื้นสี MMI ที่คาด", interactive=False, zindex=1).add_to(m)
    lats, lons = [], []
    for la, lo in _add_rupture(m, rupture):
        lats.append(la); lons.append(lo)
    stg = folium.FeatureGroup(name="สถานี (สีตาม MMI ที่วัดได้)")
    num = folium.FeatureGroup(name="ตัวเลข MMI ที่สถานี", show=False)
    prov = folium.FeatureGroup(name="ชื่อจังหวัด + MMI")
    for p in stations:
        v = p.get("mmi")
        txt = f"{v:.1f}" if v is not None else "ไม่ผ่าน QC"
        tip = f"{p['name']}: MMI {txt}" + (f" ({MMI_ROMAN[mmi_level(v)]} {MMI_TEXT[mmi_level(v)][0]})" if v is not None else "")
        folium.CircleMarker([p["lat"], p["lon"]], radius=6 if v is not None else 3.5, color="#222222",
                            weight=1.2, fill=True, fill_color=MMI_COLORS[mmi_level(v)] if v is not None else "#555555",
                            fill_opacity=1.0 if v is not None else 0.6, tooltip=tip,
                            popup=folium.Popup(p.get("popup", ""), max_width=360)).add_to(stg)
        if v is not None:
            folium.Marker([p["lat"], p["lon"]], icon=folium.DivIcon(
                icon_size=(50, 18), icon_anchor=(25, 26),
                html=("<div style='text-align:center;font:bold 13px Tahoma,sans-serif;color:#222;"
                      f"text-shadow:0 0 3px #fff,0 0 3px #fff,0 0 3px #fff'>{v:.1f}</div>"))).add_to(num)
        if v is not None and p.get("province"):
            _mmi_tag(folium, p["lat"], p["lon"], p["province"], v, 30).add_to(prov)
        lats.append(p["lat"]); lons.append(p["lon"])
    stg.add_to(m)
    num.add_to(m)
    prov.add_to(m)
    star = ("<svg width='34' height='34' viewBox='0 0 24 24'><path fill='#ffd21f' stroke='#222' stroke-width='1.2' "
            "d='M12 1.8l3.1 6.6 7.2.9-5.3 5 1.4 7.1L12 17l-6.4 3.5 1.4-7.1-5.3-5 7.2-.9z'/></svg>")
    folium.Marker([epicenter["lat"], epicenter["lon"]], tooltip="จุดเหนือศูนย์กลางจากระบบ",
                  popup=folium.Popup(epicenter.get("popup", "จุดเหนือศูนย์กลาง"), max_width=320),
                  icon=folium.DivIcon(icon_size=(34, 34), icon_anchor=(17, 17), html=star)).add_to(m)
    if epicenter.get("mmi") is not None:
        _mmi_tag(folium, epicenter["lat"], epicenter["lon"], "จุดศูนย์กลาง (คาด)", epicenter["mmi"], 42).add_to(m)
    lats.append(epicenter["lat"]); lons.append(epicenter["lon"])
    if reference:
        folium.RegularPolygonMarker([reference["lat"], reference["lon"]], number_of_sides=4, radius=9, rotation=0,
                                    color="#111111", weight=2, fill=False,
                                    tooltip=f"จุดเหนือศูนย์กลางจาก USGS: {reference['mag_type']} {reference['mag']:g}").add_to(m)
    MousePosition(position="bottomright", separator=" : ", num_digits=3, prefix="").add_to(m)
    cells = "".join(f"<td style='background:{MMI_COLORS[k]};width:30px;text-align:center;border:1px solid #888;padding:1px'>"
                    f"{MMI_ROMAN[k]}</td>" for k in range(1, 11))
    rows = ["<span style='color:#c79a00;font-size:15px'>★</span> จุดเหนือศูนย์กลางจากระบบ"]
    if reference:
        rows.append("◇ จุดเหนือศูนย์กลางจาก USGS")
    if rupture and rupture.get("trace"):
        rows.append(f"<span style='color:#7a0000;font-weight:bold'>━━</span> {rupture['label']} (ข้อมูลหลังเหตุการณ์)")
    elif rupture and rupture.get("radius_km"):
        rows.append(f"<span style='color:#7a0000'>╌╌</span> {rupture['label']}")
    elif rupture and rupture.get("note"):
        rows.append(rupture["note"])
    rows += ["● สถานี สีตาม MMI ที่วัดได้  (จุดเทาเล็ก = ไม่ผ่าน QC)",
             "ป้าย = จังหวัด + MMI ที่วัดได้จากสถานีตัวแทน",
             "สถานีอื่น: ชี้เมาส์หรือคลิกเพื่อดูชื่อสถานีและค่า MMI",
             "พื้นสี = MMI ที่คาด จาก IPE + ค่าปรับจากสถานี"]
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;top:12px;left:60px;z-index:9999;background:rgba(255,255,255,.92);padding:7px 13px;"
        "border-radius:6px;box-shadow:0 1px 5px #777;font-family:Tahoma,sans-serif'>"
        f"<div style='font-size:19px;font-weight:bold;color:#222'>{title}</div>"
        f"<div style='font-size:13px;color:#444'>{subtitle}</div></div>"))
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;bottom:34px;left:10px;z-index:9999;background:rgba(255,255,255,.93);"
        "padding:7px 10px;border-radius:6px;box-shadow:0 1px 5px #777;font:12px/1.5 Tahoma,sans-serif'>"
        "<b>มาตรา MMI (USGS)</b>"
        f"<table style='border-collapse:collapse;margin:3px 0 5px'><tr>{cells}</tr></table>"
        + "<br>".join(rows) + "</div>"))
    m.fit_bounds([[min(lats) - 0.8, min(lons) - 1.5], [max(lats) + 0.8, max(lons) + 1.5]])
    return _finish(m, out_html, collapsed=True)


def _text_on(color):
    """สีตัวอักษรที่อ่านออกบนพื้นสี #rrggbb"""
    r, g, b = (int(color[i:i + 2], 16) for i in (1, 3, 5))
    return "#111111" if 0.299 * r + 0.587 * g + 0.114 * b > 140 else "#ffffff"


def make_magnitude_map(stations, out_html, summary, epicenter=None, rupture=None, classes=(), basemap="natgeo"):
    """
    แผนที่ขนาดแผ่นดินไหว: จุดกำเนิดจากขั้น epicenter + ป้ายขนาด + สถานีสีเทา (เหตุการณ์หนึ่งมีขนาดเดียว จึงไม่แยกสีตามสถานี)
    stations : list ของ dict {"name", "lat", "lon", "used": ได้ใช้คำนวณขนาดหรือไม่}
    summary  : dict {"mag_text", "color", "word", "rows": [(หัวข้อ, ข้อความ)], "note"}
    epicenter: dict {"lat", "lon", "popup"}  None = ไม่มีตำแหน่งที่เชื่อถือได้ (เช่น teleseism) จึงไม่วาดจุดกำเนิด
    rupture  : dict {"trace", "label"} แนวรอยแตกที่รู้แล้ว   classes : (ต่ำสุด, สูงสุด, คำบอกขนาด, สี) ของคำอธิบายสี
    """
    try:
        import folium
        from folium.plugins import MousePosition
    except ImportError:
        print("(ไม่ได้ติดตั้ง folium จึงไม่วาดแผนที่จริง: pip install folium)")
        return None
    lats = [p["lat"] for p in stations]
    lons = [p["lon"] for p in stations]
    if epicenter:
        lats.append(epicenter["lat"]); lons.append(epicenter["lon"])
    center = [epicenter["lat"], epicenter["lon"]] if epicenter else [float(np.mean(lats)), float(np.mean(lons))]
    m = _new_map(center, zoom_start=6, basemap=basemap)
    for la, lo in _add_rupture(m, rupture):
        lats.append(la); lons.append(lo)
    stg = folium.FeatureGroup(name="สถานี")
    for p in stations:
        used = p["used"]
        folium.CircleMarker([p["lat"], p["lon"]], radius=6 if used else 3.5, color="#222222", weight=1.2,
                            fill=True, fill_color="#8c8c8c", fill_opacity=1.0 if used else 0.6,
                            tooltip=p["name"], popup=folium.Popup(f"<b>{p['name']}</b>", max_width=320)).add_to(stg)
    stg.add_to(m)
    color, ink = summary["color"], _text_on(summary["color"])
    if epicenter:
        star = ("<svg width='34' height='34' viewBox='0 0 24 24' style='display:block;margin:2px auto 0'>"
                "<path fill='#ffd21f' stroke='#222' stroke-width='1.2' "
                "d='M12 1.8l3.1 6.6 7.2.9-5.3 5 1.4 7.1L12 17l-6.4 3.5 1.4-7.1-5.3-5 7.2-.9z'/></svg>")
        badge = (f"<div style='width:max-content;margin:auto;background:{color};color:{ink};border:2px solid #222;"
                 "border-radius:13px;padding:1px 10px;font:bold 15px Tahoma,sans-serif;white-space:nowrap'>"
                 f"Mw {summary['mag_text']}</div>")
        folium.Marker([epicenter["lat"], epicenter["lon"]], tooltip=f"จุดกำเนิดจากระบบ: Mw {summary['mag_text']}",
                      popup=folium.Popup(epicenter.get("popup", "จุดกำเนิด"), max_width=340),
                      icon=folium.DivIcon(icon_size=(110, 62), icon_anchor=(55, 45), html=badge + star)).add_to(m)
    MousePosition(position="bottomright", separator=" : ", num_digits=3, prefix="").add_to(m)
    rows = "".join(f"<tr><td style='color:#555;padding:2px 8px 2px 0;vertical-align:top;white-space:nowrap'>{k}</td>"
                   f"<td style='padding:2px 0'>{v}</td></tr>" for k, v in summary["rows"])
    note = f"<div style='margin-top:5px;color:#a50026;font-weight:bold'>{summary['note']}</div>" if summary.get("note") else ""
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;top:12px;left:60px;z-index:9999;max-width:min(400px,calc(100vw - 120px));"
        "max-height:calc(100vh - 150px);overflow:auto;background:rgba(255,255,255,.94);padding:9px 13px;"
        "border-radius:6px;box-shadow:0 1px 5px #777;font:13px/1.45 Tahoma,sans-serif;color:#222'>"
        "<details id='mag-card' open><summary style='cursor:pointer'>"
        "<span style='font-size:15px;font-weight:bold'>ขนาดแผ่นดินไหว</span>"
        f"<div style='margin:4px 0 6px'><span style='background:{color};color:{ink};border-radius:6px;padding:1px 10px;"
        f"font-size:26px;font-weight:bold'>Mw {summary['mag_text']}</span>"
        f"<span style='margin-left:9px;font-size:15px;font-weight:bold'>{summary['word']}</span></div></summary>"
        f"<table style='border-collapse:collapse'>{rows}</table>{note}</details></div>"))
    cells = "".join(f"<td style='background:{c};color:{_text_on(c)};min-width:42px;text-align:center;"
                    f"border:1px solid #888;padding:1px 3px'>{w}<br>{lo:g}{'+' if hi >= 10 else f'–{hi:g}'}</td>"
                    for lo, hi, w, c in classes)
    legend = ["<span style='color:#c79a00;font-size:15px'>★</span> จุดกำเนิดจากขั้น epicenter" if epicenter
              else "ไม่แสดงจุดกำเนิด เพราะไม่มีตำแหน่งที่เชื่อถือได้"]
    if rupture and rupture.get("trace"):
        legend.append(f"<span style='color:#7a0000;font-weight:bold'>━━</span> {rupture['label']} (ข้อมูลหลังเหตุการณ์)")
    legend.append("<span style='color:#8c8c8c'>●</span> สถานีที่ใช้คำนวณขนาด  (จุดเล็ก = ไม่ได้ใช้คำนวณ)")
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;bottom:34px;left:10px;z-index:9999;max-width:calc(100vw - 20px);overflow:auto;"
        "background:rgba(255,255,255,.93);padding:7px 10px;border-radius:6px;box-shadow:0 1px 5px #777;"
        "font:12px/1.5 Tahoma,sans-serif'><b>กลุ่มขนาด (Mw)</b>"
        f"<table style='border-collapse:collapse;margin:3px 0 5px'><tr>{cells}</tr></table>"
        + "<br>".join(legend) + "</div>"))
    # จัดขอบเขตหลังหน้าโหลดเสร็จ และเว้นที่ให้การ์ดสรุป  จอแคบพับการ์ดไว้ก่อน
    bounds = [[float(min(lats)) - 0.8, float(min(lons)) - 1.5], [float(max(lats)) + 0.8, float(max(lons)) + 1.5]]
    mp = m.get_name()
    m.get_root().script.add_child(folium.Element(f"""
    window.addEventListener('load', function () {{
        var wide = window.innerWidth >= 800;
        if (!wide) {{ document.getElementById('mag-card').removeAttribute('open'); }}
        function fit() {{
            {mp}.invalidateSize();
            {mp}.fitBounds({bounds}, {{paddingTopLeft: wide ? [460, 20] : [20, 90], paddingBottomRight: [20, 20]}});
        }}
        fit();
        {mp}.on('resize', function () {{ if ({mp}.getZoom() < 3) {{ fit(); }} }});
    }});"""))
    return _finish(m, out_html, collapsed=True)


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