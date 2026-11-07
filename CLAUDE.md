# CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย

ระบบเตือนภัยแผ่นดินไหวล่วงหน้าที่เขียนด้วย Python (comment/เอกสารส่วนใหญ่เป็นภาษาไทย) รับข้อมูลความเร่ง 3 แกน (Z, N, E)
จากหลายสถานี → ตรวจจับคลื่น P/S → หาตำแหน่ง → คำนวณขนาด → เวลาเดินทางถึงจุดของเรา → ทำนายความรุนแรง (JMA / MMI) → ออกคำเตือน

## ลำดับการทำงาน (pipeline)

| ขั้น | ไฟล์หลัก | หน้าที่ |
|---|---|---|
| โหลดข้อมูล | `data_loader.py`, `station_db.py`, `signal_utils.py` | อ่านไฟล์ `station/<event>/<CODE>_<date>_<time>.txt`, ตาราง `stations.csv` (code,lat,lon,table,vs30), filter/haversine/azimuth |
| 1 คลื่น P | `p_wave.py` | band-pass → STA/LTA → AIC → ให้คะแนนคุณภาพ pick (weight 0-4) |
| 2 คลื่น S | `s_wave.py` | S-P, H/V ratio, แกน R/T, วิธี "jump" สำหรับสถานีไกล (regional) |
| 3 ตำแหน่ง | `epicenter.py`, `pipeline.py` | Geiger (ทรงกลม + โลก 2 ชั้น Pg/Pn/Sg/Sn/Lg), grid search, PCA back-azimuth, association, วนซ้ำจนนิ่ง |
| 4 ขนาด | `magnitude.py` | M_JMA, Pd, τc, network magnitude (ถ่วงน้ำหนัก) |
| 5 เวลาเดินทาง | `travel_time.py` | เวลาที่คลื่น S ถึงจุดผู้รับ + เวลาเตือนล่วงหน้า |
| 6 ความรุนแรง | `intensity.py`, `vs30_tool.py` | GMPE + ปรับตาม Vs30 (`global_vs30.grd`), JMA instrumental intensity, เกณฑ์ ShakeAlert |
| ประกอบ | `eew_system.py` (`EEWSystem`), `main.py`, `synthetic.py` | demo แบบข้อมูลสังเคราะห์ (M6.3 เชียงราย → มธ.รังสิต) |
| แผนที่/ภาพ | `map_plot.py` | folium → `*.html`, `*.kml`, matplotlib → `*.png` |

สคริปต์รันแต่ละขั้นกับข้อมูลจริง: `run_p_wave.py`, `run_s_wave.py`, `run_epicenter.py`, `run_network.py`, `run_magnitude.py`,
`run_travel_time.py`, `run_intensity.py`, `run_validation.py` (ผลลัพธ์ `*_detail.txt`, `*_alert.png`, `network_picks.csv` ฯลฯ ถูก commit ไว้)

## ค่าตั้งต้น
- **`config.py` เป็นที่เดียวสำหรับค่าคงที่**: Vp=6.0, Vs=3.5, Vpn=8.0, Vsn=4.6, ความหนาเปลือกโลก 30 km, พารามิเตอร์ STA/LTA ชุด local และ `P_REGIONAL`/`S_REGIONAL` (สถานี > ~300 km), เกณฑ์เตือน (Warning: JMA ≥ 4.5, Forecast: M ≥ 3.5 หรือ I ≥ 2.5)
- `pipeline.py` มี `DEFAULTS` ของขั้น network location (override ผ่าน `settings={...}` หรือชื่อเดียวกันใน `config.py`), `fix_depth=10 km` เป็นค่าตั้งต้น
- สูตรทั้งหมดเขียนแยกไว้ที่ `EEW_formulas.md` (ฉบับล่าสุด) และ `FORMULAS.md` (ฉบับเก่า) — แก้โค้ดสูตรแล้วควรอัปเดตเอกสารด้วย

## ข้อมูล
- `station/<วันที่>_M<ขนาด>_<ชื่อ>/` = ข้อมูลเหตุการณ์จริงจากสถานี TMD (เช่น `2025-03-29_M5.9_AllTMD`, `2025-07-29_M8.8_Russia_AllTMD`, `2025-3-28_Mw7.7_station`, `2026-05-18_M5.2_AllTMD`)
- `global_vs30.grd` ~610 MB (netCDF, ใช้ netCDF4 อ่าน) — **ไม่ควรอ่านทั้งไฟล์/พิมพ์ออกมา**
- `SED-EEW-SeisComP-contributions/` = โค้ด C++ ของ SeisComP/SED EEW (`apps/eew`, `libs`, `doc`) เป็นของภายนอก ใช้อ้างอิง ไม่ใช่ส่วนของระบบ Python

## การรัน
```bash
pip install -r requirements.txt   # numpy scipy matplotlib folium netCDF4
python main.py                    # demo ข้อมูลสังเคราะห์
python run_network.py             # ตำแหน่งจากเครือข่ายจริง -> network_*.png/html/kml/csv
python run_magnitude.py           # ขนาด -> magnitude_*.png/txt
python run_travel_time.py         # เวลาเดินทาง
python run_intensity.py           # ความรุนแรง -> intensity_*.png/html/kml/txt
python run_validation.py          # เทียบกับแคตตาล็อกจริง
```
มี `.venv` ในโฟลเดอร์ (Windows, Python 3.15) ไม่มีชุดทดสอบอัตโนมัติ; ตรวจผลด้วย `run_validation.py` และเทียบกับค่า "จริง" ใน `synthetic.py`

## ข้อควรระวังเวลาแก้โค้ด
- ไฟล์ใหญ่ที่รวม logic มาก: `run_intensity.py` (~86 KB), `pipeline.py` (~52 KB), `map_plot.py`, `magnitude.py` — อ่านเป็นช่วง (offset/limit) อย่าอ่านทั้งไฟล์
- ไฟล์ผลลัพธ์ (`.png/.html/.kml/_detail.txt`) ถูก track ใน git จะเปลี่ยนทุกครั้งที่รันซ้ำ
- สถานีใกล้ vs ไกล ใช้พารามิเตอร์คนละชุด (local vs regional) — แก้ threshold ต้องตรวจทั้งสองโหมด
- คอมเมนต์และ output เป็นภาษาไทย ให้คงภาษาเดิมเมื่อเพิ่มโค้ด

## Graphify
สร้าง knowledge graph ของโค้ดไว้แล้วที่ `graphify-out/` (691 nodes, 1598 edges, 29 communities; `GRAPH_REPORT.md`, `graph.json`, `graph.html`)
ใช้ https://github.com/Graphify-Labs/graphify — ก่อนค้นโค้ดกว้าง ๆ ให้ดู `graphify-out/GRAPH_REPORT.md` หรือถามผ่าน CLI:
```bash
graphify explain "detect_p_wave"
graphify path "detect_p_wave" "EEWSystem"
GRAPHIFY_BACKEND=none graphify update .   # หลังแก้โค้ด (ไม่ใช้ LLM)
```
ติดตั้งด้วย `uv tool install --python 3.12 graphifyy` (Python 3.15 ของเครื่องไม่มี wheel ของ tree-sitter) ข้อยกเว้นอยู่ใน `.graphifyignore`
God nodes: `haversine_km`, `detect_p_wave`, `azimuth_deg`, `bandpass`, `locate_event`, `detect_s_wave`, `load_station_txt`

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- When graphify-out/graph.json exists, run `graphify query "<question>"` before your first grep, search or source-file read for a codebase question, and put this rule in every subagent prompt that explores code. Read raw files only for lines the graph or an earlier read already located. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
- A committed graphify-out/graph.json is rewritten whole on every build, so a large diff to it in a pull request is expected. Do not flag it as churn or suggest skipping the graph rebuild or its hooks.
