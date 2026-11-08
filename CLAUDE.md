# CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย

ระบบเตือนภัยแผ่นดินไหวล่วงหน้าที่เขียนด้วย Python (comment/เอกสารส่วนใหญ่เป็นภาษาไทย) รับข้อมูลความเร่ง 3 แกน (Z, N, E)
จากหลายสถานี → ตรวจจับคลื่น P/S → หาตำแหน่ง → คำนวณขนาด → เวลาเดินทางถึงจุดของเรา → ทำนายความรุนแรง (JMA / MMI) → ออกคำเตือน

## ลำดับการทำงาน (pipeline)

| ขั้น | ไฟล์หลัก | หน้าที่ |
|---|---|---|
| โหลดข้อมูล | `data_loader.py`, `station_db.py`, `signal_utils.py` | อ่านไฟล์ `station/<event>/<CODE>_<date>_<time>.txt`, ตาราง `stations.csv` (code,lat,lon,table,vs30), filter/haversine/azimuth |
| 1 คลื่น P | `p_wave.py` | band-pass → STA/LTA → AIC → ให้คะแนนคุณภาพ pick (weight 0-4) |
| 2 คลื่น S | `s_wave.py` | S-P, H/V ratio, แกน R/T, วิธี "jump" สำหรับสถานีไกล (regional) |
| 3 ตำแหน่ง | `epicenter.py`, `pipeline.py` | Geiger (ทรงกลม + โลก 2 ชั้น Pg/Pn/Sg/Sn/Lg), grid search, PCA back-azimuth, association, วนซ้ำจนนิ่ง, แล้ว `teleseism_check` (ผลอยู่ใน `hypo["teleseism"]`) |
| 4 ขนาด | `magnitude.py` | M_JMA, Pd, τc, network magnitude (ถ่วงน้ำหนัก) |
| 5 เวลาเดินทาง | `travel_time.py` | เวลาที่คลื่น S ถึงจุดผู้รับ + เวลาเตือนล่วงหน้า |
| 6 ความรุนแรง | `intensity.py`, `vs30_tool.py` | GMPE + ปรับตาม Vs30 (`global_vs30.grd`), JMA instrumental intensity, เกณฑ์ ShakeAlert |
| ประกอบ | `eew_system.py` (`EEWSystem`), `main.py`, `synthetic.py` | demo แบบข้อมูลสังเคราะห์ (M6.3 เชียงราย → มธ.รังสิต) |
| แผนที่/ภาพ | `map_plot.py` | folium → `*.html`, `*.kml`, matplotlib → `*.png` |

สคริปต์รันแต่ละขั้นกับข้อมูลจริง: `run_p_wave.py`, `run_s_wave.py`, `run_epicenter.py`, `run_network.py`, `run_magnitude.py`,
`run_travel_time.py`, `run_intensity.py`, `run_validation.py` (ผลลัพธ์ `*_detail.txt`, `*_alert.png`, `network_picks.csv` ฯลฯ ถูก commit ไว้)

สคริปต์เรียกต่อกันเป็นทอด: `run_intensity.compute()` → `run_travel_time.analyse()` → `run_magnitude.compute()` → `pipeline.locate_event()`
ค่า `hypo` ชุดเดียวกันส่งต่อตลอดสาย (รวม `hypo["teleseism"]`) แก้การตัดสินใจเรื่องตำแหน่ง/เตือนที่ `pipeline.py` แล้วทุกสคริปต์ได้ผลตาม

## การเตือนและ teleseism
- ระดับเตือน `intensity.alert_level(M, MMI)`: 2 = สาธารณะ (M ≥ 5.0, MMI ≥ 3.5), 1 = แอป (M ≥ 4.5, MMI ≥ 2.5), 0 = ไม่แจ้ง; `run_intensity.simulate()` จำลองการอัปเดตทุก `UPDATE_INTERVAL`
- `pipeline.teleseism_check` (ค่า `tele_*` ใน `DEFAULTS`): เข้า ≥ 2 จาก 3 เกณฑ์ = แผ่นดินไหวไกลมาก → ไม่เตือนทุกระดับ, การ์ดไม่แสดงตำแหน่ง/ระยะจากเมือง (`run_magnitude.is_teleseism`)
  1) v_app ของหน้าคลื่น P ≥ 15 km/s  2) พบ Lg < 35% ของสถานีที่มี P  3) เวลา Lg คลาด > 20 s
  ทดสอบแล้ว: M8.8 คัมชัตกา 3/3, เหตุการณ์ในภูมิภาคทั้ง 4 ได้ 0/3 — แก้เกณฑ์แล้วต้องรัน `run_validation.py` ดูทั้ง 5 เหตุการณ์
- ทิศจาก plane-wave fit คลาดมาก (M8.8: 83° เทียบจริง ~35°) จึงไม่แสดงทิศในการ์ด

## ค่าตั้งต้น
- **`config.EVENT`** เลือกโฟลเดอร์เหตุการณ์ใน `station/` ที่สคริปต์ `run_*.py` ใช้ (ค่าตั้งต้น `2025-3-28_Mw7.7_station`) → `data_loader.DATA_DIR = station/<EVENT>` (`data_loader.event_dir()`); `"auto"`/ชื่อไฟล์/รหัสสถานีจะค้นในโฟลเดอร์นี้ก่อน แล้ว `station/` และโฟลเดอร์โค้ด
  - ค่าเลือกสถานีในแต่ละสคริปต์ (`STATIONS` ใน `run_network.py`, `SELECTED_STATIONS`/`DATA_FILE` ใน `run_magnitude.py`, `run_intensity.py`) ถ้าใส่ชื่อโฟลเดอร์เหตุการณ์จะทับ `config.EVENT` — เช็กก่อนสรุปว่ารันเหตุการณ์ไหน
  - `REF_EPICENTER` ใน `run_network.py` เป็นตำแหน่งของ Mw7.7 สะกาย ต้องเปลี่ยนเมื่อรันเหตุการณ์อื่น
  - รันเหตุการณ์อื่นโดยไม่แก้ไฟล์: `python -c "import config; config.EVENT='2026-05-18_M5.2_AllTMD'; import run_intensity as R; R.main()"` (ต้องตั้ง `config.EVENT` ก่อน import `data_loader`)
- **`config.py` เป็นที่เดียวสำหรับค่าคงที่**: Vp=6.0, Vs=3.5, Vpn=8.0, Vsn=4.6, ความหนาเปลือกโลก 35 km, พารามิเตอร์ STA/LTA ชุด local และ `P_REGIONAL`/`S_REGIONAL` (สถานี > ~300 km), เกณฑ์เตือน (Warning: JMA ≥ 4.5, Forecast: M ≥ 3.5 หรือ I ≥ 2.5)
- `pipeline.py` มี `DEFAULTS` ของขั้น network location (override ผ่าน `settings={...}` หรือชื่อเดียวกันใน `config.py`), `fix_depth=10 km` เป็นค่าตั้งต้น
- ความลึก: `pipeline.refine_depth` ลองหาความลึกจาก pick ชุดสุดท้าย รับเมื่อผ่านเกณฑ์ `depth_*` (สถานีใกล้ ≤ 150 km, ±15 km, ตำแหน่งย้าย ≤ 30 km, RMS ดีขึ้น) ไม่ผ่านใช้ 10 km; ผลอยู่ใน `hypo["depth_info"]`, การ์ดบอกว่าวัดได้หรือค่าตรึง ทั้ง 5 เหตุการณ์ปัจจุบันใช้ค่าตรึง
- `CRUST_THICKNESS_KM = 35` เป็นค่าปรับเทียบ: คลื่นหัว `(2H − h)·q` ที่ h = 10 km เท่าสูตรเดิม `2·30·q` แก้ค่านี้แล้วผลทุกเหตุการณ์เปลี่ยน ต้องรัน `run_validation.py`
- สูตรทั้งหมดเขียนแยกไว้ที่ `EEW_formulas.md` (ฉบับล่าสุด) และ `FORMULAS.md` (ฉบับเก่า) — แก้โค้ดสูตรแล้วควรอัปเดตเอกสารด้วย

## ข้อมูล
- `station/<วันที่>_M<ขนาด>_<ชื่อ>/` = ข้อมูลเหตุการณ์จริงจากสถานี TMD 5 เหตุการณ์ ค่าอ้างอิง USGS อยู่ใน `run_validation.EVENTS`:

  | โฟลเดอร์ | USGS | ผลล่าสุดของระบบ |
  |---|---|---|
  | `2025-3-28_Mw7.7_station` | Mww 7.7 สะกาย (us7000pn9s) | คลาด 51 km, Mw 7.9 |
  | `2025-03-29_M5.9_AllTMD` | **mb 5.1** เนปิดอว์ (us7000pnkj) — ชื่อโฟลเดอร์ใช้ขนาด 5.9 ของ TMD | คลาด 26 km, Mw 5.1 |
  | `2025-07-29_M8.8_Russia_AllTMD` | Mww 8.8 คัมชัตกา (us6000qw60) teleseism ~7000 km | ตรวจพบ teleseism ไม่เตือน |
  | `2025-08-21_M4.6_AllTMD` | mb 4.6 นอกชายฝั่งเมียนมา | คลาด 17 km, Mw 4.9 |
  | `2026-05-18_M5.2_AllTMD` | Mww 5.2 ชายฝั่งเมียนมา (us6000syg0) | คลาด 9 km, Mw 5.2 |

  สถานี PWSA ของ Mw7.7 มีเวลาไฟล์ช้ากว่าสถานีอื่น 40 นาที (ระบบตัดออกเองด้วย `check_time_coverage`)
  สถานีไทยอยู่ด้านเดียวของเหตุการณ์ในเมียนมา (gap > 180°) ความคลาดของตำแหน่งจึงอยู่ตามแนวระยะเป็นหลัก
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
- `.venv` ไม่มี pip และยังไม่มี `netCDF4` (ติดตั้งด้วย `uv pip install --python .venv/Scripts/python.exe netCDF4`); `global_vs30.grd` ไม่อยู่ในเครื่อง (gitignored) — ระบบยังรันได้โดยใช้ Vs30 จาก `stations.csv`
- รันด้วย `MPLBACKEND=Agg` เมื่อไม่ต้องการหน้าต่างกราฟ; `run_validation.py` ใช้ ~7 นาที, `run_intensity.py` ~2.5 นาที → รันเบื้องหลัง
- รันเสร็จแล้วไฟล์ผลลัพธ์ที่ track จะเปลี่ยน ถ้าไม่ตั้งใจ commit ให้ `git checkout --` คืนค่า
- Lint: `uvx --python 3.12 ruff check --select F *.py` (F821 `la`/`lo` ใน `pipeline.scan_locate` เป็น false positive)

## ข้อควรระวังเวลาแก้โค้ด
- ไฟล์ใหญ่ที่รวม logic มาก: `run_intensity.py` (~86 KB), `pipeline.py` (~52 KB), `map_plot.py`, `magnitude.py` — อ่านเป็นช่วง (offset/limit) อย่าอ่านทั้งไฟล์
- ไฟล์ผลลัพธ์ (`.png/.html/.kml/_detail.txt`) ถูก track ใน git จะเปลี่ยนทุกครั้งที่รันซ้ำ
- สถานีใกล้ vs ไกล ใช้พารามิเตอร์คนละชุด (local vs regional) — แก้ threshold ต้องตรวจทั้งสองโหมด
- คอมเมนต์และ output เป็นภาษาไทย ให้คงภาษาเดิมเมื่อเพิ่มโค้ด

## Graphify
สร้าง knowledge graph ของโค้ดไว้แล้วที่ `graphify-out/` (graphify 0.9.103; 705 nodes, 1628 edges, 31 communities; `GRAPH_REPORT.md`, `graph.json`, `graph.html`)
กราฟรวมเอกสาร `EEW_formulas.md` ด้วย: หัวข้อสูตรเชื่อมกับฟังก์ชันที่ใช้ (เช่น ข้อ 3.10 → `teleseism_check`) ใช้ `graphify explain "<ฟังก์ชัน>"` หาสูตรและเหตุผลประกอบได้
ใช้ https://github.com/Graphify-Labs/graphify — ก่อนค้นโค้ดกว้าง ๆ ให้ดู `graphify-out/GRAPH_REPORT.md` หรือถามผ่าน CLI:
```bash
graphify explain "detect_p_wave"
graphify path "detect_p_wave" "EEWSystem"
GRAPHIFY_BACKEND=none graphify update .   # หลังแก้โค้ด (ไม่ใช้ LLM)
```
ติดตั้งด้วย `uv tool install --python 3.12 graphifyy` (Python 3.15 ของเครื่องไม่มี wheel ของ tree-sitter) ข้อยกเว้นอยู่ใน `.graphifyignore`
God nodes (จำนวน edge): `haversine_km` (43), `azimuth_deg` (23), `detect_p_wave` (23), `bandpass` (22), `locate_event` (22), `detect_s_wave` (19), `load_station_txt` (18), `phase_travel_time`, `predict_arrival`, `locate`, `grid_search_locate`, `travel_time_2layer` — แก้ฟังก์ชันเหล่านี้กระทบหลายขั้น ใช้ `graphify affected "<ชื่อ>"` ดูผลกระทบก่อนแก้

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- When graphify-out/graph.json exists, run `graphify query "<question>"` before your first grep, search or source-file read for a codebase question, and put this rule in every subagent prompt that explores code. Read raw files only for lines the graph or an earlier read already located. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
- A committed graphify-out/graph.json is rewritten whole on every build, so a large diff to it in a pull request is expected. Do not flag it as churn or suggest skipping the graph rebuild or its hooks.
