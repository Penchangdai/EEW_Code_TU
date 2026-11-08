# Graph Report - code  (2026-11-08)

## Corpus Check
- 26 files · ~40,887 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 3 file(s) not represented in the graph (top: (none) 3)

## Summary
- 536 nodes · 1227 edges · 19 communities (17 shown, 2 thin omitted)
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 22 edges (avg confidence: 0.94)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `d877bcb4`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- run_intensity.py
- map_plot.py
- intensity.py
- pipeline.py
- numpy
- 4. magnitude.py — ขนาดแผ่นดินไหว
- run_travel_time.py
- csv
- 3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว
- build_sites
- observed
- simulate
- magnitude.py
- plot_mmi_map
- haversine_km
- _mmi_on_grid
- station_db.py
- make_map
- README.md

## God Nodes (most connected - your core abstractions)
1. `haversine_km()` - 35 edges
2. `locate_event()` - 22 edges
3. `phase_travel_time()` - 18 edges
4. `predict_arrival()` - 18 edges
5. `azimuth_deg()` - 18 edges
6. `detect_p_wave()` - 18 edges
7. `load_station_txt()` - 15 edges
8. `pick_p_window()` - 14 edges
9. `plot_mmi_map()` - 13 edges
10. `pick_s_window()` - 13 edges

## Surprising Connections (you probably didn't know these)
- `ข้อมูล` --references--> `check_time_coverage()`  [INFERRED]
  CLAUDE.md → pipeline.py
- `การเตือนและ teleseism` --references--> `simulate()`  [INFERRED]
  CLAUDE.md → run_intensity.py
- `ลำดับการทำงาน (pipeline)` --references--> `plot_mmi_map()`  [INFERRED]
  CLAUDE.md → run_intensity.py
- `ลำดับการทำงาน (pipeline)` --references--> `compute()`  [INFERRED]
  CLAUDE.md → run_intensity.py
- `ลำดับการทำงาน (pipeline)` --references--> `EEWSystem`  [INFERRED]
  CLAUDE.md → eew_system.py

## Import Cycles
- None detected.

## Communities (19 total, 2 thin omitted)

### Community 0 - "run_intensity.py"
Cohesion: 0.29
Nodes (11): compute(), detail_report(), lead_text(), main(), method_text(), _mmi_txt(), plot_alert(), plot_distance() (+3 more)

### Community 1 - "map_plot.py"
Cohesion: 0.08
Nodes (24): ลำดับการทำงาน (pipeline), _add_rupture(), _circle_coords(), _finish(), _great_circle(), _kml_color(), make_mmi_map(), make_multi_station_map() (+16 more)

### Community 2 - "intensity.py"
Cohesion: 0.11
Nodes (11): action_text(), distance_to_trace_km(), fault_distance_from_hypo(), ipe_allen2012(), mmi_class(), mmi_from_pgv(), mmi_slope(), _predict_core() (+3 more)

### Community 3 - "pipeline.py"
Cohesion: 0.07
Nodes (34): find_station_file(), resolve_data_path(), _search_folders(), 3.11 หาความลึกแบบมีเงื่อนไข — `refine_depth`, phase_name(), pick_weight(), predict_arrival(), _as_folder() (+26 more)

### Community 4 - "numpy"
Cohesion: 0.05
Nodes (40): check_units(), detect_sensor_type(), detect_station_type(), load_station_txt(), _parse_time(), azimuthal_gap(), tau_c(), aic_pick() (+32 more)

### Community 5 - "4. magnitude.py — ขนาดแผ่นดินไหว"
Cohesion: 0.06
Nodes (35): 0. signal_utils.py — ฟังก์ชันพื้นฐาน, 1.1 Characteristic Function ของ Allen (1978), 1.2 STA/LTA, 1.3 กันการ trigger บน noise burst (รถวิ่งผ่าน ประตูปิด ฯลฯ), 1.4 ปรับเวลาให้แม่นด้วย AIC (Maeda, 1985) — 2 รอบ, 1.5 คุณภาพ pick (แบบ Hypoinverse, Klein 2002), 1.6 ทิศการเคลื่อนที่แรก (First motion), 1.7 Pick คลื่น P ในหน้าต่าง (re-pick หลังรู้ตำแหน่งคร่าว ๆ) (+27 more)

### Community 6 - "run_travel_time.py"
Cohesion: 0.10
Nodes (20): make_travel_map(), main(), plot_summary(), analyse(), blind_zone_km(), fmt(), magnitude_ready_time(), main() (+12 more)

### Community 8 - "3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว"
Cohesion: 0.06
Nodes (32): CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย, Graphify, graphify, การรัน, การเตือนและ teleseism, ข้อควรระวังเวลาแก้โค้ด, ข้อมูล, ค่าตั้งต้น (+24 more)

### Community 9 - "build_sites"
Cohesion: 0.22
Nodes (5): build_sites(), site_rjb(), station_csv(), station_sensor(), station_vs30()

### Community 10 - "observed"
Cohesion: 0.33
Nodes (3): despike(), is_clipped(), observed()

### Community 11 - "simulate"
Cohesion: 0.14
Nodes (10): _base(), _candidates(), correction_model(), _fit(), flag_gross_outliers(), magnitude_at(), model_corr(), _robust_trend() (+2 more)

### Community 12 - "magnitude.py"
Cohesion: 0.06
Nodes (23): broadband_velocity(), _butter_mag(), choose_preferred(), _clip(), _distance_factor(), magnitude_from_pd(), magnitude_from_tauc(), msbb_from_velocity() (+15 more)

### Community 13 - "plot_mmi_map"
Cohesion: 0.25
Nodes (4): _halo(), _load_land(), plot_mmi_map(), _usgs_reference()

### Community 14 - "haversine_km"
Cohesion: 0.05
Nodes (28): 3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`), EEWSystem, geiger_locate_2layer(), tcalc(), observations(), phase_travel_time(), travel_time_2layer(), alert_level() (+20 more)

### Community 15 - "_mmi_on_grid"
Cohesion: 0.25
Nodes (4): mmi_field(), _mmi_on_grid(), _mmi_rgba(), _trace_dist_grid()

### Community 16 - "station_db.py"
Cohesion: 0.27
Nodes (4): get_vs30(), load_station_db(), main(), site_class()

### Community 22 - "make_map"
Cohesion: 0.18
Nodes (10): _epi_html(), make_map(), _make_map_leaflet(), make_station_map(), _mmi_rings(), obs_mmi(), _popup_html(), rupture_info() (+2 more)

## Knowledge Gaps
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว` connect `3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว` to `pipeline.py`, `4. magnitude.py — ขนาดแผ่นดินไหว`, `haversine_km`?**
  _High betweenness centrality (0.149) - this node is a cross-community bridge._
- **Should `map_plot.py` be split into smaller, more focused modules?**
  _Cohesion score 0.07878787878787878 - nodes in this community are weakly interconnected._
- **Why does `สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง` connect `4. magnitude.py — ขนาดแผ่นดินไหว` to `3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว`?**
  _High betweenness centrality (0.126) - this node is a cross-community bridge._
- **Should `intensity.py` be split into smaller, more focused modules?**
  _Cohesion score 0.11462450592885376 - nodes in this community are weakly interconnected._
- **Why does `geiger_locate_2layer()` connect `haversine_km` to `pipeline.py`, `numpy`?**
  _High betweenness centrality (0.086) - this node is a cross-community bridge._
- **Should `pipeline.py` be split into smaller, more focused modules?**
  _Cohesion score 0.06766917293233082 - nodes in this community are weakly interconnected._
- **Should `numpy` be split into smaller, more focused modules?**
  _Cohesion score 0.05350140056022409 - nodes in this community are weakly interconnected._