# Graph Report - code  (2026-11-08)

## Corpus Check
- 26 files · ~41,553 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 3 file(s) not represented in the graph (top: (none) 3)

## Summary
- 550 nodes · 1261 edges · 19 communities (18 shown, 1 thin omitted)
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 22 edges (avg confidence: 0.94)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `a1004917`
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
- make_mmi_map
- 3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว
- build_sites
- run_magnitude.py
- simulate
- magnitude.py
- plot_mmi_map
- haversine_km
- _mmi_on_grid
- _fit
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
10. `locate()` - 13 edges

## Surprising Connections (you probably didn't know these)
- `3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`)` --references--> `geiger_locate_2layer()`  [INFERRED]
  EEW_formulas.md → epicenter.py
- `3.11 หาความลึกแบบมีเงื่อนไข — `refine_depth`` --references--> `refine_depth()`  [INFERRED]
  EEW_formulas.md → pipeline.py
- `ข้อมูล` --references--> `check_time_coverage()`  [INFERRED]
  CLAUDE.md → pipeline.py
- `ลำดับการทำงาน (pipeline)` --references--> `EEWSystem`  [INFERRED]
  CLAUDE.md → eew_system.py
- `ลำดับการทำงาน (pipeline)` --references--> `make_station_mmi_map()`  [INFERRED]
  CLAUDE.md → map_plot.py

## Import Cycles
- None detected.

## Communities (19 total, 1 thin omitted)

### Community 0 - "run_intensity.py"
Cohesion: 0.42
Nodes (8): detail_report(), main(), method_text(), _mmi_txt(), plot_alert(), plot_distance(), _plt(), public_report()

### Community 1 - "map_plot.py"
Cohesion: 0.07
Nodes (22): _circle_coords(), _great_circle(), _kml_color(), make_multi_station_map(), _plain(), save_kml(), place(), style() (+14 more)

### Community 2 - "intensity.py"
Cohesion: 0.07
Nodes (16): EEWSystem, action_text(), alert_level(), distance_to_trace_km(), fault_distance_from_hypo(), ipe_allen2012(), mmi_class(), mmi_from_pgv() (+8 more)

### Community 3 - "pipeline.py"
Cohesion: 0.08
Nodes (28): find_station_file(), phase_name(), pick_weight(), predict_arrival(), _as_folder(), associate(), _cfg(), check_time_coverage() (+20 more)

### Community 4 - "numpy"
Cohesion: 0.06
Nodes (32): check_units(), detect_sensor_type(), detect_station_type(), load_station_txt(), _parse_time(), resolve_data_path(), _search_folders(), aic_pick() (+24 more)

### Community 5 - "4. magnitude.py — ขนาดแผ่นดินไหว"
Cohesion: 0.06
Nodes (35): 0. signal_utils.py — ฟังก์ชันพื้นฐาน, 1.1 Characteristic Function ของ Allen (1978), 1.2 STA/LTA, 1.3 กันการ trigger บน noise burst (รถวิ่งผ่าน ประตูปิด ฯลฯ), 1.4 ปรับเวลาให้แม่นด้วย AIC (Maeda, 1985) — 2 รอบ, 1.5 คุณภาพ pick (แบบ Hypoinverse, Klein 2002), 1.6 ทิศการเคลื่อนที่แรก (First motion), 1.7 Pick คลื่น P ในหน้าต่าง (re-pick หลังรู้ตำแหน่งคร่าว ๆ) (+27 more)

### Community 6 - "run_travel_time.py"
Cohesion: 0.11
Nodes (18): make_travel_map(), analyse(), blind_zone_km(), fmt(), magnitude_ready_time(), main(), make_map(), phase_name() (+10 more)

### Community 7 - "make_mmi_map"
Cohesion: 0.16
Nodes (9): _add_rupture(), _finish(), make_magnitude_map(), make_mmi_map(), make_station_mmi_map(), _mmi_label(), mmi_level(), _new_map() (+1 more)

### Community 8 - "3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว"
Cohesion: 0.08
Nodes (23): CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย, Graphify, graphify, การรัน, การเตือนและ teleseism, ข้อควรระวังเวลาแก้โค้ด, ข้อมูล, ค่าตั้งต้น (+15 more)

### Community 9 - "build_sites"
Cohesion: 0.22
Nodes (5): build_sites(), site_rjb(), station_csv(), station_sensor(), station_vs30()

### Community 10 - "run_magnitude.py"
Cohesion: 0.07
Nodes (34): azimuthal_gap(), despike(), is_clipped(), observed(), city_lines(), compute(), covers(), detail_report() (+26 more)

### Community 11 - "simulate"
Cohesion: 0.20
Nodes (7): _candidates(), correction_model(), flag_gross_outliers(), magnitude_at(), model_corr(), simulate(), site_corr()

### Community 12 - "magnitude.py"
Cohesion: 0.06
Nodes (25): broadband_velocity(), _butter_mag(), choose_preferred(), _clip(), _distance_factor(), magnitude_from_pd(), magnitude_from_tauc(), msbb_from_velocity() (+17 more)

### Community 13 - "plot_mmi_map"
Cohesion: 0.17
Nodes (7): ลำดับการทำงาน (pipeline), compute(), _halo(), _load_land(), plot_mmi_map(), select_rupture(), _usgs_reference()

### Community 14 - "haversine_km"
Cohesion: 0.09
Nodes (19): geiger_locate_2layer(), tcalc(), observations(), phase_travel_time(), travel_time_2layer(), scan_locate(), cost(), plot() (+11 more)

### Community 15 - "_mmi_on_grid"
Cohesion: 0.25
Nodes (4): mmi_field(), _mmi_on_grid(), _mmi_rgba(), _trace_dist_grid()

### Community 16 - "_fit"
Cohesion: 0.29
Nodes (4): _base(), _fit(), _mmi_rings(), _robust_trend()

### Community 22 - "make_map"
Cohesion: 0.20
Nodes (10): _epi_html(), lead_text(), make_map(), _make_map_leaflet(), make_station_map(), obs_mmi(), _popup_html(), rupture_info() (+2 more)

## Knowledge Gaps
- **1 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว` connect `3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว` to `4. magnitude.py — ขนาดแผ่นดินไหว`?**
  _High betweenness centrality (0.146) - this node is a cross-community bridge._
- **Should `map_plot.py` be split into smaller, more focused modules?**
  _Cohesion score 0.06938775510204082 - nodes in this community are weakly interconnected._
- **Why does `สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง` connect `4. magnitude.py — ขนาดแผ่นดินไหว` to `3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว`?**
  _High betweenness centrality (0.122) - this node is a cross-community bridge._
- **Should `intensity.py` be split into smaller, more focused modules?**
  _Cohesion score 0.06767676767676768 - nodes in this community are weakly interconnected._
- **Why does `geiger_locate_2layer()` connect `haversine_km` to `intensity.py`, `pipeline.py`, `numpy`, `3. epicenter.py + pipeline.py — หาจุดศูนย์กลางแผ่นดินไหว`, `run_magnitude.py`?**
  _High betweenness centrality (0.084) - this node is a cross-community bridge._
- **Should `pipeline.py` be split into smaller, more focused modules?**
  _Cohesion score 0.08078231292517007 - nodes in this community are weakly interconnected._
- **Should `numpy` be split into smaller, more focused modules?**
  _Cohesion score 0.05945945945945946 - nodes in this community are weakly interconnected._