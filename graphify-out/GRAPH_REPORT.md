# Graph Report - code  (2026-11-07)

## Corpus Check
- 28 files · ~59,252 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 3 file(s) not represented in the graph (top: (none) 3)

## Summary
- 705 nodes · 1628 edges · 31 communities (29 shown, 2 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 19 edges (avg confidence: 0.93)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `2b712c8b`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- simulate
- map_plot.py
- numpy
- pipeline.py
- intensity.py
- 4. magnitude.py — ขนาดแผ่นดินไหว
- run_travel_time.py
- station_weight
- run_magnitude.py
- msbb_magnitude
- config.py
- run_epicenter.py
- legacy_magnitudes
- run_network.py
- eew_system.py
- epicenter.py
- _apply_response
- _predict_core
- predict_intensity
- 4. magnitude.py — ขนาดแผ่นดินไหว
- network_magnitude
- 7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน
- run_intensity.py
- type_scores
- build_sites
- mblg_magnitude
- predict_mmi
- CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย
- observed
- mmi_field
- README.md

## God Nodes (most connected - your core abstractions)
1. `haversine_km()` - 43 edges
2. `azimuth_deg()` - 23 edges
3. `detect_p_wave()` - 23 edges
4. `locate_event()` - 22 edges
5. `bandpass()` - 22 edges
6. `detect_s_wave()` - 19 edges
7. `load_station_txt()` - 18 edges
8. `locate()` - 17 edges
9. `phase_travel_time()` - 17 edges
10. `predict_arrival()` - 17 edges

## Surprising Connections (you probably didn't know these)
- `3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`)` --references--> `geiger_locate_2layer()`  [INFERRED]
  EEW_formulas.md → epicenter.py
- `ลำดับการทำงาน (pipeline)` --references--> `EEWSystem`  [INFERRED]
  CLAUDE.md → eew_system.py
- `3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`)` --references--> `geiger_locate_2layer()`  [INFERRED]
  FORMULAS.md → epicenter.py
- `Graphify` --references--> `locate_event()`  [INFERRED]
  CLAUDE.md → pipeline.py
- `3.10 ตรวจแผ่นดินไหวไกลมาก (teleseism) — `teleseism_check`` --references--> `teleseism_check()`  [INFERRED]
  EEW_formulas.md → pipeline.py

## Import Cycles
- None detected.

## Communities (31 total, 2 thin omitted)

### Community 0 - "simulate"
Cohesion: 0.14
Nodes (10): _base(), _candidates(), correction_model(), _fit(), flag_gross_outliers(), magnitude_at(), model_corr(), _robust_trend() (+2 more)

### Community 1 - "map_plot.py"
Cohesion: 0.11
Nodes (16): _circle_coords(), _finish(), _kml_color(), make_intensity_map(), make_mmi_map(), make_station_mmi_map(), _mmi_label(), mmi_level() (+8 more)

### Community 2 - "numpy"
Cohesion: 0.06
Nodes (36): Graphify, check_units(), detect_sensor_type(), detect_station_type(), event_dir(), load_station_txt(), _parse_time(), aic_pick() (+28 more)

### Community 3 - "pipeline.py"
Cohesion: 0.08
Nodes (27): back_azimuth_pca(), phase_name(), predict_arrival(), _as_folder(), associate(), _cfg(), check_time_coverage(), delta_from_differential() (+19 more)

### Community 4 - "intensity.py"
Cohesion: 0.16
Nodes (7): action_text(), alert_level(), distance_to_trace_km(), _jma_filter(), jma_instrumental_intensity(), jma_intensity_detail(), mmi_class()

### Community 5 - "4. magnitude.py — ขนาดแผ่นดินไหว"
Cohesion: 0.04
Nodes (45): 0. signal_utils.py — ฟังก์ชันพื้นฐาน, 1.1 Characteristic Function ของ Allen (1978), 1.2 STA/LTA, 1.3 กันการ trigger บน noise burst (รถวิ่งผ่าน ประตูปิด ฯลฯ), 1.4 ปรับเวลาให้แม่นด้วย AIC (Maeda, 1985) — 2 รอบ, 1.5 คุณภาพ pick (แบบ Hypoinverse, Klein 2002), 1.6 ทิศการเคลื่อนที่แรก (First motion), 1.7 Pick คลื่น P ในหน้าต่าง (re-pick หลังรู้ตำแหน่งคร่าว ๆ) (+37 more)

### Community 6 - "run_travel_time.py"
Cohesion: 0.11
Nodes (19): make_travel_map(), analyse(), blind_zone_km(), fmt(), magnitude_ready_time(), main(), make_map(), phase_name() (+11 more)

### Community 7 - "station_weight"
Cohesion: 0.25
Nodes (4): _distance_factor(), _pick_factor(), _snr_factor(), station_weight()

### Community 8 - "run_magnitude.py"
Cohesion: 0.09
Nodes (29): find_station_file(), resolve_data_path(), _search_folders(), azimuthal_gap(), compute(), covers(), detail_report(), _epoch() (+21 more)

### Community 9 - "msbb_magnitude"
Cohesion: 0.15
Nodes (11): _check_index(), _clip(), ml_from_amplitude(), ml_magnitude(), msbb_from_velocity(), msbb_magnitude(), mwp_magnitude(), mwp_stable() (+3 more)

### Community 11 - "run_epicenter.py"
Cohesion: 0.09
Nodes (20): 3.1 เวลาเดินทาง: แบบจำลองโลก 2 ชั้น (เปลือกโลกหนา H บนเนื้อโลก), omori_distance(), single_station_location(), 3.1 เวลาเดินทาง: แบบจำลองโลก 2 ชั้น (เปลือกโลกหนา H บนเนื้อโลก), _great_circle(), make_epicenter_map(), analyze_station(), choose_consistent() (+12 more)

### Community 12 - "legacy_magnitudes"
Cohesion: 0.16
Nodes (9): jma_magnitude(), jma_magnitude_from_amplitude(), legacy_magnitudes(), max_horizontal_displacement_um(), pd_value(), reliability(), size_class(), tau_c() (+1 more)

### Community 13 - "run_network.py"
Cohesion: 0.12
Nodes (11): make_multi_station_map(), find_data_files(), main(), plot_location(), save_csv(), get_coords(), get_vs30(), load_station_db() (+3 more)

### Community 14 - "eew_system.py"
Cohesion: 0.18
Nodes (4): EEWSystem, magnitude_from_pd(), magnitude_from_tauc(), main()

### Community 15 - "epicenter.py"
Cohesion: 0.05
Nodes (34): _axial_diff_deg(), delta_from_sp(), _fit_t0_truncated(), geiger_locate(), geiger_locate_2layer(), baz_weights(), tcalc(), geiger_residuals() (+26 more)

### Community 16 - "_apply_response"
Cohesion: 0.18
Nodes (6): _apply_response(), broadband_velocity(), _butter_mag(), displacement(), _osc(), wood_anderson()

### Community 17 - "_predict_core"
Cohesion: 0.17
Nodes (7): fault_distance_from_hypo(), gmpe_weight(), ipe_allen2012(), mmi_from_gm(), observed_mmi(), pgv_bedrock_si_midorikawa(), _predict_core()

### Community 18 - "predict_intensity"
Cohesion: 0.25
Nodes (4): mmi_slope(), predict_intensity(), site_amplification(), site_mmi_term()

### Community 19 - "4. magnitude.py — ขนาดแผ่นดินไหว"
Cohesion: 0.05
Nodes (41): 0. signal_utils.py — ฟังก์ชันพื้นฐาน, 1.1 Characteristic Function ของ Allen (1978), 1.2 STA/LTA, 1.3 กันการ trigger บน noise burst (รถวิ่งผ่าน ประตูปิด ฯลฯ), 1.4 ปรับเวลาให้แม่นด้วย AIC (Maeda, 1985) — 2 รอบ, 1.5 คุณภาพ pick (แบบ Hypoinverse, Klein 2002), 1.6 ทิศการเคลื่อนที่แรก (First motion), 1.7 Pick คลื่น P ในหน้าต่าง (re-pick หลังรู้ตำแหน่งคร่าว ๆ) (+33 more)

### Community 20 - "network_magnitude"
Cohesion: 0.29
Nodes (4): network_by_type(), network_magnitude(), to_mw(), weighted_median()

### Community 21 - "7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน"
Cohesion: 0.09
Nodes (21): 7.0 ภาพรวม, 7.1 MMI ที่วัดได้ที่สถานี — GMICE ของ Worden et al. (2012) (ที่ USGS ShakeMap ใช้), 7.2 JMA Instrumental Intensity (คำนวณจากคลื่นความเร่งจริงที่สถานี), 7.3 ทำนาย MMI — วิธี ก) GMPE → GMICE, 7.4 ทำนาย MMI — วิธี ข) IPE ของ Allen, Wald & Worden (2012) (ที่ USGS PAGER ใช้) แบบระยะ hypocentral, 7.5 รวม 2 วิธี และความไม่แน่นอน, 7.6 ตรวจคุณภาพข้อมูล (QC) และปรับด้วยค่าที่วัดได้จริง, 7.7 จำลองการทำงานแบบเวลาจริง (run_intensity.py) (+13 more)

### Community 22 - "run_intensity.py"
Cohesion: 0.18
Nodes (19): compute(), detail_report(), _epi_html(), lead_text(), main(), make_map(), _make_map_leaflet(), make_station_map() (+11 more)

### Community 23 - "type_scores"
Cohesion: 0.25
Nodes (4): choose_preferred(), preferred_magnitude(), rupture_duration(), type_scores()

### Community 24 - "build_sites"
Cohesion: 0.22
Nodes (5): build_sites(), site_rjb(), station_csv(), station_sensor(), station_vs30()

### Community 25 - "mblg_magnitude"
Cohesion: 0.17
Nodes (7): mblg_from_amplitude(), mblg_magnitude(), _nth_peak(), _period_at(), wwssn_sp(), _wwssn_sp_disp_resp(), wwssn_sp_gain()

### Community 26 - "predict_mmi"
Cohesion: 0.33
Nodes (3): intensity_from_pgv(), jma_intensity_class(), predict_mmi()

### Community 27 - "CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย"
Cohesion: 0.25
Nodes (7): CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย, graphify, การรัน, ข้อควรระวังเวลาแก้โค้ด, ข้อมูล, ค่าตั้งต้น, ลำดับการทำงาน (pipeline)

### Community 29 - "observed"
Cohesion: 0.33
Nodes (3): despike(), is_clipped(), observed()

### Community 30 - "mmi_field"
Cohesion: 0.33
Nodes (3): mmi_field(), _mmi_rgba(), _trace_dist_grid()

## Knowledge Gaps
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `geiger_locate_2layer()` connect `epicenter.py` to `numpy`, `pipeline.py`, `4. magnitude.py — ขนาดแผ่นดินไหว`, `run_epicenter.py`, `4. magnitude.py — ขนาดแผ่นดินไหว`?**
  _High betweenness centrality (0.127) - this node is a cross-community bridge._
- **Should `simulate` be split into smaller, more focused modules?**
  _Cohesion score 0.14166666666666666 - nodes in this community are weakly interconnected._
- **Why does `mmi_from_pgv()` connect `7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน` to `numpy`, `intensity.py`?**
  _High betweenness centrality (0.107) - this node is a cross-community bridge._
- **Should `map_plot.py` be split into smaller, more focused modules?**
  _Cohesion score 0.11397849462365592 - nodes in this community are weakly interconnected._
- **Why does `สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง` connect `4. magnitude.py — ขนาดแผ่นดินไหว` to `7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน`?**
  _High betweenness centrality (0.092) - this node is a cross-community bridge._
- **Should `numpy` be split into smaller, more focused modules?**
  _Cohesion score 0.05980861244019139 - nodes in this community are weakly interconnected._
- **Should `pipeline.py` be split into smaller, more focused modules?**
  _Cohesion score 0.08233117483811286 - nodes in this community are weakly interconnected._