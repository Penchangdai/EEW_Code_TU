# Graph Report - code  (2026-11-08)

## Corpus Check
- 28 files · ~59,170 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 3 file(s) not represented in the graph (top: (none) 3)

## Summary
- 709 nodes · 1657 edges · 31 communities (30 shown, 1 thin omitted)
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 37 edges (avg confidence: 0.94)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `4d685095`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- p_wave.py
- map_plot.py
- s_wave.py
- pipeline.py
- intensity.py
- 4. magnitude.py — ขนาดแผ่นดินไหว
- run_travel_time.py
- station_weight
- run_magnitude.py
- msbb_magnitude
- epicenter.py
- run_epicenter.py
- legacy_magnitudes
- geiger_locate_2layer
- eew_system.py
- 1. p_wave.py — ตรวจจับคลื่น P
- _apply_response
- legacy_magnitudes
- 7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน
- 4. magnitude.py — ขนาดแผ่นดินไหว
- network_magnitude
- 7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน
- run_intensity.py
- phase_travel_time
- refine_depth
- mblg_magnitude
- สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง
- CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย
- mmi_from_pgv
- 2. s_wave.py — ตรวจจับคลื่น S และ Lg
- README.md

## God Nodes (most connected - your core abstractions)
1. `haversine_km()` - 44 edges
2. `locate_event()` - 24 edges
3. `detect_p_wave()` - 23 edges
4. `azimuth_deg()` - 23 edges
5. `bandpass()` - 22 edges
6. `phase_travel_time()` - 19 edges
7. `detect_s_wave()` - 19 edges
8. `predict_arrival()` - 18 edges
9. `travel_time_2layer()` - 18 edges
10. `load_station_txt()` - 18 edges

## Surprising Connections (you probably didn't know these)
- `3.1 เวลาเดินทาง: แบบจำลองโลก 2 ชั้น (เปลือกโลกหนา H บนเนื้อโลก)` --references--> `omori_distance()`  [INFERRED]
  EEW_formulas.md → epicenter.py
- `3.1 เวลาเดินทาง: แบบจำลองโลก 2 ชั้น (เปลือกโลกหนา H บนเนื้อโลก)` --references--> `omori_distance()`  [INFERRED]
  FORMULAS.md → epicenter.py
- `3.11 หาความลึกแบบมีเงื่อนไข — `refine_depth`` --references--> `refine_depth()`  [INFERRED]
  EEW_formulas.md → pipeline.py
- `3.10 ตรวจแผ่นดินไหวไกลมาก (teleseism) — `teleseism_check`` --references--> `teleseism_check()`  [INFERRED]
  EEW_formulas.md → pipeline.py
- `3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`)` --references--> `geiger_locate_2layer()`  [INFERRED]
  EEW_formulas.md → epicenter.py

## Import Cycles
- None detected.

## Communities (31 total, 1 thin omitted)

### Community 0 - "p_wave.py"
Cohesion: 0.15
Nodes (8): allen_cf(), detect_p_wave(), first_motion(), first_real_trigger(), is_sustained(), _noise_window(), pick_p_window(), sta_lta()

### Community 1 - "map_plot.py"
Cohesion: 0.05
Nodes (31): _circle_coords(), _finish(), _great_circle(), _kml_color(), make_epicenter_map(), make_intensity_map(), make_mmi_map(), make_multi_station_map() (+23 more)

### Community 2 - "s_wave.py"
Cohesion: 0.10
Nodes (16): aic_pick(), _cfg(), first_sustained(), pick_quality(), _cfg(), _detect_s_jump(), _detect_s_usgs(), detect_s_wave() (+8 more)

### Community 3 - "pipeline.py"
Cohesion: 0.14
Nodes (16): ข้อมูล, azimuthal_gap(), phase_name(), associate(), _cfg(), check_time_coverage(), delta_from_differential(), load_stations() (+8 more)

### Community 4 - "intensity.py"
Cohesion: 0.07
Nodes (20): action_text(), alert_level(), distance_to_trace_km(), fault_distance_from_hypo(), gmpe_weight(), intensity_from_pgv(), ipe_allen2012(), _jma_filter() (+12 more)

### Community 5 - "4. magnitude.py — ขนาดแผ่นดินไหว"
Cohesion: 0.06
Nodes (33): 0. signal_utils.py — ฟังก์ชันพื้นฐาน, 1.1 Characteristic Function ของ Allen (1978), 1.2 STA/LTA, 1.3 กันการ trigger บน noise burst (รถวิ่งผ่าน ประตูปิด ฯลฯ), 1.4 ปรับเวลาให้แม่นด้วย AIC (Maeda, 1985) — 2 รอบ, 1.5 คุณภาพ pick (แบบ Hypoinverse, Klein 2002), 1.6 ทิศการเคลื่อนที่แรก (First motion), 1.7 Pick คลื่น P ในหน้าต่าง (re-pick หลังรู้ตำแหน่งคร่าว ๆ) (+25 more)

### Community 6 - "run_travel_time.py"
Cohesion: 0.12
Nodes (18): analyse(), blind_zone_km(), fmt(), magnitude_ready_time(), main(), make_map(), phase_name(), plot() (+10 more)

### Community 7 - "station_weight"
Cohesion: 0.25
Nodes (4): _distance_factor(), _pick_factor(), _snr_factor(), station_weight()

### Community 8 - "run_magnitude.py"
Cohesion: 0.12
Nodes (22): jma_instrumental_intensity(), compute(), covers(), detail_report(), _epoch(), _fmt(), _folder_of(), is_teleseism() (+14 more)

### Community 9 - "msbb_magnitude"
Cohesion: 0.21
Nodes (7): _check_index(), _clip(), ml_from_amplitude(), ml_magnitude(), mwp_magnitude(), mwp_stable(), _result()

### Community 10 - "epicenter.py"
Cohesion: 0.10
Nodes (15): _axial_diff_deg(), _fit_t0_truncated(), geiger_locate(), geiger_residuals(), grid_search_locate(), cost(), predict(), residuals() (+7 more)

### Community 11 - "run_epicenter.py"
Cohesion: 0.06
Nodes (32): check_units(), detect_sensor_type(), detect_station_type(), find_station_file(), load_station_txt(), _parse_time(), resolve_data_path(), _search_folders() (+24 more)

### Community 12 - "legacy_magnitudes"
Cohesion: 0.13
Nodes (9): choose_preferred(), jma_magnitude(), jma_magnitude_from_amplitude(), max_horizontal_displacement_um(), preferred_magnitude(), reliability(), rupture_duration(), size_class() (+1 more)

### Community 13 - "geiger_locate_2layer"
Cohesion: 0.08
Nodes (25): 3.10 ตรวจแผ่นดินไหวไกลมาก (teleseism) — `teleseism_check`, 3.11 หาความลึกแบบมีเงื่อนไข — `refine_depth`, 3.1 เวลาเดินทาง: แบบจำลองโลก 2 ชั้น (เปลือกโลกหนา H บนเนื้อโลก), 3.2 ตรวจเวลาไฟล์ (Time coverage), 3.3 Association — เวลา P ของสองสถานีต้องเข้ากันได้ทางกายภาพ, 3.4 ตำแหน่งเริ่มต้น: scan จาก P + ทิศ + เวลา Lg + พลังงาน S (grid search), 3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`), 3.6 ขั้นตอนทั้งหมดของ pipeline (+17 more)

### Community 14 - "eew_system.py"
Cohesion: 0.07
Nodes (15): EEWSystem, single_station_location(), main(), make_station_record(), _wavelet(), alert_time_from_picks(), blind_zone_radius(), _delta_from_R() (+7 more)

### Community 15 - "1. p_wave.py — ตรวจจับคลื่น P"
Cohesion: 0.22
Nodes (9): 1.1 Characteristic Function ของ Allen (1978), 1.2 STA/LTA, 1.3 กันการ trigger บน noise burst (รถวิ่งผ่าน ประตูปิด ฯลฯ), 1.4 ปรับเวลาให้แม่นด้วย AIC (Maeda, 1985) — 2 รอบ, 1.5 คุณภาพ pick (แบบ Hypoinverse, Klein 2002), 1.6 ทิศการเคลื่อนที่แรก (First motion), 1.7 Pick คลื่น P ในหน้าต่าง (re-pick หลังรู้ตำแหน่งคร่าว ๆ), 1.8 การแปลงหน่วยข้อมูลจริง (data_loader.py) (+1 more)

### Community 16 - "_apply_response"
Cohesion: 0.14
Nodes (8): _apply_response(), broadband_velocity(), _butter_mag(), displacement(), _osc(), wood_anderson(), wwssn_sp(), _wwssn_sp_disp_resp()

### Community 17 - "legacy_magnitudes"
Cohesion: 0.16
Nodes (8): legacy_magnitudes(), magnitude_from_pd(), magnitude_from_tauc(), pd_value(), tau_c(), acc_to_disp(), acc_to_vel(), integrate()

### Community 18 - "7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน"
Cohesion: 0.22
Nodes (9): 7.0 ภาพรวม, 7.2 JMA Instrumental Intensity (คำนวณจากคลื่นความเร่งจริงที่สถานี), 7.3 ทำนาย MMI — วิธี ก) GMPE → GMICE, 7.4 ทำนาย MMI — วิธี ข) IPE ของ Allen, Wald & Worden (2012) (ที่ USGS PAGER ใช้) แบบระยะ hypocentral, 7.5 รวม 2 วิธี และความไม่แน่นอน, 7.6 ตรวจคุณภาพข้อมูล (QC) และปรับด้วยค่าที่วัดได้จริง, 7.7 จำลองการทำงานแบบเวลาจริง (run_intensity.py), 7.8 คำอธิบายระดับ MMI (ตามตาราง USGS ShakeMap)  ระดับ = ค่าที่ปัดเป็นจำนวนเต็ม (+1 more)

### Community 19 - "4. magnitude.py — ขนาดแผ่นดินไหว"
Cohesion: 0.17
Nodes (12): 4.0 จำลองเครื่องวัดมาตรฐานในโดเมนความถี่ (s = 2πif), 4.10 คำบอกขนาด (การแบ่งกลุ่มของ USGS), 4.1 ML — Local magnitude (IASPEI 2013 = Hutton & Boore, 1987), 4.2 mb_Lg — Lg magnitude (Nuttli, 1973), 4.3 Ms_BB — Broadband surface-wave magnitude (IASPEI 2013), 4.4 Mwp — Moment magnitude จากคลื่น P (Tsuboi et al., 1995; Whitmore et al., 2002), 4.5 ขนาดแบบเดิม / สำหรับเตือนภัยเร็ว (ไม่ใช่มาตรฐาน USGS ใช้เทียบ), 4.6 น้ำหนักของแต่ละสถานี (0–1) (+4 more)

### Community 20 - "network_magnitude"
Cohesion: 0.29
Nodes (4): network_by_type(), network_magnitude(), to_mw(), weighted_median()

### Community 21 - "7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน"
Cohesion: 0.22
Nodes (9): 7.0 ภาพรวม, 7.2 JMA Instrumental Intensity (คำนวณจากคลื่นความเร่งจริงที่สถานี), 7.3 ทำนาย MMI — วิธี ก) GMPE → GMICE, 7.4 ทำนาย MMI — วิธี ข) IPE ของ Allen, Wald & Worden (2012) (ที่ USGS PAGER ใช้) แบบระยะ hypocentral, 7.5 รวม 2 วิธี และความไม่แน่นอน, 7.6 ตรวจคุณภาพข้อมูล (QC) และปรับด้วยค่าที่วัดได้จริง, 7.7 จำลองการทำงานแบบเวลาจริง (run_intensity.py), 7.8 คำอธิบายระดับ MMI (ตามตาราง USGS ShakeMap)  ระดับ = ค่าที่ปัดเป็นจำนวนเต็ม (+1 more)

### Community 22 - "run_intensity.py"
Cohesion: 0.06
Nodes (41): _base(), build_sites(), _candidates(), compute(), correction_model(), despike(), detail_report(), _epi_html() (+33 more)

### Community 23 - "phase_travel_time"
Cohesion: 0.16
Nodes (10): Graphify, การรัน, phase_travel_time(), predict_arrival(), travel_time_2layer(), plot_record_section(), _s_energy_cf(), scan_locate() (+2 more)

### Community 24 - "refine_depth"
Cohesion: 0.19
Nodes (8): _locate(), _locator_dict(), _mark_notes(), _n_obs(), refine_depth(), _residuals(), review(), _review_loo()

### Community 25 - "mblg_magnitude"
Cohesion: 0.19
Nodes (9): mblg_from_amplitude(), mblg_magnitude(), msbb_from_velocity(), msbb_magnitude(), _nth_peak(), _period_at(), _snr(), station_magnitudes() (+1 more)

### Community 26 - "สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง"
Cohesion: 0.25
Nodes (7): 0. signal_utils.py — ฟังก์ชันพื้นฐาน, 5. run_magnitude.py — การ์ดแจ้งเตือน, 6. travel_time.py + run_travel_time.py — เวลาที่คลื่นเดินทางมาถึงตำแหน่งของเรา และเวลาเตือนล่วงหน้า, 8. เกณฑ์การออกคำเตือน — ShakeAlert (USGS), ข้อจำกัดที่ควรทราบ, สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง, เอกสารอ้างอิง

### Community 27 - "CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย"
Cohesion: 0.17
Nodes (9): CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย, graphify, การเตือนและ teleseism, ข้อควรระวังเวลาแก้โค้ด, ค่าตั้งต้น, ลำดับการทำงาน (pipeline), event_dir(), plane_wave_fit() (+1 more)

### Community 28 - "mmi_from_pgv"
Cohesion: 0.50
Nodes (3): 7.1 MMI ที่วัดได้ที่สถานี — GMICE ของ Worden et al. (2012) (ที่ USGS ShakeMap ใช้), 7.1 MMI ที่วัดได้ที่สถานี — GMICE ของ Worden et al. (2012) (ที่ USGS ShakeMap ใช้), mmi_from_pgv()

### Community 29 - "2. s_wave.py — ตรวจจับคลื่น S และ Lg"
Cohesion: 0.50
Nodes (4): 2.1 หมุนแกน N, E → R (radial), T (transverse), 2.2 Pick ในหน้าต่างรอบเวลาที่คาด (วิธีหลักของระบบ), 2.3 วิธีแบบไม่รู้ตำแหน่ง (blind) — ใช้เมื่อมีสถานีเดียว, 2. s_wave.py — ตรวจจับคลื่น S และ Lg

## Knowledge Gaps
- **1 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `geiger_locate_2layer()` connect `geiger_locate_2layer` to `pipeline.py`, `epicenter.py`, `run_epicenter.py`, `refine_depth`, `mblg_magnitude`?**
  _High betweenness centrality (0.107) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `locate_event()` (e.g. with `Graphify` and `ลำดับการทำงาน (pipeline)`) actually correct?**
  _`locate_event()` has 2 INFERRED edges - model-reasoned connections that need verification._
- **Should `map_plot.py` be split into smaller, more focused modules?**
  _Cohesion score 0.05443371378402107 - nodes in this community are weakly interconnected._
- **Why does `mmi_from_pgv()` connect `mmi_from_pgv` to `mblg_magnitude`, `intensity.py`?**
  _High betweenness centrality (0.099) - this node is a cross-community bridge._
- **Should `s_wave.py` be split into smaller, more focused modules?**
  _Cohesion score 0.10160427807486631 - nodes in this community are weakly interconnected._
- **Why does `สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง` connect `4. magnitude.py — ขนาดแผ่นดินไหว` to `7. intensity.py + run_intensity.py — ความรุนแรงแผ่นดินไหว (MMI แบบ USGS และ JMA) และการแจ้งเตือน`, `geiger_locate_2layer`?**
  _High betweenness centrality (0.092) - this node is a cross-community bridge._
- **Should `pipeline.py` be split into smaller, more focused modules?**
  _Cohesion score 0.14153846153846153 - nodes in this community are weakly interconnected._