# Graph Report - code  (2026-11-07)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 691 nodes · 1598 edges · 29 communities (28 shown, 1 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 11 edges (avg confidence: 0.91)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `de4bb864`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- Community 0
- Community 1
- Community 2
- Community 3
- Community 4
- Community 5
- Community 6
- Community 7
- Community 8
- Community 9
- Community 10
- Community 11
- Community 12
- Community 13
- Community 14
- Community 15
- Community 16
- Community 17
- Community 18
- Community 19
- Community 20
- Community 21
- Community 22
- Community 23
- Community 24
- Community 25
- Community 26
- Community 27
- Community 28

## God Nodes (most connected - your core abstractions)
1. `haversine_km()` - 42 edges
2. `detect_p_wave()` - 22 edges
3. `azimuth_deg()` - 22 edges
4. `bandpass()` - 21 edges
5. `locate_event()` - 20 edges
6. `detect_s_wave()` - 18 edges
7. `load_station_txt()` - 17 edges
8. `phase_travel_time()` - 17 edges
9. `predict_arrival()` - 17 edges
10. `locate()` - 17 edges

## Surprising Connections (you probably didn't know these)
- `3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`)` --references--> `geiger_locate_2layer()`  [INFERRED]
  EEW_formulas.md → epicenter.py
- `3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`)` --references--> `geiger_locate_2layer()`  [INFERRED]
  FORMULAS.md → epicenter.py
- `ลำดับการทำงาน (pipeline)` --references--> `EEWSystem`  [INFERRED]
  CLAUDE.md → eew_system.py
- `3.1 เวลาเดินทาง: แบบจำลองโลก 2 ชั้น (เปลือกโลกหนา H บนเนื้อโลก)` --references--> `omori_distance()`  [INFERRED]
  EEW_formulas.md → epicenter.py
- `3.1 เวลาเดินทาง: แบบจำลองโลก 2 ชั้น (เปลือกโลกหนา H บนเนื้อโลก)` --references--> `omori_distance()`  [INFERRED]
  FORMULAS.md → epicenter.py

## Import Cycles
- None detected.

## Communities (29 total, 1 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.06
Nodes (41): _base(), build_sites(), _candidates(), compute(), correction_model(), despike(), detail_report(), _epi_html() (+33 more)

### Community 1 - "Community 1"
Cohesion: 0.06
Nodes (30): _circle_coords(), _finish(), _great_circle(), _kml_color(), make_epicenter_map(), make_intensity_map(), make_mmi_map(), make_multi_station_map() (+22 more)

### Community 2 - "Community 2"
Cohesion: 0.06
Nodes (29): aic_pick(), allen_cf(), _cfg(), detect_p_wave(), first_motion(), first_real_trigger(), first_sustained(), is_sustained() (+21 more)

### Community 3 - "Community 3"
Cohesion: 0.09
Nodes (26): phase_name(), predict_arrival(), _as_folder(), associate(), _cfg(), check_time_coverage(), delta_from_differential(), expand_specs() (+18 more)

### Community 4 - "Community 4"
Cohesion: 0.07
Nodes (21): action_text(), alert_level(), distance_to_trace_km(), fault_distance_from_hypo(), gmpe_weight(), intensity_from_pgv(), ipe_allen2012(), _jma_filter() (+13 more)

### Community 5 - "Community 5"
Cohesion: 0.06
Nodes (33): 0. signal_utils.py — ฟังก์ชันพื้นฐาน, 1.1 Characteristic Function ของ Allen (1978), 1.2 STA/LTA, 1.3 กันการ trigger บน noise burst (รถวิ่งผ่าน ประตูปิด ฯลฯ), 1.4 ปรับเวลาให้แม่นด้วย AIC (Maeda, 1985) — 2 รอบ, 1.5 คุณภาพ pick (แบบ Hypoinverse, Klein 2002), 1.6 ทิศการเคลื่อนที่แรก (First motion), 1.7 Pick คลื่น P ในหน้าต่าง (re-pick หลังรู้ตำแหน่งคร่าว ๆ) (+25 more)

### Community 6 - "Community 6"
Cohesion: 0.11
Nodes (18): analyse(), blind_zone_km(), fmt(), magnitude_ready_time(), main(), make_map(), phase_name(), plot() (+10 more)

### Community 7 - "Community 7"
Cohesion: 0.09
Nodes (14): choose_preferred(), _distance_factor(), network_by_type(), network_magnitude(), _pick_factor(), preferred_magnitude(), reliability(), rupture_duration() (+6 more)

### Community 8 - "Community 8"
Cohesion: 0.14
Nodes (20): compute(), covers(), detail_report(), _epoch(), _fmt(), _folder_of(), locate(), location_lines() (+12 more)

### Community 9 - "Community 9"
Cohesion: 0.14
Nodes (15): _check_index(), _clip(), mblg_from_amplitude(), mblg_magnitude(), ml_from_amplitude(), ml_magnitude(), msbb_from_velocity(), msbb_magnitude() (+7 more)

### Community 10 - "Community 10"
Cohesion: 0.11
Nodes (11): CLAUDE.md — โปรเจ็ค EEW (Earthquake Early Warning) สำหรับประเทศไทย, Graphify (ยังไม่ได้ติดตั้ง), การรัน, ข้อควรระวังเวลาแก้โค้ด, ข้อมูล, ค่าตั้งต้น, ลำดับการทำงาน (pipeline), EEWSystem (+3 more)

### Community 11 - "Community 11"
Cohesion: 0.14
Nodes (11): find_station_file(), resolve_data_path(), delta_from_sp(), analyze_station(), choose_consistent(), choose_nearest(), compare_reference(), main() (+3 more)

### Community 12 - "Community 12"
Cohesion: 0.14
Nodes (11): jma_magnitude(), jma_magnitude_from_amplitude(), legacy_magnitudes(), magnitude_from_pd(), magnitude_from_tauc(), max_horizontal_displacement_um(), pd_value(), tau_c() (+3 more)

### Community 13 - "Community 13"
Cohesion: 0.09
Nodes (21): 3.1 เวลาเดินทาง: แบบจำลองโลก 2 ชั้น (เปลือกโลกหนา H บนเนื้อโลก), 3.2 ตรวจเวลาไฟล์ (Time coverage), 3.3 Association — เวลา P ของสองสถานีต้องเข้ากันได้ทางกายภาพ, 3.4 ตำแหน่งเริ่มต้น: scan จาก P + ทิศ + เวลา Lg + พลังงาน S (grid search), 3.5 Geiger's method บนผิวโลกทรงกลม + โลก 2 ชั้น (`geiger_locate_2layer`), 3.6 ขั้นตอนทั้งหมดของ pipeline, 3.7 ตรวจทาน residual (ตัด pick ที่ไม่เข้ากับตำแหน่ง), 3.8 Azimuthal gap (+13 more)

### Community 14 - "Community 14"
Cohesion: 0.13
Nodes (8): geiger_locate(), geiger_residuals(), observations(), pick_weight(), _robust_misfit(), _travel_time_and_derivs(), latlon_to_xy(), xy_to_latlon()

### Community 15 - "Community 15"
Cohesion: 0.15
Nodes (10): alert_time_from_picks(), blind_zone_radius(), _delta_from_R(), hypocentral_distance(), lg_arrival_time(), p_arrival_time(), phase_name(), s_arrival_time() (+2 more)

### Community 16 - "Community 16"
Cohesion: 0.12
Nodes (9): _apply_response(), broadband_velocity(), _butter_mag(), displacement(), _osc(), wood_anderson(), wwssn_sp(), _wwssn_sp_disp_resp() (+1 more)

### Community 17 - "Community 17"
Cohesion: 0.21
Nodes (9): geiger_locate_2layer(), baz_weights(), tcalc(), phase_travel_time(), travel_time_2layer(), scan_locate(), cost(), azimuth_deg() (+1 more)

### Community 18 - "Community 18"
Cohesion: 0.24
Nodes (8): _axial_diff_deg(), azimuthal_gap(), _fit_t0_truncated(), grid_search_locate(), cost(), predict(), residuals(), search()

### Community 19 - "Community 19"
Cohesion: 0.17
Nodes (12): 4.0 จำลองเครื่องวัดมาตรฐานในโดเมนความถี่ (s = 2πif), 4.10 คำบอกขนาด (การแบ่งกลุ่มของ USGS), 4.1 ML — Local magnitude (IASPEI 2013 = Hutton & Boore, 1987), 4.2 mb_Lg — Lg magnitude (Nuttli, 1973), 4.3 Ms_BB — Broadband surface-wave magnitude (IASPEI 2013), 4.4 Mwp — Moment magnitude จากคลื่น P (Tsuboi et al., 1995; Whitmore et al., 2002), 4.5 ขนาดแบบเดิม / สำหรับเตือนภัยเร็ว (ไม่ใช่มาตรฐาน USGS ใช้เทียบ), 4.6 น้ำหนักของแต่ละสถานี (0–1) (+4 more)

### Community 20 - "Community 20"
Cohesion: 0.29
Nodes (6): check_units(), detect_sensor_type(), detect_station_type(), load_station_txt(), _parse_time(), differentiate()

### Community 21 - "Community 21"
Cohesion: 0.22
Nodes (9): 7.0 ภาพรวม, 7.2 JMA Instrumental Intensity (คำนวณจากคลื่นความเร่งจริงที่สถานี), 7.3 ทำนาย MMI — วิธี ก) GMPE → GMICE, 7.4 ทำนาย MMI — วิธี ข) IPE ของ Allen, Wald & Worden (2012) (ที่ USGS PAGER ใช้) แบบระยะ hypocentral, 7.5 รวม 2 วิธี และความไม่แน่นอน, 7.6 ตรวจคุณภาพข้อมูล (QC) และปรับด้วยค่าที่วัดได้จริง, 7.7 จำลองการทำงานแบบเวลาจริง (run_intensity.py), 7.8 คำอธิบายระดับ MMI (ตามตาราง USGS ShakeMap)  ระดับ = ค่าที่ปัดเป็นจำนวนเต็ม (+1 more)

### Community 22 - "Community 22"
Cohesion: 0.22
Nodes (9): 1.1 Characteristic Function ของ Allen (1978), 1.2 STA/LTA, 1.3 กันการ trigger บน noise burst (รถวิ่งผ่าน ประตูปิด ฯลฯ), 1.4 ปรับเวลาให้แม่นด้วย AIC (Maeda, 1985) — 2 รอบ, 1.5 คุณภาพ pick (แบบ Hypoinverse, Klein 2002), 1.6 ทิศการเคลื่อนที่แรก (First motion), 1.7 Pick คลื่น P ในหน้าต่าง (re-pick หลังรู้ตำแหน่งคร่าว ๆ), 1.8 การแปลงหน่วยข้อมูลจริง (data_loader.py) (+1 more)

### Community 23 - "Community 23"
Cohesion: 0.22
Nodes (9): 7.0 ภาพรวม, 7.2 JMA Instrumental Intensity (คำนวณจากคลื่นความเร่งจริงที่สถานี), 7.3 ทำนาย MMI — วิธี ก) GMPE → GMICE, 7.4 ทำนาย MMI — วิธี ข) IPE ของ Allen, Wald & Worden (2012) (ที่ USGS PAGER ใช้) แบบระยะ hypocentral, 7.5 รวม 2 วิธี และความไม่แน่นอน, 7.6 ตรวจคุณภาพข้อมูล (QC) และปรับด้วยค่าที่วัดได้จริง, 7.7 จำลองการทำงานแบบเวลาจริง (run_intensity.py), 7.8 คำอธิบายระดับ MMI (ตามตาราง USGS ShakeMap)  ระดับ = ค่าที่ปัดเป็นจำนวนเต็ม (+1 more)

### Community 24 - "Community 24"
Cohesion: 0.25
Nodes (7): 0. signal_utils.py — ฟังก์ชันพื้นฐาน, 5. run_magnitude.py — การ์ดแจ้งเตือน, 6. travel_time.py + run_travel_time.py — เวลาที่คลื่นเดินทางมาถึงตำแหน่งของเรา และเวลาเตือนล่วงหน้า, 8. เกณฑ์การออกคำเตือน — ShakeAlert (USGS), ข้อจำกัดที่ควรทราบ, สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง, เอกสารอ้างอิง

### Community 25 - "Community 25"
Cohesion: 0.29
Nodes (4): delta_from_sp_lg(), locate_robust(), relabel_as_lg(), _w()

### Community 26 - "Community 26"
Cohesion: 0.50
Nodes (3): 7.1 MMI ที่วัดได้ที่สถานี — GMICE ของ Worden et al. (2012) (ที่ USGS ShakeMap ใช้), 7.1 MMI ที่วัดได้ที่สถานี — GMICE ของ Worden et al. (2012) (ที่ USGS ShakeMap ใช้), mmi_from_pgv()

### Community 28 - "Community 28"
Cohesion: 0.50
Nodes (4): 2.1 หมุนแกน N, E → R (radial), T (transverse), 2.2 Pick ในหน้าต่างรอบเวลาที่คาด (วิธีหลักของระบบ), 2.3 วิธีแบบไม่รู้ตำแหน่ง (blind) — ใช้เมื่อมีสถานีเดียว, 2. s_wave.py — ตรวจจับคลื่น S และ Lg

## Knowledge Gaps
- **1 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `geiger_locate_2layer()` connect `Community 17` to `Community 9`, `Community 3`, `Community 13`, `Community 14`?**
  _High betweenness centrality (0.133) - this node is a cross-community bridge._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.060153776571687016 - nodes in this community are weakly interconnected._
- **Why does `mmi_from_pgv()` connect `Community 26` to `Community 9`, `Community 4`?**
  _High betweenness centrality (0.111) - this node is a cross-community bridge._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.05501165501165501 - nodes in this community are weakly interconnected._
- **Why does `สูตรที่ใช้ในระบบ EEW (Earthquake Early Warning) — ฉบับปรับปรุง` connect `Community 5` to `Community 21`, `Community 13`?**
  _High betweenness centrality (0.094) - this node is a cross-community bridge._
- **Should `Community 2` be split into smaller, more focused modules?**
  _Cohesion score 0.061507936507936505 - nodes in this community are weakly interconnected._
- **Should `Community 3` be split into smaller, more focused modules?**
  _Cohesion score 0.08686868686868687 - nodes in this community are weakly interconnected._