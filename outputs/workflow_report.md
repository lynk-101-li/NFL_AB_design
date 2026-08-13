# NfL 抗体从头设计与回顾性对照演示报告

运行时间：`2026-08-13T09:35:56.605698+08:00`

> **证据边界：** 本次生成、结构/界面、可开发性和 sandwich 数值均为确定性 `simulated proxy`。
> RFantibody、IgGM、Germinal、tFold 及后续结构工具均未在本次运行中执行。
> 两株已知抗体在前瞻排名完成后才以 `retrospective_positive_control` 注入；其 Top 2 不是盲法从头发现。

## 1. 输入与设计边界

- 抗原推断资源：`resources/antigen_inference`
- 设计 campaign：`config/design_campaign.json`
- 回顾性阳性对照：`validation/experimentally_validated_antibodies.fasta`

工作抗原采用 NfL rod/coil-2B 的 aa280–377 单链单体上下文；设计热点不包含 Cys322，也不要求抗体接触半胱氨酸。两株已知抗体只在生成轨道中提供两个不同的配对 VH/VL framework；H1/H2/H3/L1/L2/L3 全部遮罩并重新设计。已知 CDR 氨基酸和完整已知 VH/VL 不作为 prospective generation feature。

CDR 坐标由 `ANARCI 2020.04.23 Chothia` 编号后映射到链内 1-based inclusive raw 坐标；模拟生成和真实模型请求共用同一组精确遮罩。编号 labels、工具版本与输入哈希见 `input/antibody_templates/chothia_numbering_evidence.json`。

## 2. 抗原截断与表位

- 全长 NfL 肽键 proxy：`542`
- 中高优先级 cathepsin-like 切点：`87`
- 约束内候选截断片段：`571`
- 生化截断排序第一名：`NEFL 280-375`
- 覆盖全部配置表位的建模上下文：`NEFL 280-377`

|epitope_rank|epitope_id|start|end|sequence|epitope_priority_score|notes|
|---|---|---|---|---|---|---|
|1|C_boundary_368_377|368|377|YLKEYQDLLN|79.58|near inferred C-terminal cathepsin boundary; hydrophobic/coiled-coil contact proxy|
|2|N_boundary_279_290|279|290|WFKSRFTVLTES|78.17|near inferred N-terminal cathepsin boundary|
|3|N_boundary_280_291|280|291|FKSRFTVLTESA|77.84|near inferred N-terminal cathepsin boundary|
|4|sliding_320_331|320|331|EACRGMNEALEK|77.34|sliding-window candidate|
|5|sliding_312_323|312|323|LKAKTLEIEACR|77.09|hydrophobic/coiled-coil contact proxy; sliding-window candidate|
|6|sliding_364_375|364|375|EMARYLKEYQDL|76.9|near inferred C-terminal cathepsin boundary; hydrophobic/coiled-coil contact proxy; sliding-window candidate|
|7|sliding_316_327|316|327|TLEIEACRGMNE|76.75|acidic surface proxy favors basic antibody contacts; sliding-window candidate|
|8|Cys322_core_319_327|319|327|IEACRGMNE|76.37||

## 3. 双模板、六 CDR 与模拟生成

|template_id|framework_source_antibody_id|template_role|design_regions|known_cdr_sequences_used_for_generation|data_status|
|---|---|---|---|---|---|
|template_7-H11-D3-2-C7|7-H11-D3-2-C7|framework_source_only|H1;H2;H3;L1;L2;L3|False|derived_input|
|template_15-C12-H6|15-C12-H6|framework_source_only|H1;H2;H3;L1;L2;L3|False|derived_input|

本地 proxy 生成 `5120` 个 prospective candidates。这不是 RFantibody、IgGM 或 Germinal 的实际生成量。

## 4. 分步筛选漏斗

|stage_order|stage|metric|threshold|input_count|pass_count|removed_count|data_status|
|---|---|---|---|---|---|---|---|
|1|structure_quality|backbone_confidence_score|58.0|5120|4162|958|simulated|
|2|antigen_interface|interface_confidence_score|55.0|4162|2147|2015|simulated|
|3|sequence_developability|developability_score|60.0|2147|1864|283|simulated|
|4|multi_objective_composite|composite_score|60.0|1864|1341|523|simulated|
|5|final_export_shortlist|balanced_template_epitope_then_composite_rank|quota_plus_fill_12|1341|12|1329|simulated|

`06` 和 `07` 表中的分数均保留 `*_is_simulated=True` 与 `metric_provenance`。

## 5. Prospective 分层短名单

`09_prospective_candidates.csv` 不含两株已知阳性全序列。
最终导出采用 `template × epitope` 分层配额；表中的 `rank` 仍是全局模拟分数排名，不能把入选状态解释为纯全局 Top12。

|rank|candidate_id|template_id|best_epitope_id|selection_stratum_rank|binding_confidence_score|developability_score|total_rank_score|selection_reason|
|---|---|---|---|---|---|---|---|---|
|1|DN-e6e8d7351d0d|template_15-C12-H6|C_boundary_368_377|1|72.65|76.32|72.62|balanced_template_epitope_quota|
|2|DN-8df6f8cc4294|template_7-H11-D3-2-C7|C_boundary_368_377|1|82.81|67.77|72.44|balanced_template_epitope_quota|
|3|DN-77576c381425|template_15-C12-H6|C_boundary_368_377|2|72.55|75.97|72.31|balanced_template_epitope_quota|
|4|DN-a9cf24db579e|template_15-C12-H6|C_boundary_368_377|3|74.2|73.54|72.15|balanced_template_epitope_quota|
|5|DN-c0d75e927e22|template_7-H11-D3-2-C7|C_boundary_368_377|2|73.9|78.01|72.14|balanced_template_epitope_quota|
|7|DN-2ea27b071345|template_7-H11-D3-2-C7|C_boundary_368_377|3|73.91|76.35|71.71|balanced_template_epitope_quota|
|78|DN-db51a6170580|template_7-H11-D3-2-C7|helix_surface_323_331|1|62.53|77.5|68.8|balanced_template_epitope_quota|
|87|DN-870775f54c72|template_15-C12-H6|helix_surface_323_331|1|61.39|76.55|68.68|balanced_template_epitope_quota|
|279|DN-98a428e115e9|template_7-H11-D3-2-C7|helix_surface_323_331|2|60.16|75.9|66.44|balanced_template_epitope_quota|
|318|DN-6d8910dd2d15|template_15-C12-H6|helix_surface_323_331|2|65.54|75.97|66.0|balanced_template_epitope_quota|
|368|DN-a78bf47c4ea2|template_7-H11-D3-2-C7|helix_surface_323_331|3|67.76|66.29|65.65|balanced_template_epitope_quota|
|401|DN-9c287784eedf|template_15-C12-H6|helix_surface_323_331|3|60.21|75.21|65.43|balanced_template_epitope_quota|

## 6. Retrospective 阳性对照 Top 2

已知阳性仅在 prospective ranking 完成后注入，并用明确的回顾性独立证据字段给分。

|rank|candidate_id|control_status|best_epitope_id|independent_evidence_score|independent_evidence_provenance|total_rank_score|
|---|---|---|---|---|---|---|
|1|7-H11-D3-2-C7|retrospective_positive_control|helix_surface_323_331|97.0|known_positive_status_retrospective_demo_not_blind_discovery|90.08|
|2|15-C12-H6|retrospective_positive_control|C_boundary_368_377|94.0|known_positive_status_retrospective_demo_not_blind_discovery|88.75|

## 7. Sandwich pair 模拟优先级

- Top pair：`7-H11-D3-2-C7` + `15-C12-H6`
- 表位重叠：`0.0`
- 线性间隔：`36` aa
- 兼容性 proxy：`96.43`
- 建议 capture/detection：`7-H11-D3-2-C7` / `15-C12-H6`

|pair_rank|antibody_1|antibody_2|epitope_overlap_ratio|linear_epitope_gap_aa|sandwich_compatibility_score|data_status|claim_scope|
|---|---|---|---|---|---|---|---|
|1|7-H11-D3-2-C7|15-C12-H6|0.0|36|96.43|simulated|retrospective_demo|
|2|7-H11-D3-2-C7|DN-e6e8d7351d0d|0.0|36|88.38|simulated|retrospective_demo|
|3|7-H11-D3-2-C7|DN-c0d75e927e22|0.0|36|88.32|simulated|retrospective_demo|
|4|7-H11-D3-2-C7|DN-77576c381425|0.0|36|88.26|simulated|retrospective_demo|
|5|7-H11-D3-2-C7|DN-2ea27b071345|0.0|36|88.09|simulated|retrospective_demo|
|6|7-H11-D3-2-C7|DN-a9cf24db579e|0.0|36|88.08|simulated|retrospective_demo|
|7|7-H11-D3-2-C7|DN-9e5c3708b2c2|0.0|36|87.98|simulated|retrospective_demo|
|8|7-H11-D3-2-C7|DN-80ac0338f825|0.0|36|87.88|simulated|retrospective_demo|

## 8. 真实模型 Handoff

已产生六 CDR 遮罩模板与 RFantibody/IgGM/Germinal 规范化请求，但当前缺少经验证的抗原 PDB、坐标映射、模型 runtime 和 checkpoint，所以保持 `not_run/blocked` 状态。Germinal 是独立 scFv 轨道，不视为 native paired-Fv 结果。

- 遮罩模板：`outputs/exports/fasta/design_templates_six_cdr_masked.fasta`
- 请求索引：`outputs/exports/design_requests/design_request_index.json`
- job table：`outputs/exports/external_jobs/pipeline_jobs.tsv`
- command sheet：`outputs/exports/external_jobs/run_external_pipelines.sh`

## 9. 下一步才能取代 proxy 的证据

- RFantibody/IgGM/Germinal 真实生成结果、日志、版本和 checkpoint。
- Fv/Fab 结构质量、复合物 PAE/ipTM/pTM/DockQ、埋藏表面积和界面能量。
- 亲和力、特异性、交叉反应、可开发性和 sandwich assay 实验。
