# NfL 抗体设计计算流程报告

生成日期：2026-07-08

## 1. 输入材料

- 故事线：`resources/project_context/storyline.txt`
- 研究计划：`resources/project_context/research_plan.txt`
- 抗原推断目录：`resources/antigen_inference`
- 已验证抗体 FASTA：`validation/experimentally_validated_antibodies.fasta`

本流程先执行 NfL 抗原截断推断，再进入表位、抗体和 sandwich pair 计算。当前工作假设是 NfL rod/coil-2B 的 aa 280-377 附近片段包含 Cys322，并形成二硫键同源二聚体。

## 2. 计算模块

1. 抗原截断推断：对 NfL 全长逐肽键进行 cathepsin-like 打分，并按 Cys322 + 22 kDa 二硫键二聚体约束枚举片段。
2. 抗原结构可靠性分层：用 rod/coil-2B 区域和截断推断结果定义 NfL aa 280-377 的优先抗原上下文。
3. NfL 特异性表位图谱：在候选片段内生成边界/Cys322/滑窗表位，并用暴露度、带电性、rod 稳定性和 PTM 风险代理指标打分。
4. 抗体序列建模准备：解析 VH/VL，做启发式 CDR 注释和 developability 体检。
5. 计算复现筛选：把两株真实抗体和 CDR 扰动阴性候选放在同一评分体系下排序。
6. Sandwich pair 兼容性：对两株真实抗体做 pair-aware 表位分配，计算表位重叠、线性距离和 clash 代理指标。
7. 外部结构工具衔接：导出 Fv/抗原复合物 FASTA、AF3-style JSON、job table 和 shell 任务模板。

## 3. 抗原截断推断

- 全长 NfL 肽键打分数：`542`。
- 中高优先级 cathepsin-like 切点数：`87`。
- 满足 Cys322 + 22 kDa 二硫键二聚体约束的候选片段数：`571`。
- 详细报告：`00_antigen_truncation_report.md`。

|fragment|length_aa|monomer_avg_mass_kDa|disulfide_homodimer_avg_mass_kDa|N_terminal_cut|C_terminal_cut|combined_boundary_score|mass_error_from_22kDa|
|---|---|---|---|---|---|---|---|
|280-375|96|11.041|22.081|W279\|F280|L375\|L376|10.2|0.081|
|282-377|96|10.993|21.984|K281\|S282|N377\|V378|8.4|0.016|
|281-376|96|11.007|22.013|F280\|K281|L376\|N377|7.3|0.013|
|280-374|95|10.928|21.854|W279\|F280|D374\|L375|8.7|0.146|
|280-376|97|11.154|22.307|W279\|F280|L376\|N377|8.9|0.307|
|281-375|95|10.894|21.786|F280\|K281|L375\|L376|8.6|0.214|
|282-376|95|10.879|21.756|K281\|S282|L376\|N377|8.6|0.244|
|282-375|94|10.766|21.53|K281\|S282|L375\|L376|9.9|0.47|

## 4. 抗原片段优先级

|antigen_rank|fragment|length_aa|disulfide_homodimer_avg_mass_kDa|combined_boundary_score|N_terminal_cut|C_terminal_cut|antigen_confidence_score|
|---|---|---|---|---|---|---|---|
|1|280-375|96|22.081|10.2|W279\|F280|L375\|L376|98.25|
|2|282-377|96|21.984|8.4|K281\|S282|N377\|V378|93.96|
|3|281-376|96|22.013|7.3|F280\|K281|L376\|N377|90.77|
|4|272-367|96|21.98|8.9|K271\|N272|R367\|Y368|87.37|
|5|272-368|97|22.307|9.8|K271\|N272|Y368\|L369|86.08|
|6|298-392|95|22.031|8.6|A297\|V298|L392\|L393|85.82|
|7|273-368|96|22.079|8.6|N272\|M273|Y368\|L369|85.57|
|8|269-365|97|22.067|8.3|L268\|A269|M365\|A366|85.38|

## 5. 候选表位窗口

|epitope_rank|epitope_id|start|end|sequence|epitope_priority_score|notes|
|---|---|---|---|---|---|---|
|1|Cys322_anchor_316_331|316|331|TLEIEACRGMNEALEK|79.76|contains Cys322 disulfide-anchor region; acidic surface proxy favors basic antibody contacts|
|2|C_boundary_368_377|368|377|YLKEYQDLLN|79.58|near inferred C-terminal cathepsin boundary; hydrophobic/coiled-coil contact proxy|
|3|Cys322_core_319_327|319|327|IEACRGMNE|79.37|contains Cys322 disulfide-anchor region|
|4|N_boundary_279_290|279|290|WFKSRFTVLTES|78.17|near inferred N-terminal cathepsin boundary|
|5|N_boundary_280_291|280|291|FKSRFTVLTESA|77.84|near inferred N-terminal cathepsin boundary|
|6|sliding_320_331|320|331|EACRGMNEALEK|77.34|contains Cys322 disulfide-anchor region; sliding-window candidate|
|7|sliding_312_323|312|323|LKAKTLEIEACR|77.09|contains Cys322 disulfide-anchor region; hydrophobic/coiled-coil contact proxy; sliding-window candidate|
|8|sliding_364_375|364|375|EMARYLKEYQDL|76.9|near inferred C-terminal cathepsin boundary; hydrophobic/coiled-coil contact proxy; sliding-window candidate|
|9|sliding_316_327|316|327|TLEIEACRGMNE|76.75|contains Cys322 disulfide-anchor region; acidic surface proxy favors basic antibody contacts; sliding-window candidate|
|10|sliding_280_291|280|291|FKSRFTVLTESA|74.84|near inferred N-terminal cathepsin boundary; sliding-window candidate|
|11|sliding_352_363|352|363|NKLENELRTTKS|73.78|sliding-window candidate|
|12|sliding_360_371|360|371|TTKSEMARYLKE|72.96|near inferred C-terminal cathepsin boundary; sliding-window candidate|

## 6. 已验证抗体序列体检

|antibody_id|HCDR3_proxy|LCDR3_proxy|n_glyco_risk_notes|hydrophobic_patch_proxy|developability_score|
|---|---|---|---|---|---|
|7-H11-D3-2-C7|TRKDY|LQLYSTPLT|none|0.222|92.0|
|15-C12-H6|ATSLLRLRDWFPY|QQTNTWPYT|VH:57:NTS:CDR;VL:41:NGS:FR|0.23|69.87|

## 7. 计算复现排序

真实抗体在包含 CDR 扰动 decoy 的候选库中的排名：

|rank|candidate_id|best_epitope_id|binding_confidence_score|developability_score|total_rank_score|top_tier|
|---|---|---|---|---|---|---|
|1|7-H11-D3-2-C7|C_boundary_368_377|74.77|92.0|81.55|yes|
|2|15-C12-H6|Cys322_anchor_316_331|79.99|69.87|81.32|yes|

Top 10 候选：

|rank|candidate_id|experimental_status|best_epitope_id|binding_confidence_score|developability_score|off_target_penalty_proxy|total_rank_score|
|---|---|---|---|---|---|---|---|
|1|7-H11-D3-2-C7|validated|C_boundary_368_377|74.77|92.0|0.0|81.55|
|2|15-C12-H6|validated|Cys322_anchor_316_331|79.99|69.87|3.0|81.32|
|3|15-C12-H6-cdr-perturb-14|in_silico_decoy|Cys322_anchor_316_331|75.04|69.75|3.6|78.25|
|4|7-H11-D3-2-C7-cdr-perturb-14|in_silico_decoy|C_boundary_368_377|69.23|92.14|0.6|78.23|
|5|15-C12-H6-cdr-perturb-01|in_silico_decoy|Cys322_anchor_316_331|73.36|69.75|4.2|76.76|
|6|15-C12-H6-cdr-perturb-07|in_silico_decoy|Cys322_anchor_316_331|71.4|69.87|3.6|76.52|
|7|7-H11-D3-2-C7-cdr-perturb-07|in_silico_decoy|C_boundary_368_377|67.58|85.0|3.6|76.34|
|8|7-H11-D3-2-C7-cdr-perturb-15|in_silico_decoy|C_boundary_368_377|66.7|92.09|1.2|76.33|
|9|7-H11-D3-2-C7-cdr-perturb-01|in_silico_decoy|C_boundary_368_377|66.1|92.14|1.2|76.05|
|10|15-C12-H6-cdr-perturb-08|in_silico_decoy|Cys322_anchor_316_331|70.63|69.83|4.2|75.46|

## 8. Sandwich Pair 结论

- 抗体 1：`7-H11-D3-2-C7`，表位 `C_boundary_368_377`。
- 抗体 2：`15-C12-H6`，表位 `Cys322_anchor_316_331`。
- 表位重叠比例：`0.0`。
- 线性表位间隔：`36` aa。
- sandwich 兼容性代理评分：`84.52`。
- 建议 capture：`7-H11-D3-2-C7`。
- 建议 detection：`15-C12-H6`。

## 9. 外部结构工具 Handoff

外部工具输入已经写入：

- `outputs/exports/fasta` for FASTA templates
- `outputs/exports/af3_json` for AF3-style JSON templates
- `outputs/exports/external_jobs/pipeline_jobs.tsv` for external pipeline jobs
- `outputs/exports/external_jobs/run_external_pipelines.sh` for editable command sheet

当前流程中的以下指标是代理指标，建议由外部结构或实验结果替换：

- Fv/Fab 模型质量：应替换为 IgFold/ABodyBuilder3 的模型质量、CDR loop 收敛性和 VH/VL packing 结果。
- 复合物可信度：应替换为 AF3/Chai-1/Boltz 的 ipTM、pTM、interface PAE、interface pLDDT、pDockQ/DockQ。
- 界面物理量：应替换为 buried surface area、Rosetta interface ΔG、shape complementarity、氢键/盐桥和 buried unsatisfied polar atoms。
- sandwich 空间兼容性：应替换为 Fab1:NfL:Fab2 三元复合物结构的实际 clash、Fab-Fab 最小距离和标记端可及性。

## 10. 方法学边界

本流程用于计算复现、候选排序和外部结构预测输入准备。候选库中的扰动序列、表位分配和多目标评分均为确定性代理计算；正式结论应结合真实结构建模、亲和力测定、交叉反应性实验和 sandwich assay 数据。
