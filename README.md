# NFL_AB_design

`NFL_AB_design` 是一个自包含的 NfL/NEFL 抗体设计课题代码包。它把课题主线整理成一条可重复运行的计算流程：

```text
抗原推断 -> 抗原处理 -> 抗体设计/复现 -> 打分过滤 -> sandwich pair 评估 -> 外部结构管线交接
```

本仓库包含项目背景、抗原截断推断资料、NfL 参考序列与注释、起始设计约束、抗体模板、实验验证抗体、可执行代码和示例输出。只看这个 repo，不依赖外部文件，也能理解课题为什么这样设计、代码如何运行、结果如何解释。

当前实现使用确定性的代理评分，用于方法复现、流程审计和外部结构工具输入准备。它不能替代真实结构预测、亲和力测定、特异性实验、可开发性实验或 sandwich assay 验证。

## Project Background

NfL，也称 NEFL，是神经丝轻链蛋白，常作为神经损伤相关检测标志物。本课题的目标不是直接从全长 NfL 随机挑选表位，而是先根据上游观察推断更可能被检测到的抗原片段，再围绕该片段设计和筛选抗体。

课题起点是一个约 22 kDa 的非还原条件 NfL 相关条带。该条带被解释为二硫键连接的同源二聚体，因此每个单体片段约为 11 kDa。人源 NfL 全长序列中只有一个 canonical cysteine，即 Cys322。如果 22 kDa 条带确实来自 NfL 二硫键二聚体，那么相关单体片段必须包含 Cys322。

这个推理把抗原空间从“全长 NfL”收缩到 rod/coil-2B 附近的 Cys322-containing fragment。结合片段质量、cathepsin-like 边界支持和结构区域可靠性，工作抗原集中到 NfL aa 280-377 附近，核心候选包括：

- `NEFL 280-375`
- `NEFL 281-376`
- `NEFL 282-377`

最终实验验证通过的两株抗体作为 validation/replay set 存放在 `validation/`，用于检查整套计算流程是否能把它们排在最前面：

- `7-H11-D3-2-C7`
- `15-C12-H6`

## End-to-End Logic

本仓库不是只对最终抗体打分，而是覆盖完整课题过程。

1. **抗原推断**：根据 22 kDa 非还原条带、11 kDa 单体假设、Cys322 必须包含、cathepsin-like 边界支持，枚举并排序 NfL 截断片段。
2. **抗原处理**：选择优先抗原片段，生成 N 端边界、Cys322 anchor、C 端边界和片段内部滑动窗口，作为后续表位分析对象。
3. **抗体设计/复现**：`input/antibody_templates/` 提供中性 Fv 模板格式；`validation/` 中的两株实验阳性抗体用于复现课题后期筛选结果，并与扰动 decoy、negative control 一起进入候选库。
4. **打分过滤**：综合抗原可信度、表位适配性、CDR 接触潜力、HCDR3 几何代理、序列可开发性、off-target penalty 和模型一致性代理，对候选抗体排序。
5. **Sandwich pair 评估**：比较两株抗体的分配表位、重叠比例、线性距离、空间冲突代理和 capture/detection 方向，判断是否值得进入 Fab1:NfL:Fab2 三元结构预测和实验验证。
6. **外部管线交接**：导出 FASTA、AF3-style JSON、任务表和可编辑命令脚本，交给 IgFold、ABodyBuilder3、AlphaFold 3、Chai-1、Boltz、Rosetta 等工具继续分析。

## Input, Resources, Validation

为避免把“最终答案”误放在 `input/`，仓库将数据分成三类。

`input/` 放起始假设和模板：

- `input/antigen_truncation/truncation_constraints.json`：抗原截断模型约束，包括 Cys322、22 kDa 非还原条带、单体质量范围和核心片段假设。
- `input/antibody_templates/template_fv_backgrounds.fasta`：前瞻性抗体设计可使用的 VH/VL 模板背景，展示输入格式。

`resources/` 放课题背景和推断资料：

- `resources/project_context/storyline.txt`：课题设计故事线。
- `resources/project_context/research_plan.txt`：研究方案。
- `resources/antigen_inference/`：NfL 截断推断过程、GenPept 注释、cathepsin-like 位点表和候选片段表。

`validation/` 放最终实验验证结果：

- `validation/experimentally_validated_antibodies.fasta`：两株实验阳性抗体的 VH/VL 序列。

FASTA header 采用：

```text
>{clone_id}|{VH_or_VL}|{sequence_id}
```

每个 clone 必须同时包含 VH 和 VL。

## Repository Layout

```text
NFL_AB_design/
  pyproject.toml
  README.md
  input/
    antigen_truncation/
      truncation_constraints.json
    antibody_templates/
      template_fv_backgrounds.fasta
      README.md
  resources/
    project_context/
      storyline.txt
      research_plan.txt
    antigen_inference/
      README.md
      NFL_22kDa_disulfide_dimer_cathepsin_truncation.md
      nfl_cathepsin_annotated_for_snapgene.gp
      nfl_all_peptide_bonds_cathepsin_scores.csv
      nfl_medium_high_cathepsin_candidate_sites.csv
      nfl_22kda_disulfide_dimer_fragment_candidates.csv
      nfl_snapgene_cathepsin_notes.md
  validation/
    experimentally_validated_antibodies.fasta
  docs/
    high_school_quickstart_zh.md
    nfl_truncation_inference_rationale.md
  config/
    external_pipelines.example.json
  scripts/
    run_nfl_ab_design.py
  src/
    nfl_ab_design/
      __init__.py
      __main__.py
      workflow.py
  outputs/
    00_antigen_truncation_all_peptide_bonds.csv
    00_antigen_truncation_medium_high_sites.csv
    00_antigen_truncation_fragment_candidates.csv
    00_antigen_truncation_report.md
    01_antigen_fragment_prioritization.csv
    02_epitope_windows.csv
    03_antibody_developability.csv
    04_candidate_library.csv
    05_candidate_ranking.csv
    06_sandwich_pair_report.md
    workflow_report.md
    exports/
```

`outputs/` 是可重复生成的示例结果。重新运行工作流会刷新这些文件。

## Method Details

### 1. 抗原截断推断

工作流从 `resources/antigen_inference/nfl_cathepsin_annotated_for_snapgene.gp` 读取 NfL 序列，并对每个 peptide bond 记录：

- P4-P4' 局部序列环境；
- 所属 NfL 区域；
- cysteine cathepsin-like score；
- cathepsin D/E-like score；
- best cleavage score；
- priority label。

随后代码配对 N 端和 C 端边界，保留满足下列条件的片段：

- 包含 Cys322；
- 单体质量落在配置范围内；
- 二硫键同源二聚体质量接近 22 kDa；
- 两端边界具有 cathepsin-like 支持；
- 与 NfL rod/coil-2B 的 aa 280-377 工作区域一致。

详细动机和推理链见：

```text
docs/nfl_truncation_inference_rationale.md
```

### 2. 抗原片段和表位窗口

优先抗原片段进入表位窗口生成模块。窗口类型包括：

- Cys322 anchor region；
- N-terminal boundary region；
- C-terminal boundary region；
- primary antigen fragment 内的 sliding windows。

每个窗口按 rod-domain placement、Cys322/boundary support、极性、带电性、疏水/coil contact proxy、low-complexity penalty 和 PTM risk 打分。

### 3. 抗体候选和可开发性检查

代码读取 VH/VL 序列，使用确定性规则标注近似 CDR 区域，并计算：

- HCDR3 长度和几何代理；
- CDR charge；
- CDR aromaticity；
- NXS/T glycosylation motif；
- deamidation-prone motif；
- oxidation-prone residue；
- hydrophobic patch proxy；
- overall developability score。

两株实验阳性抗体与 CDR-perturbed decoys、negative controls 一起进入候选库，统一排序。

### 4. 多目标排序

候选总分综合：

- antigen confidence；
- epitope priority；
- CDR charge complementarity；
- aromatic contact potential；
- HCDR3 geometry proxy；
- model-consensus proxy；
- developability score；
- off-target penalty。

这些值是透明的 proxy score，不应写成真实亲和力或真实结构能量。

### 5. Sandwich pair 兼容性

pair 模块比较两株抗体的 assigned epitope windows，输出：

- epitope overlap ratio；
- linear spacing；
- clash proxy；
- pair compatibility score；
- capture/detection orientation recommendation。

该结果用于决定是否进入三元复合物结构预测和实验 sandwich assay 确认。

## Installation

本包没有必需的第三方 Python 依赖。推荐使用 Python 3.10 或更高版本。

从仓库根目录运行：

```bash
python3 scripts/run_nfl_ab_design.py
```

也可以用模块方式运行：

```bash
PYTHONPATH=src python3 -m nfl_ab_design
```

可选 editable install：

```bash
python3 -m pip install -e .
nfl-ab-design
```

## Testing

本仓库包含最小回归测试，确认工作流能复现当前核心结论：

- primary antigen fragment 为 `280-375`；
- `7-H11-D3-2-C7` 排名第 1；
- `15-C12-H6` 排名第 2；
- 两株抗体被推荐进入 sandwich pair 结构后续分析。

运行测试：

```bash
python3 -m unittest discover -s tests
```

GitHub Actions 配置位于：

```text
.github/workflows/ci.yml
```

## Outputs

运行后会刷新 `outputs/`。

核心结果：

- `workflow_report.md`：中文总报告，按流程汇总抗原推断、表位窗口、抗体排序、sandwich pair 和外部交接。
- `00_antigen_truncation_all_peptide_bonds.csv`：NfL 每个 peptide bond 的 cathepsin-like 分数。
- `00_antigen_truncation_medium_high_sites.csv`：中高优先级切割位点。
- `00_antigen_truncation_fragment_candidates.csv`：包含 Cys322 且兼容 22 kDa 二聚体假设的候选片段。
- `00_antigen_truncation_report.md`：抗原截断推断报告。
- `01_antigen_fragment_prioritization.csv`：抗原片段优先级。
- `02_epitope_windows.csv`：候选表位窗口。
- `03_antibody_developability.csv`：抗体可开发性检查。
- `04_candidate_library.csv`：实验阳性抗体、扰动 decoy 和 negative control 候选库。
- `05_candidate_ranking.csv`：多目标候选排序。
- `06_sandwich_pair_report.md`：两株抗体的 pair compatibility 结论。

外部结构工具输入：

- `outputs/exports/fasta/antigen_fragments.fasta`
- `outputs/exports/fasta/validated_fv_chains.fasta`
- `outputs/exports/fasta/complex_*.fasta`
- `outputs/exports/fasta/sandwich_*.fasta`
- `outputs/exports/af3_json/af3_complex_*.json`
- `outputs/exports/af3_json/af3_sandwich_*.json`
- `outputs/exports/external_tool_manifest.json`

## External Pipeline Integration

外部工具配置文件：

```text
config/external_pipelines.example.json
```

默认 adapter 都是 disabled。确认本机命令、模型数据库、GPU 设置和输出路径后，将对应工具的 `enabled` 改为 `true`，并编辑 `command_template`。

支持的 input selector：

- `validated_fv_chains_fasta`：VH/VL Fv 序列，用于抗体结构建模。
- `antigen_fragments_fasta`：优先抗原片段。
- `complex_fastas`：抗体-抗原复合物 FASTA。
- `sandwich_fasta`：Fab1:NfL:Fab2 sandwich FASTA。
- `af3_json_files`：AF3-style JSON 模板。

每次运行会写出：

- `outputs/exports/external_jobs/pipeline_jobs.tsv`
- `outputs/exports/external_jobs/run_external_pipelines.sh`

`run_external_pipelines.sh` 是命令表，不是必须立即执行的脚本。disabled jobs 会被注释。正式运行前应检查工具版本、输入路径、数据库路径、GPU 设置和输出目录。

## Expected Replay Result

使用当前 validation set 时，两株实验阳性抗体在 replay screen 中排在前两名：

1. `7-H11-D3-2-C7`
2. `15-C12-H6`

pair 模块将二者分配到 NfL rod/coil-2B 区域内的非重叠窗口，并建议进入外部三元结构建模和实验 sandwich assay 确认。

## Scientific Boundary

本仓库提供的是课题逻辑的计算化复现和候选优先级排序。除非后续替换为真实结构和实验数据，结果应表述为 computational prioritization / workflow reproduction，而不是已证明的物理结合事实。

后续需要补充或替换的证据包括：

- antibody/Fv/Fab 结构模型及模型质量指标；
- antibody-antigen co-folding 或 docking 的 interface PAE、ipTM/pTM、DockQ/pDockQ、buried surface area、clash metrics；
- 可开发性实验；
- 亲和力和特异性测定；
- sandwich capture/detection 实验验证；
- 对 22 kDa 和 11 kDa 条带的 LC-MS/MS、neo-termini、Cys322-containing peptide 和 Cys322-Cys322 disulfide evidence。

## Citation and License

引用信息见：

```text
CITATION.cff
```

代码许可证见：

```text
LICENSE
```

当前许可证为 MIT。正式上传前，如需改成更严格的许可或改成个人/机构署名，请同步更新 `LICENSE`、`CITATION.cff` 和 `pyproject.toml`。

## Reproducibility

工作流使用确定性规则，并为 CDR perturbation controls 固定随机种子。输入不变时，重新运行应得到相同的排序表和报告。
