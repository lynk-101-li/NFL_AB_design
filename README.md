# NFL_AB_design

`NFL_AB_design` 是一个自包含、可审计的 NfL/NEFL 抗体计算设计示范包。当前主流程是：

```text
抗原截断推断
  -> 两个 NfL 目标表位
  -> 两个 framework-only 模板
  -> H1/H2/H3/L1/L2/L3 六 CDR 从头设计模拟
  -> 结构/界面代理筛选
  -> 可开发性代理筛选
  -> 前瞻候选排名
  -> 独立的回顾性阳性对照演示
  -> sandwich pair 优先级
  -> 外部真实模型交接
```

> **当前证据状态：模拟，不是真实模型结果。** 仓库使用确定性 proxy 生成和 proxy 打分来演示设计漏斗。仓库已把 RFantibody、IgGM、Germinal、tFold、IgFold、ImmuneBuilder/ABodyBuilder2、AlphaFold 3、Chai-1 和 Boltz-2 的官方源码固定为 submodule，但当前运行未执行这些模型，也没有产生其预测结果。两株已知实验阳性抗体如果在回顾表中位居 Top 2，只能解释为 **retrospective positive-control demonstration**，不是盲法发现。

完整的模拟实验体系、各阶段输入输出、对照和进入真实模型前的门禁见 `docs/simulated_experimental_system.md`。

## Project Background

NfL（NEFL）是神经丝轻链蛋白，常用于神经损伤相关检测研究。本课题的上游观察是：NfL 相关片段在非还原条件下约为 25–35 kDa，DTT 还原后约为 6–12 kDa。约四倍的表观质量变化优先支持四聚体解释，但仍不能排除“二硫键二聚体＋未知结合伙伴”。

人 NfL 通常每条链只有一个关键半胱氨酸 Cys322。当前优先生化模型写作 `4M → (M–S–S–M)₂`：两个 Cys322–Cys322 二硫键分别形成两个二聚体，两个二聚体再以反平行、错位方式通过卷曲螺旋界面结合成四聚体。这里的非共价界面不能简化成“疏水作用”，而应包括卷曲螺旋疏水核心、静电作用以及几何/构象互补。DTT 响应说明二硫键对复合物稳定性重要，但不能单独证明四聚体化学计量，也不能证明二聚体间界面纯粹由疏水作用形成。

还原后 6–12 kDa 的范围、Cys322 位置、cathepsin-like 边界支持和结构区域共同将抗原空间收缩到 rod/coil-2B 附近。工作抗原集中在 NfL aa 280–377，主要片段包括 `280-375`、`281-376` 和 `282-377`。非还原 SDS-PAGE 的表观迁移不再作为精确二聚体理论质量约束。详细推理见 `docs/nfl_truncation_inference_rationale.md`。

流程保留 `280-375` 作为生化截断排序第一名，但结构/生成 handoff 使用 `280-377` 建模上下文，以完整覆盖配置中的 368–377 表位；二者在 manifest 中分别记录，不能混称同一个实验片段。

上游条带解释不再作为抗体接触约束。当前 campaign 只使用单链单体结构，并把第一设计表位收窄到 323–331；热点仅为 Met325/Leu329。Cys322 虽自然存在于抗原坐标中，但不在表位窗口、热点、名称或直接接触要求中。范围说明与机器可读状态分别见 `docs/antigen_conformation_strategy.md` 和 `config/antigen_conformation_tracks.json`。

## Prospective Design Campaign

主配置位于 `config/design_campaign.json`。默认 campaign 使用两个目标表位：

| 表位 ID | NfL 坐标（1-based, inclusive） | 序列 | 设计角色 |
|---|---:|---|---|
| `helix_surface_323_331` | 323–331 | `RGMNEALEK` | Met325/Leu329 helical face |
| `C_boundary_368_377` | 368–377 | `YLKEYQDLLN` | C-terminal boundary |

最终 12 个前瞻候选采用分层短名单：在四个 `template × epitope` 组合都有合格候选时，每层保留 3 个；各层内部仍按模拟综合分排序，并同时保留全局 rank。这样不会让某个表位因代理分数尺度偏差而完全挤出实验候选池。

### 两模板与六 CDR

`7-H11-D3-2-C7` 和 `15-C12-H6` 只为生成阶段提供两个不同的 VH/VL **framework source**。程序在生成请求中将 H1、H2、H3、L1、L2、L3 全部遮罩，并将这六个区域都定义为待设计区域。

因此，前瞻生成阶段：

- 可以使用两株抗体的 framework residues；
- 不读取已知 CDR 氨基酸身份作为生成特征；
- 不使用已知完整 VH/VL 序列作为生成特征；
- 不允许已知阳性序列进入 prospective candidate table。

CDR 范围在 campaign 配置中以链内 1-based inclusive 原始序列坐标显式记录，并已锁定为 `ANARCI 2020.04.23 Chothia` 编号结果。两个模板的区间不强制相同；例如 `52A`、`31A` 和 `100A–100C` 等插入号保留在编号证据中，而模拟生成与 IgGM 遮罩使用对应的连续 raw 坐标。完整工具版本、输入哈希、原始 Chothia labels 和 gap 说明见 `input/antibody_templates/chothia_numbering_evidence.json`。

### 模拟规模与筛选阈值

默认随机种子是 `20260812`。每个“模板 × 表位”组合模拟 1280 个设计，即 `2 × 2 × 1280 = 5120` 个 prospective candidates。这是用于检验漏斗、字段和排序稳定性的确定性 proxy library，不代表 RFantibody、IgGM 或 Germinal 必须各自完成 5120 次昂贵结构推理。真实模型采用“低成本大库 → 分层结构复核 → 少量共同终筛”的预算。

默认漏斗参数为：

| 参数 | 默认值 | 用途 |
|---|---:|---|
| `structure_min` | 58 | 结构 proxy 最低分 |
| `interface_min` | 55 | 界面 proxy 最低分 |
| `developability_min` | 60 | 可开发性 proxy 最低分 |
| `composite_min` | 60 | 综合 proxy 最低分 |
| `selection_count` | 12 | 最多导出的前瞻候选数 |

所有模拟分数都应保留 `data_status=simulated`、`*_is_simulated` 和 machine-readable provenance。这些数值不是 pLDDT、PAE、ipTM、DockQ、结合自由能或真实亲和力。

## Prospective 与 Retrospective 的隔离

已知实验阳性序列存放在 `validation/experimentally_validated_antibodies.fasta`，但它们的完整 VH/VL 只在 prospective generation、screening 和 ranking 完成后，才以 `retrospective_positive_control` 状态注入独立的回顾性演示。

| 轨道 | 允许的输入 | 用途 | 能否声称盲法发现 |
|---|---|---|---|
| Prospective | 遮罩六 CDR 的 framework-only 模板 + 目标表位 | 模拟从头设计和筛选 | 否：当前是 proxy simulation |
| Retrospective | 前瞻排名完成后注入的已知阳性全序列 | 阳性对照、管线演示和回顾性排名 | 否：Top 2 不是 de novo discovery |

两株已知阳性在 retrospective table 中排名 Top 2，只说明当前代理评分能将阳性对照排到前面。这不能证明管线在不知道答案时能找到它们，也不能证明其表位、结构或 sandwich 几何已被真实模型确认。

## Pipeline Stages and Outputs

`outputs/` 中的主表按下列阶段组织：

| 阶段 | 输出 | 内容 |
|---:|---|---|
| 00 | `00_antigen_truncation_*` | 全肽键cathepsin-like 分数、候选边界和截断报告 |
| 01 | `01_antigen_fragment_prioritization.csv` | 包含 Cys322 的抗原片段优先级 |
| 02 | `02_epitope_windows.csv` | NfL 表位窗口与代理分数 |
| 03 | `03_template_frameworks.csv` | 两个 framework-only、六 CDR 遮罩的模板 |
| 04 | `04_backbone_generation.csv` | 确定性模拟的骨架生成记录 |
| 05 | `05_sequence_candidates.csv` | 六 CDR 候选序列及生成来源 |
| 06 | `06_structure_interface_screen.csv` | 结构和界面 proxy screen，显式标注 simulated |
| 07 | `07_developability_screen.csv` | 可开发性 proxy screen，显式标注 simulated |
| 08 | `08_screening_funnel.csv` | 每个阈值前后的输入、通过和剔除计数 |
| 09 | `09_prospective_candidates.csv` | 不含已知阳性全序列的前瞻候选排名 |
| 10 | `10_retrospective_demo_candidates.csv` | 单独的回顾性阳性对照演示 |
| 11 | `11_sandwich_pair_ranking.csv` | 候选配对的非重叠表位/冲突 proxy 排名 |

`outputs/workflow_report.md` 是总报告。`outputs/intermediate/source_manifest.json` 和 `outputs/exports/external_tool_manifest.json` 记录输入、输出、运行状态与来源；需要日期时应以实际运行清单中的字段为准，不要从文档推断运行日期。

## Repository Layout

```text
NFL_AB_design/
  README.md
  pyproject.toml
  config/
    design_campaign.json
    external_pipelines.example.json
  input/
    antigen_truncation/truncation_constraints.json
    antibody_templates/
  resources/
    project_context/
    antigen_inference/
  validation/
    experimentally_validated_antibodies.fasta
  docs/
    high_school_quickstart_zh.md
    nfl_truncation_inference_rationale.md
  scripts/run_nfl_ab_design.py
  src/nfl_ab_design/
  tests/
  outputs/
```

`input/` 放起始约束和格式示例；`resources/` 放课题背景与 NfL 截断推断资料；`validation/` 中的完整已知序列是回顾性对照，不是前瞻候选库。

## Installation and Run

本包的本地 proxy workflow 需要 Python 3.10 或更高版本及 NumPy；editable install 会自动安装声明的依赖。

```bash
python3 scripts/run_nfl_ab_design.py
```

也可以用模块方式：

```bash
PYTHONPATH=src python3 -m nfl_ab_design
```

可选 editable install：

```bash
python3 -m pip install -e .
nfl-ab-design
```

相同输入、campaign 配置和种子应生成相同的 proxy 表。改变阈值会改变 funnel 通过数和导出候选，因此比较两次运行时必须同时保存配置与 manifest。

如果目标是在 GPU 服务器上完成 RFantibody、IgGM、Germinal 的完整部署，请从 [`deploy/autodl/README.md`](deploy/autodl/README.md) 开始。部署脚本固定上游 commit，但不会自动接受第三方许可、伪造人工审核或绕过 runtime attestation。

## External Model Handoff

`config/external_pipelines.example.json` 包含默认 `enabled: false` 的命令模板：

- RFantibody 六 CDR 生成 adapter；
- IgGM 六 CDR 生成 adapter；
- Germinal 表位条件 scFv 生成 adapter（独立单链轨道，不等同于 native paired Fv）；
- tFold、IgFold、ImmuneBuilder/ABodyBuilder2 用于抗体结构预测与交叉复核；
- AlphaFold 3、Chai-1、Boltz-2 用于抗体—抗原复合物交叉预测。

配置中的 `<PATH_TO_...>` 是故意保留的明确占位符。不得在未替换占位符、未核对安装版本 schema、未提供 NfL antigen PDB 和模型 checkpoints 时启用任务。

主要 selector 包括：

- `design_request_index`：RFantibody/IgGM/Germinal 规范化生成请求索引；
- `design_request_files`：按 engine 记录的生成请求 JSON；
- `masked_template_fasta`：六 CDR 遮罩模板；
- `candidate_fv_chains_fasta`：筛选后 prospective VH/VL 候选；
- `complex_fastas` 与 `af3_json_files`：候选-抗原复合物交接输入。

运行 proxy workflow 只会生成请求、FASTA/JSON 模板、任务表和默认注释掉的命令表。这些是 **handoff artifacts**，不是已执行模型的证据。真实执行后，应将工具名、版本、checkpoint/database、硬件、命令、原始输出和运行状态追加到 manifest。

三个真实后端已有严格的、只生成计划而不执行模型的 Python adapter：

- `src/nfl_ab_design/adapters/rfantibody.py`：要求 target PDB、两份 Chothia HLT PDB、显式坐标映射和人工选择的 PDB hotspots，生成 `rfdiffusion → proteinmpnn → rf2` 四组 job；
- `src/nfl_ab_design/adapters/iggm.py`：要求 target PDB 链序列和 full→local 1-based 映射，输出 IgGM 所需的 H/L/A 三链 FASTA 与四组 job；
- `src/nfl_ab_design/adapters/germinal.py`：要求 target PDB 和两份 VH-linker-VL scFv PDB，输出 Germinal target YAML 与四组 Hydra job，并明确记录其单链几何限制。

`config/model_components.json` 固定九个模型的 submodule、commit、许可证边界和环境；`config/real_model_backends.json` 保留三种生成器的执行合同。学生服务器从零部署、权重放置和逐模型运行命令见 `deploy/autodl/README.md` 和 `docs/real_model_installation.md`。

克隆时应初始化 submodule：

```bash
git clone --recurse-submodules https://github.com/lynk-101-li/NFL_AB_design.git
cd NFL_AB_design
python3 scripts/verify_model_components.py
```

补齐并人工复核 target manifest 后，将 `execution_state` 设为 `reviewed_ready_for_handoff`，填写带时区的 reviewer/timestamp 并确认 review contract，再用统一编译器生成三个模型的作业交接：

```bash
python3 scripts/prepare_real_model_jobs.py \
  --target-manifest config/target_structure_manifest.json \
  --profile smoke \
  --job-scope canary \
  --iggm-repo-dir <PINNED_IGGM_CHECKOUT> \
  --germinal-repo-dir <PINNED_GERMINAL_CHECKOUT> \
  --germinal-af-params-dir <AF_MULTIMER_PARAMS>
```

该命令会验证全部 `2 templates × 2 epitopes` 计划，但 smoke 默认只授权三个引擎中同一个 template×epitope canary；它仍不执行模型。真实 handoff 和结果必须写入 `real_runs/`，不能写入可再生的 `outputs/`。proxy workflow 只清理明确列出的可再生 artifact，不会删除 `real_runs/` 或误放在 `outputs/external_results/` 的文件。

### 真实运行前的两个阻断条件

1. **抗原构象路线**：当前 campaign 唯一采用三模型都能审计的 NfL aa280–377 单链单体，以 `helix_surface_323_331` 和 `C_boundary_368_377` 为设计目标。任何二聚体或其他多链状态均不属于本 campaign。
2. **学生服务器资源与许可**：每台服务器必须独立运行只读 preflight。建议系统盘 80GB、安装前至少 35GB 可用，并另备至少 50GB 文件存储。PyRosetta、AlphaFold 参数及其他受限资产必须由学生按各自许可合法获取；仓库不提供这些文件。

## Testing

```bash
python3 -m unittest discover -s tests
```

测试覆盖的核心 contract 包括：两个独立 framework source、六 CDR 均为设计区、固定种子可重现、prospective tables 不含已知阳性、回顾性对照单独标注、模拟指标具有可机读 provenance，以及 funnel 计数内部一致。

## Scientific Boundary

该仓库适合用于设计策略演示、软件 contract 测试、漏斗审计和外部工具输入准备。它不能代替：

- 抗体/Fv/Fab 的真实结构建模与模型质量指标；
- 抗体-抗原 co-folding/docking 的 PAE、ipTM/pTM、DockQ/pDockQ、埋藏表面积和 clash metrics；
- 亲和力、特异性、交叉反应和可开发性实验；
- sandwich capture/detection 实验；
- 对非还原 25–35 kDa 与 DTT 后 6–12 kDa 条带的 LC-MS/MS、neo-termini、Cys322-containing peptide 和 Cys322–Cys322 disulfide evidence；
- C322S、还原/非还原 SDS-PAGE、SEC-MALS 或原生质谱对“四聚体”与“二聚体＋未知伙伴”模型的判别。

在完成这些验证之前，结果应表述为 **simulated computational prioritization**，而不是已证明的物理结合或盲法抗体发现。

## Citation and License

引用信息见 `CITATION.cff`，项目代码许可证见 `LICENSE`。AlphaFold DB 坐标以及外部模型、权重和参数仍受各自许可约束，详见 `THIRD_PARTY_NOTICES.md`。
