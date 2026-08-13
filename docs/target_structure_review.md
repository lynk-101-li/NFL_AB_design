# NEFL 280–377 单链单体候选结构审核

## 当前状态

**`blocked_pending_human_review`：候选结构已生成和机器校验，但尚未经人工审核，不能进入真实模型 handoff。**

当前 campaign 只采用单链单体，并将第一表位定义为 `helix_surface_323_331`。热点仅为 Met325/Leu329；Cys322 不属于表位窗口、热点、名称或直接接触约束。多链状态不建模、不生成、不排名，详见 `docs/antigen_conformation_strategy.md`。

该候选使用 AlphaFold Protein Structure Database 中人 NEFL/P07196 的 v6 全长单体预测，裁剪 A 链 280–377 残基。AlphaFold DB 文件自身明确声明这是 **theoretical modelling only**；它不是实验解析结构，也不能证明该片段在溶液中的实际构象、表位可及性或聚集状态。

## 来源和可复现性

- 来源：[AlphaFold DB AF-P07196-F1 model v6](https://alphafold.ebi.ac.uk/files/AF-P07196-F1-model_v6.pdb)
- 全长来源 SHA-256：`37912aac5cefd85b177e754a7c55c10c0f50166baf7f30b012151492eae300b1`
- 裁剪 PDB SHA-256：`80686267d6a93eda5829bab187406c14caf0365ba7bee21a4eb0d2ed082fc73b`
- 坐标范围：UniProt P07196 1-based inclusive `280–377`
- PDB 链：`A`
- 局部长度：98 aa；PDB 中保留全长残基编号 280–377
- 序列：`FKSRFTVLTESAAKNTDAVRAAKDEVSESRRLLKAKTLEIEACRGMNEALEKQLQELEDKQNADISAMQDTINKLENELRTTKSEMARYLKEYQDLLN`

脚本会在写出任何结果前检查：来源 SHA-256、P07196 DBREF、单个 MODEL、SEQRES 543 aa、A 链 ATOM 残基 1–543 连续完整、SEQRES/ATOM 序列一致、每个残基的标准氨基酸重原子拓扑、目标序列以及两个表位序列。任意一项不匹配均 fail closed。

重新生成：

```bash
python3 scripts/prepare_target_structure.py --overwrite
```

## pLDDT 概要

PDB B-factor 列在 AlphaFold DB 文件中存储 pLDDT，不是实验 B-factor。

| 区域 | 残基数 | 平均 | 中位数 | 最小 | 最大 | ≥90 的残基数 |
|---|---:|---:|---:|---:|---:|---:|
| NEFL 280–377 | 98 | 95.15 | 95.44 | 87.88 | 98.12 | 94 |
| Helix surface 323–331 | 9 | 92.40 | 93.12 | 87.88 | 94.50 | 8 |
| C-boundary 368–377 | 10 | 96.79 | 96.78 | 95.81 | 97.81 | 10 |

这些数值支持“模型对局部预测几何较有信心”，不支持“表位已暴露”、“抗体可结合”或“单体是生理构象”。

## 建议热点与几何证据

下列残基只是 **proposed pending human review**。不应将整个 323–331 或 368–377 窗口盲目作为热点。

| 表位 | 建议残基 | pLDDT | 侧链探针原子 | CA→表位 CA 质心 (Å) | 探针→最近非局部重原子 (Å) | 10 Å 非局部 CA 邻居数 | 6 Å 探针接触残基数 |
|---|---:|---:|---|---:|---:|---:|---:|
| Helix surface | Met325 | 93.88 | CE | 3.981 | 5.640 | 8 | 1 |
| Helix surface | Leu329 | 91.00 | CD1 | 4.018 | 3.917 | 8 | 2 |
| C-boundary | Tyr368 | 96.75 | OH | 7.016 | 5.424 | 8 | 1 |
| C-boundary | Tyr372 | 97.31 | OH | 2.293 | 5.714 | 8 | 1 |
| C-boundary | Leu375 | 97.12 | CD2 | 4.366 | 5.020 | 8 | 2 |

原子完整性检查通过。建议热点的 CA 成对距离为：

- Met325/Leu329：6.191 Å。
- Tyr368/Tyr372：6.024 Å；Tyr368/Leu375：10.300 Å；Tyr372/Leu375：5.021 Å。

### SASA 与螺旋表面独立复核

已用 BioPython 1.87 `Bio.PDB.ShrakeRupley` 在**全长 543 aa 模型上下文**中重算 SASA：溶剂探针半径 1.4 Å，每原子 960 个球面点，相对 SASA 用 Wilke 残基最大可及面积归一化。完整可机读证据见 `input/structures/NEFL_P07196_AFDB_v6_280-377_sasa_review.json`。

| 表位 | 热点 | 全长上下文绝对 SASA (Å²) | Wilke 相对 SASA | 侧链 SASA (Å²) | pLDDT |
|---|---:|---:|---:|---:|---:|
| Helix surface | Met325 | 116.67 | 0.521 | 113.39 | 93.88 |
| Helix surface | Leu329 | 94.31 | 0.469 | 87.30 | 91.00 |
| C-boundary | Tyr368 | 132.71 | 0.505 | 124.99 | 96.75 |
| C-boundary | Tyr372 | 152.34 | 0.579 | 147.90 | 97.31 |
| C-boundary | Leu375 | 113.25 | 0.563 | 108.69 | 97.12 |

复核支持当前两组建议热点 `325/329` 与 `368/372/375`。`323–331` 的 φ/ψ 分别落在 `-70.2…-59.4°` 和 `-47.8…-35.0°`；`368–377` 分别落在 `-65.5…-61.3°` 和 `-47.7…-37.3°`，两个窗口都是连续 α 螺旋。Met325/Leu329 是间距 6.191 Å 的 i/i+4 同面组合，不需要引入半胱氨酸接触约束。

`368–377` 触及裁剪结构的人工 C 端，只在 crop 中计算会高估末端暴露。例如 Leu375 的绝对/相对 SASA 从全长上下文的 `113.25/0.563` 升至 crop 的 `142.84/0.711`；Leu376 和 Asn377 的虚高更明显。现有方案避开 376/377，而 375 在全长上下文中仍有较高可及性，因此结论不变。真实设计/对接宜为 377 后补充原生下游结构缓冲，或防止模型利用人工 C 端截面。

以上仍只是基于 AlphaFold DB 理论单体坐标的确定性计算建议，**不是实验表位、可及性、聚集状态或抗体结合验证**。单体模型也不能表征生理多聚体/神经丝组装中的界面遮挡。

## 人工审核清单

审核人应亲自完成下列项目，不应由脚本自动代填：

- [ ] 确认来源为 AlphaFold DB P07196 v6，并理解其为理论预测。
- [ ] 确认裁剪链 A、序列、全长编号 280–377 和 local 1–98 映射。
- [ ] 可视化检查 323–331 与 368–377，确认侧链方向、断链、clash 和截短边界影响。
- [ ] 审阅 SASA 复核 JSON、人工 C 端边界差异和侧链可视化；不把静态单体 SASA 当作实验表位验证。
- [ ] 审核 `325/329` 与 `368/372/375`，可接受、删除或替换，但不选择整个表位窗口。
- [ ] 确认本次 canary 只回答单链单体输入可运行性，不回答任何多链构象问题。
- [ ] 确认两份 RFantibody HLT PDB 和两份 Germinal scFv PDB 是彼此独立的模板几何。
- [ ] 在正式 `config/target_structure_manifest.json` 中记录 reviewer、带时区的 ISO-8601 时间、最终热点和 `contracts_acknowledged: true`。

### 审核记录（待填）

- Reviewer：
- 日期时间（含时区）：
- 审核结论：
- `helix_surface_323_331` 最终热点：
- `C_boundary_368_377` 最终热点：
- 独立 SASA/可及性工具和版本：
- 备注：

## 产物索引

- 全长锁定来源：`input/structures/AF-P07196-F1-model_v6.pdb`
- 裁剪 PDB：`input/structures/NEFL_P07196_AFDB_v6_280-377_chainA.pdb`
- FASTA：`input/structures/NEFL_P07196_AFDB_v6_280-377.fasta`
- 机器证据：`input/structures/NEFL_P07196_AFDB_v6_280-377_evidence.json`
- SASA 与螺旋表面复核：`input/structures/NEFL_P07196_AFDB_v6_280-377_sasa_review.json`
- 阻断状态候选 manifest：`input/structures/target_structure_manifest.candidate.blocked.json`
- 生成脚本：`scripts/prepare_target_structure.py`

`input/structures/target_structure_manifest.candidate.blocked.json` 故意保留空的 `selected_hotspots_by_epitope`、空的 reviewer 字段以及未补齐的抗体模板路径。不应将其改名后直接运行；必须先完成上述审核，再创建正式 manifest。

## 组装待审候选（仍不可执行）

当结构候选、SASA 复核和两份中性模板坐标均已生成后，可运行：

```bash
python3 scripts/assemble_review_candidate_manifest.py
```

该命令会严格复核输入 schema、blocked 状态、文件存在性、SHA-256 绑定、两份独立模板以及 SASA 建议热点，然后生成 `input/structures/target_structure_manifest.review_candidate.blocked.json`。该文件虽已填入计算建议热点和坐标路径，仍保留 `execution_state: blocked_pending_human_review`、空 reviewer 和 `contracts_acknowledged: false`，因此会被 `prepare_real_model_jobs.py` 拒绝。

默认拒绝覆盖；只有确认要重新组装 blocked 候选时才使用 `--overwrite`。脚本无论如何都不会创建或覆盖 `config/target_structure_manifest.json`。完成人工审核后，应根据输出中的 `promotion_instructions` **单独创建**正式 manifest。
