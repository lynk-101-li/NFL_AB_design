# NfL 抗体六 CDR 设计模拟：高中生快速上手

这份指南帮你理解课题问题、运行双模板六 CDR 设计模拟、阅读筛选漏斗，并分清“模拟分数”、“回顾性对照”和“真实模型/实验结果”。

## 1. 先记住最重要的边界

当前代码是一个确定性的 **proxy simulation**：它用可重复的规则模拟生成、结构/界面筛选和可开发性筛选。

当前仓库运行：

- 没有实际运行 RFantibody、IgGM 或 Germinal；
- 没有实际运行 tFold、AlphaFold 3、Chai-1、Boltz 或 Rosetta；
- 没有产生真实的 pLDDT、PAE、ipTM、DockQ 或能量；
- 没有用实验证明新候选一定结合 NfL。

因此，输出中的 simulated score 只用来学习流程、比较候选和检查代码，不能当成真实结构预测结果。

## 2. 课题要解决什么

NfL，也叫 NEFL，是一种神经丝蛋白。课题希望获得能识别 NfL 的抗体，并最终组成 sandwich assay。

Sandwich assay 可以想成两只手同时抓住一个目标：

- capture antibody 先抓住 NfL；
- detection antibody 再抓住 NfL 的另一处，并提供检测信号。

两株抗体识别的位置需尽量不重叠，三维结构中也不应严重冲突。

## 3. 为什么先分析 NfL 片段

上游实验观察到 NfL 相关片段在非还原条件下约为 25–35 kDa，DTT 后约为 6–12 kDa。当前优先假设是四条链组成复合物：Cys322 分别连接成两个二硫键二聚体，两个二聚体再通过卷曲螺旋的疏水核心、静电作用和形状互补形成反平行四聚体。二聚体＋未知结合伙伴仍不能排除。

NfL 全长序列中只有一个 canonical cysteine：`Cys322`。所以，如果二硫键二聚体假设正确，相关片段就必须包含 Cys322。

代码根据下列线索给 NfL 截断片段排序：

- 包含 Cys322；
- 单体质量接近 11 kDa；
- 单链理论质量与 DTT 后约 6–12 kDa 的宽范围相容；
- 非还原/还原表观质量变化支持高阶复合物，但不被当作精确理论质量约束；
- 两端有 cathepsin-like 切割支持；
- 位于 NfL rod/coil-2B 附近。

当前工作片段集中在 aa 280–377，其中 `280-375` 是主要候选。这仍然是需要质谱等证据检验的生物学假设。

## 4. 这次“从头设计”是怎么做的

### 两个目标表位

表位是抗体在抗原上识别的局部区域。默认 campaign 针对：

| 表位 | NfL 坐标 | 序列 |
|---|---:|---|
| Helix surface | 323–331 | `RGMNEALEK` |
| C-terminal boundary | 368–377 | `YLKEYQDLLN` |

第一表位的设计热点仅为 `Met325/Leu329`。Cys322 不作为抗体接触热点；上游条带假设与本轮抗体设计约束分开处理。

最后的 12 个候选不是简单取一个全局分数榜的前 12，而是先保证两个模板和两个表位的四种组合都有代表，再在每组内部按模拟综合分选优。

### 两个 framework-only 模板

代码从两株已知抗体获取两个不同的 VH/VL framework source：

- `7-H11-D3-2-C7`；
- `15-C12-H6`。

这不等于把已知抗体整体拿来做候选。在生成请求里，六个 CDR 都被遮罩：

```text
H1  H2  H3  L1  L2  L3
```

生成阶段只使用 framework residues：不使用已知 CDR 氨基酸身份，也不把已知完整 VH/VL 当成生成特征。六个 CDR 都是待设计区域。

这些 CDR 不再使用启发式大致范围，而是使用 `ANARCI 2020.04.23 Chothia` 编号后映射回每条输入链的 1-based 原始坐标。两个模板的 CDR 长度可以不同；程序只遮罩和替换配置中的精确位点，其余 framework 位点必须与源序列一致。可审计的编号 labels 和输入哈希在 `input/antibody_templates/chothia_numbering_evidence.json`。

### 生成规模

默认配置是：

```text
2 个模板 × 2 个表位 × 每组合 1280 个模拟设计 = 5120 个前瞻候选
seed = 20260812
```

`seed` 是随机种子，它让同样的配置能重现同样的模拟结果。5120 是本地 proxy library 的候选数，不是 RFantibody/IgGM/Germinal 各自必须真实运行的数量；真实结构模型会分层缩小计算量。

## 5. 筛选漏斗怎么读

漏斗就是逐层删掉不符合条件的候选。默认参数在 `config/design_campaign.json`：

| 关卡 | 默认阈值 | 简单理解 |
|---|---:|---|
| Structure | 58 | 模拟结构是否像一个可用候选 |
| Interface | 55 | 模拟界面是否适合目标表位 |
| Developability | 60 | 序列风险是否可接受 |
| Composite | 60 | 综合 proxy 是否达标 |
| Final selection | 12 | 最多导出多少个前瞻候选 |

表中的 58、55、60 都是此代码的 proxy 阈值，不是生物学通用标准。不要把它们和真实的结构模型质量或结合常数等同。

## 6. 为什么有两张排名表

这是本课题最容易误解的部分。

### Prospective table

`09_prospective_candidates.csv` 只包含从遮罩六 CDR 模板生成的模拟候选。已知阳性的完整序列不允许出现在这条轨道中。

### Retrospective table

`10_retrospective_demo_candidates.csv` 在前瞻生成和排名结束后，才加入两株已知阳性全序列，并把它们标记为：

```text
retrospective_positive_control
```

已知的 Top 2 是：

1. `7-H11-D3-2-C7`
2. `15-C12-H6`

这只是回顾性阳性对照演示，不是盲法发现。换句话说，代码在这一步已经“知道它们是阳性对照”，所以不能说代码从 5120 个新候选里盲法找回了它们。

## 7. 如何运行

需要 Python 3.10 或更高版本。在终端进入仓库根目录：

```bash
cd NFL_AB_design
python3 --version
python3 scripts/run_nfl_ab_design.py
```

也可以运行：

```bash
PYTHONPATH=src python3 -m nfl_ab_design
```

运行后先打开：

```text
outputs/workflow_report.md
```

每次比较结果时，还要保存 campaign 配置和 manifest。如果需要查运行日期，读实际运行 manifest 的时间字段，不要根据文档或 seed 猜测日期。

## 8. 主要输出是什么

| 文件 | 要回答的问题 |
|---|---|
| `00_antigen_truncation_report.md` | 为什么重点考虑 DTT 后 6–12 kDa 范围内且包含 Cys322 的片段？ |
| `01_antigen_fragment_prioritization.csv` | 哪些 NfL 片段更符合工作假设？ |
| `02_epitope_windows.csv` | 哪些局部表位更值得设计？ |
| `03_template_frameworks.csv` | 六 CDR 是否都已遮罩，两模板来自何处？ |
| `04_backbone_generation.csv` | 每个模板-表位组合生成了哪些模拟记录？ |
| `05_sequence_candidates.csv` | 六个 CDR 各自的模拟序列是什么？ |
| `06_structure_interface_screen.csv` | 哪些候选通过结构/界面 proxy？ |
| `07_developability_screen.csv` | 哪些候选通过可开发性 proxy？ |
| `08_screening_funnel.csv` | 每一关进入、通过、被剔除多少个？ |
| `09_prospective_candidates.csv` | 哪些新模拟候选排在前面？ |
| `10_retrospective_demo_candidates.csv` | 已知阳性对照的回顾性演示怎样？ |
| `11_sandwich_pair_ranking.csv` | 哪些候选对的表位更不重叠？ |

阅读 06、07 时，要查看 `data_status`、`*_is_simulated` 和 `metric_provenance`，它们用来防止把模拟分数误写成真实模型值。

## 9. 外部真实模型如何接入

代码会准备规范化请求、FASTA/JSON 和任务表，但不会默认运行外部工具。配置位于：

```text
config/external_pipelines.example.json
```

其中包含默认关闭的：

- RFantibody 六 CDR 生成 adapter；
- IgGM 六 CDR 生成 adapter；
- Germinal 表位条件 scFv 生成 adapter；它把 VH 和 VL 用 linker 连成一条链，不能冒充原生配对 Fv 几何；
- 可选 tFold 结构预测/复核 adapter；
- 其他可选建模、复合物预测和界面分析 adapter。

命令中的 `<PATH_TO_...>` 是占位符，必须按实际安装修改。还需要 NfL antigen PDB、模型 checkpoint/database、对应版本的 adapter 和可用硬件。

看到以下文件只能说明“已准备交接输入”：

```text
outputs/exports/design_requests/
outputs/exports/fasta/
outputs/exports/af3_json/
outputs/exports/external_jobs/pipeline_jobs.tsv
outputs/exports/external_jobs/run_external_pipelines.sh
```

它们不是外部模型已运行的证明。只有真实工具输出、日志、版本/checkpoint 记录和完整 manifest 同时存在，才可以声称某个模型已经运行。

三个真实 adapter 都要求抗原的坐标 PDB，并要求把本文中的 NfL 全长坐标明确映射到 PDB 残基。RFantibody 还要求 HLT 模板 PDB，Germinal 要求 scFv 模板 PDB；如果这些输入缺失，程序会直接停止，不会猜一个看似合理的坐标。

补齐并人工复核 `config/target_structure_manifest.json` 后，只有在 reviewer、带时区时间戳和 review contract 都记录完成、状态改为 `reviewed_ready_for_handoff` 时，才可用 `scripts/prepare_real_model_jobs.py` 统一编译三种模型的作业。`--profile smoke --job-scope canary` 会校验全部四个模板-表位组合，但每个引擎只授权同一个 canary；这仍然不会运行模型。真实交接和结果放在 `real_runs/`，不要放入可再生的 `outputs/`。

## 10. 运行测试

```bash
python3 -m unittest discover -s tests
```

测试会检查：

- 确实有两个不同的 framework source；
- H1/H2/H3/L1/L2/L3 六个 CDR 都是设计区；
- 相同 seed 的模拟可重现；
- prospective tables 不含已知阳性全序列；
- retrospective controls 被单独标记；
- simulated metrics 有可机读的来源标记；
- funnel 的计数前后一致。

## 11. 一节课可以怎么讲

1. NfL 和 sandwich assay 是什么。
2. 非还原 25–35 kDa、DTT 后 6–12 kDa 和 Cys322 如何形成“四聚体优先、二聚体＋伙伴未排除”的抗原假设。
3. framework 和 CDR 分别做什么。
4. 为什么生成时必须遮罩已知的六 CDR。
5. 运行代码，看 5120 个模拟候选如何通过漏斗。
6. 对比 prospective 和 retrospective 两张表。
7. 解释为什么已知 Top 2 不能写成盲法发现。
8. 讨论下一步需要的真实结构模型和实验。

## 12. 最快操作清单

```bash
cd NFL_AB_design
python3 --version
python3 scripts/run_nfl_ab_design.py
python3 -m unittest discover -s tests
```

然后依次看：

```text
outputs/workflow_report.md
outputs/08_screening_funnel.csv
outputs/09_prospective_candidates.csv
outputs/10_retrospective_demo_candidates.csv
outputs/11_sandwich_pair_ranking.csv
```

最后问自己三个问题：这个数值是 simulated 还是 measured？这个候选是 prospective 还是 retrospective control？这个结论是否还需要真实模型和实验验证？
