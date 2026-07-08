
# NfL 抗体设计课题快速上手指南

这份指南面向第一次接触这个课题的高中生，目标是让你理解课题在做什么、代码包怎么运行、结果怎么看，以及后续可以搭配哪些软件继续做结构分析。

## 1. 这个课题在研究什么

本课题围绕一种神经损伤标志物蛋白 NfL，也叫 NEFL。NfL 可以出现在神经损伤相关样本中，因此常被用于诊断检测研究。

我们的目标是设计并筛选能识别 NfL 的抗体，尤其是能组成 sandwich assay 的一对抗体。Sandwich assay 可以理解为“两只手同时抓住同一个目标蛋白”：一只抗体负责捕获 NfL，另一只抗体负责检测信号。如果两只抗体识别的位置互不重叠，就更有机会形成好用的检测抗体对。

这个课题的最终实验阳性抗体是：

- `7-H11-D3-2-C7`
- `15-C12-H6`

代码包的作用是把课题中的设计逻辑变成一套可重复运行的计算流程：先推断 NfL 的可能抗原片段，再寻找候选表位，然后对抗体序列做质量检查、排序，并准备外部结构预测软件需要的输入文件。

## 2. 课题设计主线

### 第一步：推断抗原片段

实验中观察到一个约 22 kDa 的 NfL 相关条带，并判断它可能是二硫键连接的二聚体。二聚体是两个相同片段连在一起，所以每个单体大约是 11 kDa。

NfL 全长蛋白里只有一个半胱氨酸 `Cys322`，它可以形成二硫键。因此，合理的 NfL 片段必须包含 `Cys322`。

代码会对 NfL 全长的每个肽键做 cathepsin-like 切割倾向打分，然后筛选出：

- 包含 `Cys322`；
- 单体质量接近 11 kDa；
- 二聚体质量接近 22 kDa；
- N 端和 C 端边界有蛋白酶切割支持；
- 位于 NfL rod/coil-2B 区域附近。

当前最重要的候选抗原片段是：

```text
NfL aa 280-375
```

### 第二步：寻找抗体可能识别的表位

表位就是抗体在抗原上识别的位置。代码会重点检查：

- `Cys322` 附近；
- 片段 N 端边界附近；
- 片段 C 端边界附近；
- `280-375` 区间内的滑动窗口。

这些窗口会按极性、带电性、结构区域、潜在修饰风险等指标排序。

### 第三步：检查抗体序列

抗体有重链 `VH` 和轻链 `VL`。代码会读取两株抗体的 VH/VL 序列，做基础检查：

- CDR 区域的大致位置；
- HCDR3 长度；
- 是否有潜在糖基化 motif，例如 `NXS/T`；
- 是否有氧化、脱酰胺等风险；
- CDR 区域的带电性和芳香族氨基酸比例。

### 第四步：候选库排序

代码会把两个真实阳性抗体和一些确定性扰动出来的阴性候选放到同一个候选库里，统一打分。

当前结果中，两个实验阳性抗体排在前两名：

```text
Rank 1: 7-H11-D3-2-C7
Rank 2: 15-C12-H6
```

这说明这套计算流程可以复现课题后期筛选得到的有效抗体对。

### 第五步：评估 sandwich pair

代码会判断两株抗体的预测表位是否重叠。如果两个表位距离足够远、重叠低，就更适合组成 sandwich assay。

当前计算结论是：

- `7-H11-D3-2-C7` 更偏向识别 C 端边界附近；
- `15-C12-H6` 更偏向识别 `Cys322` anchor 附近；
- 两者表位不重叠；
- 推荐 `7-H11-D3-2-C7` 做 capture，`15-C12-H6` 做 detection。

## 3. 代码包在哪里

代码包目录是：

```text
NFL_AB_design/
```

最重要的文件：

```text
NFL_AB_design/README.md
NFL_AB_design/scripts/run_nfl_ab_design.py
NFL_AB_design/src/nfl_ab_design/workflow.py
NFL_AB_design/input/antigen_truncation/truncation_constraints.json
NFL_AB_design/input/antibody_templates/template_fv_backgrounds.fasta
NFL_AB_design/resources/project_context/storyline.txt
NFL_AB_design/resources/project_context/research_plan.txt
NFL_AB_design/resources/antigen_inference/nfl_cathepsin_annotated_for_snapgene.gp
NFL_AB_design/validation/experimentally_validated_antibodies.fasta
NFL_AB_design/docs/nfl_truncation_inference_rationale.md
NFL_AB_design/config/external_pipelines.example.json
```

这里要特别区分：

- `input/` 放起始约束和模板，不放最终实验答案；
- `resources/` 放课题故事线、研究方案、NfL 截断分析过程和参考数据；
- `validation/` 放最终实验验证通过的两个抗体，用来检查计算流程是否能把它们排到最前面；
- `docs/` 放方法解释文档，适合讲课时使用。

运行后结果会写到：

```text
NFL_AB_design/outputs/
```

## 4. 推荐安装的软件

最小上手只需要：

- macOS / Windows / Linux 任意系统；
- Python 3.10 或更高版本；
- 一个终端软件；
- VS Code、Cursor 或其他代码编辑器；
- Excel、Numbers 或 LibreOffice，用来看 `.csv` 表格。

推荐搭配的软件：

- SnapGene：查看蛋白序列和注释；
- PyMOL 或 UCSF ChimeraX：查看蛋白结构；
- Excel / Numbers / LibreOffice：查看输出表格；
- IgFold 或 ABodyBuilder3：建模抗体 Fv/Fab；
- AlphaFold 3、Chai-1 或 Boltz：预测抗体-抗原复合物；
- Rosetta：做界面能量和结构精修分析。

高中生快速理解阶段，不需要一开始就安装所有结构预测软件。先能跑通 Python 流程、看懂输出表格即可。

## 5. 终端快速上手

先打开终端，进入 `NFL_AB_design` 仓库根目录。如果是从 GitHub 拉到本地，通常是：

```bash
git clone <你的仓库地址>
cd NFL_AB_design
```

如果已经在包含 `NFL_AB_design/` 的上一级目录，也可以直接：

```bash
cd NFL_AB_design
```

检查 Python 版本：

```bash
python3 --version
```

运行代码包：

```bash
python3 scripts/run_nfl_ab_design.py
```

如果运行成功，终端会看到类似输出：

```text
NfL antibody workflow complete.
Primary antigen fragment: NEFL 280-375
Rank  1: 7-H11-D3-2-C7
Rank  2: 15-C12-H6
```

也可以用模块方式运行：

```bash
PYTHONPATH=src python3 -m nfl_ab_design
```

## 6. 运行后应该看哪些结果

### 总报告

```text
NFL_AB_design/outputs/workflow_report.md
```

这是最适合先看的文件，里面按顺序总结了：

- 抗原截断推断；
- 表位窗口；
- 抗体序列检查；
- 候选排序；
- sandwich pair 结论；
- 外部结构工具输入文件。

### 抗原截断推断

```text
NFL_AB_design/outputs/00_antigen_truncation_report.md
NFL_AB_design/outputs/00_antigen_truncation_all_peptide_bonds.csv
NFL_AB_design/outputs/00_antigen_truncation_medium_high_sites.csv
NFL_AB_design/outputs/00_antigen_truncation_fragment_candidates.csv
```

这些文件回答：NfL 哪些位置可能被 cathepsin-like 蛋白酶切开？哪些片段最符合 22 kDa 二硫键二聚体假设？

### 抗原片段优先级

```text
NFL_AB_design/outputs/01_antigen_fragment_prioritization.csv
```

这个文件回答：最终哪些 NfL 片段最值得作为抗原重点分析？

### 表位窗口

```text
NFL_AB_design/outputs/02_epitope_windows.csv
```

这个文件回答：抗体可能识别 NfL 的哪些局部区域？

### 抗体序列检查

```text
NFL_AB_design/outputs/03_antibody_developability.csv
```

这个文件回答：两株抗体序列有没有明显的可开发性风险？

### 候选排序

```text
NFL_AB_design/outputs/05_candidate_ranking.csv
```

这个文件回答：在同一套计算评分下，哪些抗体候选排在前面？

### Sandwich pair 报告

```text
NFL_AB_design/outputs/06_sandwich_pair_report.md
```

这个文件回答：两个抗体是否适合组成检测抗体对？

## 7. 如何阅读 CSV 表格

CSV 是表格文件，可以用 Excel、Numbers 或 LibreOffice 打开。

常见字段含义：

- `fragment`：NfL 片段范围，例如 `280-375`；
- `monomer_avg_mass_kDa`：单体理论质量；
- `disulfide_homodimer_avg_mass_kDa`：二硫键二聚体理论质量；
- `N_terminal_cut`：N 端切割边界；
- `C_terminal_cut`：C 端切割边界；
- `combined_boundary_score`：两个边界的切割支持分；
- `epitope_priority_score`：表位优先级；
- `developability_score`：抗体可开发性评分；
- `total_rank_score`：候选抗体总排序分。

初学者可以先按这些列排序：

- `antigen_confidence_score`
- `epitope_priority_score`
- `developability_score`
- `total_rank_score`

分数越高，说明越值得优先关注。

## 8. 外部结构软件怎么接

代码包已经帮你准备了外部结构工具输入文件。

FASTA 输入：

```text
NFL_AB_design/outputs/exports/fasta/
```

AF3-style JSON 输入：

```text
NFL_AB_design/outputs/exports/af3_json/
```

外部任务表：

```text
NFL_AB_design/outputs/exports/external_jobs/pipeline_jobs.tsv
```

可编辑命令脚本：

```text
NFL_AB_design/outputs/exports/external_jobs/run_external_pipelines.sh
```

外部工具配置文件：

```text
NFL_AB_design/config/external_pipelines.example.json
```

默认情况下，外部任务都是关闭的，因为不同电脑上安装的软件命令不一样。等你确认本机已经安装好 IgFold、Chai-1、Boltz、Rosetta 等工具后，再修改配置文件里的命令和 `enabled` 字段。

## 9. 一节课可以怎么讲

可以按下面顺序讲：

1. 什么是 NfL，为什么它可以作为神经损伤标志物。
2. 什么是抗体，VH/VL 和 CDR 是什么。
3. 什么是 sandwich assay，为什么需要两个非重叠表位抗体。
4. 为什么 22 kDa 条带提示可能有 11 kDa 单体片段。
5. 为什么 Cys322 对二硫键二聚体推断很关键。
6. 运行代码，观察 `280-375` 片段如何被选出来。
7. 查看两个抗体为什么排在第 1 和第 2。
8. 讨论为什么计算结果还需要结构预测和实验验证。

## 10. 常见问题

### 这套代码是不是直接证明抗体一定有效？

不是。代码给出的是计算优先级和复现逻辑。真正有效需要实验验证。

### 为什么最后还要做结构预测？

因为表位是否真的不重叠、两个 Fab 是否会空间冲突，需要三维结构进一步确认。

### 高中生需要理解所有打分公式吗？

不需要。一开始只要理解每个分数代表什么方向：质量是否接近、边界是否合理、表位是否合适、抗体序列是否健康、两个抗体是否能同时结合。

### 最重要的结论是什么？

本流程从 NfL 片段推断开始，最终把两个实验阳性抗体排在最高优先级，并支持它们作为 sandwich 检测抗体对继续做结构和实验验证。

## 11. 最快操作清单

```bash
git clone <你的仓库地址>
cd NFL_AB_design
python3 --version
python3 scripts/run_nfl_ab_design.py
open outputs/workflow_report.md
```

如果 `open` 命令不能用，就直接在文件管理器或编辑器中打开：

```text
outputs/workflow_report.md
```
