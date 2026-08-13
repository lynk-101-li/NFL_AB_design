# NfL 单体抗体设计：模拟实验体系

## 目标与证据边界

本体系模拟“表位定义 → 从头生成 → 结构/界面筛选 → 可开发性筛选 → 多目标短名单 → 回顾性阳性对照 → 三模型真实 handoff”的完整计算实验。所有生成分数和漏斗结果都是确定性 proxy；RFantibody、IgGM、Germinal 尚未执行，候选也没有获得实验结合或可开发性验证。

设计只使用 NfL aa280–377 单链单体。两个表位为：

- `helix_surface_323_331`：序列 `RGMNEALEK`，待人工确认热点 `Met325/Leu329`；不包含 Cys322，也不要求半胱氨酸接触。
- `C_boundary_368_377`：序列 `YLKEYQDLLN`，待人工确认热点 `Tyr368/Tyr372/Leu375`。

任何二聚体或其他多链构象均不属于本 campaign。

## 模拟设计矩阵

| 维度 | 设置 |
|---|---|
| framework 来源 | `template_7-H11-D3-2-C7`、`template_15-C12-H6` |
| 设计区域 | H1/H2/H3/L1/L2/L3 全部重设计 |
| 目标表位 | 2 个单体表面窗口 |
| 每个 template×epitope | 1280 个模拟候选 |
| 初始候选总数 | 5120 |
| 已知抗体使用方式 | 仅提供 framework；已知 CDR 不进入生成特征 |
| 随机性 | 固定 seed `20260812` 的确定性模拟 |

## 漏斗与产物

1. `04_backbone_generation.csv` / `05_sequence_candidates.csv`：5120 个 prospective 模拟设计。
2. `06_structure_interface_screen.csv`：模拟 backbone、VH/VL packing、抗原界面、shape、clash 和模型分歧指标。
3. `07_developability_screen.csv`：模拟溶解性、聚集风险、化学 liability、免疫原性 proxy 和综合可开发性。
4. `08_screening_funnel.csv`：按结构、界面、可开发性和多目标阈值逐级过滤。
5. `09_prospective_candidates.csv`：所有通过多目标阈值的候选及全局 rank；`selected_for_export=True` 标出 12 个分层短名单。
6. `10_retrospective_demo_candidates.csv`：在 prospective 排名完成后才注入两株已知阳性。它们的 Top2 是回顾性对照演示，不是盲法发现。

最终 12 个采用 `template × epitope` 平衡策略：当四层均有合格候选时，每层选 3 个。每层内部按综合分选优；全局 rank 继续保留，以公开多样性配额造成的位次差异。

## 对照体系

- 前瞻阴性边界：prospective 表和导出 FASTA 中不得出现两株已知抗体的完整 VH/VL。
- 回顾性阳性对照：`7-H11-D3-2-C7`、`15-C12-H6` 只在最终独立表中注入并标记 `retrospective_positive_control`。
- 模拟证据标签：所有 proxy 指标必须保留 `*_is_simulated`、`data_status=simulated` 和机器可读 provenance。
- 结构输入对照：两株 template 分别使用独立 HLT PDB 和独立 scFv PDB，不复用同一坐标模板。

## 真实模型桥接（当前不执行）

同一套规范化请求分别映射到：

- RFantibody：paired-Fv HLT 模板、六 CDR 设计、RF2 复核；
- IgGM：H/L 六 CDR mask + 单链抗原 PDB；
- Germinal：独立 scFv 几何轨道。

`scripts/prepare_real_model_jobs.py` 只编译和校验 handoff；`scripts/execute_real_model_jobs.py` 只有在正式 target manifest、人工审核、runtime attestation、checkpoint hash 和显式 `--execute` 全部满足后才允许串行执行。当前正式 `config/target_structure_manifest.json` 不存在，因此真实模型保持 fail closed。

## 进入真实 canary 前的人工门禁

- 审核 AlphaFold DB 单体结构、323–331 和 368–377 的表面几何与裁剪边界风险；
- 最终确认热点 `325/329` 与 `368/372/375`；
- 审核两份 RFantibody HLT 和两份 Germinal scFv 模板；
- 填写 reviewer、带时区时间和 contracts acknowledgement；
- 在 AutoDL 上补齐 IgGM/Germinal 独立环境与 checkpoint 哈希；
- 先运行每个引擎一个共同 canary，再决定是否扩展到四个 template×epitope 作业。
