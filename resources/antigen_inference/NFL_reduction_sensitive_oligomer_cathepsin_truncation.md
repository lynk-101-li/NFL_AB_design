# NfL 还原敏感高阶复合物片段的 cathepsin 截断分析

日期：2026-07-08

## 1. 分析目的

已观察到 NfL 相关片段在非还原条件下约为 25–35 kDa，DTT 后约为 6–12 kDa。当前优先模型是两个 Cys322–Cys322 二硫键二聚体通过反平行、错位的卷曲螺旋界面形成四聚体；二聚体＋未知结合伙伴仍不能排除。目标是从 cathepsin 切割偏好出发，寻找可能产生单链片段的 NfL rod 区截断边界。

## 2. 输入与关键约束

使用序列：canonical human NEFL/NfL，UniProt P07196，长度 543 aa。

结构域标注采用：

| 区域 | 氨基酸位置 |
|---|---:|
| Head | 2-92 |
| IF rod | 90-400 |
| Coil 2B | 281-396 |
| Tail | 401-543 |

关键实验约束：

1. 非还原条件下片段约 25–35 kDa，DTT 后约 6–12 kDa。
2. 约四倍的表观变化支持但不能单独证明四聚体；备选解释为二硫键二聚体＋未知结合伙伴。
3. canonical human NfL 中只有一个半胱氨酸：Cys322。
4. 因此，能够进入当前二硫键模型的 NfL 单体片段必须包含 Cys322。
5. Cys322 位于 rod/coil 2B 区域，因此候选片段应集中在 coil 2B 附近。

## 3. 方法过程与结果

### 步骤 1：判断 SnapGene 是否可自动显示蛋白酶切位点

方法：

- 查询 SnapGene 功能边界。
- 确认 SnapGene 可导入蛋白序列、GenPept/UniProt 等蛋白记录。
- 确认 SnapGene 的酶切显示主要针对 DNA/RNA restriction enzyme sites，不是 cathepsin/protease 切割预测。

结果：

- SnapGene 不能直接像显示限制性内切酶位点那样自动预测 cathepsin 切点。
- 解决方案是生成带 feature 注释的 GenPept 蛋白文件，在 SnapGene 中作为 protein features 查看。

输出：

- `nfl_cathepsin_annotated_for_snapgene.gp`

该文件可直接用 SnapGene 打开，查看：

- NfL 结构域；
- Cys322；
- 还原敏感高阶复合物的 Cys322-containing 候选片段；
- coil 2B 附近中高优先级 cathepsin-like 切点。

### 步骤 2：建立 cathepsin-like 切点打分规则

方法：

对 NfL 全长 542 个肽键逐一打分。打分分为两类：

1. Cysteine cathepsins，包括 Cathepsin B/L/S/K/V-like：
   - 重点考虑 P2 位点的疏水或芳香氨基酸；
   - P1 碱性、极性或疏水残基作为辅助因素；
   - P1' 为 Pro 时降权。

2. Aspartyl cathepsins，包括 Cathepsin D/E-like：
   - 重点考虑 P1/P1' 的疏水或芳香残基；
   - Leu/Phe/Tyr/Trp 等残基附近切割优先级较高。

结果：

- 全部 542 个肽键均完成打分。
- 其中 34 个肽键被标记为 medium/high cathepsin-like candidate。
- coil 2B 附近的重点候选切点包括：

| 切点 | 上下文 | 主要解释 |
|---|---|---|
| W279\|F280 | AEEW\|FKSR | CatD/E-like，高疏水/芳香边界 |
| F280\|K281 | EEWF\|KSRF | CatD-like 或后续修剪边界 |
| K281\|S282 | EWFK\|SRFT | Cysteine cathepsin-like，P2=F，P1=K |
| Y368\|L369 | MARY\|LKEY | CatD/E-like，疏水/芳香边界 |
| L375\|L376 | YQDL\|LNVK | CatD/E-like，L-L 边界 |
| L376\|N377 | QDLL\|NVKM | 可作为相邻 C 端边界 |
| N377\|V378 | DLLN\|VKMA | 产生 282-377 候选片段 |
| K379\|M380 | LNVK\|MALD | Cysteine cathepsin-like，P2=V，P1=K |
| R390\|K391 | AAYR\|KLLE | Cysteine cathepsin-like，P2=Y，P1=R |
| L392\|L393 | YRKL\|LEGE | CatD/E-like，L-L 边界 |

输出：

- `nfl_all_peptide_bonds_cathepsin_scores.csv`
- `nfl_medium_high_cathepsin_candidate_sites.csv`

### 步骤 3：按还原条带质量范围与 Cys322 约束筛选片段

方法：

枚举所有满足以下条件的片段：

1. 片段包含 Cys322；
2. 片段位于 rod/coil 2B 附近；
3. 单体理论分子量约 10.0-12.6 kDa；
4. 单体理论质量与 DTT 后 6–12 kDa 的宽范围相容；
5. N 端和 C 端边界均由 cathepsin-like 打分支持。

结果：

优先候选片段如下：

表中的“二硫键二聚体理论质量”只表示优先四聚体模型中的共价中间体；对应两二聚体四聚体的序列理论质量约为该列的两倍。SDS-PAGE 表观迁移受构象、SDS 结合和复合物解离影响，不能要求其等于序列理论质量。

| 优先级 | 候选片段 | 单体理论质量 | 二硫键二聚体理论质量 | N 端切点 | C 端切点 | 解释 |
|---:|---|---:|---:|---|---|---|
| 1 | NfL 280-375 | 11.041 kDa | 22.081 kDa | W279\|F280 | L375\|L376 | 单体质量落在还原条带范围，两个边界均符合 CatD/E-like 疏水切割 |
| 2 | NfL 281-376 | 11.007 kDa | 22.013 kDa | F280\|K281 | L376\|N377 | 单体质量落在还原条带范围，边界位于 280/281 与 376/377 附近 |
| 3 | NfL 282-377 | 10.993 kDa | 21.985 kDa | K281\|S282 | N377\|V378 | 单体质量落在还原条带范围，N 端边界支持 cysteine cathepsin-like |
| 4 | NfL 282-376 | 10.879 kDa | 21.756 kDa | K281\|S282 | L376\|N377 | 同样支持 K281\|S282，C 端略短 |
| 5 | NfL 280-376 | 11.155 kDa | 22.307 kDa | W279\|F280 | L376\|N377 | 与 280-375 相邻，可作为 C 端修剪变体 |

输出：

- `nfl_reduction_sensitive_fragment_candidates.csv`

## 4. 当前结论

当前优先工作假设是：

> 非还原约 25–35 kDa 的 NfL 相关复合物由 aa 280–377 附近、包含 Cys322 的 rod/coil-2B 片段构成。两个 Cys322–Cys322 二硫键二聚体再通过卷曲螺旋疏水核心、静电作用和几何/构象互补形成反平行、错位四聚体。该模型优先于但不能排除“二聚体＋未知结合伙伴”。

从 cathepsin 角度看：

1. Cathepsin D/E-like 切割可解释 `W279|F280`、`L375|L376`、`L376|N377` 等疏水边界。
2. Cysteine cathepsin B/L/S/K/V-like 切割或修剪可解释 `K281|S282`，因为该位点具有 `EWFK|SRFT` 上下文，其中 P2=F、P1=K。
3. 上述候选片段都包含唯一 Cys322，因此能够形成模型中的 Cys322-Cys322 二硫键二聚体中间体。

因此，抗体设计上优先关注的区域应覆盖：

- Cys322 附近：`TLEIEACRGMNEALEK`
- 片段 N 端附近：aa 279-286，`AEEWFKSRFTVL`
- 片段 C 端附近：aa 368-377，`YLKEYQDLLN`
- 候选片段主体：aa 280-377

## 5. 后续实验验证建议

### 5.1 质谱验证片段身份

建议对非还原 25–35 kDa 条带和 DTT 后 6–12 kDa 条带分别做 LC-MS/MS，并用 C322S、SEC-MALS 或原生质谱判别化学计量。

重点验证：

1. 是否能检测到 Cys322 所在肽段：
   - `TLEIEACR`
   - `GMNEALEK`

2. 非还原条件下是否能检测到 Cys322-Cys322 二硫键同源二聚肽：
   - `TLEIEACR` - `TLEIEACR`

3. 还原烷基化后是否出现烷基化 Cys322 肽段：
   - `TLEIEA(Cam)R` 或 NEM 标记形式

### 5.2 验证真实 N/C 端

建议质谱搜索时加入：

- no-enzyme search；
- semi-tryptic search；
- multiple protease digestion，例如 trypsin、Glu-C、Lys-C、chymotrypsin。

重点寻找 neo-termini：

- N 端：F280、K281、S282；
- C 端：L375、L376、N377、V378。

### 5.3 体外 cathepsin 消化验证

建议使用重组 full-length NfL 或 aa 253-396 / aa 281-396 片段，分别加入：

- Cathepsin D；
- Cathepsin B；
- Cathepsin L；
- Cathepsin S；
- Cathepsin K。

条件：

- pH 4.5-5.5；
- 时间梯度；
- 非还原/还原 Western blot 或胶内质谱；
- 抑制剂对照：pepstatin A 用于 CatD/E，E64 用于 cysteine cathepsins，CA-074 用于 CatB 偏向验证。

判断标准：

如果某一 cathepsin 或 cathepsin 组合能重复产生非还原约 25–35 kDa、还原约 6–12 kDa 的片段，并且质谱端点落在 279–282 与 375–377 附近，则支持截断部分；四聚体化学计量仍需 SEC-MALS 或原生质谱确认。

### 5.4 排除二硫键人为形成

样本处理时建议立即加入 NEM 或 IAA 封闭游离 Cys，并比较：

- 未封闭样本；
- 采样后立即封闭样本；
- 还原后再烷基化样本。

如果立即封闭后非还原复合物明显减少，说明二硫键可能部分来自体外氧化；如果仍稳定存在，则支持体内或样本原位存在的 Cys322-Cys322 二聚体中间体，但仍不能单独证明最终为四聚体。

## 6. 文件说明

| 文件 | 用途 |
|---|---|
| `NFL_reduction_sensitive_oligomer_cathepsin_truncation.md` | 本报告 |
| `nfl_cathepsin_annotated_for_snapgene.gp` | SnapGene 可打开的 GenPept 蛋白注释文件 |
| `nfl_reduction_sensitive_fragment_candidates.csv` | Cys322-containing 单链片段候选，不把非还原条带当作精确二聚体质量约束 |
| `nfl_all_peptide_bonds_cathepsin_scores.csv` | NfL 全部肽键 cathepsin-like 打分表 |
| `nfl_medium_high_cathepsin_candidate_sites.csv` | 中高优先级切点筛选表 |
| `nfl_snapgene_cathepsin_notes.md` | SnapGene 使用说明 |
