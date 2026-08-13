# 单体螺旋表面唯一构象范围

## 当前决策

本 campaign **只运行单体轨道**：NfL 280–377 单链单体。项目不建立二聚体结构、不生成二聚体候选、不运行二聚体模型，也不进行跨构象排名。

第一表位为 `helix_surface_323_331`，序列 `RGMNEALEK`。单体轨道仅使用 `Met325/Leu329` 这组 i/i+4 热点限定同一段 α 螺旋的紧凑暴露面。Cys322 不属于表位窗口、热点、名称或直接接触约束。

上游片段在非还原条件下约为 25–35 kDa，DTT 后约为 6–12 kDa。当前优先解释是两个 Cys322–Cys322 二硫键二聚体通过反平行、错位卷曲螺旋界面形成四聚体；二聚体＋未知结合伙伴仍不能排除。该生化解释不参与本轮抗体生成或排名。

## 单体轨道

- 抗原：NfL 280–377 单链单体候选结构。
- 第一表位热点：`325/329`。
- 排除的直接接触约束：`322`。
- 目的：运行 RFantibody、IgGM 和 Germinal 的共同单链 handoff，并在同一单体构象范围内生成和筛选候选。

## 结果解释边界

所有候选、漏斗和排名只对当前 AFDB 单体输入成立。它们不能确认非还原复合物是四聚体还是“二聚体＋未知伙伴”，也不能直接外推任何多聚体中的表位可及性或局部几何。多链状态不会被复制、压平或伪装成单链输入。

## 生化模型与证据边界

- 优先模型：`4M → (M–S–S–M)₂`。两个二硫键分别稳定两个二聚体，两个二聚体再形成反平行、错位四聚体。
- 二聚体间非共价界面：卷曲螺旋疏水核心、静电作用和几何/构象互补共同贡献，不能只写成“疏水作用”。
- DTT 能证明的范围：二硫键对复合物稳定性重要。
- DTT 不能单独证明的范围：四聚体化学计量、二聚体间界面组成，以及是否存在未知结合伙伴。
- 判别实验：C322S；还原/非还原 SDS-PAGE；SEC-MALS 或原生质谱；必要时辅以交联质谱。

中间丝的一般装配框架为平行同向卷曲螺旋二聚体，再形成反平行、错位四聚体并进入更高阶组装；NfL 碎片是否严格遵循同一界面仍需实验确认。

参考背景：Herrmann 与 Aebi 对中间丝结构/装配的综述，以及中间丝四聚体结构工作均支持“平行卷曲螺旋二聚体 → 反平行错位四聚体”的一般路径：

- Herrmann H, Aebi U. *Intermediate Filaments: Structure and Assembly*. Cold Spring Harb Perspect Biol. 2016. https://doi.org/10.1101/cshperspect.a018242
- Eldirany SA et al. *Human keratin 1/10-1B tetramer structures reveal a knob-pocket mechanism in intermediate filament assembly*. EMBO J. 2019. https://doi.org/10.15252/embj.2018100741

机器可读范围见 `config/antigen_conformation_tracks.json`。
