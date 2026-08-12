# 抗体模板坐标审核

## 状态

**`blocked_pending_human_review`**：本工具只把已有坐标整理成 RFantibody HLT 和 Germinal scFv 输入，不运行结构模型、网络请求或外部进程。生成的 manifest 只是阻断状态片段，不是正式 target manifest。

## 中性 CDR seed

这些序列仅为坐标初始 seed，不是候选抗体或模型预测。它们只替换由锁定的 ANARCI 2020.04.23 Chothia 结果确定的 CDR raw positions；长度不变，每段均只含 20 种标准氨基酸，且与对应已知阳性 CDR 不同。

| 模板 | H1 | H2 | H3 | L1 | L2 | L3 |
|---|---|---|---|---|---|---|
| `template_7-H11-D3-2-C7` | `ASGTQSS` | `GSTASQ` | `GST` | `QASSSGTSTLA` | `GASTSAS` | `QQSGTSPRT` |
| `template_15-C12-H6` | `GASSTQSA` | `STAGQ` | `SGTQADSGTAY` | `RASQSGTSTLA` | `SASTGAS` | `QQAGTSPAT` |

完整 neutral VH/VL 会记入 `antibody_template_evidence.json`。不允许以已知阳性完整 VH/VL 或已知 CDR 生成坐标输入。

## 精确 Chothia 合同

- 完整 raw position→Chothia label 映射来自锁定的 `nfl_H.csv` 和 `nfl_KL.csv`，两个文件哈希必须与 `chothia_numbering_evidence.json` 一致。
- HLT 中的 `REMARK PDBinfo-LABEL` 是 H 后接 L 的绝对 pose index，只标记本 campaign 的精确 CDR：模板 7 长度为 `7/6/3/11/7/9`，模板 15 长度为 `8/5/11/11/7/9`。
- RFantibody 官方 `chothia2HLT.py` 的宽区间 REMARK 会被丢弃和重建，因为宽区间会把邻近 framework 也标为 loop。该规范化会在 evidence 中显式记录，不会声称是官方 converter 原样输出。
- 输入 PDB 必须仅含 H/L 两条链，序列必须精确等于 neutral seed，每个 Chothia 编号和插入码必须与原始 ANARCI 表一致。

## 几何和原子检查

- 所有残基必须是标准氨基酸，具有完整规范重原子，不允许多 model、重复原子或未处理 altloc。
- 每一个链内相邻残基的 peptide `C–N` 距离必须在 `1.15–1.55 Å`；未精修 ABodyBuilder2 坐标若越界会直接拒绝，不会因 manifest 未审核而放行。
- scFv 的 VH 和 VL 分别使用排除六 CDR 后的 framework CA 对应作 Kabsch 刚体对齐；对齐对数、RMSD、linker 两端 peptide `C–N` 距离均写入 evidence。
- 按 Germinal adapter 的明确切片合同（前 `len(VH)`、后 `len(VL)`），锁定官方 `pdbs/scfv.pdb` 中间间隔序列是 `AGGGGSGGGGSGGGS`。它不是标准 `(GGGGS)3`；工具原样保留这 15 aa 的残基身份和坐标，不会删除、突变或猜测域边界。evidence 会显式记录“官方 observed linker、非标准 `(GGGGS)3`、无坐标/残基修改”。

## 用法

```bash
python3 scripts/prepare_antibody_template_inputs.py \
  --paired-fv template_7-H11-D3-2-C7=/path/to/template7_refined_chothia_HL.pdb \
  --paired-fv template_15-C12-H6=/path/to/template15_refined_chothia_HL.pdb \
  --official-broad-hlt template_7-H11-D3-2-C7=/path/to/template7_official_broad_HLT.pdb \
  --official-broad-hlt template_15-C12-H6=/path/to/template15_official_broad_HLT.pdb \
  --prepared-scfv template_7-H11-D3-2-C7=/path/to/template7_refined_scfv.pdb \
  --prepared-scfv template_15-C12-H6=/path/to/template15_refined_scfv.pdb \
  --neutral-heavy-anarci-csv /path/to/neutral_chothia_H.csv \
  --neutral-light-anarci-csv /path/to/neutral_chothia_KL.csv \
  --upstream-provenance /path/to/template_coordinate_generation_provenance.json \
  --generic-scfv /path/to/germinal_generic_scfv.pdb
```

默认输出到 `input/template_structures/`。已有文件时默认拒绝覆盖；审核后确实需要重生成时才使用 `--overwrite`。
`--official-broad-hlt` 是可选的审计证据：工具会验证它与 paired-Fv 原子坐标的最大差不超过 `0.002 Å`，记录宽区间 REMARK 计数和 SHA-256，然后丢弃这些宽标签。

`--prepared-scfv` 必须两份同时提供，不接受只有一份。这个入口用于已在外部做过受约束精修的 domain-graft scFv；工具不信任“已精修”声明，仍独立强制检查单链 A、连续编号、精确 `neutral VH + AGGGGSGGGGSGGGS + neutral VL` 序列、重原子、所有 peptide `C–N` 距离和严重 clash。通过后 evidence 标记 `scfv_coordinate_source=externally_refined_domain_graft`并记录输入路径/哈希；默认仍由本地 graft 生成并 fail closed。

neutral ANARCI 两份 CSV 也必须同时提供。它们用于证明替换 neutral seed 后的整链 Chothia label/插入码布局没有变化；原始已知抗体 ANARCI 表仍是 CDR raw range 的锁定证据。

`--upstream-provenance` 是可选 JSON object，schema 必须为 `nfl_ab_design.template_coordinate_generation_provenance.v1`。其 `input_artifact_sha256` 必须精确包含 `paired_fv_pdbs`、`prepared_scfv_pdbs`、`generic_scfv_pdb` 和 `official_broad_hlt_pdbs`，并与本次 CLI 实际文件哈希逐项一致。工具只用它绑定文件来源，不信任或证实其中的模型/精修结论。evidence 中的 `preparation_cli_model_or_network_execution_performed: false` 仅描述本地准备 CLI，不否认上游曾运行结构模型。

## 人工审核清单

- [ ] 确认两份 paired-Fv 均来自 neutral VH/VL，而非已知阳性 CDR。
- [ ] 核对 H/L 链、完整序列、Chothia 插入码和六组 exact HLT REMARK。
- [ ] 可视化检查两份 paired-Fv 的 CDR loop、VH/VL packing、clash 和异常键长。
- [ ] 审核 scFv 的 framework CA 对齐 RMSD、linker 构象及两个域–linker接口。
- [ ] 确认两份 HLT 和两份 scFv 文件的 SHA-256 彼此独立且与 evidence 一致。
- [ ] 在正式 target manifest 中手工记录 reviewer、含时区时间和 `contracts_acknowledged: true`；不直接改名或使用本工具的 blocked fragment。
