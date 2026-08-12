# RFantibody、IgGM 与 Germinal 安装说明

本页记录本项目固定的上游版本及学生部署边界。上游仓库可能更新；本项目以 `config/real_model_backends.json` 的 commit 为准。

## RFantibody

- 仓库：<https://github.com/RosettaCommons/RFantibody>
- 固定 revision：`8fe311415754e0276d1a39c87c57e69c88927a2d`
- 许可：MIT
- 环境：Python 3.10、PyTorch 2.3、CUDA 11.8
- 官方安装：下载权重后运行 `uv sync`

最小验证：

```bash
cd /root/apps/RFantibody
git rev-parse HEAD
uv run python -c 'import torch; assert torch.cuda.is_available(); print(torch.cuda.get_device_name())'
uv run rfdiffusion --help
uv run proteinmpnn --help
uv run rf2 --help
```

RFantibody 需要真实 target PDB、两份各自独立的 H/L/T 模板、精确 CDR REMARK 和人工热点；FASTA 不能替代坐标输入。

## IgGM

- 仓库：<https://github.com/TencentAI4S/IgGM>
- 固定 revision：`06abc563b3fc8c7ea020543add16b69b6f8a1c8d`
- 许可：MIT
- 上游锁定环境：Python 3.10.14、PyTorch 2.0.1、CUDA 11.7
- 权重：约 3.3GB，可由上游代码自动下载或手动放入 `checkpoints/`

最小验证：

```bash
cd /root/apps/IgGM
git rev-parse HEAD
conda run -n iggm python -c 'import torch; assert torch.cuda.is_available(); print(torch.cuda.get_device_name())'
conda run -n iggm python design.py --help
```

输入 FASTA 必须包含 `>H`、`>L`、`>A`；`A` 的序列和顺序必须与 PDB 链完全一致，表位位置是 PDB 抗原序列的局部 1-based 坐标。

## Germinal

- 仓库：<https://github.com/SantiagoMille/germinal>
- 固定 revision：`1e1c1a5b79884ae45abae030c9df90d9423a990a`
- 仓库代码许可：Apache-2.0
- 额外许可：PyRosetta 需要单独许可；其他模型和参数也遵循各自条款
- 推荐：CUDA 12+、40GB 以上显存；较大任务可能需要 60GB 以上

学生必须自行合法获取 PyRosetta 和 AlphaFold-Multimer 参数。项目脚本不分发这些资产。最小验证：

```bash
cd /root/apps/germinal
git rev-parse HEAD
conda run -n germinal python validate_install.py
```

本项目把 Germinal 作为独立的单链 scFv 几何轨道。它的输出不能直接声称保持 paired-Fv 的原生 VH/VL 几何。

## 许可与引用

本仓库的 MIT 许可证只覆盖本仓库代码，不会重许可第三方仓库、模型权重或参数。用于教学、论文或商业工作前，学生必须阅读并遵守每个上游项目的许可证、模型条款与引用要求。
