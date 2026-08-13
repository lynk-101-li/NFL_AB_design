# 学生服务器部署指南

本目录用于学生在自己的 Linux GPU 服务器上从零部署本项目。九个模型的固定源码通过 Git submodule 随项目交付；仓库不假定环境已经安装，也不包含模型权重、AlphaFold 参数、数据库、PyRosetta、服务器凭据或真实运行结果。

## 1. 推荐硬件与存储

- Ubuntu 22.04 或兼容 Linux；
- NVIDIA GPU，驱动可运行 `nvidia-smi`；
- 48GB 显存适合多数最小 canary；Germinal、大复合物和部分 AF3 配置可能需要更大显存；
- 系统盘建议 80GB，总安装前至少 35GB 可用；
- 持久化文件盘建议至少 100GB，用于权重、参数和结果；
- 同一 GPU 一次只运行一个模型任务。

如果使用 AutoDL，推荐目录：

```text
/root/workspace/NFL_AB_design/third_party/  固定的上游源码
/root/miniconda3/envs/                      每个模型的独立环境
/root/autodl-fs/shared/models/           可复用模型权重
/root/autodl-fs/topics/nfl-antibody-design/runs/  正式运行
```

不要把 Conda 环境放到网络文件存储。不要把缓存、环境、权重或真实结果提交到 GitHub。

## 2. 克隆项目并运行本地模拟

```bash
git clone --recurse-submodules https://github.com/lynk-101-li/NFL_AB_design.git
cd NFL_AB_design
python3 scripts/verify_model_components.py
python3 -m pip install -e .
python3 scripts/run_nfl_ab_design.py
python3 -m unittest discover -s tests -v
```

这里运行的是确定性模拟，不会调用九个真实模型。

## 3. 服务器预检

```bash
bash deploy/autodl/preflight.sh
```

预检只读，不安装软件。返回 `blocked_*` 时先处理硬件、磁盘或命令缺失。环境安装会产生大量网络流量和磁盘占用，不应绕过门禁。

## 4. 初始化固定源码并按需安装

先阅读 [模型部署与运行教程](../../docs/real_model_installation.md)，然后初始化并校验全部源码：

```bash
bash deploy/autodl/bootstrap_models.sh sources
```

`config/model_components.json` 是九个组件的机器可读清单。安装脚本按需创建隔离环境；它不会替学生接受许可证，也不会运行本项目真实任务。

也可逐个安装：

```bash
bash deploy/autodl/bootstrap_models.sh rfantibody
bash deploy/autodl/bootstrap_models.sh iggm
bash deploy/autodl/bootstrap_models.sh immunebuilder
bash deploy/autodl/bootstrap_models.sh chai1
bash deploy/autodl/bootstrap_models.sh boltz2

# 以下组件必须先阅读相应条款，再设置教程指定的 ACK 变量
bash deploy/autodl/bootstrap_models.sh germinal
bash deploy/autodl/bootstrap_models.sh tfold
bash deploy/autodl/bootstrap_models.sh igfold
bash deploy/autodl/bootstrap_models.sh alphafold3
```

## 5. 审核输入并编译 canary

仓库中的 `input/structures/target_structure_manifest.review_candidate.blocked.json` 是待审候选，不可直接执行。学生必须完成 [目标结构审核](../../docs/target_structure_review.md)，再从 `config/target_structure_manifest.example.json` 创建本地文件：

```text
config/target_structure_manifest.json
```

该文件默认被 `.gitignore` 排除。填写 reviewer、带时区审核时间、完整坐标映射、最终热点及 `contracts_acknowledged: true` 后，编译共同 canary：

```bash
python3 scripts/prepare_real_model_jobs.py \
  --target-manifest config/target_structure_manifest.json \
  --profile smoke \
  --job-scope canary \
  --canary-template-id template_7-H11-D3-2-C7 \
  --canary-epitope-id C_boundary_368_377 \
  --iggm-repo-dir "$PWD/third_party/IgGM" \
  --germinal-repo-dir "$PWD/third_party/germinal" \
  --germinal-af-params-dir "${NFL_GERMINAL_AF_PARAMS:?set this path}"
```

编译器验证全部 `2 templates × 2 epitopes`，但 smoke 只把每个引擎的同一个 canary 放入 `execution_jobs`。它不执行模型。

## 6. Runtime attestation 与执行

学生应实际验证 GPU import、上游 revision 和 checkpoint 文件，再按 `config/runtime_attestation.example.json` 创建本地 `runtime_attestation.json`。不得填写虚构哈希。

先 dry-run：

```bash
python3 scripts/execute_real_model_jobs.py \
  --handoff-manifest /absolute/path/to/unified_handoff_manifest.json \
  --runtime-attestation /absolute/path/to/runtime_attestation.json
```

审阅命令、路径和输入哈希后才执行：

```bash
python3 scripts/execute_real_model_jobs.py \
  --handoff-manifest /absolute/path/to/unified_handoff_manifest.json \
  --runtime-attestation /absolute/path/to/runtime_attestation.json \
  --execute
```

任何命令失败都会停止；修复后使用 `--execute --resume`。真实输出只能写入 `real_runs/` 或服务器的持久化 run 目录，不能覆盖 `outputs/` 中的模拟结果。

## 7. 结构预测复核

RFantibody、IgGM、Germinal 的真实生成候选进入同一个候选池后，先用 IgFold/ImmuneBuilder 检查 Fv，再用 tFold、AF3、Chai-1、Boltz-2 对复合物进行交叉预测。每种模型的安装、命令和输入边界见 [完整教程](../../docs/real_model_installation.md)。模型之间不要求分数数值相等，重点检查表位、结合方向、CDR 接触和界面几何是否收敛。

## 8. 学生验收标准

- 本地测试全部通过；
- 九个 submodule 的实际 commit 与配置一致；
- GPU smoke 通过；
- 所有 checkpoint/参数都有实际 SHA-256；
- 正式 manifest 已人工审核；
- dry-run 只显示三个共同 canary job；
- 执行报告明确区分 planned/running/succeeded/failed；
- canary 成功后才允许 `--profile full --job-scope all`。
