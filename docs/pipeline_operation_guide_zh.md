# NfL 抗体设计 pipeline 完整操作指南（学生版）

## 1. 这条 pipeline 做什么

本项目从 NfL/NEFL 抗原片段和候选表位出发，准备两个 framework-only 抗体模板，然后分别使用 RFantibody、IgGM 和 Germinal 生成抗体。Germinal 内部还会经过 AlphaFold-Multimer 式优化、Chai 复折叠、PyRosetta 界面过滤和 AbMPNN 重设计。

项目同时包含两条不同的证据轨道：

- **Proxy workflow**：本地 CPU 可重复运行，用模拟指标演示候选漏斗。这些分数不是 pLDDT、PAE、DockQ 或亲和力。
- **Real-model workflow**：在 GPU 服务器上运行真实模型，需要独立环境、checkpoint、runtime attestation 和明确授权。

> 技术运行成功不等于获得可实验候选；计算候选也不等于已证明结合。

## 2. 学生开始前必须知道的四个概念

### 2.1 Run、attempt 和 wave

- **Run**：一个稳定的科学目标，例如“对 NfL 368–377 表位设计抗体”。
- **Attempt**：实现、输入、参数或资产发生改变时的新尝试。不要覆盖旧 attempt。
- **Wave**：一个 attempt 内有界的小批次，用来控制 GPU 成本和早停。

### 2.2 Manifest

`unified_handoff_manifest.json` 是“要运行什么”的不可变合同，其中 `engines[].execution_jobs` 才是权威执行队列。`native_plan` 和 `jobs` 只用于审计，不得整表提交。

### 2.3 Runtime attestation

`runtime_attestation.json` 说明“现在这台服务器是否真的能跑”，必须绑定：

- handoff ID 和 manifest SHA-256；
- 实际安装的上游 revision；
- checkpoint/参数/权重 manifest SHA-256；
- 明确的 `ready=true`。

### 2.4 Dry-run 与 `--execute`

不加 `--execute` 时，executor 只验证输入、哈希、路径、roster 和 staging 计划，不启动模型。真实 GPU 运行必须有当前 attempt/wave 的明确授权。

## 3. 软件和目录准备

### 3.1 本地 proxy 环境

```bash
cd NFL_AB_design
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -e '.[visualization]'
python3 -m unittest discover -s tests
```

### 3.2 AutoDL 推荐目录

```text
/root/apps/                       固定源码/worktree
/root/workspace/NFL_AB_design/    项目代码
/root/miniconda3/envs/            隔离环境
/root/.cache/                     本地编译缓存
/root/autodl-fs/shared/models/    持久化权重
/root/autodl-fs/topics/.../runs/  正式 run/attempt/result
```

先读：

```bash
sed -n '1,120p' /root/workspace/AUTODL_CONTEXT.md
/root/workspace/pipeline/status.sh
nvidia-smi
screen -ls
```

## 4. 可视化统一用法

先复制配置：

```bash
cp config/pipeline_visualization.example.json config/pipeline_visualization.json
```

将服务器的真实路径填入 `paths` 字段。运行全部图：

```bash
python3 scripts/pipeline_visualization/run_all.py \
  --config config/pipeline_visualization.json \
  --output-dir visualizations/pipeline
```

- 练习/尚未运行的阶段：缺失真实路径时会生成 `NOT RUN` 占位图。
- 正式验收：加 `--strict`，任何必要数据缺失都会非零退出。
- 默认生成 180 dpi PNG，不修改原始数据。

## 5. Pipeline 逐步操作

## Step 01：抗原片段和表位

**目的**：从生化边界、序列和结构中选出要建模的 NfL 片段和热点。

**输入**：

- `input/antigen_truncation/truncation_constraints.json`
- `config/design_campaign.json`
- `input/structures/NEFL_P07196_AFDB_v6_280-377_chainA.pdb`

**操作**：

```bash
python3 scripts/run_nfl_ab_design.py
```

**重点结果**：

- aa280–377 是结构建模上下文；
- 设计表位为 `helix_surface_323_331` 和 `C_boundary_368_377`；
- 368–377 当前热点为 368/372/375。

**图**：

```bash
python3 scripts/pipeline_visualization/01_antigen.py --config config/pipeline_visualization.json
```

图中 proxy 分数只能用于本漏斗内部比较。

## Step 02：抗体 framework-only 模板

**目的**：保留 VH/VL framework，将 H1/H2/H3/L1/L2/L3 全部设为待设计区，防止已知阳性 CDR 泄漏。

**输入**：`input/antibody_templates/` 中的 FASTA、Chothia 编号证据和两份结构模板。

**操作**：

```bash
python3 scripts/prepare_antibody_template_inputs.py --help
python3 scripts/run_nfl_ab_design.py
```

**验收**：

- 两个 template ID 对应两份不同坐标；
- 六个 CDR 都在 `design_regions`；
- `known_cdr_sequences_used_for_generation=False`。

**图**：

```bash
python3 scripts/pipeline_visualization/02_templates.py --config config/pipeline_visualization.json
```

## Step 03：目标结构和坐标映射

**目的**：将完整 NEFL 坐标、target PDB 残基编号和各模型的链内编号对齐。

**操作**：

```bash
python3 scripts/prepare_target_structure.py --help
cp config/target_structure_manifest.example.json config/target_structure_manifest.json
```

完成人工复核后，填写 reviewer、带时区 timestamp 和 review contract，再将 `execution_state` 改为 `reviewed_ready_for_handoff`。

**必查**：

- target 只有授权的链；
- 热点残基类型与序列一致；
- 原始 `A368/A372/A375` 映射到 Chai/PyRosetta 链内 `A89/A93/A96`。

**图**：

```bash
python3 scripts/pipeline_visualization/03_target.py --config config/pipeline_visualization.json
```

## Step 04：本地 proxy 筛选漏斗

**目的**：在不占用 GPU 的情况下验证数据契约、分层筛选和候选导出逻辑。

```bash
python3 scripts/run_nfl_ab_design.py
python3 -m unittest discover -s tests
```

**输出**：`outputs/00_*` 至 `outputs/11_*`。

**不得声称**：proxy pass 不代表真实结构、结合或可开发性 pass。

**图**：

```bash
python3 scripts/pipeline_visualization/04_proxy_funnel.py --config config/pipeline_visualization.json
```

## Step 05：编译 handoff 和 runtime attestation

**目的**：把已复核的输入编译成不可变的实际执行队列。

```bash
python3 scripts/prepare_real_model_jobs.py \
  --target-manifest config/target_structure_manifest.json \
  --output-dir /root/autodl-fs/topics/nfl-antibody-design/runs/<RUN>/attempts/<ATTEMPT>/handoff \
  --profile smoke \
  --job-scope canary \
  --canary-template-id template_7-H11-D3-2-C7 \
  --canary-epitope-id C_boundary_368_377 \
  --iggm-repo-dir /root/workspace/NFL_AB_design/third_party/IgGM \
  --iggm-python /root/miniconda3/envs/iggm/bin/python \
  --germinal-repo-dir /root/apps/<PINNED_GERMINAL_WORKTREE> \
  --germinal-python /root/workspace/pipeline/germinal_python_chai.sh \
  --germinal-af-params-dir /root/autodl-fs/shared/models/alphafold_multimer/params
```

对于 Germinal-only pilot，再加：

```bash
--germinal-profile pilot --execution-engine Germinal
```

使用 `config/runtime_attestation.example.json` 为开始，根据 manifest 和实际资产 manifest 填写 attestation。不得复用旧 handoff 的哈希。

**Dry-run**：

```bash
python3 scripts/execute_real_model_jobs.py \
  --handoff-manifest <HANDOFF>/unified_handoff_manifest.json \
  --runtime-attestation <HANDOFF>/runtime_attestation.json
```

**图**：

```bash
python3 scripts/pipeline_visualization/05_handoff.py --config config/pipeline_visualization.json
```

### Step 06–10 的执行原则（必读）

Step 06–10 不是五个互相独立的 shell 任务。它们是 Step 05 所编译的 `engines[].execution_jobs` 内的模型阶段，由同一 executor 按 manifest 顺序执行。这样设计是为了防止学生：

- 手工复制过期命令；
- 绕过 runtime attestation 和输入哈希；
- 意外执行未选中的 2×2 全部 job；
- 在单 GPU 上并发启动多个大模型。

下面的命令用于所有 Step 06–10。首先设置当前 handoff：

```bash
export HANDOFF=/root/autodl-fs/topics/nfl-antibody-design/runs/<RUN>/attempts/<ATTEMPT>/handoff
export MANIFEST="$HANDOFF/unified_handoff_manifest.json"
export ATTESTATION="$HANDOFF/runtime_attestation.json"
```

查看真正会执行的命令（只读）：

```bash
python3 scripts/show_real_model_execution_plan.py \
  --manifest "$MANIFEST"
```

执行 dry-run：

```bash
python3 scripts/execute_real_model_jobs.py \
  --handoff-manifest "$MANIFEST" \
  --runtime-attestation "$ATTESTATION"
```

只有在当前 attempt/wave 得到明确真实运行授权后，才启动：

```bash
export SESSION=nfl_<ATTEMPT>_<WAVE>
screen -L -Logfile "$HANDOFF/execute.screen.log" -dmS "$SESSION" \
  conda run --no-capture-output -n iggm \
  python3 scripts/execute_real_model_jobs.py \
    --handoff-manifest "$MANIFEST" \
    --runtime-attestation "$ATTESTATION" \
    --execute
```

如果只想运行一个 engine，不得手工删改 manifest。应在新 attempt 中重新执行 Step 05，加入下列选项之一：

```bash
--execution-engine RFantibody
--execution-engine IgGM
--execution-engine Germinal
```

## Step 06：RFantibody

**目的**：生成抗体骨架和六 CDR 序列，再用 RF2 复核结构。

**内部次序**：

```text
prepare output
  → RFdiffusion
  → RFdiffusion score export
  → ProteinMPNN
  → RF2
  → RF2 score export
  → final PDB extraction
```

**执行命令**：

1. 全三引擎 handoff：使用上一节的统一 executor，RFantibody 会作为第一个 engine 执行。
2. RFantibody-only：用 Step 05 原命令编译新 handoff，加 `--execution-engine RFantibody`，然后执行：

```bash
python3 scripts/show_real_model_execution_plan.py \
  --manifest "$MANIFEST" --engine RFantibody

python3 scripts/execute_real_model_jobs.py \
  --handoff-manifest "$MANIFEST" \
  --runtime-attestation "$ATTESTATION"          # dry-run

# 授权后使用上一节的 screen + --execute 命令
```

`show_real_model_execution_plan.py` 会显示每个 RFdiffusion/ProteinMPNN/RF2 stage 的 cwd 和 argv，但绝不执行。

**验收**：最终 PDB 非空、链正确、Quiver/score companion 齐全，且所有命令零退出。

**图**：

```bash
python3 scripts/pipeline_visualization/06_rfantibody.py --config config/pipeline_visualization.json
```

## Step 07：IgGM

**目的**：在 H/L/A 三链输入上生成配对 VH/VL 抗体结构。

**执行命令**：

- 全三引擎 handoff：由统一 executor 在 RFantibody 成功后串行执行 IgGM。
- IgGM-only：在新 handoff 编译命令中加 `--execution-engine IgGM`。

```bash
python3 scripts/show_real_model_execution_plan.py \
  --manifest "$MANIFEST" --engine IgGM

python3 scripts/execute_real_model_jobs.py \
  --handoff-manifest "$MANIFEST" \
  --runtime-attestation "$ATTESTATION"          # dry-run

# 授权后使用通用 screen + --execute 命令
```

**验收**：每个 job 至少有一份 PDB 和一份 FASTA；残基数和链 ID 符合 manifest。

**图**：

```bash
python3 scripts/pipeline_visualization/07_iggm.py --config config/pipeline_visualization.json
```

## Step 08：Germinal 起始几何

**目的**：将 scFv CDR 朝向热点，避免上游固定平移对任意 PDB 坐标系的敏感性。

**操作命令**：这一步在 `prepare_real_model_jobs.py` 编译期间生成，不消耗 GPU。Germinal-only pilot 必须在 Step 05 命令中加：

```bash
--germinal-profile pilot --execution-engine Germinal
```

生成后检查授权的 Germinal job 和 seed 哈希：

```bash
python3 scripts/show_real_model_execution_plan.py \
  --manifest "$MANIFEST" --engine Germinal

python3 - <<'PY'
import json, os
manifest = json.load(open(os.environ["MANIFEST"]))
job = next(e for e in manifest["engines"] if e["engine"] == "Germinal")["execution_jobs"][0]
for artifact in job["input_artifacts"]:
    if artifact["role"] == "prepositioned_starting_complex":
        print(artifact["path"], artifact["sha256"])
PY
```

**门禁**：

- target–binder 重原子不得小于 2.5 Å；
- 热点–CDR 最近距离不大于 12 Å；
- seed PDB 路径和 SHA-256 写入 handoff。

**图**：

```bash
python3 scripts/pipeline_visualization/08_germinal_seed.py --config config/pipeline_visualization.json
```

## Step 09：Germinal 优化

**目的**：用梯度优化寻找抗体–抗原界面。

**执行命令**：Step 09 和 Step 10 属于同一个 Germinal job，不得在中间手工拆分或更换输出。

```bash
python3 scripts/show_real_model_execution_plan.py \
  --manifest "$MANIFEST" --engine Germinal

python3 scripts/execute_real_model_jobs.py \
  --handoff-manifest "$MANIFEST" \
  --runtime-attestation "$ATTESTATION"          # dry-run

# 当前 Germinal attempt/wave 获得授权后：
screen -L -Logfile "$HANDOFF/germinal.execute.screen.log" -dmS "$SESSION" \
  conda run --no-capture-output -n iggm \
  python3 scripts/execute_real_model_jobs.py \
    --handoff-manifest "$MANIFEST" \
    --runtime-attestation "$ATTESTATION" \
    --execute
```

Pilot 当前使用：

```text
logits_steps=60
softmax_steps=35
search_steps=10
max_trajectories=4
num_seqs=8
max_mpnn_sequences=4
```

**初始置信度门禁**：

- pLDDT > 0.80
- iPTM > 0.68
- iPAE < 0.27

不得为了获得候选而临时放宽门槛。

**图**：

```bash
python3 scripts/pipeline_visualization/09_germinal_optimization.py --config config/pipeline_visualization.json
```

## Step 10：Chai、PyRosetta 和 AbMPNN

**目的**：

1. Chai 独立复折叠抗体–抗原复合物；
2. PyRosetta 计算碰撞、热点接触、界面几何和能量项；
3. 通过初筛的结构进入 AbMPNN 重设计；
4. 对 redesign candidates 重新执行严格最终过滤。

**执行命令**：无需再启动第二个命令；Step 09 的 Germinal executor 会在每条轨迹达到门槛后自动执行 Step 10。不要同时手工运行 Chai 或 AbMPNN。

仅在修复后验证既有轨迹、且不重新设计序列时，可以使用 filter-only 诊断命令：

```bash
cd /root/apps/<ATTESTED_GERMINAL_WORKTREE>
PYTHONPATH="$PWD" /root/workspace/pipeline/germinal_python_chai.sh \
  /root/workspace/NFL_AB_design/scripts/validate_germinal_existing_trajectory.py \
  --final-config <PREVIOUS_RESULT>/final_config.yaml \
  --trajectory-pdb <PREVIOUS_RESULT>/trajectories/structures/<DESIGN>.pdb \
  --output-dir <NEW_ATTEMPT>/filter_only/<DESIGN> \
  --external-target-hotspots A89,A93,A96
```

filter-only 是修复诊断，不是候选晋级捷径；输出必须进入新 attempt，不得覆盖原轨迹。

**编号规则**：

- hallucination 使用原 PDB 编号 `A368/A372/A375`；
- Chai restraint 使用链内 `Y93`；
- Chai/PyRosetta 过滤使用 `A89/A93/A96`。

**pDockQ2**：Chai 必须提供完整 PAE；当且仅当复合物有两条链时，允许使用 binder-chain 本地 pDockQ2。多链复合物仍 fail-closed。

**图**：

```bash
python3 scripts/pipeline_visualization/10_filters.py --config config/pipeline_visualization.json
```

## Step 11：执行、监控和终态

真实执行：

```bash
screen -L -Logfile <ATTEMPT>/execute.screen.log -dmS <SESSION> \
  conda run --no-capture-output -n iggm \
  python3 scripts/execute_real_model_jobs.py \
    --handoff-manifest <HANDOFF>/unified_handoff_manifest.json \
    --runtime-attestation <HANDOFF>/runtime_attestation.json \
    --execute
```

查看：

```bash
/root/workspace/pipeline/status.sh
screen -ls
nvidia-smi
python3 -m json.tool <HANDOFF>/execution/execution_report.json
```

完成只能在以下条件同时满足时声明：

- execution report 为 `succeeded`；
- 计划的 canonical outputs 存在且可解析；
- status 为 `completed`；
- screen/model process 正常结束；
- GPU 显存已释放；
- 日志无未处理 traceback。

**图**：

```bash
python3 scripts/pipeline_visualization/11_summary.py --config config/pipeline_visualization.json
```

## 6. 常见错误处理

### Dry-run 成功，`--execute` 失败

Dry-run 只证明合同和路径正确。真实运行还可能因 checkpoint、CUDA、输出 schema、模型内部编号或运行时下载失败。查看 execution report 的首个失败 command 及其 stdout/stderr 日志。

### status 显示 running，但 screen/GPU 不存在

这是 stale status。先读 execution report，不要自动重启；根据报告结果将 status 收口为 failed 或 completed。

### Germinal 没有 accepted candidate

先区分：

- 技术失败：异常、缺资产、编号错误或指标为 `None`；
- 科学拒绝：所有指标存在，但未达预定阈值。

只修复第一类，不要将第二类当作 bug 放宽门槛。

### 修复后是否可以覆盖旧输出

不可以。源码、输入、参数、权重或坐标映射改变时创建新 attempt，并保留旧 execution report 和已完成输出。

## 7. 每次运行的学生记录模板

```text
Run ID:
Attempt:
Wave:
Scientific objective:
Template / epitope:
Source commit(s):
Environment path(s):
Checkpoint manifest SHA-256:
Handoff ID / SHA-256:
Runtime attestation SHA-256:
Dry-run SHA-256:
Screen session:
Start / finish time:
Execution report SHA-256:
Canonical outputs and SHA-256:
Technical outcome:
Scientific filter outcome:
Automatic candidate promotion: false
Reviewer:
```

## 8. 学生最终交付清单

- 已填写和复核的 target manifest；
- handoff manifest 与 runtime attestation；
- dry-run preview；
- execution report 及命令日志；
- canonical PDB/FASTA/CSV 与 SHA-256；
- `visualizations/pipeline/` 中 11 张阶段图；
- 明确分开“技术成功”、“科学过滤通过”和“实验验证”；
- 不得将未公开序列或结构上传到公共服务。
