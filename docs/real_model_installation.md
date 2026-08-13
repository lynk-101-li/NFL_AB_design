# 模型源码、部署与运行教程

本项目把九个上游模型以 **Git submodule** 固定到精确 commit。这样学生 clone 后能取得同一份源码，同时不会把第三方代码复制进本项目历史。源码存在不等于依赖、权重或数据库已经安装，更不等于模型已经运行。

## 1. 克隆与源码校验

推荐一次性克隆：

```bash
git clone --recurse-submodules https://github.com/lynk-101-li/NFL_AB_design.git
cd NFL_AB_design
python3 scripts/verify_model_components.py
```

如果已经普通 clone：

```bash
git submodule sync --recursive
git submodule update --init --recursive
python3 scripts/verify_model_components.py
```

固定版本、许可证状态和用途见 `config/model_components.json`。更新模型必须同时更新 gitlink、注册表、测试和发行说明，不能只执行 `git pull`。

## 2. 目录与环境原则

- 九个源码目录位于 `third_party/`；不要直接修改。
- 每个模型使用独立 Conda/uv 环境，避免 PyTorch、JAX 和 CUDA 依赖互相覆盖。
- 权重、AlphaFold 数据库、PyRosetta 和缓存放在服务器文件盘，不提交到 Git。
- 真实运行放在 `real_runs/`；可再生的 proxy 输出放在 `outputs/`。
- 单卡服务器一次只执行一个 GPU job，先 canary，再扩量。

AutoDL 推荐布局：

```text
/root/workspace/NFL_AB_design/                    # 本仓库和 submodule
/root/miniconda3/envs/                            # 独立软件环境（系统盘）
/root/autodl-fs/shared/models/                     # 权重与数据库（文件盘）
/root/autodl-fs/topics/nfl-antibody-design/runs/  # 正式运行和结果
/root/autodl-tmp/nfl-antibody-design/             # 可删除 scratch
```

先检查资源：

```bash
bash deploy/autodl/preflight.sh
bash deploy/autodl/bootstrap_models.sh sources
```

## 3. 抗体从头设计模型

### RFantibody

```bash
bash deploy/autodl/bootstrap_models.sh rfantibody
cd third_party/RFantibody
uv run rfdiffusion --help
uv run proteinmpnn --help
uv run rf2 --help
```

真实输入不是普通 FASTA，而是 target PDB、每个模板各自的 Chothia H/L/T PDB、全长坐标到 PDB 残基映射和人工热点。项目通过 `scripts/prepare_real_model_jobs.py` 生成 `rfdiffusion → ProteinMPNN → RF2` 的 shell-free argv，不应手工拼接坐标。

### IgGM

```bash
bash deploy/autodl/bootstrap_models.sh iggm
conda run -n iggm python third_party/IgGM/design.py --help
```

输入 FASTA 必须包含 `>H`、`>L`、`>A`；H/L 的 `X` 是设计区，A 必须与 target PDB 抗原链的序列和顺序完全一致。权重首次运行可能下载，运行前后记录文件 SHA-256。

### Germinal

Germinal 作为独立 scFv 轨道，不代表原生 paired-Fv 几何。先阅读 PyRosetta、ColabDesign、结构模型和参数的独立条款：

```bash
export NFL_ACK_PYROSETTA_AND_DEPENDENCY_TERMS=1
bash deploy/autodl/bootstrap_models.sh germinal
conda run -n germinal python third_party/germinal/validate_install.py
```

环境变量只是“已阅读”的本地门禁，不授予许可证。PyRosetta 和 AlphaFold-Multimer 参数须由学生合法取得并放在 Git 之外。

### 统一 canary 作业

只有人工审核后的正式 `config/target_structure_manifest.json` 才能通过：

```bash
python3 scripts/prepare_real_model_jobs.py \
  --target-manifest config/target_structure_manifest.json \
  --iggm-repo-dir third_party/IgGM \
  --germinal-repo-dir third_party/germinal \
  --germinal-af-params-dir /root/autodl-fs/shared/models/alphafold_multimer \
  --profile smoke \
  --job-scope canary \
  --canary-template-id template_7-H11-D3-2-C7 \
  --canary-epitope-id C_boundary_368_377
```

该命令只编译和校验 handoff，不运行模型。随后按 `deploy/autodl/README.md` 生成 runtime attestation，并先 dry-run：

```bash
python3 scripts/execute_real_model_jobs.py \
  --handoff-manifest real_runs/handoffs/手工选择对应目录/unified_handoff_manifest.json \
  --runtime-attestation /绝对路径/runtime_attestation.json
```

确认命令、输入哈希、GPU和权重哈希后，才追加 `--execute`。

## 4. 抗体结构预测

### tFold

tFold 当前许可证是 PolyForm Noncommercial 1.0.0，先确认用途符合条款：

```bash
export NFL_ACK_TFOLD_NONCOMMERCIAL_TERMS=1
bash deploy/autodl/bootstrap_models.sh tfold
conda run -n tfold python third_party/tfold/projects/tfold_ab/predict.py \
  --fasta real_runs/structure/one_candidate_HL.fasta \
  --output real_runs/structure/tfold_ab/one_candidate.pdb
```

`one_candidate_HL.fasta` 必须只含该候选的 H、L 两条链；不要把含 12 个候选的汇总 FASTA 直接交给单任务入口。

抗体—抗原复合物走 tFold-Ag，除复合物 FASTA 外还需要按上游格式准备 MSA：

```bash
conda run -n tfold python third_party/tfold/projects/tfold_ag/predict.py \
  --fasta real_runs/structure/one_complex.fasta \
  --msa real_runs/structure/one_complex.a3m \
  --output real_runs/structure/tfold_ag
```

### IgFold

IgFold 使用 JHU Academic Software License。确认符合学术/非商业条款后：

```bash
export NFL_ACK_IGFOLD_ACADEMIC_TERMS=1
bash deploy/autodl/bootstrap_models.sh igfold
```

单个候选的最小 Python 入口：

```python
from igfold import IgFoldRunner
runner = IgFoldRunner()
runner.fold(
    "real_runs/structure/igfold/candidate.pdb",
    sequences={"H": "VH_SEQUENCE", "L": "VL_SEQUENCE"},
    do_refine=False,
    do_renum=False,
)
```

### ImmuneBuilder / ABodyBuilder2

仓库采用公开的 ABodyBuilder2，不再使用不准确的 “ABodyBuilder3” 名称：

```bash
bash deploy/autodl/bootstrap_models.sh immunebuilder
```

```python
from ImmuneBuilder import ABodyBuilder2
predictor = ABodyBuilder2()
model = predictor.predict({"H": "VH_SEQUENCE", "L": "VL_SEQUENCE"})
model.save("real_runs/structure/immunebuilder/candidate.pdb")
```

IgFold 和 ImmuneBuilder 主要用于 Fv 折叠、CDR 构象和 VH/VL packing 的交叉检查，不提供抗原结合证据。

## 5. 抗体—抗原复合物预测

### AlphaFold 3

代码是 Apache-2.0，但模型参数、输出及数据库有独立条款，参数不得随本仓库分发。遵循上游容器安装说明：

```bash
export NFL_ACK_AF3_MODEL_TERMS=1
bash deploy/autodl/bootstrap_models.sh alphafold3
```

仓库已经为选中候选导出 AF3 JSON 到 `outputs/exports/af3_inputs/`。官方容器命令示例：

```bash
docker run --rm --gpus all \
  -v /root/autodl-fs/shared/models/alphafold3:/root/af_input:ro \
  -v "$PWD/outputs/exports/af3_json":/root/af_json:ro \
  -v /root/autodl-fs/topics/nfl-antibody-design/runs/af3:/root/af_output \
  alphafold3 \
  python run_alphafold.py \
  --json_path=/root/af_json/CANDIDATE.json \
  --model_dir=/root/af_input/models \
  --db_dir=/root/af_input/public_databases \
  --output_dir=/root/af_output
```

镜像名和挂载结构必须与学生实际构建的官方镜像一致。

### Chai-1

```bash
bash deploy/autodl/bootstrap_models.sh chai1
conda run -n chai1 chai-lab fold \
  outputs/exports/fasta/complex_CANDIDATE.fasta \
  real_runs/structure/chai1/CANDIDATE
```

### Boltz-2

```bash
bash deploy/autodl/bootstrap_models.sh boltz2
conda run -n boltz2 boltz predict real_runs/structure/boltz2/CANDIDATE.yaml \
  --use_msa_server \
  --out_dir real_runs/structure/boltz2/CANDIDATE
```

Boltz 输入使用其 YAML schema，不应把多链 FASTA 路径直接当作 YAML。按上游示例把 VH、VL、抗原分别声明为 protein entity，并保留各自链 ID。

## 6. 结果判读与留痕

建议每个候选至少保存：

- 输入序列/PDB、命令 argv、环境锁、GPU、源码 commit；
- 每个 checkpoint/database 的路径与 SHA-256；
- 原始 PDB/mmCIF、模型置信度和未加工日志；
- Fv 内部 CDR/VH-VL packing 指标；
- 复合物的 interface PAE、ipTM/pTM、界面接触与 clash；
- 多模型是否在表位、方向和 CDR 接触上收敛。

任何真实模型失败都必须保留失败日志。proxy 分数、模型预测和实验数据是三类不同证据，不可互相改名。

## 7. 许可证边界

本仓库 MIT 许可证仅覆盖本项目代码。Git submodule 仍由各自上游许可证管理；权重、数据库和参数还可能有独立条款。尤其是 tFold、IgFold、AlphaFold 3、PyRosetta，在使用前必须阅读固定源码内的许可证和上游当前条款。
