#!/usr/bin/env bash
set -euo pipefail

mode="${1:-help}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
conda_exe="${NFL_CONDA_EXE:-$(command -v conda || true)}"

usage() {
  cat <<'EOF'
Usage: bash deploy/autodl/bootstrap_models.sh MODE

MODE:
  sources        initialize and verify all pinned source submodules
  design-open    install RFantibody and IgGM
  structure-open install ImmuneBuilder, Chai-1 and Boltz-2
  rfantibody | iggm | germinal | tfold | igfold
  immunebuilder | alphafold3 | chai1 | boltz2

Restricted components require an explicit acknowledgement environment variable:
  germinal:   NFL_ACK_PYROSETTA_AND_DEPENDENCY_TERMS=1
  tfold:      NFL_ACK_TFOLD_NONCOMMERCIAL_TERMS=1
  igfold:     NFL_ACK_IGFOLD_ACADEMIC_TERMS=1
  alphafold3: NFL_ACK_AF3_MODEL_TERMS=1

An acknowledgement records that the student reviewed upstream terms. It does
not grant a license. Model weights, databases and PyRosetta are never committed.
EOF
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || {
    printf 'ERROR: required command not found: %s\n' "$1" >&2
    exit 2
  }
}

require_conda() {
  [[ -n "$conda_exe" ]] || { printf 'ERROR: conda not found\n' >&2; exit 2; }
}

require_ack() {
  local variable="$1" component="$2"
  [[ "${!variable:-}" == 1 ]] || {
    printf 'ERROR: review upstream terms for %s, then set %s=1 for this command.\n' "$component" "$variable" >&2
    exit 3
  }
}

ensure_sources() {
  require_command git
  git -C "$repo_root" submodule sync --recursive
  git -C "$repo_root" submodule update --init --recursive
  python3 "$repo_root/scripts/verify_model_components.py"
}

create_env() {
  local name="$1" python_version="$2"
  require_conda
  if "$conda_exe" env list | awk '{print $1}' | grep -Fxq "$name"; then
    printf 'INFO: environment %s already exists; refusing to mutate it automatically.\n' "$name"
    return 1
  fi
  "$conda_exe" create -y -n "$name" "python=$python_version" pip
}

install_rfantibody() {
  ensure_sources
  require_command uv
  (
    cd "$repo_root/third_party/RFantibody"
    bash include/download_weights.sh
    uv sync --frozen
    uv run rfdiffusion --help >/dev/null
    uv run proteinmpnn --help >/dev/null
    uv run rf2 --help >/dev/null
  )
  printf 'READY source_cli_weights\trfantibody\n'
}

install_iggm() {
  ensure_sources
  require_conda
  if create_env iggm 3.10; then
    "$conda_exe" env update -n iggm -f "$repo_root/third_party/IgGM/environment.yaml"
  fi
  "$conda_exe" run -n iggm python "$repo_root/third_party/IgGM/design.py" --help >/dev/null
  printf 'PARTIAL source_cli\tiggm\tcheckpoints download on first model run\n'
}

install_germinal() {
  require_ack NFL_ACK_PYROSETTA_AND_DEPENDENCY_TERMS Germinal
  ensure_sources
  require_conda
  if create_env germinal 3.10; then
    "$conda_exe" run -n germinal python -m pip install uv
    "$conda_exe" run -n germinal uv pip install \
      pandas matplotlib numpy biopython scipy seaborn tqdm ffmpeg py3dmol \
      chex dm-haiku dm-tree joblib ml-collections immutabledict optax cvxopt \
      mdtraj colabfold ipsae
    "$conda_exe" run -n germinal uv pip install -e "$repo_root/third_party/germinal/colabdesign"
    "$conda_exe" run -n germinal uv pip install \
      iglm 'torchvision==0.21.*' 'chai-lab==0.6.1' 'torch==2.6.*' \
      'torchaudio==2.6.*' 'torchtyping==0.1.5' 'torch_geometric==2.6.*'
    "$conda_exe" run -n germinal uv pip install -e "$repo_root/third_party/germinal"
    "$conda_exe" run -n germinal uv pip install \
      'jax==0.5.3' 'dm-haiku==0.0.13' hydra-core omegaconf
    "$conda_exe" run -n germinal uv pip install \
      'jax[cuda12_pip]==0.5.3' \
      -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    "$conda_exe" run -n germinal uv pip install ablang2 --no-deps
    "$conda_exe" run -n germinal uv pip install rotary_embedding_torch --no-deps
  fi
  printf 'PARTIAL source_environment\tgerminal\n'
  printf 'ACTION REQUIRED: legally obtain PyRosetta and required AlphaFold-Multimer parameters, then run validate_install.py.\n'
}

install_tfold() {
  require_ack NFL_ACK_TFOLD_NONCOMMERCIAL_TERMS tFold
  ensure_sources
  require_conda
  if create_env tfold 3.10; then
    "$conda_exe" env update -n tfold -f "$repo_root/third_party/tfold/environment.yaml"
  fi
  printf 'PARTIAL source_environment\ttfold\tweights not attested\n'
}

install_igfold() {
  require_ack NFL_ACK_IGFOLD_ACADEMIC_TERMS IgFold
  ensure_sources
  require_conda
  if create_env igfold 3.10; then
    "$conda_exe" run -n igfold python -m pip install -e "$repo_root/third_party/IgFold"
  fi
  "$conda_exe" run -n igfold python -c 'import igfold' >/dev/null
  printf 'PARTIAL source_environment\tigfold\tmodel assets not attested\n'
}

install_immunebuilder() {
  ensure_sources
  require_conda
  if create_env immunebuilder 3.10; then
    "$conda_exe" run -n immunebuilder python -m pip install -e "$repo_root/third_party/ImmuneBuilder"
  fi
  "$conda_exe" run -n immunebuilder python -c 'from ImmuneBuilder import ABodyBuilder2' >/dev/null
  printf 'PARTIAL source_environment\timmunebuilder\tmodel assets download on first run\n'
}

install_alphafold3() {
  require_ack NFL_ACK_AF3_MODEL_TERMS "AlphaFold 3"
  ensure_sources
  require_command docker
  printf 'SOURCE VERIFIED\talphafold3\n'
  printf 'ACTION REQUIRED: follow third_party/alphafold3/docs/installation.md to build the official container.\n'
  printf 'ACTION REQUIRED: obtain parameters and databases under their separate terms; keep them outside Git.\n'
}

install_chai1() {
  ensure_sources
  require_conda
  if create_env chai1 3.10; then
    "$conda_exe" run -n chai1 python -m pip install -e "$repo_root/third_party/chai-lab"
  fi
  "$conda_exe" run -n chai1 chai-lab --help >/dev/null
  printf 'PARTIAL source_cli\tchai1\tweights download on first run\n'
}

install_boltz2() {
  ensure_sources
  require_conda
  if create_env boltz2 3.10; then
    "$conda_exe" run -n boltz2 python -m pip install -e "$repo_root/third_party/boltz"
  fi
  "$conda_exe" run -n boltz2 boltz --help >/dev/null
  printf 'PARTIAL source_cli\tboltz2\tweights download on first run\n'
}

case "$mode" in
  help|-h|--help) usage ;;
  sources) ensure_sources ;;
  design-open) install_rfantibody; install_iggm ;;
  structure-open) install_immunebuilder; install_chai1; install_boltz2 ;;
  rfantibody) install_rfantibody ;;
  iggm) install_iggm ;;
  germinal) install_germinal ;;
  tfold) install_tfold ;;
  igfold) install_igfold ;;
  immunebuilder) install_immunebuilder ;;
  alphafold3) install_alphafold3 ;;
  chai1) install_chai1 ;;
  boltz2) install_boltz2 ;;
  *) usage >&2; exit 64 ;;
esac
