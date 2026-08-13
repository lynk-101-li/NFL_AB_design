#!/usr/bin/env bash
set -euo pipefail

mode="${1:-all}"
case "$mode" in
  all|rfantibody|iggm|germinal) ;;
  *) printf 'Usage: %s {all|rfantibody|iggm|germinal}\n' "$0" >&2; exit 64 ;;
esac

apps_root="${NFL_APPS_ROOT:-/root/apps}"
conda_exe="${NFL_CONDA_EXE:-}"
if [[ -z "$conda_exe" ]]; then
  conda_exe="$(command -v conda || true)"
fi

rf_revision="8fe311415754e0276d1a39c87c57e69c88927a2d"
iggm_revision="06abc563b3fc8c7ea020543add16b69b6f8a1c8d"
germinal_revision="1e1c1a5b79884ae45abae030c9df90d9423a990a"

require_command() {
  command -v "$1" >/dev/null 2>&1 || {
    printf 'ERROR: required command not found: %s\n' "$1" >&2
    exit 2
  }
}

clone_pinned() {
  local url="$1" destination="$2" revision="$3"
  if [[ -e "$destination" && ! -d "$destination/.git" ]]; then
    printf 'ERROR: existing path is not a git checkout: %s\n' "$destination" >&2
    exit 2
  fi
  if [[ ! -d "$destination/.git" ]]; then
    git clone "$url" "$destination"
  fi
  git -C "$destination" fetch --tags origin
  git -C "$destination" checkout --detach "$revision"
  [[ "$(git -C "$destination" rev-parse HEAD)" == "$revision" ]] || {
    printf 'ERROR: revision verification failed: %s\n' "$destination" >&2
    exit 2
  }
}

install_rfantibody() {
  local destination="$apps_root/RFantibody"
  require_command git
  require_command uv
  clone_pinned https://github.com/RosettaCommons/RFantibody.git "$destination" "$rf_revision"
  (
    cd "$destination"
    bash include/download_weights.sh
    uv sync --frozen
    uv run rfdiffusion --help >/dev/null
    uv run proteinmpnn --help >/dev/null
    uv run rf2 --help >/dev/null
  )
  printf 'READY source_and_cli\tRFantibody\t%s\n' "$rf_revision"
}

install_iggm() {
  local destination="$apps_root/IgGM"
  require_command git
  [[ -n "$conda_exe" ]] || { printf 'ERROR: conda not found\n' >&2; exit 2; }
  clone_pinned https://github.com/TencentAI4S/IgGM.git "$destination" "$iggm_revision"
  if "$conda_exe" env list | awk '{print $1}' | grep -Fxq iggm; then
    printf 'INFO: conda environment iggm already exists; refusing to mutate it automatically.\n'
  else
    "$conda_exe" env create -n iggm -f "$destination/environment.yaml"
    "$conda_exe" run -n iggm python -m pip install \
      pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
      -f https://data.pyg.org/whl/torch-2.0.1+cu117.html
  fi
  "$conda_exe" run -n iggm python "$destination/design.py" --help >/dev/null
  printf 'READY source_and_cli\tIgGM\t%s\n' "$iggm_revision"
  printf 'NOTICE: IgGM checkpoints may download on first model run; record their SHA-256 before attestation.\n'
}

install_germinal() {
  local destination="$apps_root/germinal"
  require_command git
  [[ -n "$conda_exe" ]] || { printf 'ERROR: conda not found\n' >&2; exit 2; }
  clone_pinned https://github.com/SantiagoMille/germinal.git "$destination" "$germinal_revision"
  if "$conda_exe" env list | awk '{print $1}' | grep -Fxq germinal; then
    printf 'INFO: conda environment germinal already exists; refusing to mutate it automatically.\n'
  else
    "$conda_exe" create -y -n germinal python=3.10 pip
    "$conda_exe" run -n germinal python -m pip install uv
    "$conda_exe" run -n germinal uv pip install \
      pandas matplotlib numpy biopython scipy seaborn tqdm ffmpeg py3dmol \
      chex dm-haiku dm-tree joblib ml-collections immutabledict optax cvxopt \
      mdtraj colabfold ipsae
    "$conda_exe" run -n germinal uv pip install -e "$destination/colabdesign"
    "$conda_exe" run -n germinal uv pip install \
      iglm 'torchvision==0.21.*' 'chai-lab==0.6.1' 'torch==2.6.*' \
      'torchaudio==2.6.*' 'torchtyping==0.1.5' 'torch_geometric==2.6.*'
    "$conda_exe" run -n germinal uv pip install -e "$destination"
    "$conda_exe" run -n germinal uv pip install \
      'jax==0.5.3' 'dm-haiku==0.0.13' hydra-core omegaconf
    "$conda_exe" run -n germinal uv pip install \
      'jax[cuda12_pip]==0.5.3' \
      -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    "$conda_exe" run -n germinal uv pip install ablang2 --no-deps
    "$conda_exe" run -n germinal uv pip install rotary_embedding_torch --no-deps
  fi
  printf 'PARTIAL open_dependencies_only\tGerminal\t%s\n' "$germinal_revision"
  printf 'ACTION REQUIRED: obtain PyRosetta under its license and place AlphaFold-Multimer parameters in a reviewed directory.\n'
  printf 'ACTION REQUIRED: run conda run -n germinal python %s/validate_install.py after restricted assets are installed.\n' "$destination"
}

mkdir -p "$apps_root"
if [[ "$mode" == all || "$mode" == rfantibody ]]; then install_rfantibody; fi
if [[ "$mode" == all || "$mode" == iggm ]]; then install_iggm; fi
if [[ "$mode" == all || "$mode" == germinal ]]; then install_germinal; fi
