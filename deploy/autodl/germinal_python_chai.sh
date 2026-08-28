#!/usr/bin/env bash
set -euo pipefail

germinal_python="${NFL_GERMINAL_PYTHON:-/root/miniconda3/envs/germinal/bin/python}"
chai_downloads_dir="${CHAI_DOWNLOADS_DIR:-/autodl-fs/data/shared/models/chai1}"

required_assets=(
  "conformers_v1.apkl"
  "models_v2/feature_embedding.pt"
  "models_v2/bond_loss_input_proj.pt"
  "models_v2/token_embedder.pt"
  "models_v2/trunk.pt"
  "models_v2/diffusion_module.pt"
  "models_v2/confidence_head.pt"
  "esm/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt"
)

if [[ ! -x "$germinal_python" ]]; then
  echo "ERROR: Germinal Python is not executable: $germinal_python" >&2
  exit 2
fi
for relative in "${required_assets[@]}"; do
  if [[ ! -s "$chai_downloads_dir/$relative" ]]; then
    echo "ERROR: Chai runtime asset is missing: $chai_downloads_dir/$relative" >&2
    exit 2
  fi
done

export CHAI_DOWNLOADS_DIR="$chai_downloads_dir"
exec "$germinal_python" "$@"
