#!/usr/bin/env bash
set -euo pipefail

# Read-only student-server preflight. It performs no install or directory write.
required_free_gb="${NFL_REQUIRED_FREE_SYSTEM_GB:-35}"
required_file_free_gb="${NFL_REQUIRED_FREE_FILE_GB:-50}"
file_root="${NFL_FILE_STORAGE_ROOT:-/root/autodl-fs}"

for command_name in git curl tar sha256sum; do
  if ! command -v "$command_name" >/dev/null 2>&1; then
    printf 'missing_command\t%s\n' "$command_name"
    printf 'state\tblocked_missing_prerequisites\n'
    exit 2
  fi
done

if ! command -v nvidia-smi >/dev/null 2>&1; then
  printf 'state\tblocked_nvidia_smi_unavailable\n'
  exit 2
fi

free_kb="$(df -Pk / | awk 'NR==2 {print $4}')"
free_gb="$((free_kb / 1024 / 1024))"

printf 'gpu\t'
nvidia-smi --query-gpu=name,memory.total,memory.used,utilization.gpu --format=csv,noheader
printf 'system_disk_free_gb\t%s\n' "$free_gb"
printf 'required_free_system_gb\t%s\n' "$required_free_gb"
if [[ ! -d "$file_root" ]]; then
  printf 'state\tblocked_file_storage_missing\n'
  printf 'file_storage_root\t%s\n' "$file_root"
  exit 2
fi
file_free_kb="$(df -Pk "$file_root" | awk 'NR==2 {print $4}')"
file_free_gb="$((file_free_kb / 1024 / 1024))"
printf 'file_storage_root\t%s\n' "$file_root"
printf 'file_storage_free_gb\t%s\n' "$file_free_gb"
printf 'required_free_file_gb\t%s\n' "$required_file_free_gb"

if (( free_gb < required_free_gb )); then
  printf 'state\tblocked_insufficient_system_disk\n'
  exit 2
fi

if (( file_free_gb < required_file_free_gb )); then
  printf 'state\tblocked_insufficient_file_storage\n'
  exit 2
fi

printf 'state\tready_for_environment_bootstrap\n'
