#!/usr/bin/env bash
set -euo pipefail

context_path="${NFL_AUTODL_CONTEXT:-/root/workspace/AUTODL_CONTEXT.md}"
status_path="${NFL_STATUS_PATH:-}"
if [[ -z "$status_path" && -f "$context_path" ]]; then
  formal_run_root="$(
    awk -F': ' '$1 == "formal_run_root" {print $2; exit}' "$context_path"
  )"
  if [[ -n "$formal_run_root" ]]; then
    status_path="$formal_run_root/status.json"
  fi
fi

if [[ -z "$status_path" || ! -f "$status_path" ]]; then
  echo "active_run=none"
  echo "assessment=failed"
  echo "reason=status_file_missing:$status_path"
  exit 2
fi

readarray -t status_fields < <(
  python3 - "$status_path" <<'PY'
import json
import sys

value = json.load(open(sys.argv[1], encoding="utf-8"))
for key in (
    "run_id",
    "attempt",
    "stage",
    "state",
    "updated_at",
    "heartbeat_at",
    "screen_session",
    "planned_outputs",
    "completed_outputs",
):
    raw = value.get(key)
    print("none" if raw is None else raw)
PY
)

run_id="${status_fields[0]}"
attempt="${status_fields[1]}"
stage="${status_fields[2]}"
state="${status_fields[3]}"
updated_at="${status_fields[4]}"
heartbeat_at="${status_fields[5]}"
screen_session="${status_fields[6]}"
planned_outputs="${status_fields[7]}"
completed_outputs="${status_fields[8]}"

screen_state="none"
if [[ "$screen_session" != "none" ]]; then
  if screen -ls 2>/dev/null | grep -Fq ".${screen_session}"; then
    screen_state="alive"
  else
    screen_state="missing"
  fi
fi

gpu_process="$(
  nvidia-smi \
    --query-compute-apps=pid,process_name,used_memory \
    --format=csv,noheader,nounits 2>/dev/null | head -n 1 || true
)"
gpu_process="${gpu_process:-none}"
system_disk_available="$(df -h --output=avail / | tail -n 1 | tr -d ' ')"
persistent_disk_available="$(
  df -h --output=avail /autodl-fs/data | tail -n 1 | tr -d ' '
)"

assessment="normal"
if [[ "$state" == "failed" ]]; then
  assessment="failed"
elif [[ "$state" == "running" && "$screen_state" != "alive" ]]; then
  assessment="stale"
elif [[ "$state" == "running" && "$gpu_process" == "none" ]]; then
  assessment="stale"
fi

echo "active_run=$run_id"
echo "attempt=$attempt"
echo "stage=$stage"
echo "state=$state"
echo "status_updated_at=$updated_at"
echo "heartbeat_at=$heartbeat_at"
echo "screen=$screen_state"
echo "gpu_process=$gpu_process"
echo "outputs=$completed_outputs/$planned_outputs"
echo "system_disk_available=$system_disk_available"
echo "persistent_disk_available=$persistent_disk_available"
echo "assessment=$assessment"
