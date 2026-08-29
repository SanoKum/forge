#!/usr/bin/env bash
# run ディレクトリを S3 へ回収する。インスタンス stop/terminate 前に実行する。
#   usage: sync_results_s3.sh <run_dir> [<run_dir> ...]
#   env:   FORGE_S3_BUCKET (default: s3://forge-runs)
# S3 側のキーはリポジトリルートからの相対パスにする (case/08.bump/run_0001_... の形で追える)。
set -euo pipefail

bucket="${FORGE_S3_BUCKET:-s3://forge-runs}"

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <run_dir> [<run_dir> ...]" >&2
  exit 1
fi

for run_dir in "$@"; do
  run_dir="$(cd "$run_dir" && pwd)"
  repo_root="$(git -C "$run_dir" rev-parse --show-toplevel)"
  rel="${run_dir#"$repo_root"/}"
  echo "== sync $rel -> $bucket/$rel =="
  aws s3 sync "$run_dir" "$bucket/$rel" --no-progress
done
echo "DONE"
