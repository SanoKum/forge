#!/usr/bin/env bash
# run_case.sh と同じ (実行後に check_convergence を強制) だが、forge バイナリを引数で指定する (別ビルドの検証用)。
#   run_case_bin.sh <forge_binary> <run_dir>
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
BIN="$1"; RUNDIR="$2"
cd "$RUNDIR" || { echo "run_case_bin: bad run dir $RUNDIR"; exit 1; }
"$BIN" > forge_run.log 2>&1; rc=$?
python3 "$ROOT/solver_density_cuda/tools/check_convergence.py" . > CONVERGENCE_VERDICT.txt 2>&1 || true
python3 "$ROOT/solver_density_cuda/tools/plot_implicit_residuals.py" --input residual_history.csv --output residual_history.png > /dev/null 2>&1 || true
exit $rc
