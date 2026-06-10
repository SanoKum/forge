#!/usr/bin/env bash
# 08.bump 回帰検証ランナー。
#   低M/高M × 陽解法/陰解法 を新しい run_* で流し、y≈0.3 プロファイルと収束度合いを
#   committed baseline と比較して PASS/FAIL を出す。
#
# usage: run_verification.sh [quick|full]
#   quick (既定): 陰解法のみ (loM imp + hiM imp、~70s)。収束・解の主要回帰を素早く確認。
#   full        : 陽解法+陰解法の4本 (loM exp は ~4分かかる)。
#
# 前提: solver_density_cuda を native ビルド済み (.build-native/relwithdebinfo/forge)。
#   stale チェックは AGENTS.md「実行前のバイナリ鮮度チェック」に従い各自で行うこと。
set -uo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUMP="$(cd -- "$HERE/.." && pwd)"
FORGE="${FORGE_BIN:-/home/sano/work/forge/solver_density_cuda/.build-native/relwithdebinfo/forge}"
MODE="${1:-quick}"
WORK="$BUMP/verify/_work"

if [[ ! -x "$FORGE" ]]; then echo "ERROR: forge not found: $FORGE"; exit 2; fi
rm -rf "$WORK"; mkdir -p "$WORK"

# matrix: case method cfl nOut Pt mesh  (full は4本、quick は imp 2本)
if [[ "$MODE" == "full" ]]; then
  JOBS=(
    "loM exp 1.0 40000 120193 bump.h5"
    "loM imp 10  8000  120193 bump.h5"
    "hiM exp 1.0 8000  100000 bump_4pct.h5"
    "hiM imp 10  3000  100000 bump_4pct.h5"
  )
else
  JOBS=(
    "loM imp 10  8000  120193 bump.h5"
    "hiM imp 10  3000  100000 bump_4pct.h5"
  )
fi

for j in "${JOBS[@]}"; do
  read CASE METHOD CFL NOUT PT MESH <<< "$j"
  RUN="$WORK/run_${CASE}_${METHOD}"
  echo ">>> $CASE $METHOD cfl=$CFL nOut=$NOUT -> $RUN"
  bash "$HERE/gen_config.sh" "$RUN" "$CASE" "$METHOD" "$CFL" "$NOUT"
  ( cd "$RUN" && "$FORGE" > forge_run.log 2>&1 ); rc=$?
  echo "    forge exit=$rc"
  # 最終 res_*.h5
  LAST=$(ls "$RUN"/res_*.h5 2>/dev/null | awk -F'[_.]' '{print $(NF-1), $0}' | sort -n | tail -1 | cut -d' ' -f2)
  python3 "$HERE/extract_profile.py" --mesh "$RUN/$MESH" --result "$LAST" --pt "$PT" \
      --output "$RUN/profile_y03.csv" --y 0.3
  python3 /home/sano/work/forge/solver_density_cuda/tools/plot_implicit_residuals.py \
      --input "$RUN/residual_history.csv" --output "$RUN/residual_history.png" \
      --equations ro,roUx,roUy,roUz,roe --phases outer_begin >/dev/null 2>&1 || true
done

echo "=================================================================="
python3 "$HERE/compare.py" --workdir "$WORK"
