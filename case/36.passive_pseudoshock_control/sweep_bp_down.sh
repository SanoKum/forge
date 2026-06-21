#!/usr/bin/env bash
# 背圧 down-sweep: Ps を段階的に下げ、各点を前段から継続 (cfl5/implicitRelax1/detectNaN1/80k)。
# run_case.sh 経由なので各 run に CONVERGENCE_VERDICT.txt が残る。発散 (res_nan) したら停止。
set -uo pipefail
ROOT=/home/sano/work/forge
CD=$ROOT/case/36.passive_pseudoshock_control
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
cd "$CD"

mkbcond(){ # $1=file $2=Ps[Pa]
cat > "$1" <<EOF
inlet:
  physID: 1
  kind: inlet_uniformVelocity
  outputHDFflg: 0
  ints:
  floats:
    Ux: 458.6
    Uy: 0.0
    Uz: 0.0
    ro: 11.733
    Ps: 617800.0
    k: 31.547
    omega: 3.0018e7
outlet:
  physID: 2
  kind: outlet_statPress
  outputHDFflg: 0
  ints:
  floats:
    Ps: $2
    Pt: $2
    Tt: 288.15
wall_slip:
  physID: 3
  kind: slip
  outputHDFflg: 0
  ints:
  floats:
wall:
  physID: 4
  kind: wall
  outputHDFflg: 0
  ints:
  floats:
EOF
}

prev="run_0027_sst_porous_bp1p9_cfl5_cont"
tags=(run_0028_sst_porous_bp1p8 run_0029_sst_porous_bp1p7 run_0030_sst_porous_bp1p6)
pss=(1800000.0 1700000.0 1600000.0)
for i in 0 1 2; do
  R=${tags[$i]}; ps=${pss[$i]}
  mkdir -p "$R"
  cp "$prev/solverConfig.yaml" "$prev/probe.yaml" "$prev/passive_porous.h5" "$R/"
  python3 restart_field.py "$prev/res_80000.h5" "$R/passive_porous.h5" >/dev/null
  mkbcond "$R/bcondConfig.yaml" "$ps"
  echo "=== sweep $R  Ps=$ps  (continue from $prev) ==="
  bash "$ROOT/solver_density_cuda/tools/run_case.sh" "$CD/$R"
  if ls "$R"/res_nan_*.h5 >/dev/null 2>&1; then echo "DIVERGED at $R (res_nan) — stop sweep"; break; fi
  prev="$R"
done
echo "SWEEP DONE"
