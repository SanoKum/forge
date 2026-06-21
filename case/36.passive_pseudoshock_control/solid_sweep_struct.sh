#!/usr/bin/env bash
# 単純 wall (キャビティ・多孔板なし) 構造メッシュで背圧 sweep: Ps 1.80..1.90 を 0.02 刻み (6 点)。
# fine-sweep 完了 & GPU 空きを待って開始。初点は porous Ps1.8 (run_0034) から cross-mesh interp、
# 以降は同一 solid メッシュで restart_field。porous との比較基準。
set -uo pipefail
ROOT=/home/sano/work/forge
CD=$ROOT/case/36.passive_pseudoshock_control
T=$ROOT/solver_density_cuda/tools
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
cd "$CD"

until grep -q "FINE-SWEEP DONE" fine_sweep_struct.out 2>/dev/null; do sleep 60; done
until [ -z "$(pgrep -x forge)" ]; do sleep 30; done
echo "[solid_sweep] fine-sweep done, GPU free -> start solid sweep"

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

tags=(run_0044_solid_struct_bp1p80 run_0045_solid_struct_bp1p82 run_0046_solid_struct_bp1p84 \
      run_0047_solid_struct_bp1p86 run_0048_solid_struct_bp1p88 run_0049_solid_struct_bp1p90)
pss=(1800000.0 1820000.0 1840000.0 1860000.0 1880000.0 1900000.0)
prev=""
for i in 0 1 2 3 4 5; do
  R=${tags[$i]}; ps=${pss[$i]}
  if [ "$i" -eq 0 ]; then
    # run_0044 は既に mesh+config+bcond あり。初期値は porous Ps1.8 から interp。
    python3 "$T/interp_field.py" run_0034_struct_porous_bp1p8/res_80000.h5 "$R/passive_solid.h5"
  else
    mkdir -p "$R"
    cp run_0044_solid_struct_bp1p80/solverConfig.yaml run_0044_solid_struct_bp1p80/probe.yaml \
       run_0044_solid_struct_bp1p80/passive_solid.h5 "$R/"
    python3 restart_field.py "$prev/res_80000.h5" "$R/passive_solid.h5" >/dev/null
    mkbcond "$R/bcondConfig.yaml" "$ps"
  fi
  echo "=== solid-sweep $R Ps=$ps ==="
  bash "$T/run_case.sh" "$CD/$R"
  if ls "$R"/res_nan_*.h5 >/dev/null 2>&1; then echo "DIVERGED at $R — stop"; break; fi
  prev="$R"
done
echo "SOLID-SWEEP DONE"
