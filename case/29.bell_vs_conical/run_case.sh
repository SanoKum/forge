#!/bin/bash
# case 29: 1 ノズルを convert + 等エントロピー IC 貼付 + forge 実行する。
# 使い方: ./run_case.sh <run_dir> <nozzle:bell|conical>
set -e
RUN=$1; NOZ=$2
FORGE=/home/sano/work/forge/solver_density_cuda/build-verify
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE/$RUN"
cp "$HERE/mesh/$NOZ.msh" nozzle.msh
"$FORGE/convertGmshToForge" nozzle.msh nozzle.h5 > convert.log 2>&1
python3 "$HERE/mesh/set_isentropic_ic.py" nozzle.h5 "$HERE/mesh/contour_$NOZ.csv" \
    --pt 4.0e6 --tt 1500 --rt 10.0
"$FORGE/forge" > forge_run.log 2>&1
echo "[$RUN] forge done (exit $?)"
