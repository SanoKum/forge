#!/usr/bin/env bash
# soft 段 (1次, cfl 0.5, 3000 step) → 本段 (runner_axismach.run_staged と同手順)
set -e
T=/home/sano/work/forge-axisnode/solver_density_cuda/tools
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
cp solverConfig.yaml solverConfig_main.yaml
sed -e 's/cfl: 4.0, cfl_pseudo: 4.0/cfl: 0.5, cfl_pseudo: 0.5/' -e 's/convMethod: 1/convMethod: 0/' -e 's/nStepOuter: [0-9]*/nStepOuter: 3000/' -e 's/outStepInterval: [0-9]*/outStepInterval: 3000/' solverConfig_main.yaml > solverConfig.yaml
$T/run_case.sh . > soft_stdout.log 2>&1
python3 $T/interp_field.py res_3000.h5 nozzle.h5
mkdir -p soft && mv res_* residual_history.* CONVERGENCE_VERDICT.txt forge_run.log soft/ 2>/dev/null || true
cp solverConfig_main.yaml solverConfig.yaml
$T/run_case.sh . > run_case_stdout.log 2>&1
