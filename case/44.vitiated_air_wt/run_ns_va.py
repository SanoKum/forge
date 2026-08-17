"""推奨点 R2/L_U6/L_c8 の NS (A12/A13) 投入ドライバ: coarse 中継 → y+~1 本計算 (v1 相関 δ*) → [v2 CFD 抽出 δ*]。
使い方: cd design && ./.venv-opt/bin/python ../case/44.vitiated_air_wt/run_ns_va.py <stage> <run_dir> [--ic-from DIR] [--dstar-csv CSV] [--prepare-only]
  stage: coarse | main
"""
import argparse, json, sys, subprocess
from pathlib import Path
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
CASE = Path("/home/sano/work/forge/case/44.vitiated_air_wt")
ap = argparse.ArgumentParser(); ap.add_argument("stage", choices=("coarse", "main")); ap.add_argument("run_dir")
ap.add_argument("--ic-from"); ap.add_argument("--dstar-csv"); ap.add_argument("--dstar-blend", default="6,9", help="相関→CSV の smoothstep 区間 [r_t] (CSV 全域採用なら -1,-0.5)")
ap.add_argument("--prepare-only", action="store_true")
a = ap.parse_args()
prob = CASE / ("problem_va_R2_LU6_Lc8_ns_coarse.yaml" if a.stage == "coarse" else "problem_va_R2_LU6_Lc8_ns.yaml")
blend = tuple(float(v) for v in a.dstar_blend.split(","))
info = prepare_ns(prob, a.run_dir, ic_from=a.ic_from, dstar_csv=a.dstar_csv, dstar_blend=blend)
print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "mesh", "nStepOuter", "cfl_main")}, indent=1), flush=True)
print(Path(a.run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
if a.prepare_only: sys.exit(0)
rc = run_staged_ns(a.run_dir)
print("forge rc", rc, flush=True)
m = collect(prob, a.run_dir); print(json.dumps(m, indent=1), flush=True)
subprocess.run([sys.executable, "/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py", a.run_dir], check=False)
