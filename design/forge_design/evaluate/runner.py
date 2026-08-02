"""バッチ評価: 問題定義 YAML + dv → run ディレクトリ準備 → forge 実行 → 判定・抽出。

パイプライン (Phase 0, thruster_bell):
  1. NozzleWall 構築 + validate (幾何フィルタ)
  2. 構造化 quad メッシュ → nozzle.msh (msh4.1 直書き)
  3. solverConfig/bcondConfig/probe を生成
  4. convertGmshToForge → nozzle.h5 (wall_dist・幾何量は変換器が計算)
  5. check_mesh_quality.py ゲート (FAIL → 評価失敗)
  6. 準 1D 等エントロピー IC 貼り付け
  7. run_case.sh (forge 実行 + 収束判定 + 残差 png)  ※hook 準拠の起動経路
  8. CONVERGENCE_VERDICT / metrics.json
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..geometry.wall import NozzleWall, WallParams
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from .ic import paste_isentropic_ic

FORGE_TOOLS = Path("/home/sano/work/forge/solver_density_cuda/tools")
FORGE_BUILD = Path("/home/sano/work/forge/solver_density_cuda/build")
_ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")


def build_wall(p: Problem) -> tuple:
    g = p.geometry
    theta_a = np.deg2rad(dv_value(p, "theta_a_deg"))
    L_div = dv_value(p, "L_over_rt")
    wp = WallParams(
        r_inlet=float(g.get("r_inlet", 2.5)),
        L_pipe=float(g.get("L_pipe", 0.5)),
        L_contract=float(g.get("L_contract", 3.0)),
        Ru=float(g.get("Ru", 1.5)),
        Rd=float(g.get("Rd", 0.4)),
        phi_u=np.deg2rad(float(g.get("phi_u_deg", 25.0))),
        theta_a=theta_a,
        L_div=L_div,
        r_exit=float(np.sqrt(p.spec["eps"])),
        theta_e=np.deg2rad(float(g.get("theta_e_deg", 8.0))),
    )
    wall = NozzleWall(wp)
    msgs = wall.validate()
    if msgs:
        raise ValueError("幾何フィルタ不合格: " + "; ".join(msgs))
    return wall, wp


def _solver_config(p: Problem, nsteps: int, out_interval: int) -> str:
    sst = p.evaluate.get("turbulence", "sst") == "sst"
    turb = ('turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2}'
            if sst else 'turbulence: {model: "none"}')
    return f"""mesh: {{meshFormat: "hdf5", meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, isAxisymmetric: 1, thermalMethod: 0, viscMethod: 1,
           ro: 1.2, visc: 1.8e-5, thermCond: 0.0257, cp: {p.cp}, gamma: {p.gamma}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 4.0,
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, lowMachPrecond: 0}}
  outStepStart: 0
  outStepInterval: {out_interval}
  timeIntegration: 11
  nStepInner: 20
space: {{convMethod: 1, limiter: 2}}
{turb}
initial: "uniform_p101325_u10"
"""


def _bcond_config(p: Problem) -> str:
    Pt = float(p.spec["Pt"])
    Tt = float(p.spec["Tt"])
    pa = float(p.spec["p_ambient"])
    return f"""inlet:  {{physID: 1, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {{Pt: {Pt}, Tt: {Tt}, k: 1.0, omega: 18000.0}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {{Ps: {pa}, Pt: {pa}, Tt: 300.0}}}}
wall:   {{physID: 3, kind: wall,             outputHDFflg: 1, ints: , floats: }}
axis:   {{physID: 4, kind: axis,             outputHDFflg: 0, ints: , floats: }}
"""


PROBE_STUB = "outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n"


def prepare(problem_path, run_dir, nsteps=None) -> dict:
    """メッシュ・config・IC まで準備する (forge は起動しない)。"""
    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    wall, wp = build_wall(p)
    scale = float(p.spec["r_throat"])

    mp = Mesh2DParams(
        ni=int(p.mesh.get("ni", 241)),
        nj=int(p.mesh.get("nj", 81)),
        wall_first_frac=float(p.mesh.get("wall_first_frac", 2.0e-3)),
        throat_refine=float(p.mesh.get("throat_refine", 3.0)),
        scale=scale,
    )
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    # 壁座標も記録 (IC・後処理・帰還エンジンの入力)
    xs = np.linspace(wall.x_in, wall.x_e, 800)
    np.savetxt(run_dir / "contour.csv", np.c_[xs * scale, wall.r(xs) * scale],
               delimiter=",", header="x_m,r_m", comments="")

    n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
    (run_dir / "solverConfig.yaml").write_text(_solver_config(p, n, int(p.evaluate.get("outStepInterval", max(n // 3, 1)))))
    (run_dir / "bcondConfig.yaml").write_text(_bcond_config(p))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)

    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"), "nozzle.h5"],
                       cwd=run_dir, env=_ENV, capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(q.stdout + q.stderr)
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")

    exit1d = paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                                 float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp)
    info = {"problem": str(problem_path), "run_dir": str(run_dir), "nsteps": n,
            "scale_m": scale, "exit_1d": exit1d,
            "mesh": {"ni": mp.ni, "nj": mp.nj, "cells": int(quads.shape[0])}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_forge(run_dir) -> int:
    """run_case.sh 経由で forge を実行 (収束判定・残差 png も生成される)。"""
    r = subprocess.run([str(FORGE_TOOLS / "run_case.sh"), str(Path(run_dir).resolve())],
                       env=_ENV, capture_output=True, text=True)
    (Path(run_dir) / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode


def collect(problem_path, run_dir) -> dict:
    """収束 VERDICT + 最終 res からメトリクスを抽出し metrics.json に保存。"""
    from ..metrics.extract import thrust_metrics

    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    verdict = (run_dir / "CONVERGENCE_VERDICT.txt").read_text() if (run_dir / "CONVERGENCE_VERDICT.txt").exists() else "(no verdict)"
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    out = {"convergence_verdict": verdict.strip().splitlines()[-2:] if verdict else []}
    if res:
        out["res_file"] = res[-1].name
        out.update(thrust_metrics(run_dir / "nozzle.h5", res[-1], p.gamma,
                                  float(p.spec["Pt"]), float(p.spec["p_ambient"]),
                                  float(p.spec["r_throat"])))
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None):
    import argparse

    ap = argparse.ArgumentParser(description="forge_design バッチ評価 (Phase 0)")
    ap.add_argument("problem")
    ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    a = ap.parse_args(argv)
    info = prepare(a.problem, a.run_dir, a.steps)
    print(json.dumps(info, indent=1))
    if a.prepare_only:
        return 0
    rc = run_forge(a.run_dir)
    out = collect(a.problem, a.run_dir)
    print(json.dumps(out, indent=1))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
