#!/usr/bin/env python3
"""W0 判別実験 (R=5 A/B) の cold node run 準備・起動。

walldriven plan W0 / 調査 §8。run_0034_node_base (R=2, node, cold 単段 24000 step,
cfl4, convMethod 1) と同一 CFD 条件で、R=5 の Sauer アンカー pass-1 仮壁を評価する。
`runner_wt.prepare` は euler=True で cell config を書くため、node 版の準備をここで
明示的に行う (run_0034 と同じ構成)。

使い方: design/.venv-opt/bin/python run_r5_cold.py [--prepare-only]
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE.parents[1] / "design"))

from forge_design.evaluate.ic import paste_isentropic_ic  # noqa: E402
from forge_design.evaluate.runner import (  # noqa: E402
    FORGE_BUILD, FORGE_TOOLS, PROBE_STUB, _ENV, run_forge)
from forge_design.evaluate.runner_wt import (  # noqa: E402
    _bcond, _config_euler, _config_euler_node, design_chain)
from forge_design.meshing.mesh2d import (  # noqa: E402
    Mesh2DParams, generate_axisym_mesh, write_msh41_2d)
from forge_design.probdef import load_problem  # noqa: E402

PROBLEM = CASE / "problem_m4_r5.yaml"
RUN = CASE / "run_0039_r5_cold"


def main() -> int:
    p = load_problem(PROBLEM)
    RUN.mkdir(parents=True, exist_ok=False)
    d = design_chain(p, viscous=False)
    wall = d["wall"]
    scale = float(p.spec["r_throat"])
    mp = Mesh2DParams(ni=int(p.mesh["ni"]), nj=int(p.mesh["nj"]),
                      wall_first_frac=float(p.mesh["wall_first_frac"]),
                      throat_refine=float(p.mesh["throat_refine"]), scale=scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(RUN / "nozzle.msh", coords, quads, bedges)

    # 設計記録 (J = 設計壁始点 = wall_inv 先頭。C⁻ 追跡・うねり窓の入力)
    wi = d["wall_inv"]
    np.savetxt(RUN / "wall_design.csv", np.c_[wi * [scale, scale, 1.0, 1.0]],
               delimiter=",", header="x_m,r_m,theta_rad,M_wall", comments="")
    xs = np.linspace(d["x0"], d["xd"], 400)
    np.savetxt(RUN / "target_axis_M.csv",
               np.c_[xs * scale, [d["target"](x) for x in xs]], delimiter=",",
               header="x_m,M_target", comments="")
    info = {"R": d["R"], "anchor": "sauer_pass1", "discretization": "node",
            "x0": d["x0"], "xd": d["xd"], "Md": d["Md"],
            "x_j": float(wi[0, 0]), "r_j": float(wi[0, 1]),
            "x_reach_MOC": d["x_reach"], "mdot_ratio_moc": d["mdot_ratio_moc"],
            "nStepOuter": int(p.evaluate["nStepOuter"]),
            "control_runs": ["run_0031_baseline_cold", "run_0034_node_base"]}
    (RUN / "baseline_info.json").write_text(json.dumps(info, indent=1))

    n = int(p.evaluate["nStepOuter"])
    out_int = int(p.evaluate["outStepInterval"])
    (RUN / "bcondConfig.yaml").write_text(_bcond(p, euler=True))
    (RUN / "probe.yaml").write_text(PROBE_STUB)

    # 品質検査は cell 変換の一時コピーで (品質ツールは node CONNE 非対応)
    (RUN / "solverConfig.yaml").write_text(_config_euler(p, n, out_int, 4.0, 1))
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle_qc.h5"],
                   cwd=RUN, env=_ENV, check=True, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                        "nozzle_qc.h5"], cwd=RUN, env=_ENV,
                       capture_output=True, text=True)
    (RUN / "MESH_QUALITY.txt").write_text(
        "# cell 変換コピーで検査 (品質は primal の性質)\n" + q.stdout + q.stderr)
    (RUN / "nozzle_qc.h5").unlink()
    print(q.stdout.strip().splitlines()[-1] if q.stdout else q.stderr)
    if q.returncode != 0:
        raise RuntimeError("メッシュ品質 FAIL — 投入中止")

    # node 変換 (node config を書いてから変換 — 必須) + 等エントロピー IC
    (RUN / "solverConfig.yaml").write_text(_config_euler_node(p, n, out_int, 4.0, 1))
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                   cwd=RUN, env=_ENV, check=True, capture_output=True, text=True)
    paste_isentropic_ic(RUN / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp)
    print(json.dumps(info, indent=1))
    if "--prepare-only" in sys.argv:
        return 0
    rc = run_forge(RUN)
    print(f"forge exit={rc}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
