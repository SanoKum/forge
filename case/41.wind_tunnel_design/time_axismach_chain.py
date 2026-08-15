#!/usr/bin/env python3
"""axis-Mach チェーン各工程の wall-clock 内訳を実測する (ユーザ要望 2026-08-15)。

design → mesh → forge → 帰還抽出 の全段を個別に計測し、`axismach_timing.json` に
書き出す。CFD は実際に 1 本回す (soft 3000 + 本段 12000 step)。
計測用 run は `run_9901_timing_probe` (使い捨て・完了後に残すが破棄可)。
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))

from forge_design.evaluate.runner import FORGE_BUILD, FORGE_TOOLS, _ENV  # noqa: E402
from forge_design.geometry.axis_law import QuinticHermiteAxisLaw  # noqa: E402
from forge_design.geometry.moc_inverse import inverse_design  # noqa: E402
from forge_design.geometry.transonic import HallThroat  # noqa: E402
from forge_design.geometry.wall_axismach import (AxisMachCFDWall,  # noqa: E402
                                                 area_ratio_isentropic, wall_qa)
from forge_design.meshing.mesh2d import (Mesh2DParams,  # noqa: E402
                                         generate_axisym_mesh, write_msh41_2d)

G, MD, R, SCALE = 1.4, 4.0, 2.0, 0.01
X_REACH, ANCHOR, L_C = 1.7839, (2.07549, 0.64351, -0.04266), 7.8465
RUN = Path(__file__).parent / "run_9901_timing_probe"
SEED = Path(__file__).parent / "run_0047_axismach_p3"

T: dict[str, float] = {}


class step:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t = time.perf_counter()
        return self

    def __exit__(self, *a):
        T[self.name] = time.perf_counter() - self.t
        print(f"  {self.name:38s} {T[self.name]:8.2f} s", flush=True)


def main() -> int:
    if RUN.exists():
        shutil.rmtree(RUN)
    RUN.mkdir(parents=True)
    print("=== 設計段 (CPU, design/) ===", flush=True)
    with step("1. Hall 遷音速解 (構築+starting_line+anchor)"):
        ht = HallThroat(R=R, gamma=G)
        x0 = ht.x_axis_of_mach(1.05)
        ht.starting_line(M_start=1.05, n=121)
        ht.axis_anchor(x0)
    with step("2. 軸 Mach law: 許容 L_c 窓探索"):
        QuinticHermiteAxisLaw.admissible_Lc_range(X_REACH, *ANCHOR, MD)
    with step("3. 軸 Mach law: 構築 + ゲート判定"):
        law = QuinticHermiteAxisLaw(X_REACH, L_C, *ANCHOR, MD)
        law.gates()
    rF = float(np.sqrt(area_ratio_isentropic(MD, G)))
    x_end = law.x_E + 2.3 * rF * float(np.sqrt(MD * MD - 1.0))
    with step("4. 逆 MOC (n_axis=2000, 場充填+壁流線)"):
        res = inverse_design(ht, lambda x: float(law(x)), x_axis_end=x_end,
                             n_axis=2000, n_start=121, gamma=G, dx_wall=0.005,
                             th_wall0=float(np.arctan(x0 / R)), M_start=1.05)
    with step("5. 壁 QA + CFD 壁組立 (spline)"):
        wall_qa(res["wall"], MD, law.x_E, G)
        wall = AxisMachCFDWall(res["wall"][:, :2], R=R, r_U=2.5, L_U=3.5, L_pipe=0.5)
        wall.validate()
    print("=== メッシュ段 ===", flush=True)
    mp = Mesh2DParams(ni=321, nj=65, wall_first_frac=5e-3, throat_refine=3.0, scale=SCALE)
    with step("6. TFI メッシュ生成 + msh4.1 書き出し"):
        coords, quads, bedges = generate_axisym_mesh(wall, mp)
        write_msh41_2d(RUN / "nozzle.msh", coords, quads, bedges)
    for f in ("solverConfig.yaml", "bcondConfig.yaml", "probe.yaml"):
        shutil.copy(SEED / f, RUN / f)
    cfg_main = (SEED / "solverConfig.yaml").read_text()
    with step("7. convertGmshToForge (cell, 品質検査用)"):
        cfg_cell = cfg_main.replace('discretization: "node"', 'discretization: "cell"')
        (RUN / "solverConfig.yaml").write_text(cfg_cell)
        subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle_qc.h5"],
                       cwd=RUN, env=_ENV, check=True, capture_output=True, text=True)
    with step("8. check_mesh_quality.py"):
        subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                        "nozzle_qc.h5"], cwd=RUN, env=_ENV, capture_output=True, text=True)
    (RUN / "nozzle_qc.h5").unlink()
    with step("9. convertGmshToForge (node, median-dual)"):
        (RUN / "solverConfig.yaml").write_text(cfg_main)
        subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                       cwd=RUN, env=_ENV, check=True, capture_output=True, text=True)
    with step("10. 準1D 等エントロピー IC 貼り"):
        from forge_design.evaluate.ic import paste_isentropic_ic
        paste_isentropic_ic(RUN / "nozzle.h5", wall, SCALE, 1.0e6, 800.0, G, 1004.5)
    with step("11. interp_field (warm start 移植)"):
        src = sorted(SEED.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(src), str(RUN / "nozzle.h5")], env=_ENV, check=True,
                       capture_output=True, text=True)
    print("=== CFD 段 (GPU) ===", flush=True)
    from forge_design.evaluate.runner import run_forge
    import re
    cfg_soft = cfg_main.replace("cfl: 4.0, cfl_pseudo: 4.0", "cfl: 0.5, cfl_pseudo: 0.5")
    cfg_soft = cfg_soft.replace("convMethod: 1", "convMethod: 0")
    cfg_soft = re.sub(r"nStepOuter: \d+", "nStepOuter: 3000", cfg_soft)
    cfg_soft = re.sub(r"outStepInterval: \d+", "outStepInterval: 3000", cfg_soft)
    (RUN / "solverConfig.yaml").write_text(cfg_soft)
    with step("12. forge soft 段 (1次+cfl0.5, 3000 step)"):
        run_forge(RUN)
    with step("13. interp_field (soft→本段の場移植)"):
        res_f = sorted(RUN.glob("res_[0-9]*.h5"),
                       key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(res_f[-1]), str(RUN / "nozzle.h5")], env=_ENV, check=True,
                       capture_output=True, text=True)
        for f in RUN.glob("res_*"):
            f.unlink()
    (RUN / "solverConfig.yaml").write_text(cfg_main)
    with step("14. forge 本段 (2次+cfl4, 12000 step)"):
        run_forge(RUN)
    print("=== 評価・帰還段 ===", flush=True)
    with step("15. check_convergence.py"):
        subprocess.run([sys.executable, str(FORGE_TOOLS / "check_convergence.py"), str(RUN)],
                       env=_ENV, capture_output=True, text=True)
    with step("16. check_quasisteady.py"):
        subprocess.run([sys.executable, str(FORGE_TOOLS / "check_quasisteady.py"), str(RUN)],
                       env=_ENV, capture_output=True, text=True)
    with step("17. 軸 M 抽出 + 出口一様性 (コア C⁺ トレース込み)"):
        from forge_design.metrics.extract import axis_mach, exit_uniformity
        r_last = sorted(RUN.glob("res_[0-9]*.h5"),
                        key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        axis_mach(RUN / "nozzle.h5", r_last, axis_band=0.9 * SCALE * 0.1)
        exit_uniformity(RUN / "nozzle.h5", r_last, MD, x_d=law.x_E * SCALE, gamma=G)
    with step("18. x_reach,CFD 抽出 (C⁻ トレース + 感度一式)"):
        from forge_design.geometry.cminus_cfd import extract_x_reach
        extract_x_reach(RUN, x0, 1.0 + x0 * x0 / (2.0 * R))
    with step("19. 軸 M 平滑化スプライン + アンカー評価"):
        from forge_design.evaluate.runner_axismach import axis_curve_node
        spl, _, _ = axis_curve_node(RUN, SCALE, lam=1e-5)
        [float(spl.derivative(k)(1.79)) if k else float(spl(1.79)) for k in (0, 1, 2)]

    tot = sum(T.values())
    groups = {
        "設計 (CPU)": sum(T[k] for k in T if k[0] in "12345"and k.split(".")[0] in
                          {"1", "2", "3", "4", "5"}),
        "メッシュ+IC": sum(T[k] for k in T if k.split(".")[0] in {"6", "7", "8", "9", "10", "11"}),
        "CFD (GPU)": sum(T[k] for k in T if k.split(".")[0] in {"12", "13", "14"}),
        "評価+帰還": sum(T[k] for k in T if k.split(".")[0] in {"15", "16", "17", "18", "19"}),
    }
    print("\n=== グループ内訳 ===")
    for k, v in groups.items():
        print(f"  {k:14s} {v:8.2f} s  ({v/tot*100:4.1f}%)")
    print(f"  {'合計 (1 pass)':14s} {tot:8.2f} s  ({tot/60:.1f} 分)")
    out = {"steps_sec": T, "groups_sec": groups, "total_sec": tot,
           "note": "1 pass = 設計→メッシュ→CFD→評価→次パス用アンカー抽出"}
    (Path(__file__).parent / "axismach_timing.json").write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
