"""① 風洞 (モード F) の評価: 逆設計 → メッシュ → forge → 軸 M 達成度 (Phase 3)。

v1 チェーン: dv (Bézier 自由 CP + 軸長 L_ax + 円弧 R) → Sauer starting line →
逆 MOC マーチ → 壁 (Euler 用は非粘性壁 / NS 用は +δ* 経験式) → ModeFWall →
TFI メッシュ → forge (Euler=cell/slip/visc0 | RANS-SST=node) → 軸 M vs 目標。

使い方:
  design/.venv-opt/bin/python -m forge_design.evaluate.runner_wt \
      case/41.wind_tunnel_design/problem_m4.yaml run_dir [--euler] [--steps N]
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..geometry.bezier import MachBezier
from ..geometry.moc_inverse import inverse_design
from ..geometry.transonic import SauerThroat
from ..geometry.wall_modef import ModeFWall
from ..feedback.deltastar import deltastar_offset
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from .ic import paste_isentropic_ic
from .runner import FORGE_BUILD, FORGE_TOOLS, PROBE_STUB, _ENV, run_forge


def design_chain(p: Problem, viscous: bool) -> dict:
    """dv → 目標 Bézier → 逆設計壁 (+δ*) → ModeFWall。決定的。"""
    g = p.gamma
    Md = float(p.spec["M_design"])
    R = float(p.geometry.get("R", 2.0))
    st = SauerThroat(R=R, gamma=g)
    x0, rr0, MM0, tt0 = st.starting_line(M_start=1.05, n=int(p.geometry.get("n_start", 41)))
    h = 1e-4
    M0 = float(st.mach(x0, 0.0))
    dM0 = float((st.mach(x0 + h, 0.0) - st.mach(x0 - h, 0.0)) / (2 * h))
    xd = dv_value(p, "L_ax")
    cps = [dv_value(p, f"mc_cp{i+1}") for i in range(3)]
    bz = MachBezier.from_constraints(x0, xd, start=(M0, dM0), free_cp=cps,
                                     end=(Md, 0.0, 0.0))

    def target(x: float) -> float:
        return Md if x >= xd else float(np.atleast_1d(bz(float(x)))[0])

    eps = (1.0 / Md) * ((2.0 + (g - 1.0) * Md * Md) / (g + 1.0)) ** (
        (g + 1.0) / (2.0 * (g - 1.0)))
    x_end = xd + 2.3 * np.sqrt(eps) / np.tan(np.arcsin(1.0 / Md))
    res = inverse_design(st, target, x_axis_end=float(x_end),
                         n_axis=int(p.geometry.get("n_axis_inv", 500)),
                         n_start=int(p.geometry.get("n_start", 41)), gamma=g,
                         th_wall0=float(np.arcsin(min(x0 / R, 1.0))))
    wall_inv = res["wall"]
    pts = (deltastar_offset(wall_inv, float(p.spec["r_throat"]),
                            float(p.spec["Pt"]), float(p.spec["Tt"]), g, p.cp)
           if viscous else wall_inv[:, :2])
    wall = ModeFWall(pts, R=R,
                     r_inlet=float(p.geometry.get("r_inlet", 2.5)),
                     L_pipe=float(p.geometry.get("L_pipe", 0.5)),
                     L_contract=float(p.geometry.get("L_contract", 3.0)),
                     phi_u=np.deg2rad(float(p.geometry.get("phi_u_deg", 25.0))))
    msgs = wall.validate()
    if msgs:
        raise ValueError("モード F 壁フィルタ不合格: " + "; ".join(msgs))
    return {"wall": wall, "wall_inv": wall_inv, "target": target, "x0": float(x0),
            "xd": float(xd), "Md": Md,
            "mdot_ratio_moc": res["mdot_exit"] / res["mdot_start"]}


def _config_euler(p: Problem, nsteps: int, out_int: int, cfl: float, conv: int) -> str:
    return f"""mesh: {{meshFormat: "hdf5", discretization: "cell", isAxisymmetric: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 0,
           ro: 1.2, visc: 0.0, thermCond: 0.0, cp: {p.cp}, gamma: {p.gamma}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl},
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: {conv}, limiter: 2}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
"""


def _bcond(p: Problem, euler: bool) -> str:
    Pt, Tt = float(p.spec["Pt"]), float(p.spec["Tt"])
    pa = float(p.spec.get("p_ambient", 1000.0))
    wall_kind = "slip" if euler else "wall"
    return f"""inlet:  {{physID: 1, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {{Pt: {Pt}, Tt: {Tt}, k: 1.0, omega: 18000.0}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {{Ps: {pa}, Pt: {pa}, Tt: 300.0}}}}
wall:   {{physID: 3, kind: {wall_kind},             outputHDFflg: 1, ints: , floats: }}
axis:   {{physID: 4, kind: axis,             outputHDFflg: 0, ints: , floats: }}
"""


def prepare(problem_path, run_dir, euler: bool = True, nsteps=None) -> dict:
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym":
        raise ValueError("runner_wt は wind_tunnel_axisym 専用")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    d = design_chain(p, viscous=not euler)
    wall = d["wall"]
    scale = float(p.spec["r_throat"])
    mp = Mesh2DParams(ni=int(p.mesh.get("ni", 321)), nj=int(p.mesh.get("nj", 65)),
                      wall_first_frac=float(p.mesh.get("wall_first_frac", 5.0e-3)),
                      throat_refine=float(p.mesh.get("throat_refine", 3.0)),
                      scale=scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    # 目標軸分布と設計壁の記録 (帰還 v2 とΔM 評価の入力)
    xs = np.linspace(d["x0"], d["xd"], 400)
    np.savetxt(run_dir / "target_axis_M.csv",
               np.c_[xs * scale, [d["target"](x) for x in xs]], delimiter=",",
               header="x_m,M_target", comments="")
    np.savetxt(run_dir / "wall_design.csv",
               np.c_[d["wall_inv"] * [scale, scale, 1.0, 1.0]], delimiter=",",
               header="x_m,r_m,theta_rad,M_wall", comments="")
    n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
    (run_dir / "solverConfig.yaml").write_text(
        _config_euler(p, n, int(p.evaluate.get("outStepInterval", max(n // 3, 1))),
                      4.0, 1))
    (run_dir / "bcondConfig.yaml").write_text(_bcond(p, euler))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                        "nozzle.h5"], cwd=run_dir, env=_ENV, capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(q.stdout + q.stderr)
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp)
    info = {"x0": d["x0"], "xd": d["xd"], "Md": d["Md"],
            "mdot_ratio_moc": d["mdot_ratio_moc"], "scale_m": scale,
            "mesh": {"ni": mp.ni, "nj": mp.nj}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_staged(problem_path, run_dir, nsteps=None) -> None:
    """soft 段 (1次+cfl0.5, 3000 step) → 本段 — bell と同じ 2 段起動を config
    書換えで同一 run dir 内実行 (soft 完走後に本 config へ戻し restart なしで
    冷間から回すのでなく、soft の最終場を IC に採用する)。"""
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    cfg_soft = cfg_main.replace("cfl: 4.0, cfl_pseudo: 4.0", "cfl: 0.5, cfl_pseudo: 0.5")
    cfg_soft = cfg_soft.replace("convMethod: 1", "convMethod: 0")
    import re
    cfg_soft = re.sub(r"nStepOuter: \d+", "nStepOuter: 3000", cfg_soft)
    cfg_soft = re.sub(r"outStepInterval: \d+", "outStepInterval: 3000", cfg_soft)
    (run_dir / "solverConfig.yaml").write_text(cfg_soft)
    run_forge(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if not res or int("".join(c for c in res[-1].stem if c.isdigit())) < 3000:
        raise RuntimeError("soft 段が失敗")
    # soft 最終場を nozzle.h5 の IC に移植 (同一メッシュ → interp は index 一致)
    subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                    str(res[-1]), str(run_dir / "nozzle.h5")],
                   env=_ENV, check=True, capture_output=True, text=True)
    for f in run_dir.glob("res_*"):
        f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    run_forge(run_dir)


def collect(problem_path, run_dir) -> dict:
    """達成軸 M vs 目標の ΔM (固定格子)・質量流量比を評価。"""
    from ..metrics.extract import axis_mach

    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    tgt = np.loadtxt(run_dir / "target_axis_M.csv", delimiter=",", skiprows=1)
    x_ach, M_ach = axis_mach(run_dir / "nozzle.h5", res[-1],
                             axis_band=0.9 * float(p.spec["r_throat"]) * 0.1)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    scale = info["scale_m"]
    # 評価窓: M>1.05 の少し先 〜 目標終端手前 (§4.7(c) のマスク方針)
    xs = np.linspace(info["x0"] * scale + 0.5 * scale, info["xd"] * scale - 0.3 * scale, 200)
    Mt = np.interp(xs, tgt[:, 0], tgt[:, 1])
    Ma = np.interp(xs, x_ach, M_ach)
    dM = Ma - Mt
    out = {"res_file": res[-1].name,
           "dM_max": float(np.max(np.abs(dM))),
           "dM_max_rel_Md": float(np.max(np.abs(dM)) / info["Md"]),
           "dM_rms": float(np.sqrt(np.mean(dM ** 2))),
           "M_axis_exit": float(np.interp(info["xd"] * scale, x_ach, M_ach)),
           "mdot_ratio_moc": info["mdot_ratio_moc"]}
    np.savetxt(run_dir / "achieved_vs_target.csv", np.c_[xs, Mt, Ma],
               delimiter=",", header="x_m,M_target,M_achieved", comments="")
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="モード F 風洞評価 (Phase 3 v1)")
    ap.add_argument("problem")
    ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    a = ap.parse_args(argv)
    info = prepare(a.problem, a.run_dir, euler=True, nsteps=a.steps)
    print(json.dumps(info, indent=1))
    if a.prepare_only:
        return 0
    run_staged(a.problem, a.run_dir, a.steps)
    print(json.dumps(collect(a.problem, a.run_dir), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
