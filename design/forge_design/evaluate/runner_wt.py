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
from ..geometry.cminus_map import build_cminus_map
from ..geometry.moc_inverse import inverse_design
from ..geometry.transonic import SauerThroat
from ..geometry.wall_modef import ModeFWall
from ..feedback.deltastar import deltastar_offset
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from .ic import paste_isentropic_ic
from .runner import FORGE_BUILD, FORGE_TOOLS, PROBE_STUB, _ENV, run_forge


def design_chain(p: Problem, viscous: bool) -> dict:
    r"""dv → 目標 Bézier → 逆設計壁 (+δ*) → ModeFWall。決定的。

    **B7 (plan §9.2)**: 目標曲線は $x_{\mathrm{reach}}$ (設計壁始点 $x_0$ から出て
    最初に成功する C⁻ の軸着地点 — 壁の帰還がそもそも届かない境界) で分割する。
    $x_0\le x<x_{\mathrm{reach}}$ は実測 (CFD アンカー時は `axis_segment`、Sauer
    時は解析式) をそのまま目標にし、自由 CP は $x\ge x_{\mathrm{reach}}$ だけを
    担当する。Bézier のグローバル台 (Bernstein 基底は局所支持を持たない) が遠方
    CP の形を到達不能域まで押し戻す問題を解消する (実測: 分割前は
    $x/r_t\approx1.38$ で CP2 が既に 11% の重みを持ち $\Delta M\approx+0.05$ の
    主残差を作っていた)。$x_{\mathrm{reach}}$ 自体は下流形状にほぼ依存しない
    ため、**まず旧来の単一 Bézier ("ブートストラップ") で仮設計して MOC マップの
    被覆から $x_{\mathrm{reach}}$ を得てから、本目標を組んで MOC をやり直す**。
    """
    g = p.gamma
    Md = float(p.spec["M_design"])
    R = float(p.geometry.get("R", 2.0))
    M_start = float(p.geometry.get("M_start", 1.05))   # 物理定数でなく数値上の選択
    # **B6 CFD アンカー** (plan §9.2): `geometry.throat_anchor_run` があれば
    # Sauer をその run の CFD 抽出値で置き換える (アンカーは抽出元 run に凍結)
    anchor_run = p.geometry.get("throat_anchor_run")
    if anchor_run:
        from ..feedback.cfd_anchor import from_run
        st = from_run(anchor_run, float(p.spec["r_throat"]), R, g)
    else:
        st = SauerThroat(R=R, gamma=g)
    n_start = int(p.geometry.get("n_start", 41))
    x0, rr0, MM0, tt0 = st.starting_line(M_start=M_start, n=n_start)
    h = 1e-4
    M0 = float(st.mach(x0, 0.0))
    dM0 = float((st.mach(x0 + h, 0.0) - st.mach(x0 - h, 0.0)) / (2 * h))
    d2M0 = float((st.mach(x0 + h, 0.0) - 2.0 * M0 + st.mach(x0 - h, 0.0)) / h ** 2)
    # **B1: ブートストラップ段の中心線接続次数のみに効く** (plan §9.2)。B7 導入後は
    # 本目標側の x_reach 接続が常に C2 (M,M',M'') を課すため、ここは x_reach の
    # 推定に使う仮設計の質にしか影響しない。
    order = int(p.geometry.get("centerline_start_order", 1))
    if order not in (1, 2):
        raise ValueError("geometry.centerline_start_order は 1 か 2")
    start_bc = (M0, dM0) if order == 1 else (M0, dM0, d2M0)
    xd = dv_value(p, "L_ax")
    cps = [dv_value(p, f"mc_cp{i+1}") for i in range(3)]
    eps = (1.0 / Md) * ((2.0 + (g - 1.0) * Md * Md) / (g + 1.0)) ** (
        (g + 1.0) / (2.0 * (g - 1.0)))
    x_end = xd + 2.3 * np.sqrt(eps) / np.tan(np.arcsin(1.0 / Md))
    n_axis_inv = int(p.geometry.get("n_axis_inv", 500))
    th_wall0 = float(np.arcsin(min(x0 / R, 1.0)))

    # === B8 (2026-08-16): **設計 = B6 グローバル Bézier / 評価 = B7 分割帳簿** ===
    # 経緯: B7 は x_reach で設計目標を再アンカーした (Hermite + 再アンカー Bézier)
    # が、意図非依存のうねり指標 (実流 M(x) の 3 次カーブからの片振幅) で測ると
    # スロート近傍の**実流が後退**していた: Sauer ±0.073 → B6 ±0.023 → B7 ±0.062。
    # B6 のグローバル Bézier は (CP ジグザグにも関わらず) 逆設計を通じ実質的に
    # 接合波をキャンセルする壁を作っており、設計目標の構成としては B6 が上。
    # B7 の帳簿 (到達不能域 ±0.009) は「凍結アンカー = Sauer 壁の波だらけの場」を
    # 参照していたため後退を隠した — アンカーは B6 壁の場 (run_0018) へ更新する
    # (problem YAML の throat_anchor_run)。
    # 残す B7 資産: x_reach の概念 (マスク w=0 境界 + 評価分割点)、実測こぶ曲線に
    # よる到達不能域の評価、被覆外/隙間ステーションの除外、規約整合抽出。
    # 設計はブートストラップ 1 回で完結する (2 回目の逆設計は不要になった)。
    bz = MachBezier.from_constraints(x0, xd, start=start_bc, free_cp=cps,
                                     end=(Md, 0.0, 0.0))

    def target(x: float) -> float:          # 設計用 = B6 グローバル (滑らか)
        return Md if x >= xd else float(np.atleast_1d(bz(float(x)))[0])

    res = inverse_design(st, target, x_axis_end=float(x_end), n_axis=n_axis_inv,
                         n_start=n_start, gamma=g, th_wall0=th_wall0, M_start=M_start)
    cmap = build_cminus_map(res["pts"], res["wall"], g)
    if len(cmap) < 4:
        raise ValueError("C⁻ マップ点が不足 (x_reach 未確定)")
    x_reach = float(cmap[:, 1].min())

    # 評価用: 到達不能域 [x0,x_reach) は実測こぶ曲線 (帯規約 — axis_mach と同じ
    # 母集団で抽出。規約が違うだけで系統残差 0.04 になる)、以降は設計 Bézier。
    # 「狙い通りか」を到達可能性に即して測る。x_reach の段差は「届く波の残余 =
    # 帰還が打ち消すべき量」として残る。
    if hasattr(st, "axis_segment_curve"):   # CFDThroat
        seg = st.axis_segment_curve(x0, x_reach)

        def target_eval(x: float) -> float:
            x = float(x)
            if x < x_reach:
                return float(seg(np.float64(x)))
            return target(x)
    else:                                    # SauerThroat: 分割は自然に消える
        target_eval = target

    wall_inv = res["wall"]
    pts = (deltastar_offset(wall_inv, float(p.spec["r_throat"]),
                            float(p.spec["Pt"]), float(p.spec["Tt"]), g, p.cp)
           if viscous else wall_inv[:, :2])
    wall = ModeFWall(pts, R=R,
                     r_inlet=float(p.geometry.get("r_inlet", 2.5)),
                     L_pipe=float(p.geometry.get("L_pipe", 0.5)),
                     L_contract=float(p.geometry.get("L_contract", 3.0)),
                     phi_u=np.deg2rad(float(p.geometry.get("phi_u_deg", 25.0))),
                     blend_len=float(p.geometry.get("blend_len", 0.0)))
    msgs = wall.validate()
    if msgs:
        raise ValueError("モード F 壁フィルタ不合格: " + "; ".join(msgs))
    return {"wall": wall, "wall_inv": wall_inv, "target": target,
            "target_eval": target_eval, "x0": float(x0),
            "xd": float(xd), "Md": Md, "pts": res["pts"], "R": R,
            "centerline_start_order": order, "sauer_d2M0": d2M0,
            "x_reach": x_reach, "target_d2M0": float(np.atleast_1d(bz.deriv(x0, 2))[0]),
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
