r"""axis-Mach チェーン A4: Hall + 5次 Hermite 軸 Mach → 逆 MOC → node Euler 評価。

plan: plans/accepted/tooling-nozzle-axismach-chain.md §6 A4。既存経路 (runner_wt の
モード F / runner_walldriven) には触れない。

チェーン: HallThroat (初期値線 + 軸アンカー) → QuinticHermiteAxisLaw
(自由度 L_c のみ) → inverse_design (E→F 閉包込み) → wall_qa →
AxisMachCFDWall (直管 + U→T Hermite + 設計壁 spline。旧・縦線構成のみ
[T, x0] に骨接放物線が入る) →
TFI メッシュ → node Euler (段階起動)。

CFD-in-the-loop アンカー更新 (A5) は problem YAML の geometry キーで受ける:
  x_reach_cfd:     壁始点発 C⁻ の軸着地点 (node 基準 run から)
  x_reach_anchor:  [M, M', M''] (同 run の軸平滑化フィットから)
  axis_segment_run: [x0, x_reach) の実測軸 M を target_moc に渡す元 run
これらが無ければ初回 (Hall アンカー、x_A = x0)。
初期値線は geometry.start_line で選ぶ ('throat_char' = スロート特性線 [A8] /
'vertical' = M_start の縦線 [旧構成])。壁の決め方は geometry.wall_mode
('flux' = 断面の質量流束閉包 [A9] / 'streamline' = 流線積分 [旧構成])。

使い方:
  design/.venv-opt/bin/python -m forge_design.evaluate.runner_axismach \
      case/41.wind_tunnel_design/problem_m4_axismach.yaml run_dir [--prepare-only]
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..geometry.axis_law import QuinticHermiteAxisLaw
from ..geometry.moc_inverse import inverse_design
from ..geometry.transonic import HallThroat
from ..geometry.wall_axismach import (AxisMachCFDWall, area_ratio_isentropic,
                                      wall_qa)
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from .ic import paste_isentropic_ic
from .runner import FORGE_BUILD, FORGE_TOOLS, PROBE_STUB, _ENV, run_forge
from .runner_wt import _bcond, _config_euler, _config_euler_node


def axis_curve_node(run_dir, scale: float, lam: float = 1e-5):
    r"""node run の軸ノード列から平滑化スプライン M(x) を作る (§15 の実装)。

    node は軸上に真の DOF を持つ (帯平均不要)。末尾 3 スナップ平均 → r=0 ノード
    抽出 → `make_smoothing_spline`。M, M', M'' は**同一スプライン**から評価する
    (生の 2 階差分禁止)。戻り値: (spline, x_min, x_max) — x は r_t 単位。"""
    import h5py
    from scipy.interpolate import make_smoothing_spline
    rd = Path(run_dir)
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    res = [f for f in res if int("".join(c for c in f.stem if c.isdigit())) > 0]
    if not res:
        raise FileNotFoundError(f"{rd} に res_*.h5 が無い")
    with h5py.File(rd / "nozzle.h5") as nz:
        nc = nz["/MESH/COORD"][:].reshape(-1, 3)
    Ms = []
    for rf in res[-3:]:
        with h5py.File(rf) as f:
            Ux, Uy, son = f["/VALUE/Ux"][:], f["/VALUE/Uy"][:], f["/VALUE/sonic"][:]
        if len(Ux) != len(nc):
            raise ValueError("node run でない (VALUE 長 != 節点数)")
        Ms.append(np.hypot(Ux, Uy) / np.maximum(son, 1e-9))
    M = np.mean(Ms, axis=0)
    ax = nc[:, 1] < 1e-12 * max(float(nc[:, 1].max()), 1.0)
    x = nc[ax, 0] / scale
    o = np.argsort(x)
    x, Ma = x[o], M[ax][o]
    spl = make_smoothing_spline(x, Ma, lam=lam)
    return spl, float(x.min()), float(x.max())


def design_chain(p: Problem) -> dict:
    """Hall (+CFD アンカー) → Hermite law → 逆 MOC → 壁 QA → CFD 壁。決定的。"""
    g = p.gamma
    Md = float(p.spec["M_design"])
    R = float(p.geometry.get("R", 2.0))
    M_start = float(p.geometry.get("M_start", 1.05))
    n_start = int(p.geometry.get("n_start", 41))
    ht = HallThroat(R=R, gamma=g)
    # 初期値線: 'throat_char' = スロート壁点発 C⁻ (CONTUR 流、壁が T から MOC 出力に
    # なる) / 'vertical' = M_start の縦線 (旧構成、回帰対照)。x0 = その軸着地点。
    start_line = str(p.geometry.get("start_line", "vertical"))
    if start_line == "throat_char":
        x0 = float(ht.throat_characteristic(n=n_start)[0][0])
    elif start_line == "vertical":
        x0 = ht.x_axis_of_mach(M_start)
    else:
        raise ValueError("geometry.start_line は 'throat_char' か 'vertical'")

    # --- アンカー: 初回 = Hall 解析値 (x_A = x0)。反復 = CFD 実測 (x_A = x_reach) ---
    x_reach_cfd = p.geometry.get("x_reach_cfd")
    seg = None
    if x_reach_cfd is not None:
        x_A = float(x_reach_cfd)
        anc = p.geometry.get("x_reach_anchor")
        if anc is None or len(anc) != 3:
            raise ValueError("x_reach_cfd 指定時は x_reach_anchor: [M,M',M''] が必須")
        M_A, Mp_A, Mpp_A = float(anc[0]), float(anc[1]), float(anc[2])
        seg_run = p.geometry.get("axis_segment_run")
        if seg_run:
            seg, _, _ = axis_curve_node(seg_run, float(p.spec["r_throat"]))
    else:
        x_A = x0
        M_A, Mp_A, Mpp_A = ht.axis_anchor(x0)

    # --- L_c: 許容窓を必ず計算し、explicit は窓内検査・max は窓上限×0.98 ---
    lo, hi = QuinticHermiteAxisLaw.admissible_Lc_range(x_A, M_A, Mp_A, Mpp_A, Md)
    mode = str(p.geometry.get("Lc_mode", "explicit"))
    if mode == "max":
        L_c = hi * 0.98
    elif mode == "explicit":
        L_c = float(dv_value(p, "L_c"))
        if not (lo <= L_c <= hi):
            raise ValueError(f"L_c = {L_c:.4g} が許容窓 ({lo:.4g}, {hi:.4g}) 外 "
                             "(M'≥0 単調ゲート)")
    else:
        raise ValueError("geometry.Lc_mode は 'explicit' か 'max'")
    law = QuinticHermiteAxisLaw(x_A, L_c, M_A, Mp_A, Mpp_A, Md)
    gates = law.gates()
    if gates["violations"]:
        raise ValueError("軸 Mach law ゲート不合格: " + "; ".join(gates["violations"]))

    # --- target: [x0, x_A) は実測 (反復時) / [x_A, x_E] law / 以降 M_d ---
    def target_moc(x: float) -> float:
        x = float(x)
        if x < x_A:
            if seg is not None:
                return float(seg(np.float64(x)))
            return float(ht.mach(x, 0.0)) if x_reach_cfd is None else float(law(x_A))
        return float(law(x))

    rF_pred = float(np.sqrt(area_ratio_isentropic(Md, g)))
    x_end = law.x_E + float(p.geometry.get("x_end_margin", 2.3)) \
        * rF_pred * float(np.sqrt(Md * Md - 1.0))
    res = inverse_design(ht, target_moc, x_axis_end=float(x_end),
                         n_axis=int(p.geometry.get("n_axis_inv", 500)),
                         n_start=n_start, gamma=g,
                         dx_wall=float(p.geometry.get("dx_wall", 0.02)),
                         th_wall0=float(np.arctan(x0 / R)), M_start=M_start,
                         exit_mode=str(p.geometry.get("exit_mode", "characteristic")),
                         x_E=law.x_E, M_d=Md, start_line=start_line,
                         wall_mode=str(p.geometry.get("wall_mode", "streamline")),
                         blend_width=float(p.geometry.get("wall_blend_width", 1.0)))
    qa = wall_qa(res["wall"], Md, law.x_E, g)
    if qa["violations"]:
        raise ValueError("壁 QA 不合格: " + "; ".join(qa["violations"]))
    wall = AxisMachCFDWall(res["wall"][:, :2], R=R,
                           r_U=float(p.geometry.get("r_inlet", 2.5)),
                           L_U=float(p.geometry.get("L_U", 3.5)),
                           L_pipe=float(p.geometry.get("L_pipe", 0.5)))
    msgs = wall.validate()
    if msgs:
        raise ValueError("axis-Mach 壁フィルタ不合格: " + "; ".join(msgs))
    return {"wall": wall, "wall_inv": res["wall"], "law": law, "qa": qa,
            "exit": {k: v for k, v in res["exit"].items() if k != "term_path"},
            "x0": float(x0), "x_A": float(x_A), "x_E": float(law.x_E),
            "L_c": float(L_c), "Lc_window": (float(lo), float(hi)),
            "anchor": (float(M_A), float(Mp_A), float(Mpp_A)),
            "anchor_source": ("cfd" if x_reach_cfd is not None else "hall"),
            "start_line": start_line, "wall_mode": res["wall_mode"],
            "Md": Md, "R": R, "gates": gates,
            "mdot_ratio_moc": float(res["mdot_exit"] / res["mdot_start"]),
            "cd_series": float(ht.cd_series())}


def prepare(problem_path, run_dir, nsteps=None, ic_from=None) -> dict:
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym_axismach":
        raise ValueError("runner_axismach は wind_tunnel_axisym_axismach 専用")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    d = design_chain(p)
    wall = d["wall"]
    scale = float(p.spec["r_throat"])
    mp = Mesh2DParams(ni=int(p.mesh.get("ni", 321)), nj=int(p.mesh.get("nj", 65)),
                      wall_first_frac=float(p.mesh.get("wall_first_frac", 5.0e-3)),
                      throat_refine=float(p.mesh.get("throat_refine", 3.0)),
                      scale=scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    # 記録: 目標軸分布 (x0 → x_E) と設計壁
    law = d["law"]
    xs = np.linspace(d["x0"], d["x_E"], 400)
    tgt = [float(law(x)) if x >= d["x_A"] else float("nan") for x in xs]
    np.savetxt(run_dir / "target_axis_M.csv", np.c_[xs * scale, tgt],
               delimiter=",", header="x_m,M_target", comments="")
    np.savetxt(run_dir / "wall_design.csv",
               np.c_[d["wall_inv"] * [scale, scale, 1.0, 1.0]], delimiter=",",
               header="x_m,r_m,theta_rad,M_wall", comments="")
    n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 3, 1)))
    (run_dir / "bcondConfig.yaml").write_text(_bcond(p, euler=True))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    # 品質検査は cell 変換の一時コピー (品質ツールは node CONNE 非対応)
    (run_dir / "solverConfig.yaml").write_text(_config_euler(p, n, out_int, 4.0, 1))
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle_qc.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                        "nozzle_qc.h5"], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(
        "# cell 変換コピーで検査 (品質は primal の性質)\n" + q.stdout + q.stderr)
    (run_dir / "nozzle_qc.h5").unlink()
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    # node 変換 (node config を書いてから変換 — 必須) + 等エントロピー IC
    (run_dir / "solverConfig.yaml").write_text(_config_euler_node(p, n, out_int, 4.0, 1))
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp)
    if ic_from is not None:
        src = sorted(Path(ic_from).glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(src), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
    info = {"chain": "axismach", "discretization": "node",
            "x0": d["x0"], "x_A": d["x_A"], "x_E": d["x_E"], "L_c": d["L_c"],
            "Lc_window": list(d["Lc_window"]), "anchor": list(d["anchor"]),
            "anchor_source": d["anchor_source"], "start_line": d["start_line"],
            "wall_mode": d["wall_mode"],
            "Md": d["Md"], "R": d["R"],
            "qa": {k: v for k, v in d["qa"].items() if k != "violations"},
            "exit": d["exit"],
            "mdot_ratio_moc": d["mdot_ratio_moc"], "cd_series": d["cd_series"],
            "nStepOuter": n, "scale_m": scale, "ic_from": str(ic_from) if ic_from else None,
            "mesh": {"ni": mp.ni, "nj": mp.nj, "wall_first_frac": mp.wall_first_frac}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_staged(run_dir) -> int:
    """soft 段 (1次+cfl0.5, 3000 step) → 本段 (walldriven W3 と同方式・段階起動必須)。"""
    import re
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    cfg_soft = cfg_main.replace("cfl: 4.0, cfl_pseudo: 4.0", "cfl: 0.5, cfl_pseudo: 0.5")
    cfg_soft = cfg_soft.replace("convMethod: 1", "convMethod: 0")
    cfg_soft = re.sub(r"nStepOuter: \d+", "nStepOuter: 3000", cfg_soft)
    cfg_soft = re.sub(r"outStepInterval: \d+", "outStepInterval: 3000", cfg_soft)
    (run_dir / "solverConfig.yaml").write_text(cfg_soft)
    rc = run_forge(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if rc != 0 or not res or int("".join(c for c in res[-1].stem if c.isdigit())) < 3000:
        raise RuntimeError("soft 段が失敗 (発散切り分けは res_nan_*.h5 を見る)")
    subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                    str(res[-1]), str(run_dir / "nozzle.h5")],
                   env=_ENV, check=True, capture_output=True, text=True)
    for f in run_dir.glob("res_*"):
        f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    return run_forge(run_dir)


def collect(problem_path, run_dir) -> dict:
    """軸 M vs 目標 (x_A ≤ x ≤ x_E−0.3 マスク) + 出口一様性 + overshoot。"""
    from ..metrics.extract import axis_mach, exit_uniformity

    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    info = json.loads((run_dir / "prepare_info.json").read_text())
    scale = info["scale_m"]
    Md = float(info["Md"])
    x_ach, M_ach = axis_mach(run_dir / "nozzle.h5", res[-1],
                             axis_band=0.9 * scale * 0.1)
    tgt = np.loadtxt(run_dir / "target_axis_M.csv", delimiter=",", skiprows=1)
    ok = np.isfinite(tgt[:, 1])
    xs = np.linspace(max(info["x_A"], tgt[ok, 0].min() / scale) * scale + 1e-9,
                     (info["x_E"] - 0.3) * scale, 200)
    Mt = np.interp(xs, tgt[ok, 0], tgt[ok, 1])
    Ma = np.interp(xs, x_ach, M_ach)
    dM = Ma - Mt
    # overshoot: x_E 以降〜出口手前の軸 M の最大値
    m_dn = x_ach > info["x_E"] * scale
    M_max_dn = float(np.max(M_ach[m_dn])) if m_dn.any() else float("nan")
    try:
        uni = exit_uniformity(run_dir / "nozzle.h5", res[-1], Md,
                              x_d=info["x_E"] * scale, gamma=p.gamma)
    except Exception as e:  # noqa: BLE001 — 一様性は診断 (核心は軸 M)
        uni = {"error": str(e)}
    out = {"res_file": res[-1].name,
           "dM_max": float(np.max(np.abs(dM))),
           "dM_max_rel_Md": float(np.max(np.abs(dM)) / Md),
           "dM_rms": float(np.sqrt(np.mean(dM ** 2))),
           "M_axis_exit": float(np.interp(info["x_E"] * scale, x_ach, M_ach)),
           "M_axis_max_downstream": M_max_dn,
           "overshoot_rel": float((M_max_dn - Md) / Md) if np.isfinite(M_max_dn) else None,
           "exit_uniformity": uni,
           "mdot_ratio_moc": info["mdot_ratio_moc"]}
    np.savetxt(run_dir / "achieved_vs_target.csv", np.c_[xs, Mt, Ma],
               delimiter=",", header="x_m,M_target,M_achieved", comments="")
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="axis-Mach チェーン評価 (node Euler)")
    ap.add_argument("problem")
    ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    a = ap.parse_args(argv)
    info = prepare(a.problem, a.run_dir, nsteps=a.steps)
    print(json.dumps(info, indent=1))
    if a.prepare_only:
        return 0
    rc = run_staged(a.run_dir)
    print(f"forge exit={rc}")
    print(json.dumps(collect(a.problem, a.run_dir), indent=1))
    return rc


if __name__ == "__main__":
    sys.exit(main())
