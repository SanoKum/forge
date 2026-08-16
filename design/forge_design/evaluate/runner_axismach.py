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
from .runner_wt import (_bcond, _config_euler, _config_euler_node,
                        _config_sst_node)


def axis_curve_node(run_dir, scale: float, lam: float = 1e-5, mode: str = "evenfit"):
    r"""node run の軸列から平滑化スプライン M(x) を作る (§15 の実装)。

    mode="evenfit" (既定): 軸ノードを除外し r>0 の 4 点で M=a₀+a₂r² を最小二乗した a₀ を軸値に
    使う (`axis_extract.axis_curve_evenfit`)。旧生産設定 (`nodeAxisDirichlet: 1`, 撤去済) の軸ノードは第一内点の
    コピーで ½M_rr r₁² のバイアス (case/41 で ~0.08% M_d) を持っていたため。現行 (軸 DOF) でも外挿の方が精度が良い。
    mode="node": 従来の軸ノード直読。末尾 3 スナップ平均 → `make_smoothing_spline`。
    M, M', M'' は**同一スプライン**から評価する (生の 2 階差分禁止)。
    戻り値: (spline, x_min, x_max) — x は r_t 単位。"""
    if mode == "evenfit":
        from .axis_extract import axis_curve_evenfit
        return axis_curve_evenfit(run_dir, scale, lam=lam)
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


def _gam_or_gas(p: Problem):
    """MOC/幾何関数に渡す γ: cpg なら float、semiperfect ならガスモデル
    (cfd_gas: cpg のときは設計も CPG に落とす — 熱力学の一致)。"""
    if str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        return p.gamma
    return p.gas_model if p.is_semiperfect else p.gamma


def _apply_gas_to_config(cfg: str, p: Problem, run_dir) -> str:
    """semi-perfect のとき forge config を **単一擬似種 TP** (thermalMethod 2) に書き換え、
    NASA-9 混合擬似種を run_dir/species_db.yaml へ出す。cpg なら無変更。
    設計 (MOC) と同一係数の熱力学で CFD が回る (forge 内蔵 DB と同じ CEA 値)。"""
    # 切り分け用: evaluate.axisym_method で CPG でも SU2 流軸対称に切替可
    if int(p.evaluate.get("axisym_method", 0)) == 1:
        cfg = cfg.replace("isAxisymmetric: 1", "isAxisymmetric: 1, axisymMethod: 1", 1)
    if not p.is_semiperfect or str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        # cfd_gas: cpg = 設計は semi-perfect のまま CFD だけ CPG(γ*, cp 参照値) で回す
        # (TP × node 軸対称の forge 側発散 [case/42 run_0001] の回避。相対比較には十分)
        return cfg
    import yaml as _yaml
    from ..gas.semiperfect import mixture_pseudo_species
    Y = dict(p.raw["gas"]["species"])
    db = mixture_pseudo_species(Y, "MIX")
    (Path(run_dir) / "species_db.yaml").write_text(_yaml.safe_dump(db, sort_keys=False))
    # thermalMethod 0 → 2、species/speciesDBFile を physProp に追加 (cp/gamma は参照値のまま
    # 残すが TP では NASA-9 が優先される)
    cfg = cfg.replace("thermalMethod: 0", "thermalMethod: 2", 1)
    # [2026-08-16] nodeAxisDirichlet は撤去済み (node は軸ノードを DOF として解く整合セットが常時 ON、
    # TP スカラーの軸ピン問題は消滅)。axisRFloor は evaluate.axis_r_floor 指定時のみ従来どおり付ける。
    floor = float(p.evaluate.get("axis_r_floor", 0.0))
    if floor > 0:
        cfg = cfg.replace("isAxisymmetric: 1", f"isAxisymmetric: 1, axisRFloor: {floor}", 1)
    cfg = cfg.replace("cp: %s, gamma: %s}" % (p.cp, p.gamma),
                      "cp: %s, gamma: %s,\n           species: [MIX], speciesDBFile: \"species_db.yaml\"}"
                      % (p.cp, p.gamma), 1)
    if "species: [MIX]" not in cfg:
        raise RuntimeError("_apply_gas_to_config: physProp の書き換えに失敗 (テンプレート変更?)")
    return cfg



def design_chain(p: Problem) -> dict:
    """Hall (+CFD アンカー) → Hermite law → 逆 MOC → 壁 QA → CFD 壁。決定的。

    ガス: `p.gas_model` (cpg なら float γ と等価 / semiperfect なら NASA-9 テーブル)。
    Hall 遷音速級数は定数 γ 前提なので**スロートの局所 γ* を渡す。MOC の ν↔M・
    質量流束・面積比はガスモデル経由 (`moc_kernel._is_gas` 規約)。"""
    gas = p.gas_model
    # cfd_gas: cpg のときは**設計も** CPG(γ 参照値) で作る。設計だけ semi-perfect にすると
    # 壁 (A/A*=13.1) と CFD (CPG γ=1.309 なら A/A*=15.3) の熱力学が食い違い、出口 M が
    # 3.86 で止まる (case/42 run_0003–0011 で実測、破棄)。設計と CFD の熱力学は必ず一致。
    if str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        from ..gas import GasCPG
        gas = GasCPG(p.gamma, p.cp)
    use_gas = getattr(gas, "kind", "cpg") == "semiperfect"
    g = gas if use_gas else p.gamma                   # MOC/面積比に渡す「γ or ガス」
    g_hall = float(gas.gamma_throat(float(p.spec["Tt"])))
    Md = float(p.spec["M_design"])
    R = float(p.geometry.get("R", 2.0))
    M_start = float(p.geometry.get("M_start", 1.05))
    n_start = int(p.geometry.get("n_start", 41))
    ht = HallThroat(R=R, gamma=g_hall)
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
    qa = wall_qa(res["wall"], Md, law.x_E, g, R=R)
    if qa["violations"]:
        raise ValueError("壁 QA 不合格: " + "; ".join(qa["violations"]))
    # 壁表現 (A14): 'interp' = 補間 5 次 B-spline (現行・比較基準) /
    # 'lsq' = 制約付き最小二乗 B-spline (n_cp は誤差ゲートで自動、または wall_ncp)
    wall_repr = str(p.geometry.get("wall_repr", "interp"))
    wkw = dict(R=R, r_U=float(p.geometry.get("r_inlet", 2.5)),
               L_U=float(p.geometry.get("L_U", 3.5)),
               L_pipe=float(p.geometry.get("L_pipe", 0.5)))
    if wall_repr == "lsq":
        from ..geometry.wall_axismach import LSQBsplineCFDWall
        wall = LSQBsplineCFDWall(res["wall"], n_cp=p.geometry.get("wall_ncp"),
                                 tol_dr=float(p.geometry.get("wall_lsq_tol_dr", 5e-4)),
                                 tol_dtheta_deg=float(p.geometry.get("wall_lsq_tol_dtheta", 0.05)),
                                 **wkw)
    elif wall_repr == "interp":
        wall = AxisMachCFDWall(res["wall"], **wkw)
    else:
        raise ValueError("geometry.wall_repr は 'interp' か 'lsq'")
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
            "wall_repr": wall_repr,
            "gas": (gas.summary() if hasattr(gas, "summary")
                    else {"kind": "cpg", "gamma": p.gamma, "cp": p.cp}),
            "gamma_hall": g_hall,
            "wall_fit": getattr(wall, "fit_diag", None),
            "Md": Md, "R": R, "gates": gates,
            # cplus 閉包では構成的に 1 になる循環指標なので出さない (A9 の教訓)
            "mdot_ratio_moc": (None if not np.isfinite(res["mdot_exit"])
                               else float(res["mdot_exit"] / res["mdot_start"])),
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
    (run_dir / "solverConfig.yaml").write_text(
        _apply_gas_to_config(_config_euler(p, n, out_int, 4.0, 1), p, run_dir))
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
    (run_dir / "solverConfig.yaml").write_text(
        _apply_gas_to_config(_config_euler_node(p, n, out_int, 4.0, 1), p, run_dir))
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp,
                        gas=(None if str(p.evaluate.get('cfd_gas', 'same')) == 'cpg' else p.gas_model))
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
            "wall_mode": d["wall_mode"], "wall_repr": d["wall_repr"],
            "gas": d["gas"], "gamma_hall": d["gamma_hall"],
            "wall_fit": d["wall_fit"],
            "Md": d["Md"], "R": d["R"],
            "qa": {k: v for k, v in d["qa"].items() if k != "violations"},
            "exit": d["exit"],
            "mdot_ratio_moc": d["mdot_ratio_moc"], "cd_series": d["cd_series"],
            "nStepOuter": n, "scale_m": scale, "ic_from": str(ic_from) if ic_from else None,
            "mesh": {"ni": mp.ni, "nj": mp.nj, "wall_first_frac": mp.wall_first_frac}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_staged(run_dir, cfl_main: float | None = None, mid_stage: bool = False) -> int:
    """soft 段 (1次+cfl0.5, 3000 step) → [mid 段 (2次+cfl1, 3000 step)] → 本段。
    walldriven W3 と同方式・段階起動必須。`cfl_main` で本段 CFL を上書き
    (semi-perfect TP は cfl4 で本段 step ~60 に爆発 — case/42 実測。TP は cfl≤2 が実績
    [[cutler-cpg-vs-tp-dplur-sst]])。`mid_stage=True` で 2 次化と CFL 上げを分離する。"""
    import re
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    if cfl_main is not None:
        cfg_main = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+",
                          f"cfl: {cfl_main}, cfl_pseudo: {cfl_main}", cfg_main)

    def _stage(cfg, nsteps, label):
        cfg = re.sub(r"nStepOuter: \d+", f"nStepOuter: {nsteps}", cfg)
        cfg = re.sub(r"outStepInterval: \d+", f"outStepInterval: {nsteps}", cfg)
        (run_dir / "solverConfig.yaml").write_text(cfg)
        rc = run_forge(run_dir)
        res = sorted(run_dir.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res or int("".join(c for c in res[-1].stem if c.isdigit())) < nsteps:
            raise RuntimeError(f"{label} 段が失敗 (発散切り分けは res_nan_*.h5 を見る)")
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(res[-1]), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
        for f in run_dir.glob("res_*"):
            f.unlink()

    soft = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 0.5, cfl_pseudo: 0.5", cfg_main)
    soft = soft.replace("convMethod: 1", "convMethod: 0")
    _stage(soft, 3000, "soft")
    if mid_stage:
        mid = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 1.0, cfl_pseudo: 1.0", cfg_main)
        _stage(mid, 3000, "mid")
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
    pp = load_problem(a.problem)
    rc = run_staged(a.run_dir, cfl_main=pp.evaluate.get("cfl_main"),
                    mid_stage=bool(pp.evaluate.get("mid_stage", pp.is_semiperfect)))
    print(f"forge exit={rc}")
    print(json.dumps(collect(a.problem, a.run_dir), indent=1))
    return rc


if __name__ == "__main__":
    sys.exit(main())


# --- A12: 粘性 δ* 補正 (RANS 経路) ------------------------------------------------
def prepare_ns(problem_path, run_dir, nsteps=None, ic_from=None,
               dstar_csv=None) -> dict:
    r"""**物理壁 (inviscid + δ*) の RANS run** を準備する (A12)。

    plan: plans/active/tooling-nozzle-axismach-viscous-deltastar.md。

    - `dstar_csv` なし → v1: `feedback.deltastar.deltastar_offset` (相関) で法線オフセット
    - `dstar_csv` あり → v2: CSV (x_rt, dstar_rt) を補間して法線オフセット
      (CFD 抽出 δ* の固定点反復)
    - メッシュ/段階起動レシピは B8 系 NS v1 (run_0028-0030) で確立したものを流用。
      coarse 中継 (y+~50) は YAML の mesh/evaluate 設定だけの違いで同じ関数で作る
    """
    from ..feedback.deltastar import _sutherland
    from ..geometry.wall_axismach import PhysicalNozzleWall
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym_axismach":
        raise ValueError("runner_axismach は wind_tunnel_axisym_axismach 専用")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    d = design_chain(p)
    scale = float(p.spec["r_throat"])
    wall_inv = d["wall_inv"]                      # (n,4) [x,r,th,M] r_t 単位
    # 物理壁 (A13): 上流履歴込み δ* + 真のスロート探索 + 上流 Hermite 再生成。
    # v2 (dstar_csv) は「下流は CFD 抽出、x<6 は履歴相関、[6,9] smoothstep」の合成。
    dstar_x = None
    if dstar_csv is not None:
        tbl = np.loadtxt(dstar_csv, delimiter=",", skiprows=1)
        base = PhysicalNozzleWall(d["wall"], wall_inv, scale, float(p.spec["Pt"]),
                                  float(p.spec["Tt"]), _gam_or_gas(p), p.cp)

        def dstar_x(x, _tbl=tbl, _corr=base._dstar_hist):
            x = np.asarray(x, dtype=float)
            w = np.clip((x - 6.0) / 3.0, 0.0, 1.0)
            w = w * w * (3.0 - 2.0 * w)
            csv = np.interp(np.clip(x, _tbl[0, 0], _tbl[-1, 0]), _tbl[:, 0], _tbl[:, 1])
            return (1.0 - w) * _corr(x) + w * csv
    wall = PhysicalNozzleWall(d["wall"], wall_inv, scale, float(p.spec["Pt"]),
                              float(p.spec["Tt"]), _gam_or_gas(p), p.cp, dstar_x=dstar_x)
    dstar_src = "correlation_hist_v1" if dstar_csv is None else str(dstar_csv)
    msgs = wall.validate()
    if msgs:
        raise ValueError("物理壁フィルタ不合格: " + "; ".join(msgs))
    mp = Mesh2DParams(ni=int(p.mesh.get("ni", 561)), nj=int(p.mesh.get("nj", 97)),
                      wall_first_frac=float(p.mesh.get("wall_first_frac", 4.5e-5)),
                      throat_refine=float(p.mesh.get("throat_refine", 3.0)),
                      scale=scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    law = d["law"]
    xs = np.linspace(d["x0"], d["x_E"], 400)
    tgt = [float(law(x)) if x >= d["x_A"] else float("nan") for x in xs]
    np.savetxt(run_dir / "target_axis_M.csv", np.c_[xs * scale, tgt],
               delimiter=",", header="x_m,M_target", comments="")
    np.savetxt(run_dir / "wall_design.csv",
               np.c_[d["wall_inv"] * [scale, scale, 1.0, 1.0]], delimiter=",",
               header="x_m,r_m,theta_rad,M_wall", comments="")
    xs_p = np.linspace(wall.x_throat, wall.x_e, 1500)
    np.savetxt(run_dir / "wall_physical.csv", np.c_[xs_p * scale, wall.r(xs_p) * scale],
               delimiter=",", header="x_m,r_m", comments="")
    n = nsteps or int(p.evaluate.get("nStepOuter", 48000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    cfl_main = float(p.evaluate.get("cfl_main", 1.0))
    (run_dir / "bcondConfig.yaml").write_text(_bcond(p, euler=False))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    # 品質は cell 変換コピーで検査 (品質ツールは node CONNE 非対応)
    (run_dir / "solverConfig.yaml").write_text(
        _apply_gas_to_config(_config_euler(p, n, out_int, 4.0, 1), p, run_dir))
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
    # node/SST 変換 (config を先に書く — wall_dist は no-slip 壁で作られる)
    (run_dir / "solverConfig.yaml").write_text(
        _apply_gas_to_config(_config_sst_node(p, n, out_int, cfl_main), p, run_dir))
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp,
                        gas=(None if str(p.evaluate.get('cfd_gas', 'same')) == 'cpg' else p.gas_model))
    if ic_from is not None:
        src = sorted(Path(ic_from).glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(src), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
        import h5py as _h5
        # ソースが Euler (k/omega~0) のときだけ入口相当を貼る (runner_wt と同じ規約)
        with _h5.File(src) as fs:
            om_src_max = float(np.max(fs["/VALUE/omega"][:])) if "/VALUE/omega" in fs else 0.0
        if om_src_max < 1.0:
            with _h5.File(run_dir / "nozzle.h5", "r+") as f:
                ro = f["/VALUE/ro"][:]
                f["/VALUE/roK"][:] = ro * 1.0
                f["/VALUE/roOmega"][:] = ro * 18000.0
        # 近壁 omega の粘性底層フロア omega = 6 nu / (beta1 y^2) (run_0030 の教訓)
        g = p.gamma
        R_gas = p.cp * (g - 1.0) / g
        with _h5.File(run_dir / "nozzle.h5", "r+") as f:
            ro = f["/VALUE/ro"][:]
            roUx, roUy = f["/VALUE/roUx"][:], f["/VALUE/roUy"][:]
            roe = f["/VALUE/roe"][:]
            wd = f["/VALUE/wall_dist"][:]
            T = np.maximum((roe - 0.5 * (roUx ** 2 + roUy ** 2) / np.maximum(ro, 1e-12))
                           * (g - 1.0) / (np.maximum(ro, 1e-12) * R_gas), 50.0)
            nu = _sutherland(T) / np.maximum(ro, 1e-12)
            om_floor = 6.0 * nu / (0.075 * np.maximum(wd, 1e-9) ** 2)
            f["/VALUE/roOmega"][:] = np.maximum(f["/VALUE/roOmega"][:], ro * om_floor)
    info = {"chain": "axismach_ns", "discretization": "node", "viscous": True,
            "dstar_source": dstar_src, "mdot_ratio_moc": None,
            "throat_physical": {"x": wall.x_throat, "r": wall.r_throat,
                                "kappa": wall.kappa_throat,
                                "dstar_throat": float(wall._dstar_hist(0.0))},
            "x0": d["x0"], "x_A": d["x_A"], "x_E": d["x_E"], "L_c": d["L_c"],
            "anchor": list(d["anchor"]), "anchor_source": d["anchor_source"],
            "start_line": d["start_line"], "wall_mode": d["wall_mode"],
            "Md": d["Md"], "R": d["R"],
            "qa": {k: v for k, v in d["qa"].items() if k != "violations"},
            "nStepOuter": n, "cfl_main": cfl_main, "scale_m": scale,
            "ic_from": str(ic_from) if ic_from else None,
            "mesh": {"ni": mp.ni, "nj": mp.nj, "wall_first_frac": mp.wall_first_frac}}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def run_staged_ns(run_dir) -> int:
    """NS の 3 段起動 (run_0030 レシピ): soft (1次 cfl0.5 ni10) → mid (1次 cfl1)
    → 本段 (2次 cfl_main)。各段の最終場を IC に引き継ぐ。"""
    import re
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()

    def _stage(cfg, nsteps):
        cfg = re.sub(r"nStepOuter: \d+", f"nStepOuter: {nsteps}", cfg)
        cfg = re.sub(r"outStepInterval: \d+", f"outStepInterval: {nsteps}", cfg)
        (run_dir / "solverConfig.yaml").write_text(cfg)
        rc = run_forge(run_dir)
        res = sorted(run_dir.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res or int("".join(c for c in res[-1].stem if c.isdigit())) < nsteps:
            raise RuntimeError(f"段階起動が失敗 (rc={rc}, res={res[-1].name if res else None})")
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(res[-1]), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
        for f in run_dir.glob("res_*"):
            f.unlink()

    soft = cfg_main
    soft = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 0.5, cfl_pseudo: 0.5", soft)
    soft = soft.replace("convMethod: 1", "convMethod: 0")
    soft = soft.replace("nStepInner: 5", "nStepInner: 10")
    _stage(soft, 3000)
    mid = cfg_main
    mid = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 1.0, cfl_pseudo: 1.0", mid)
    mid = mid.replace("convMethod: 1", "convMethod: 0")
    mid = mid.replace("nStepInner: 5", "nStepInner: 10")
    _stage(mid, 3000)
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    return run_forge(run_dir)
