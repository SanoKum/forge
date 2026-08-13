"""v2: Euler 帰還ループ (親計画 §4.7 の 2 段目)。

1 パス = forge Euler 収束計算 → 軸 ΔM 評価 → **v1 特性線網の凍結 C⁻ マップ**
(壁点 → その波が着地する軸点) で残差を壁へ帰還 (PM 換算 Δθ_w = ω·Δν) →
壁再構築 (θ_w に平滑化增分を加え r=∫tanθ で再積分, 接合点アンカー) →
同一トポロジ再メッシュ → 前パス場から warm restart。

マスク (§4.7(c)): 壁足が円弧接合近傍に落ちる station (帰還で直せない —
接合波スパイクを含む) と目標終端近傍はランプで重みゼロ。
収束: masked ‖ΔM‖∞ ≤ tol (既定 0.5% Md) or パス上限。

使い方:
  design/.venv-opt/bin/python -m forge_design.feedback.euler_loop \
      case/41.wind_tunnel_design/problem_m4.yaml <work_dir> [--passes 8]
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator

from ..geometry.moc_kernel import pm_nu, pm_mach
from ..geometry.wall_modef import ModeFWall
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import load_problem
from .. import evaluate
from ..evaluate.ic import paste_isentropic_ic
from ..evaluate.runner import FORGE_BUILD, FORGE_TOOLS, PROBE_STUB, _ENV, run_forge
from ..evaluate import runner_wt


def build_cminus_map_cfd(run_dir, wall_pts, scale: float, gamma: float,
                         m_floor: float = 1.05) -> np.ndarray:
    """**v3**: そのパスの CFD 最終場の中で C⁻ を壁→軸へ積分してマップを作る。

    注: 「収束場」とは呼ばない — 収束判定は `evaluate.health` が
    `check_convergence.py` の VERDICT として別途記録し、読み替えはしない (A3)。

    v2 の凍結マップ (v1 の MOC 設計場= 完全気体・等エントロピー・単一組成) と違い、
    局所の流れ角 θ と マッハ角 μ を**実場から**読むため、粘性による有効形状の変化・
    組成分布・物性の温度依存が自動的に入る (ユーザ指摘 2026-08-14)。
    ガード: トレースが $M<m\_floor$ 域に入った壁点は棄却 (§4.7(c))。
    戻り: (n,2) [x_w, x_axis] (無次元 r*)。"""
    run_dir = Path(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    with h5py.File(run_dir / "nozzle.h5") as nz:
        cc = nz["CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res) as f:
        Ux, Uy, son = f["VALUE/Ux"][:], f["VALUE/Uy"][:], f["VALUE/sonic"][:]
    x, r = cc[:, 0] / scale, cc[:, 1] / scale          # → r* 単位
    M = np.hypot(Ux, Uy) / np.maximum(son, 1e-9)
    th = np.arctan2(Uy, Ux)
    itp_M = LinearNDInterpolator(np.c_[x, r], M)
    itp_th = LinearNDInterpolator(np.c_[x, r], th)
    # **cell モード限定の回避策**: cell 中心は軸上に無いので補間の凸包が r_min で
    # 切れ、軸到達前に必ず NaN になり全点棄却される (2026-08-14 実測)。そこまで
    # 積分して最後の C⁻ 勾配で r=0 へ外挿する (θ→0 で勾配は緩やかなので誤差は小)。
    # **node (median-dual) 化したらこの外挿は不要** — node は軸上に DOF を持つので
    # 凸包が r=0 まで届く (ユーザ方針: 今後は node 主体)。r_axis=0 で自然に動く。
    r_axis = float(r.min()) * 1.5
    rows = []
    for xw, rw in wall_pts[::4, :2]:
        xx, rr = float(xw), float(rw) * 0.985          # 壁のわずか内側から出発
        ok, slope = True, None
        for _ in range(4000):
            if rr <= r_axis:
                break
            Mv, tv = float(itp_M(xx, rr)), float(itp_th(xx, rr))
            if not (np.isfinite(Mv) and np.isfinite(tv)) or Mv < m_floor:
                ok = False
                break
            mu = np.arcsin(1.0 / max(Mv, 1.0 + 1e-9))
            dr = -min(0.02, rr - r_axis)
            slope = np.tan(tv - mu)
            if abs(slope) < 1e-9:
                ok = False
                break
            xx += dr / slope
            rr += dr
        if ok and slope is not None and rr <= r_axis * 1.01:
            rows.append((float(xw), xx - rr / slope))   # r=0 へ外挿
    return np.asarray(rows)


def build_cminus_map(pts, wall_inv, gamma: float) -> np.ndarray:
    """**v2**: v1 の MOC 設計場で各壁点から C⁻ を軸まで前進積分し (x_w → x_axis)
    マップを作る (以後のパスで凍結して使う)。"""
    xy = np.array([[p.x, p.r] for p in pts])
    th = np.array([p.th for p in pts])
    nu = np.array([p.nu for p in pts])
    itp_th = LinearNDInterpolator(xy, th)
    itp_nu = LinearNDInterpolator(xy, nu)
    rows = []
    for xw, rw in wall_inv[::4, :2]:
        x, r = float(xw), float(rw)
        ok = True
        for _ in range(4000):
            t = itp_th(x, r)
            v = itp_nu(x, r)
            if not (np.isfinite(t) and np.isfinite(v)):
                ok = False
                break
            M = pm_mach(max(float(v), 1e-9))
            mu = np.arcsin(1.0 / max(M, 1.0 + 1e-9))
            dr = -min(0.02, r)
            x += dr / np.tan(float(t) - mu)
            r += dr
            if r <= 1e-9:
                break
        if ok and r <= 1e-6:
            rows.append((float(xw), x))
    return np.asarray(rows)


def _rebuild_wall(wall_pts: np.ndarray, dtheta_fn, x_anchor: float) -> np.ndarray:
    """現壁 (n,2)[x,r] の壁角に Δθ(x) を加えて r=∫tanθ で再構築 (接合点固定)。"""
    x, r = wall_pts[:, 0], wall_pts[:, 1]
    th = np.arctan(np.gradient(r, x))
    th_new = th + dtheta_fn(x)
    r_new = np.empty_like(r)
    i0 = int(np.argmin(np.abs(x - x_anchor)))
    r_new[i0] = r[i0]
    for i in range(i0 + 1, len(x)):
        r_new[i] = r_new[i - 1] + 0.5 * (np.tan(th_new[i - 1]) + np.tan(th_new[i])) * (x[i] - x[i - 1])
    for i in range(i0 - 1, -1, -1):
        r_new[i] = r_new[i + 1] - 0.5 * (np.tan(th_new[i]) + np.tan(th_new[i + 1])) * (x[i + 1] - x[i])
    return np.c_[x, r_new]


def run_loop(problem_path, work_dir, n_pass: int = 8, omega: float = 0.4,
             tol_rel: float = 0.005, clip_deg: float = 0.5, nsteps=None,
             map_mode: str = "frozen", gate_mode: str = "formal") -> dict:
    """map_mode: "frozen" = v2 (v1 MOC 設計場の凍結マップ) /
    "cfd" = **v3** (毎パスその回の CFD 収束場から C⁻ を引き直す)。"""
    p = load_problem(problem_path)
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    d = runner_wt.design_chain(p, viscous=False)
    Md, xd, x0 = d["Md"], d["xd"], d["x0"]
    scale = float(p.spec["r_throat"])
    cmap = build_cminus_map(d["pts"], d["wall_inv"], p.gamma)
    x_w0 = float(d["wall_inv"][0, 0])

    def make_ax2w(cm):
        o = np.argsort(cm[:, 1])
        return lambda xi: np.interp(xi, cm[o, 1], cm[o, 0])

    ax2w = make_ax2w(cmap)   # v3 では各パスの CFD 場で毎回置き換える

    tgt = np.array([[x, d["target"](x)] for x in np.linspace(x0, xd, 500)])
    wall_pts = d["wall_inv"][:, :2].copy()
    prev_res = None
    hist = []
    # trust-region 型の決定的下り探索: 悪化したら最良壁に巻き戻し ω 半減
    best = {"wall": wall_pts.copy(), "dM_inf": np.inf, "resid": None,
            "ax2w": ax2w}   # C1: マップも最良状態と一緒に保持し reject 時に戻す
    for k in range(1, n_pass + 1):
        rd = work / f"pass_{k:02d}"
        if not (rd / "metrics_pass.json").exists():
            rd.mkdir(parents=True, exist_ok=True)
            wall = ModeFWall(wall_pts, R=d["R"],
                             r_inlet=float(p.geometry.get("r_inlet", 2.5)),
                             L_pipe=float(p.geometry.get("L_pipe", 0.5)),
                             L_contract=float(p.geometry.get("L_contract", 3.0)),
                             phi_u=np.deg2rad(float(p.geometry.get("phi_u_deg", 25.0))),
                     blend_len=float(p.geometry.get("blend_len", 0.0)))
            mp = Mesh2DParams(ni=int(p.mesh.get("ni", 321)), nj=int(p.mesh.get("nj", 65)),
                              wall_first_frac=float(p.mesh.get("wall_first_frac", 5.0e-3)),
                              throat_refine=float(p.mesh.get("throat_refine", 3.0)),
                              scale=scale)
            coords, quads, bedges = generate_axisym_mesh(wall, mp)
            write_msh41_2d(rd / "nozzle.msh", coords, quads, bedges)
            n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
            (rd / "solverConfig.yaml").write_text(
                runner_wt._config_euler(p, n, max(n // 3, 1), 4.0, 1))
            (rd / "bcondConfig.yaml").write_text(runner_wt._bcond(p, euler=True))
            (rd / "probe.yaml").write_text(PROBE_STUB)
            subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                           cwd=rd, env=_ENV, check=True, capture_output=True, text=True)
            paste_isentropic_ic(rd / "nozzle.h5", wall, scale,
                                float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp)
            if prev_res is not None:
                subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                                str(prev_res), str(rd / "nozzle.h5")],
                               env=_ENV, check=True, capture_output=True, text=True)
                run_forge(rd)
            else:
                runner_wt.run_staged(problem_path, rd, nsteps)
        res = sorted(rd.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        prev_res = res[-1]
        if map_mode == "cfd":
            cm = build_cminus_map_cfd(rd, wall_pts, scale, p.gamma)
            if len(cm) >= 10:
                ax2w = make_ax2w(cm)
                np.savetxt(rd / "cminus_map.csv", cm, delimiter=",",
                           header="x_wall,x_axis_landing", comments="")
            else:
                print(f"[pass {k}] CFD マップ生成失敗 ({len(cm)} 点) — 前回の表を継続",
                      flush=True)
        # 軸 ΔM 評価
        from ..metrics.extract import axis_mach
        x_ach, M_ach = axis_mach(rd / "nozzle.h5", res[-1], axis_band=0.09 * scale)
        xi = np.linspace(x0 + 0.3, xd - 0.2, 240)
        Mt = np.interp(xi, tgt[:, 0], tgt[:, 1])
        Ma = np.interp(xi * scale, x_ach, M_ach)
        dM = Mt - Ma
        # マスク/重みランプ: 壁足が接合近傍 or 目標終端近傍
        xw = ax2w(xi)
        w = np.clip((xw - (x_w0 + 0.3)) / 1.0, 0.0, 1.0) * np.clip((xd - 0.5 - xi) / 1.0, 0.0, 1.0)
        w *= (Mt > 1.15)
        dM_inf = float(np.max(np.abs(dM * (w > 0.5))))
        rejected = bool(dM_inf > best["dM_inf"] * 1.02)
        # 出口一様性は**診断として**毎パス記録する (安価)。役割の区別 (2026-08-14):
        #   本ループの目的 = ΔM→0 = 「宣言した目標 M_c(x) を実際に達成する壁を作る」。
        #     これは出口品質のためでなく**評価の忠実性**のため — dv が「達成された
        #     目標分布」を意味して初めて dv→メトリクス写像が決定的になり、
        #     サロゲートが劣化しない (親計画 §4.7(c))。
        #   ε_M/ε_θ = **最適化ループ (§5 ステップ 7) の目的** — Bézier CP・L・R を
        #     振って最小化する量であり、1 評価内で追うものではない。したがって
        #     ΔM と ε が単調に一致しなくても (実測: pass9 ε_M=0.020% < pass13 0.039%)
        #     本ループの欠陥ではない。
        from ..metrics.extract import exit_uniformity
        from ..evaluate.health import gate
        # A3: **正式ゲート** — PASS 以外では壁もマップも更新しない (未収束誤差で
        # 帰還を汚さない)。探索目的なら gate_mode="exploratory" を明示し、正式
        # 評価とは別の run ディレクトリに分離すること。
        hl = gate(rd, M_design=Md, x_d=xd * scale, mode=gate_mode)
        try:
            eu = exit_uniformity(rd / "nozzle.h5", res[-1], Md,
                                 x_d=xd * scale)   # 有効菱形コア (実寸)
        except Exception as e:  # 測定面が取れない等
            eu = {"error": f"{type(e).__name__}: {e}"}
        hist.append({"pass": k, "dM_inf_masked": dM_inf, "omega": omega,
                     "rejected": rejected,
                     "dM_rms": float(np.sqrt(np.mean((dM * w) ** 2))),
                     "eps_M_rms": eu.get("eps_M_rms"),
                     "eps_theta_max_deg": eu.get("eps_theta_max_deg"),
                     "convergence_verdict": hl["convergence_verdict"],
                     "quasisteady_verdict": hl["quasisteady_verdict"],
                     "eps_series_verdict": hl.get("eps_steady", {}).get("verdict"),
                     "gate_ok": hl["ok"], "gate_reasons": hl["gate_reasons"]})
        print(f"[pass {k}] masked ‖ΔM‖∞ = {dM_inf:.4f} ({dM_inf/Md*100:.2f}% Md)"
              f"  εM={('%.4f%%' % (eu['eps_M_rms']*100)) if eu.get('eps_M_rms') is not None else 'n/a'}"
              f"  εθ={('%.3f°' % eu['eps_theta_max_deg']) if eu.get('eps_theta_max_deg') is not None else 'n/a'}"
              f"  [{hl['convergence_verdict']}/{hl['quasisteady_verdict']}]"
              f"{' [reject→ω/2]' if rejected else ''}", flush=True)
        np.savetxt(rd / "achieved_vs_target.csv", np.c_[xi, Mt, Ma, w],
                   delimiter=",", header="x_rt,M_target,M_achieved,weight", comments="")
        (rd / "metrics_pass.json").write_text(json.dumps(hist[-1]))
        if not hl["ok"]:
            print(f"[pass {k}] **中止**: ゲート不合格 ({'; '.join(hl['gate_reasons'])}) "
                  f"— 未収束/非定常の場から壁とマップを作らない", flush=True)
            break
        if rejected:
            omega *= 0.5           # 過ゲイン → 最良壁から小さい ω でやり直し
            xi_c, Mt_c, Ma_c, w_c = best["resid"]
            ax2w = best["ax2w"]    # C1: 棄却場のマップで最良場の残差を戻さない
        else:
            best = {"wall": wall_pts.copy(), "dM_inf": dM_inf,
                    "resid": (xi, Mt, Ma, w), "ax2w": ax2w}
            if dM_inf <= tol_rel * Md:
                break
            xi_c, Mt_c, Ma_c, w_c = xi, Mt, Ma, w
        # 帰還: Δν → Δθ_w (壁座標へ写像し平滑化・クリップ)
        dnu = np.array([pm_nu(max(float(m1), 1.001)) - pm_nu(max(float(m2), 1.001))
                        for m1, m2 in zip(Mt_c, Ma_c)])
        dth_i = np.clip(omega * dnu * w_c, -np.deg2rad(clip_deg), np.deg2rad(clip_deg))
        xw_s, o2 = np.unique(ax2w(xi_c), return_index=True)
        dth_w = dth_i[o2]
        for _ in range(3):
            dth_w = np.convolve(dth_w, np.ones(9) / 9.0, mode="same")

        def dtheta_fn(x):
            return np.interp(x, xw_s, dth_w, left=0.0, right=dth_w[-1])

        wall_pts = _rebuild_wall(best["wall"], dtheta_fn, x_anchor=x_w0)
        np.savetxt(rd / "wall_next.csv", wall_pts, delimiter=",",
                   header="x_rt,r_rt", comments="")
    out = {"history": hist, "passes": len(hist),
           "converged": bool(hist and hist[-1]["dM_inf_masked"] <= tol_rel * Md),
           "tol": tol_rel * Md}
    (work / "loop_summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1), flush=True)
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="v2 Euler 帰還ループ")
    ap.add_argument("problem")
    ap.add_argument("work_dir")
    ap.add_argument("--passes", type=int, default=8)
    ap.add_argument("--omega", type=float, default=0.4)
    ap.add_argument("--gate-mode", choices=("formal", "exploratory"), default="formal",
                    help="formal = PASS 以外は停止 (正式) / exploratory = 発散以外は継続")
    ap.add_argument("--map-mode", choices=("frozen", "cfd"), default="frozen",
                    help="frozen = v2 (v1 MOC 場の凍結マップ) / cfd = v3 (毎パス実場)")
    ap.add_argument("--steps", type=int, default=None)
    a = ap.parse_args(argv)
    run_loop(a.problem, a.work_dir, n_pass=a.passes, omega=a.omega, nsteps=a.steps,
             map_mode=a.map_mode, gate_mode=a.gate_mode)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
