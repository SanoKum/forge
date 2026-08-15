#!/usr/bin/env python3
r"""axis-Mach チェーンの出口位置 $x_F$ を終端特性線から確定する (ユーザ指摘 2026-08-15)。

現行 `inverse_design` は壁流線を「$\theta$ が最大値を経て 0.05° を切った点」で切って
いる (角度しきい値)。本スクリプトは E $(x_E,0)$ から出る**終端 C⁺ 特性線**を MOC 場の
中で追跡し、壁流線との交点として $x_F$ を決め、3 つの定義を比較する:

  (a) lip-cut    : 現行 (θ ≤ 0.05°)
  (b) char-trace : 終端 C⁺ を MOC 場で追跡 → 壁との交点 (**正しい定義**)
  (c) straight   : 一様域の直線 Mach line 近似 $x_F=x_E+r_F\sqrt{M_d^2-1}$

さらに軸 target 延長 `x_end_margin` に対する $x_F$ の感度を測る (延長が足りていれば
$x_F$ は margin に依存しないはず = 「長さは MOC が決めている」ことの検証)。

出力: `axismach_exit_analysis.json` と `axismach_field.npz` (可視化用の場・壁・特性線)。
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))

from forge_design.geometry.axis_law import QuinticHermiteAxisLaw  # noqa: E402
from forge_design.geometry.moc_inverse import InverseMOC, _Pt, pm_nu  # noqa: E402
from forge_design.geometry.transonic import HallThroat  # noqa: E402
from forge_design.geometry.wall_axismach import area_ratio_isentropic  # noqa: E402

G, MD, R = 1.4, 4.0, 2.0
# pass3 (run_0047) の設定
X_REACH, ANCHOR, L_C = 1.7839, (2.07549, 0.64351, -0.04266), 7.8465


def build_field(n_axis: int, n_start: int, dx_wall: float, margin: float,
                seg_run: str | None = None):
    """pass3 相当の設計を再構築し、(pts, 未カット壁流線, law, ht) を返す。"""
    ht = HallThroat(R=R, gamma=G)
    x0 = ht.x_axis_of_mach(1.05)
    law = QuinticHermiteAxisLaw(X_REACH, L_C, *ANCHOR, MD)
    seg = None
    if seg_run:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
        from forge_design.evaluate.runner_axismach import axis_curve_node
        seg, _, _ = axis_curve_node(seg_run, 0.01, lam=1e-5)

    def target(x):
        x = float(x)
        if x < X_REACH:
            return float(seg(np.float64(x))) if seg is not None else float(ht.mach(x, 0.0))
        return float(law(x))

    rF_pred = float(np.sqrt(area_ratio_isentropic(MD, G)))
    x_end = law.x_E + margin * rF_pred * float(np.sqrt(MD * MD - 1.0))
    inv = InverseMOC(gamma=G, delta=1.0)
    xs0, rr, MM, tt = ht.starting_line(M_start=1.05, n=n_start)
    ax = np.linspace(xs0, x_end, n_axis)[1:][::-1]
    init = [_Pt(float(x), 0.0, 0.0, float(pm_nu(float(target(x)), G)), G) for x in ax]
    init += [_Pt(float(xs0), float(rr[i]), float(tt[i]), float(pm_nu(float(MM[i]), G)), G)
             for i in range(n_start)]
    init[-1].th = float(np.arctan(xs0 / R))
    pts = inv.fill(init)
    wall = inv.wall_streamline(pts, xs0, float(rr[-1]), dx=dx_wall)  # 未カット
    return pts, wall, law, ht, float(x_end)


def lip_cut_index(wall) -> int:
    """現行 `inverse_design` と同じ θ≤0.05° カットの index。"""
    i_pk = int(np.argmax(wall[:, 2]))
    after = wall[i_pk:, 2] <= np.deg2rad(0.05)
    return i_pk + int(np.argmax(after)) if np.any(after) else len(wall) - 1


def trace_cplus(pts, x_start: float, r_start: float, wall, ds: float = 2e-3,
                n_max: int = 400000) -> dict:
    r"""$(x_s,r_s)$ から C⁺ ($dr/dx=\tan(\theta+\mu)$) を MOC 場で追跡し壁と交わる点を返す。

    場の外に出たら**一様出口状態で凍結して直線延長**する (終端特性線の下流側は定義上
    一様域なので、これは近似ではなく構成そのもの)。戻り値に軌跡も含める。
    """
    from scipy.interpolate import LinearNDInterpolator
    from forge_design.geometry.moc_inverse import pm_mach
    xy = np.array([[p.x, p.r] for p in pts])
    itp_th = LinearNDInterpolator(xy, np.array([p.th for p in pts]))
    itp_nu = LinearNDInterpolator(xy, np.array([p.nu for p in pts]))
    wx, wr = wall[:, 0], wall[:, 1]
    mu_d = float(np.arcsin(1.0 / MD))
    x, r = float(x_start), float(r_start)
    path = [(x, r)]
    n_frozen = 0
    for _ in range(n_max):
        th = float(itp_th(x, r))
        nu = float(itp_nu(x, r))
        if np.isfinite(th) and np.isfinite(nu):
            M = pm_mach(max(nu, 1e-9), G)
            slope = np.tan(th + np.arcsin(1.0 / max(M, 1.0 + 1e-12)))
        else:
            slope = np.tan(mu_d)      # 一様域 (θ=0, M=Md) — 凍結延長
            n_frozen += 1
        x_n, r_n = x + ds, r + ds * slope
        rw_n = float(np.interp(x_n, wx, wr))
        if r_n >= rw_n:               # 壁を越えた — 線形補間で交点
            rw = float(np.interp(x, wx, wr))
            t = (rw - r) / max((r_n - r) - (rw_n - rw), 1e-30)
            xf = x + t * ds
            path.append((xf, float(np.interp(xf, wx, wr))))
            return {"x_F": xf, "r_F": float(np.interp(xf, wx, wr)),
                    "hit_wall": True, "n_frozen": n_frozen,
                    "path": np.asarray(path)}
        x, r = x_n, r_n
        path.append((x, r))
        if x > wx[-1]:
            break
    return {"x_F": float("nan"), "r_F": float("nan"), "hit_wall": False,
            "n_frozen": n_frozen, "path": np.asarray(path)}


def main() -> int:
    out: dict = {"config": {"Md": MD, "R": R, "x_reach": X_REACH,
                            "anchor": list(ANCHOR), "L_c": L_C}}
    t0 = time.time()
    # --- 生産解像度で本解析 ---
    pts, wall, law, ht, x_end = build_field(2000, 121, 0.005, 2.3)
    i_lip = lip_cut_index(wall)
    rF_1d = float(np.sqrt(area_ratio_isentropic(MD, G)))
    tc = trace_cplus(pts, law.x_E, 0.0, wall)
    out["x_E"] = float(law.x_E)
    out["x_axis_end"] = x_end
    out["wall_streamline_end"] = float(wall[-1, 0])
    out["definitions"] = {
        "lip_cut": {"x_F": float(wall[i_lip, 0]), "r_F": float(wall[i_lip, 1]),
                    "theta_deg": float(np.rad2deg(wall[i_lip, 2])),
                    "M_wall": float(wall[i_lip, 3])},
        "char_trace": {"x_F": tc["x_F"], "r_F": tc["r_F"],
                       "hit_wall": tc["hit_wall"], "n_frozen_steps": tc["n_frozen"],
                       "theta_deg": float(np.rad2deg(np.interp(tc["x_F"], wall[:, 0], wall[:, 2]))),
                       "M_wall": float(np.interp(tc["x_F"], wall[:, 0], wall[:, 3]))},
        "straight_machline": {"x_F": float(law.x_E + rF_1d * np.sqrt(MD * MD - 1.0)),
                              "r_F_1d": rF_1d},
    }
    print(json.dumps(out["definitions"], indent=1))
    print(f"[{time.time()-t0:.0f}s] main field done")

    # --- 感度: 軸 target 延長 margin を変えても x_F が動かないか (低解像で可) ---
    sens = {}
    for margin in (1.6, 2.3, 3.0):
        p2, w2, l2, _, xe2 = build_field(500, 41, 0.02, margin)
        t2 = trace_cplus(p2, l2.x_E, 0.0, w2)
        sens[f"margin={margin}"] = {
            "x_axis_end": xe2, "x_F_char": t2["x_F"],
            "x_F_lip": float(w2[lip_cut_index(w2), 0]),
            "wall_end": float(w2[-1, 0])}
        print(f"[{time.time()-t0:.0f}s] margin={margin}: {sens[f'margin={margin}']}")
    out["x_end_margin_sensitivity"] = sens

    here = Path(__file__).parent
    (here / "axismach_exit_analysis.json").write_text(json.dumps(out, indent=1))
    np.savez_compressed(
        here / "axismach_field.npz",
        pts=np.array([[p.x, p.r, p.th, p.M] for p in pts]),
        wall=wall, i_lip=i_lip, term_path=tc["path"],
        x_E=law.x_E, x_F_char=tc["x_F"],
        axis_x=np.linspace(ht.x_axis_of_mach(1.05), law.x_E, 600),
        axis_M=np.array([float(law(x)) if x >= X_REACH else float("nan")
                         for x in np.linspace(ht.x_axis_of_mach(1.05), law.x_E, 600)]))
    print(f"saved. total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
