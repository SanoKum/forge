#!/usr/bin/env python3
"""wave-cancellation MOC (Phase W4, moc_cancel.py) の解析検証。

**検証範囲 (2026-08-15, 全面改訂)**: `DataLine`/`build_field` (内点充填) に加え、
壁生成 (`_wall_cancel`/`march_wall_from_dataline`/`cancellation_contour`) も
今回の訂正 (壁の閉じ条件に目標出口一様状態 $M_{\rm target}$ を明示的に使う —
moc_cancel.py モジュール docstring の「閉じ方」節) 後は検証済み:

  (i) 放射源流厳密解: `build_field` の各点が厳密解と個別に一致 (<0.1%)。
  (ii) 格子収束: データ線点数を変えても `build_field`/壁とも座標が安定。
  (iii) 壁専用単位過程 `_wall_cancel` の退化検証: 平面単純波 (δ=0、$J^-$ 厳密
       一定) で、隣接内点 A が既に目標と同じ $J^-$ を持つ場合、壁が A を
       厳密に再現する (機械精度)。
  (iv) march 全体の物理的整合性: 放射源流データ線 + 任意の $M_{\rm target}$ で
       march した壁が流線条件 $dr/dx\approx\tan\theta$ を満たし、θ が単調
       減衰し、格子収束する。θ=0 交差の安全策 (振動域回避) も確認する。

(v) 自己一貫性 (B8 場との照合) は実 CFD run を要するため
    `case/41.wind_tunnel_design/` 側で別途確認する。
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.evaluate.ic import invert_area_ratio  # noqa: E402
from forge_design.geometry.moc_cancel import (  # noqa: E402
    DataLine, _CancelKernel, build_field, march_wall_from_dataline)
from forge_design.geometry.moc_kernel import (  # noqa: E402
    _Pt, _mu, pm_mach, pm_nu, trace_cminus)

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


g = 1.4


def _fM(x, r):
    R = max(np.hypot(x, r), 1.0001)
    return float(invert_area_ratio(np.array([R ** 2]), np.array([True]), g)[0])


def _fth(x, r):
    return float(np.arctan2(r, x))


def make_dataline(thc_deg, x_w, n):
    """半角 thc の円錐上の点から出る C⁻ を厳密場でトレースし DataLine 化する。"""
    thc = np.deg2rad(thc_deg)
    r_w = x_w * np.tan(thc)
    trace = trace_cminus(_fM, _fth, x_w, r_w, n, g)  # [(x,r,th,nu)] 軸→壁順
    xs = np.array([p[0] for p in trace])
    rs = np.array([p[1] for p in trace])
    ths = np.array([p[2] for p in trace])
    Ms = np.array([pm_mach(p[3], g) for p in trace])
    return DataLine(x=xs, r=rs, theta=ths, M=Ms, gamma=g, p0_diag={},
                    x_j=x_w, r_j=r_w, x_axis=float(xs[0]))


# --- DataLine 構成の健全性 ------------------------------------------------------
thc_deg = 10.0
dl = make_dataline(thc_deg, x_w=1.10, n=41)
check("データ線 軸端 θ=0", abs(dl.theta[0]) < 1e-9)
check("データ線 壁端 θ=thc",
      abs(dl.theta[-1] - np.deg2rad(thc_deg)) < 1e-6)
check("データ線 軸端は壁端より下流 (C⁻ の浅い角度で軸へ着地)",
      dl.x[0] > dl.x[-1])

# --- (i) build_field: 放射源流厳密解照合 (各点個別に厳密解と一致) ---------------
all_pts = build_field(dl, n_rounds=25)
check("build_field が十分な点数を生成", len(all_pts) > 3 * len(dl.x))

xs_p = np.array([p.x for p in all_pts])
rs_p = np.array([p.r for p in all_pts])
th_p = np.array([p.th for p in all_pts])
M_p = np.array([p.M for p in all_pts])
M_exact = np.array([_fM(x, r) for x, r in zip(xs_p, rs_p)])
th_exact = np.array([_fth(x, r) for x, r in zip(xs_p, rs_p)])
M_err = np.max(np.abs(M_p - M_exact) / M_exact)
th_err = np.max(np.abs(th_p - th_exact))
check(f"build_field: 全点が厳密解と個別に一致 (M max rel err={M_err:.2e} <0.1%)",
      M_err < 1e-3)
check(f"build_field: 全点が厳密解と個別に一致 (θ max|Δ|={np.rad2deg(th_err):.4f}° <0.5°)",
      th_err < np.deg2rad(0.5))
check("build_field: 全点が r≥0", bool(np.all(rs_p >= -1e-9)))

# --- (ii) 格子収束: build_field の厳密解一致度が劣化しない ----------------------
dl2 = make_dataline(thc_deg, x_w=1.10, n=81)
all_pts2 = build_field(dl2, n_rounds=25)
xs2 = np.array([p.x for p in all_pts2])
rs2 = np.array([p.r for p in all_pts2])
M2 = np.array([p.M for p in all_pts2])
M_exact2 = np.array([_fM(x, r) for x, r in zip(xs2, rs2)])
M_err2 = np.max(np.abs(M2 - M_exact2) / M_exact2)
check(f"格子収束: n=41→81 でも build_field の厳密解一致度が劣化しない "
      f"(n=41: {M_err:.2e} → n=81: {M_err2:.2e})", M_err2 < 1e-3)

# --- (iii) _wall_cancel 退化検証: 平面単純波 (δ=0, J⁻ 厳密一定) -----------------
M_target = 1.95
K0 = pm_nu(M_target, g)


def pt_on_simplewave(theta, t):
    nu = K0 - theta
    M = pm_mach(nu, g)
    mu = _mu(M)
    beta = theta - mu
    return _Pt(t * np.cos(beta), t * np.sin(beta), theta, nu, g)


kmoc_planar = _CancelKernel(gamma=g, delta=0.0, n_corr=3)
A = pt_on_simplewave(np.deg2rad(3.0), 2.0)
W_prev = pt_on_simplewave(np.deg2rad(8.0), 2.3)
W = kmoc_planar._wall_cancel(A, W_prev, K0)
check(f"_wall_cancel 退化検証 (平面単純波, A が既に目標 J⁻ を持つ場合): "
      f"θ_W={np.rad2deg(W.th):.6f}° (期待 3.000000°)",
      abs(np.rad2deg(W.th) - 3.0) < 1e-6)
check(f"_wall_cancel 退化検証: M_W={W.M:.6f} (期待 {A.M:.6f})",
      abs(W.M - A.M) < 1e-6)

# --- (iv) march 全体: 放射源流データ線 + M_target で物理的整合性 ----------------
res = march_wall_from_dataline(dl, M_target=M_target, theta_tol_deg=0.2,
                               eps_M_rel=0.01, max_iter=400)
wall = res["wall"]
check("march: theta_tol で正常終了 (max_iter で打ち切りではない)",
      res["terminated_by"] == "theta_tol")
check("march: 壁が下流へ単調 (x 単調増加)",
      bool(np.all(np.diff(wall[:, 0]) > 0)))
check("march: θ が単調非増加 (振動なし)",
      bool(np.all(np.diff(wall[:, 2]) <= 1e-6)))
check(f"march: 終端 M が目標に到達 (M={wall[-1,3]:.4f}, target={M_target})",
      abs(wall[-1, 3] - M_target) < 0.05 * M_target)

dx = np.diff(wall[:, 0])
dr = np.diff(wall[:, 1])
th_avg = 0.5 * (wall[:-1, 2] + wall[1:, 2])
mask = np.abs(dx) > 1e-9
stream_err = np.max(np.abs(dr[mask] / dx[mask] - np.tan(th_avg[mask])))
check(f"march: 壁が流線条件 dr/dx≈tanθ を満たす (max|Δ|={stream_err:.2e})",
      stream_err < 1e-4)

# 格子収束 (march 全体)
res2 = march_wall_from_dataline(dl2, M_target=M_target, theta_tol_deg=0.2,
                                eps_M_rel=0.01, max_iter=800)
w2 = res2["wall"]
check("march 格子収束: n=41→81 で終端壁座標の変化 <1%",
      abs(w2[-1, 0] - wall[-1, 0]) / wall[-1, 0] < 0.01
      and abs(w2[-1, 1] - wall[-1, 1]) / wall[-1, 1] < 0.01)

# θ=0 交差の安全策: theta_tol を意図的に厳しく (ステップ幅未満) 指定しても
# 振動域に入らず綺麗に θ=0, M=M_target で終わることを確認する
res3 = march_wall_from_dataline(dl, M_target=M_target, theta_tol_deg=0.05,
                                eps_M_rel=0.01, max_iter=400)
check("θ=0 交差の安全策: 厳しすぎる theta_tol でも振動域に入らず正常終了",
      res3["terminated_by"] == "theta_tol"
      and abs(res3["wall"][-1, 2]) < 1e-9
      and abs(res3["wall"][-1, 3] - M_target) < 1e-4)

print("-" * 60)
print(f"{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
