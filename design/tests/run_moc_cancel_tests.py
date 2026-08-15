#!/usr/bin/env python3
"""wave-cancellation MOC (Phase W4, moc_cancel.py) の解析検証。

**検証範囲 (2026-08-15)**: `DataLine` 構成と `build_field` (interior()+軸反射)
のみ。壁生成 (`march_wall_from_dataline`/`cancellation_contour`/`_wall_cancel`) は
放射源流厳密解照合で有意な乖離が残り**未検証** — moc_cancel.py モジュール
docstring の「状態」節に 3 回の実装試行の記録がある。ここではその検証済み
部分のみを恒久テストにし、壁生成は別途(`plans/.../` W4 節)で継続調査する。

plan §7 (i) の「放射源流厳密解 <0.1%」は `build_field` の各点個別評価で満たす。
(ii) 格子収束も `build_field` の点数で確認する。
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.evaluate.ic import invert_area_ratio  # noqa: E402
from forge_design.geometry.moc_cancel import DataLine, build_field  # noqa: E402
from forge_design.geometry.moc_kernel import pm_mach, trace_cminus  # noqa: E402

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

# --- build_field: 放射源流厳密解照合 (各点個別に厳密解と一致) --------------------
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

# --- 格子収束: データ線の離散点数を変えても build_field の被覆域の厳密解一致度は
#     劣化しない (点数を増やしても同水準の精度を保つか) --------------------------
dl2 = make_dataline(thc_deg, x_w=1.10, n=81)
all_pts2 = build_field(dl2, n_rounds=25)
xs2 = np.array([p.x for p in all_pts2])
rs2 = np.array([p.r for p in all_pts2])
th2 = np.array([p.th for p in all_pts2])
M2 = np.array([p.M for p in all_pts2])
M_exact2 = np.array([_fM(x, r) for x, r in zip(xs2, rs2)])
M_err2 = np.max(np.abs(M2 - M_exact2) / M_exact2)
check(f"格子収束: n=41→81 でも build_field の厳密解一致度が劣化しない "
      f"(n=41: {M_err:.2e} → n=81: {M_err2:.2e})", M_err2 < 1e-3)

print("-" * 60)
print(f"{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
print()
print("注意: 壁生成 (march_wall_from_dataline/cancellation_contour) は未検証。")
print("moc_cancel.py モジュール docstring の「状態」節と")
print("plans/active/tooling-nozzle-walldriven-chain.md の W4 節を参照。")
sys.exit(1 if FAIL else 0)
