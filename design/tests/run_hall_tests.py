#!/usr/bin/env python3
"""axis-Mach チェーン A2: HallThroat (Hall–KL 軸対称遷音速) の unit test。

検証対象 (plans/active/tooling-nozzle-axismach-chain.md §5.2):
  1. スロート壁 u(0,1) = CONTUR 本文 Eq. (8) (級数の機械精度一致)
  2. v(0,1) = 0 (壁で v=0 — V 係数表の厳密恒等式)
  3. du/dx(0,1) = Eq. (9) (数値微分 vs 解析級数)
  4. R→∞ で Sauer 一次解へ退化 (M 場の差が O(1/S²) スケーリング)
  5. C_D: 数値積分 (級数場+厳密等エントロピー密度) vs Eq. (14) 級数
  6. axis_anchor の解析微分 vs 中心差分 (機械精度級)
  7. x_axis_of_mach の自己整合・軸ソニック点がスロート下流
  8. (診断) 壁 θ vs 骨接放物線接線 — Sauer の既知ギャップ (B0/B5) が縮むこと
"""
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.geometry.transonic import HallThroat, SauerThroat  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


# --- 1-3. CONTUR 本文式との一致 (γ, R を数点) -----------------------------------
for g, R in ((1.4, 2.0), (1.4, 5.0), (1.25, 2.0), (1.67, 3.0)):
    ht = HallThroat(R=R, gamma=g)
    S = R + 1.0
    u_eq8 = (1.0 + 1.0 / (4.0 * S) - (14.0 * g - 57.0) / (288.0 * S * S)
             + (2364.0 * g * g - 3915.0 * g + 14337.0) / (82944.0 * S ** 3))
    check(f"γ={g} R={R}: u(0,1) = Eq.(8)  ({float(ht.u(0.0, 1.0)):.8f})",
          abs(float(ht.u(0.0, 1.0)) - u_eq8) < 1e-13)
    check(f"γ={g} R={R}: v(0,1) = 0", abs(float(ht.v(0.0, 1.0))) < 1e-13)
    h = 1e-6
    dudx = (float(ht.u(h, 1.0)) - float(ht.u(-h, 1.0))) / (2 * h)
    eq9 = ht.lam * (1.0 + 7.0 / (8.0 * S)
                    - (64.0 * g * g + 117.0 * g - 1026.0) / (1152.0 * S * S))
    check(f"γ={g} R={R}: du/dx(0,1) = Eq.(9)", abs(dudx - eq9) < 1e-8)

# --- 4. Sauer 退化 --------------------------------------------------------------
def max_dM(R):
    ht = HallThroat(R=R)
    st = SauerThroat(R=R)
    xs = np.linspace(-0.3, 0.3, 13)
    rs = np.linspace(0.0, 1.0, 9)
    X, Rr = np.meshgrid(xs, rs)
    return float(np.max(np.abs(ht.mach(X, Rr) - st.mach(X, Rr))))

d50, d500 = max_dM(50.0), max_dM(500.0)
check(f"R=500 で Sauer と一致 |ΔM|={d500:.2e} < 3e-5", d500 < 3e-5)
# Hall−Sauer 差の主項は Sauer 原点シフト (xs_throat ~ R^{-1/2}) 由来の O(S^{-3/2})
# → 1 decade あたり ~31.6 倍で縮む (実測 ~28)
check(f"退化レート d(R=50)/d(R=500) = {d50/d500:.1f} ≈ 31.6 (O(S^-3/2))",
      10.0 < d50 / d500 < 100.0)

# --- 5. C_D --------------------------------------------------------------------
for R, tol in ((2.0, 5e-4), (5.0, 1e-4)):
    ht = HallThroat(R=R)
    cd_n, cd_s = ht.cd_estimate(n=2000), ht.cd_series()
    check(f"R={R}: C_D 数値 {cd_n:.6f} vs 級数 {cd_s:.6f} (|Δ|<{tol})",
          abs(cd_n - cd_s) < tol)

# --- 6. axis_anchor 解析微分 vs 中心差分 -----------------------------------------
ht = HallThroat(R=2.0)
for x in (0.1, 0.26, 0.5):
    M, Mp, Mpp = ht.axis_anchor(x)
    h = 1e-5
    Mn = float(ht.mach(x, 0.0))
    Mpn = (float(ht.mach(x + h, 0.0)) - float(ht.mach(x - h, 0.0))) / (2 * h)
    Mppn = (float(ht.mach(x + h, 0.0)) - 2.0 * Mn + float(ht.mach(x - h, 0.0))) / h ** 2
    check(f"x={x}: anchor M 一致", abs(M - Mn) < 1e-12)
    check(f"x={x}: anchor M' 一致 ({Mp:.6f})", abs(Mp - Mpn) < 1e-6)
    check(f"x={x}: anchor M'' 一致 ({Mpp:.6f})", abs(Mpp - Mppn) < 1e-4)

# --- 7. x_axis_of_mach ----------------------------------------------------------
x0 = ht.x_axis_of_mach(1.05)
check(f"x0(M=1.05) = {x0:.4f}: 自己整合", abs(float(ht.mach(x0, 0.0)) - 1.05) < 1e-12)
check(f"軸ソニック点 {ht.x_sonic_axis:.4f} > 0 (スロート下流)", ht.x_sonic_axis > 0.0)
st100, ht100 = SauerThroat(R=100.0), HallThroat(R=100.0)
check("R=100: 軸ソニック点が Sauer と一致 (1%)",
      abs(ht100.x_sonic_axis - st100.x_sonic_axis) / st100.x_sonic_axis < 0.01)

# --- 8. 壁 θ vs 骨接放物線接線 (診断 — B0/B5 の Sauer ギャップが縮むか) ----------
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    st2 = SauerThroat(R=2.0)
th_wall_geo = np.rad2deg(np.arctan(x0 / 2.0))
th_hall = float(np.rad2deg(ht.theta(x0, 1.0 + x0 * x0 / 4.0)))
th_sauer = float(np.rad2deg(st2.theta(x0, 1.0 + x0 * x0 / 4.0)))
gap_h, gap_s = abs(th_hall - th_wall_geo), abs(th_sauer - th_wall_geo)
print(f"info 壁足 θ [deg]: 幾何 {th_wall_geo:.3f} / Hall {th_hall:.3f} / Sauer {th_sauer:.3f}")
check(f"Hall の壁 θ ギャップ {gap_h:.3f}° < Sauer {gap_s:.3f}°", gap_h < gap_s)

# --- 9. スロート特性線 (A8: 初期値線の CONTUR 化) --------------------------------
for R in (2.0, 5.0):
    htc = HallThroat(R=R)
    xs, rs, Ms, ts = htc.throat_characteristic(n=81)
    check(f"R={R}: 特性線が軸→壁 (r 単調増 0→1)",
          abs(rs[0]) < 1e-12 and abs(rs[-1] - 1.0) < 1e-12
          and bool(np.all(np.diff(rs) > 0)))
    check(f"R={R}: 壁足 = 幾何スロート (x={xs[-1]:.2e}, θ={np.degrees(ts[-1]):.2e}°)",
          abs(xs[-1]) < 1e-12 and abs(ts[-1]) < 1e-12)
    check(f"R={R}: 壁足は超音速 M={Ms[-1]:.4f} > 1", Ms[-1] > 1.0)
    check(f"R={R}: 軸着地 x_A={xs[0]:.4f} > 0 かつ M が Hall 軸値と一致",
          xs[0] > 0.0 and abs(Ms[0] - float(htc.mach(xs[0], 0.0))) < 1e-12)
    check(f"R={R}: 軸着地は縦線 (M=1.05) の軸位置より下流",
          xs[0] > htc.x_axis_of_mach(1.05))
    # C⁻ の定義式 dr/dx = tan(θ-μ) を線上で検証 (中点で評価、端は resample 誤差が乗る)
    xm, rm = 0.5 * (xs[1:] + xs[:-1]), 0.5 * (rs[1:] + rs[:-1])
    slope_fd = np.diff(rs) / np.diff(xs)
    Mm, tm = htc.mach(xm, rm), htc.theta(xm, rm)
    slope_ch = np.tan(tm - np.arcsin(1.0 / Mm))
    err = float(np.max(np.abs(slope_fd[2:-2] / slope_ch[2:-2] - 1.0)))
    check(f"R={R}: 線上で dr/dx = tan(θ-μ) (相対誤差 {err:.2e} < 5e-3)", err < 5e-3)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
