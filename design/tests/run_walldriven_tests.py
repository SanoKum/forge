#!/usr/bin/env python3
"""wall-driven チェーン Phase W1 の解析 unit test (pytest 非依存・素の assert)。

検証対象 (plans/active/tooling-nozzle-walldriven-chain.md §7):
  1. U→T / T→D の端点条件 (1e-12 級)
  2. 原案の閉形式 (U→T 勾配, T→D r'', r_D) と実装の一致
  3. 単調条件 μ ≤ 20 / 変曲条件 L_D ≤ 3 R_t tanθ_D の境界 ±0.1% 判別
  4. 条件違反時に実際に非単調半径・余分な変曲点が発生し validate() が検出すること
  5. 解析 κ・dκ/ds と数値微分の一致 (実装式の独立検証)
  6. geometry option (§10): λ 窓 [5/3, 5/2]・λ=2 の 4 次退化
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.geometry.wall_walldriven import (  # noqa: E402
    ContractionToConeQuintic, ThroatExpansionPoly, UpstreamThroatPoly,
    WallDrivenThroatRegion)

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


# --- 1. 端点条件 ---------------------------------------------------------------
r_U, R_t, L_U = 2.5, 2.0, 3.5
up = UpstreamThroatPoly(r_U=r_U, R_t=R_t, L_U=L_U)
check("U→T r(-L_U)=r_U", abs(up.r(-L_U) - r_U) < 1e-12)
check("U→T r'(-L_U)=0", abs(up.r(-L_U, 1)) < 1e-12)
check("U→T r''(-L_U)=0", abs(up.r(-L_U, 2)) < 1e-12)
check("U→T r(0)=r_t", abs(up.r(0.0) - 1.0) < 1e-12)
check("U→T r'(0)=0", abs(up.r(0.0, 1)) < 1e-12)
check("U→T r''(0)=1/R_t", abs(up.r(0.0, 2) - 1.0 / R_t) < 1e-12)

theta_D = np.deg2rad(12.0)
L_D = 0.9 * 3.0 * R_t * np.tan(theta_D)
dn = ThroatExpansionPoly(R_t=R_t, L_D=L_D, theta_D=theta_D)
check("T→D r(0)=r_t", abs(dn.r(0.0) - 1.0) < 1e-12)
check("T→D r'(0)=0", abs(dn.r(0.0, 1)) < 1e-12)
check("T→D r''(0)=1/R_t", abs(dn.r(0.0, 2) - 1.0 / R_t) < 1e-12)
check("T→D r'(L_D)=tanθ_D", abs(dn.r(L_D, 1) - np.tan(theta_D)) < 1e-12)
check("T→D r''(L_D)=0", abs(dn.r(L_D, 2)) < 1e-12)
check("T→D r(L_D)=r_D 閉形式", abs(dn.r(L_D) - dn.r_D) < 1e-12)

# --- 2. 閉形式と実装の一致 (乱数 50 試行) ---------------------------------------
rng = np.random.default_rng(41)
ok_grad, ok_rpp, ok_rD = True, True, True
for _ in range(50):
    ru = rng.uniform(1.2, 8.0)
    Rt = rng.uniform(0.5, 10.0)
    dr_ = ru - 1.0
    Lu = rng.uniform(0.1, 1.0) * np.sqrt(20 * Rt * dr_)
    u = UpstreamThroatPoly(r_U=ru, R_t=Rt, L_U=Lu)
    xi = np.linspace(0, 1, 501)
    x = -Lu + xi * Lu
    mu = u.mu
    cf = dr_ / (2 * Lu) * xi**2 * (xi - 1) * (5 * mu * xi - 3 * mu - 60 * xi + 60)
    ok_grad &= bool(np.max(np.abs(u.r(x, 1) - cf))
                    <= 1e-9 * max(1.0, np.max(np.abs(cf))))
    th = rng.uniform(2, 40) * np.pi / 180
    Ld = rng.uniform(0.05, 1.0) * 3 * Rt * np.tan(th)
    d = ThroatExpansionPoly(R_t=Rt, L_D=Ld, theta_D=th)
    xi2 = np.linspace(0, 1, 501)
    rpp_cf = (1 - xi2) / Rt * (1 + xi2 * (6 * Rt * np.tan(th) / Ld - 3))
    ok_rpp &= bool(np.max(np.abs(d.r(xi2 * Ld, 2) - rpp_cf))
                   <= 1e-9 * max(1.0, np.max(np.abs(rpp_cf))))
    # r_D = r_t + ∫q dx (高解像 Simpson)
    from scipy.integrate import simpson
    xs = np.linspace(0, Ld, 4001)
    ok_rD &= bool(abs(d.r_D - (1.0 + simpson(d.r(xs, 1), x=xs))) < 1e-9 * d.r_D)
check("U→T 勾配閉形式 (50 乱数試行)", ok_grad)
check("T→D r'' 閉形式 (50 乱数試行)", ok_rpp)
check("T→D r_D = ∫q + r_t (50 乱数試行)", ok_rD)

# --- 3./4. 境界判別と違反検出 ---------------------------------------------------
dr_ = 1.5
for fac, expect_ok in ((0.999, True), (1.001, False)):
    Lu = np.sqrt(20.0 * fac * R_t * dr_)
    u = UpstreamThroatPoly(r_U=1.0 + dr_, R_t=R_t, L_U=Lu)
    v = u.validate(n=200001)
    check(f"U→T μ=20×{fac}: validate {'空' if expect_ok else '違反検出'}",
          (len(v) == 0) == expect_ok)
    if not expect_ok:
        x = np.linspace(-Lu, 0, 200001)
        check("U→T μ>20 で実際に dr/dx>0 (非単調半径) が発生",
              bool(np.any(u.r(x, 1) > 1e-14)))

for fac, expect_ok in ((0.999, True), (1.05, False)):
    Ld = 3.0 * R_t * np.tan(theta_D) * fac
    d = ThroatExpansionPoly(R_t=R_t, L_D=Ld, theta_D=theta_D)
    v = d.validate(n=200001)
    check(f"T→D L_D=上限×{fac}: validate {'空' if expect_ok else '違反検出'}",
          (len(v) == 0) == expect_ok)
    if not expect_ok:
        x = np.linspace(0, Ld, 200001)
        rpp = d.r(x, 2)
        # 内部で符号が負→正へ戻る = D 以前に余分な変曲点
        check("T→D L_D>上限で実際に r''<0 の内部区間 (余分な変曲点) が発生",
              bool(np.any(rpp < -1e-14)) and rpp[0] > 0)

# --- 5. 合成壁: κ・dκ/ds の解析式 vs 数値微分 ----------------------------------
w = WallDrivenThroatRegion(r_U=r_U, R_t=R_t, L_U=L_U, L_D=L_D, theta_D=theta_D)
check("合成壁 validate 空", w.validate() == [])
check("合成壁 κ(T)=1/R_t (両側連続)",
      abs(w.kappa(-1e-12) - 1.0 / R_t) < 1e-9
      and abs(w.kappa(+1e-12) - 1.0 / R_t) < 1e-9)
# κ の数値微分 (弧長) と解析 dκ/ds — 区分内部のみ (T の κ' 跳びを跨がない)
for lo, hi, tag in ((-L_U * 0.95, -1e-3, "U→T"), (1e-3, L_D * 0.95, "T→D")):
    x = np.linspace(lo, hi, 20001)
    kap = w.kappa(x)
    s = np.concatenate([[0.0], np.cumsum(
        np.sqrt(1.0 + w.r(0.5 * (x[1:] + x[:-1]), 1) ** 2) * np.diff(x))])
    dk_num = np.gradient(kap, s)
    dk_ana = w.dkappa_ds(x)
    err = np.max(np.abs(dk_num[5:-5] - dk_ana[5:-5])) / max(
        1e-30, np.max(np.abs(dk_ana)))
    check(f"dκ/ds 解析式 vs 数値微分 ({tag})", err < 1e-4)
check("κ' 跳び診断が有限値", np.isfinite(w.kappa_prime_jump_at_throat))
diag = w.diagnostics()
check("diagnostics キー揃い",
      set(diag) >= {"mu", "L_D_over_max", "r_D", "kappa_prime_jump_at_throat",
                    "max_abs_dkappa_ds", "violations"})

# --- 6. geometry option §10 -----------------------------------------------------
th_a = np.deg2rad(15.0)
dr2 = 1.2
for lam, expect_ok in ((1.6, False), (5.0 / 3.0, True), (2.0, True),
                       (2.5, True), (2.6, False)):
    L = lam * dr2 / np.tan(th_a)
    c = ContractionToConeQuintic(r_in=2.2, dr=dr2, theta_a=th_a, L=L)
    check(f"§10 λ={lam:.3g}: validate {'空' if expect_ok else '違反検出'}",
          (len(c.validate(n=100001)) == 0) == expect_ok)
c2 = ContractionToConeQuintic(r_in=2.2, dr=dr2, theta_a=th_a,
                              L=2.0 * dr2 / np.tan(th_a))
xi = np.linspace(0, 1, 101)
quart = 2.2 - 2.0 * dr2 * xi**3 + dr2 * xi**4
check("§10 λ=2 で 4 次式 r_in - 2Δr s³ + Δr s⁴ に退化",
      np.max(np.abs(c2.r(xi * c2.L) - quart)) < 1e-10)
check("§10 端点 (値・勾配・曲率)",
      abs(c2.r(0.0) - 2.2) < 1e-12 and abs(c2.r(0.0, 1)) < 1e-12
      and abs(c2.r(0.0, 2)) < 1e-12
      and abs(c2.r(c2.L, 1) + np.tan(th_a)) < 1e-12
      and abs(c2.r(c2.L, 2)) < 1e-12)

print("-" * 60)
print(f"{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
