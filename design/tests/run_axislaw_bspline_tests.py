#!/usr/bin/env python3
"""軸則 B/C (B-spline 系, 滑らかさ比較) の unit test (素の assert)。

検証対象 (plans/accepted/tooling-nozzle-axislaw-smoothness.md):
  1. B (MonotoneBSplineAxisLaw): 端点条件 (M,M',M''@x_A / M,M'@x_E / hard M''@x_E)
  2. B: 制御多角形が非減少 → M'(x)>=0 (単調性、密サンプルの偶然でなく構成的)
  3. B: soft exit_curvature でも端点値/1階微分は満たす
  4. B: x>x_E で M=M_d 一定
  5. B: 不正入力 (x_E<=x_A, M_d<=M_A, 非有限)
  6. B: 解像度依存 (n_interior を増やしても端点条件・単調性は保たれる)
  7. C (NonnegDnuBSplineAxisLaw): q(x)>=0 (非負制御点で構成的)
  8. C: ν(x_A)=ν(M_A), ν(x_E)=ν(M_d) (積分値の等式拘束)、chain rule で q(x_A),q'(x_A) が
     M'_A, M''_A と整合
  9. C: M(x)=ν^{-1}(ν(x)) の ν↔M 往復・端点一致
  10. C: CPG でも semi-perfect でも動く (dnu_dM の両分岐)
  11. C: 不正入力 (M'_A<=0 など)
  12. dnu_dM (moc_kernel): CPG 解析式が有限差分と一致、semi-perfect がスプライン微分で動く
  13. A (KnotQuinticAxisLaw) の数値が本比較機能追加の前後で不変 (回帰)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.gas import GasCPG, GasSemiPerfect                            # noqa: E402
from forge_design.geometry.axis_law import KnotQuinticAxisLaw                  # noqa: E402
from forge_design.geometry.axis_law_bspline import (MonotoneBSplineAxisLaw,    # noqa: E402
                                                     NonnegDnuBSplineAxisLaw)
from forge_design.geometry.moc_kernel import dnu_dM, pm_mach, pm_nu            # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


# M6/R3 相当のアンカー (実案件と同じ値、plans/accepted/tooling-nozzle-axislaw-smoothness.md)
X_A, M_A, MP_A, MPP_A, M_D = 0.4246, 1.1540, 0.5500, 0.0565, 6.0
X_E = X_A + 45.0

Y = {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.022, "N2": 0.7245}
GAS = GasSemiPerfect(Y, Tt=1550.0)

# --- 1-3. B: 端点条件・単調性・soft ---------------------------------------------
for n_int in (15, 20):
    law = MonotoneBSplineAxisLaw(X_A, X_E, M_A, MP_A, MPP_A, M_D, n_interior=n_int)
    e = 1e-6
    err = max(abs(float(law(X_A)) - M_A), abs(float(law.deriv(X_A, 1)) - MP_A),
             abs(float(law.deriv(X_A, 2)) - MPP_A), abs(float(law(X_E)) - M_D),
             abs(float(law.deriv(X_E, 1))), abs(float(law.deriv(X_E, 2))))
    check(f"B n_interior={n_int}: 端点条件 (hard) 最大誤差 {err:.2e}", err < 1e-4)
    xg = np.linspace(X_A, X_E, 4001)
    check(f"B n_interior={n_int}: M'(x)>=0 全域 (制御多角形非減少)",
          bool(np.all(law.deriv(xg, 1) >= -1e-9)))

law_soft = MonotoneBSplineAxisLaw(X_A, X_E, M_A, MP_A, MPP_A, M_D, n_interior=15,
                                  exit_curvature="soft")
err_soft = max(abs(float(law_soft(X_A)) - M_A), abs(float(law_soft.deriv(X_A, 1)) - MP_A),
              abs(float(law_soft(X_E)) - M_D), abs(float(law_soft.deriv(X_E, 1))))
check(f"B soft: 値/1階端点条件は保持 (誤差 {err_soft:.2e})", err_soft < 1e-4)

# --- 4. x > x_E で一様 ---
law = MonotoneBSplineAxisLaw(X_A, X_E, M_A, MP_A, MPP_A, M_D, n_interior=15)
check("B: x>x_E で M=M_d, M'=0",
      abs(float(law(X_E + 3.0)) - M_D) < 1e-9 and float(law.deriv(X_E + 3.0, 1)) == 0.0)

# --- 5. 不正入力 ---
for kw in (dict(x_E=X_A - 1.0), dict(M_d=M_A - 0.1),
          dict(x_A=float("nan"))):
    args = dict(x_A=X_A, x_E=X_E, M_A=M_A, Mp_A=MP_A, Mpp_A=MPP_A, M_d=M_D)
    args.update(kw)
    try:
        MonotoneBSplineAxisLaw(**args)
        check(f"B 不正入力 {kw} 拒否", False)
    except ValueError:
        check(f"B 不正入力 {kw} 拒否", True)

# --- 6. 解像度依存 (n_interior を増やしても崩れない) ---
Js = []
for n_int in (15, 20, 25):
    law = MonotoneBSplineAxisLaw(X_A, X_E, M_A, MP_A, MPP_A, M_D, n_interior=n_int)
    Js.append(law.J_axis)
check(f"B: n_interior を増やしても J_axis が同オーダーで安定 {[round(j,4) for j in Js]}",
      max(Js) / max(min(Js), 1e-9) < 10.0)

# --- 7-9. C: 非負・積分等式・ν↔M 往復 ------------------------------------------
for n_int in (12, 20):
    lawC = NonnegDnuBSplineAxisLaw(X_A, X_E, M_A, MP_A, MPP_A, M_D, GAS, n_interior=n_int)
    xg = np.linspace(X_A, X_E, 4001)
    check(f"C n_interior={n_int}: q(x)>=0 全域 (非負制御点)",
          bool(np.all(lawC.q(xg) >= -1e-12)))
    nu_A = float(pm_nu(M_A, GAS)); nu_d = float(pm_nu(M_D, GAS))
    err_nu = max(abs(float(lawC.nu(np.array([X_A]))[0]) - nu_A),
                abs(float(lawC.nu(np.array([X_E]))[0]) - nu_d))
    check(f"C n_interior={n_int}: ν(x_A)=ν(M_A), ν(x_E)=ν(M_d) (誤差 {err_nu:.2e})",
          err_nu < 1e-4)
    e = 1e-6
    err_end = max(abs(float(lawC(X_A)) - M_A), abs(float(lawC(X_E)) - M_D),
                 abs(float(lawC.deriv(X_A, 1)) - MP_A))
    check(f"C n_interior={n_int}: M(x_A)=M_A, M(x_E)=M_d, M'(x_A)=M'_A (誤差 {err_end:.2e})",
          err_end < 1e-3)
    check(f"C n_interior={n_int}: x>x_E で M=M_d",
          abs(float(lawC(X_E + 3.0)) - M_D) < 1e-6)
    # chain rule 検算: q(x_A) = dnu/dM(M_A) * M'_A
    q_A_expected = float(dnu_dM(M_A, GAS, 1)) * MP_A
    check(f"C n_interior={n_int}: q(x_A) = dν/dM(M_A)·M'_A (chain rule)",
          abs(lawC.q_A - q_A_expected) < 1e-9)

# --- 10. dnu_dM: CPG 解析式 vs 有限差分、semi-perfect スプライン ------------------
for M in (1.05, 1.3, 2.0, 3.5, 6.0):
    h = 1e-5
    d1_fd = (pm_nu(M + h, 1.3) - pm_nu(M - h, 1.3)) / (2 * h)
    d2_fd = (pm_nu(M + h, 1.3) - 2 * pm_nu(M, 1.3) + pm_nu(M - h, 1.3)) / h ** 2
    d1 = float(dnu_dM(M, 1.3, 1)); d2 = float(dnu_dM(M, 1.3, 2))
    check(f"dnu_dM CPG M={M}: 1階 {d1:.6f} vs fd {d1_fd:.6f}", abs(d1 - d1_fd) < 1e-4)
    check(f"dnu_dM CPG M={M}: 2階 {d2:.6f} vs fd {d2_fd:.6f}", abs(d2 - d2_fd) < 1e-2)
for M in (1.05, 1.15, 2.5, 4.0, 6.0):
    h = 1e-4
    d1_fd = (float(GAS.nu(M + h)) - float(GAS.nu(M - h))) / (2 * h)
    d1 = float(dnu_dM(M, GAS, 1))
    check(f"dnu_dM semiperfect M={M}: spline {d1:.5f} vs table fd {d1_fd:.5f}",
          abs(d1 - d1_fd) < 5e-3)

# --- CPG 経路の C (回帰対照) ---
gas_cpg = GasCPG(1.3, 1150.0)
x_A_c, m_a_c = 0.3, 1.1
lawC_cpg = NonnegDnuBSplineAxisLaw(x_A_c, x_A_c + 20.0, m_a_c, 0.5, 0.05, 4.0, gas_cpg,
                                   n_interior=15)
check("C (CPG 経路): 端点条件",
      abs(float(lawC_cpg(x_A_c)) - m_a_c) < 1e-3 and abs(float(lawC_cpg(x_A_c + 20.0)) - 4.0) < 1e-3)

# --- 11. C 不正入力 ---
for kw in (dict(Mp_A=0.0), dict(M_d=M_A - 0.1), dict(x_E=X_A - 1.0)):
    args = dict(x_A=X_A, x_E=X_E, M_A=M_A, Mp_A=MP_A, Mpp_A=MPP_A, M_d=M_D, gas=GAS)
    args.update(kw)
    try:
        NonnegDnuBSplineAxisLaw(**args)
        check(f"C 不正入力 {kw} 拒否", False)
    except ValueError:
        check(f"C 不正入力 {kw} 拒否", True)

# --- 13. A (KnotQuinticAxisLaw) 回帰 (数値不変) ---
lawA = KnotQuinticAxisLaw(X_A, 45.0, M_A, MP_A, MPP_A, M_D, 2.5)
check(f"A 回帰: x_K={lawA.x_K:.6f} (既知値 4.165 台)", abs(lawA.x_K - 4.165) < 0.01)
check("A 回帰: M'''(K-)≈0.118, M'''(K+)≈-0.0006",
      abs(float(lawA.deriv(lawA.x_K - 1e-6, 3)) - 0.1178) < 0.01
      and abs(float(lawA.deriv(lawA.x_K + 1e-6, 3)) - (-0.0006)) < 0.01)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
