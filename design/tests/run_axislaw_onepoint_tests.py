#!/usr/bin/env python3
"""軸則 D (`OnePointC4AxisLaw`: 端点アンカー + 内部補間点 1 点の C⁴ 区分 5 次) の unit test。

検証対象 (plans/accepted/tooling-nozzle-axislaw-onepoint.md):
  1. A/E 端点 6 条件 (M,M',M''@x_A / M_d,0,0@x_E) を機械精度で満たす (乱数 100 試行)
  2. P を正確に通る (M(x_P)=M_P)
  3. P で 0〜4 階導関数が連続 (有限差分に頼らず、左右の Bernstein 評価で比較)
  4. x>x_E で M=M_d, M'=M''=M'''=0
  5. M'(x)>=0 gate: 単調可能領域の外では monotone_ok=False、内では True
     (事前計算表 L_c=45: ℓ_P=3.74 → M_P∈[2.92,3.11] 単調)
  6. 不可能な ξ_P, η_P (範囲外) の拒否、不正入力 (非有限, L_c<=0, M_d<=M_A)
  7. 無次元座標 ↔ 物理座標の変換 (x_P=x_A+ξ L_c, M_P=M_A+η ΔM, ell_P=ξ L_c)
  8. 4×4 係数系の条件数が有限・拘束残差が機械精度
  9. CPG / semi-perfect のどちらのアンカーでも構成できる (則自体はガス非依存)
  10. 導関数の解析値が有限差分と一致 (Bernstein 微分公式の検算)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.gas import GasCPG, GasSemiPerfect                            # noqa: E402
from forge_design.geometry.axis_law_onepoint import OnePointC4AxisLaw          # noqa: E402
from forge_design.geometry.transonic import HallThroat                        # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


rng = np.random.default_rng(7)
X_A, M_A, MP_A, MPP_A, M_D = 0.4246, 1.1540, 0.5500, 0.0565, 6.0

# --- 1-3, 8. 端点・P・C⁴・条件数 (乱数 100 試行) ------------------------------
err_end = err_P = err_jump = 0.0
cond_max = 0.0
for _ in range(100):
    L_c = rng.uniform(20.0, 80.0)
    xi = rng.uniform(0.03, 0.30)
    eta = rng.uniform(0.1, 0.9)
    law = OnePointC4AxisLaw(X_A, L_c, M_A, MP_A, MPP_A, M_D, xi, eta)
    err_end = max(err_end, abs(float(law(X_A)) - M_A), abs(float(law.deriv(X_A, 1)) - MP_A),
                  abs(float(law.deriv(X_A, 2)) - MPP_A), abs(float(law(law.x_E)) - M_D),
                  abs(float(law.deriv(law.x_E, 1))), abs(float(law.deriv(law.x_E, 2))))
    err_P = max(err_P, abs(float(law(law.x_P)) - law.M_P))
    # 左右の区間多項式を直接評価 (x_P で t1=1, t2=0)
    from forge_design.geometry.axis_law_onepoint import _bernstein_eval
    for k in range(5):
        l = float(_bernstein_eval(law._b, 1.0, k, law.L_1))
        r = float(_bernstein_eval(law._c, 0.0, k, law.L_2))
        err_jump = max(err_jump, abs(l - r) / max(1.0, abs(l), abs(r)))
    cond_max = max(cond_max, law.cond_number)
    err_end = max(err_end, law.constraint_residual)
check(f"端点 6 条件 + 拘束残差 最大誤差 {err_end:.2e}", err_end < 1e-9)
check(f"P を通る: max|M(x_P)-M_P| = {err_P:.2e}", err_P < 1e-9)
check(f"P で 0〜4 階連続 (相対 jump 最大 {err_jump:.2e})", err_jump < 1e-8)
check(f"4×4 系の条件数 最大 {cond_max:.3g} (有限)", np.isfinite(cond_max) and cond_max < 1e9)

# --- 4. x > x_E ---
law = OnePointC4AxisLaw(X_A, 45.0, M_A, MP_A, MPP_A, M_D, 3.74 / 45.0, (3.0 - M_A) / (M_D - M_A))
check("x>x_E で M=M_d, M'=M''=M'''=0",
      abs(float(law(law.x_E + 5.0)) - M_D) < 1e-12
      and float(law.deriv(law.x_E + 5.0, 1)) == 0.0
      and float(law.deriv(law.x_E + 5.0, 2)) == 0.0
      and float(law.deriv(law.x_E + 5.0, 3)) == 0.0)

# --- 5. 単調 gate (事前計算表と整合) ---
xi = 3.74 / 45.0
for MP, expect in ((2.80, False), (3.00, True), (3.20, False)):
    eta = (MP - M_A) / (M_D - M_A)
    g = OnePointC4AxisLaw(X_A, 45.0, M_A, MP_A, MPP_A, M_D, xi, eta).gates()
    check(f"単調 gate: ℓ_P=3.74, M_P={MP} → monotone_ok={g['monotone_ok']} (期待 {expect})",
          g["monotone_ok"] == expect)
check("monotone_feasible ヘルパが gates と一致",
      OnePointC4AxisLaw.monotone_feasible(X_A, 45.0, M_A, MP_A, MPP_A, M_D, xi, (3.0 - M_A) / (M_D - M_A))
      and not OnePointC4AxisLaw.monotone_feasible(X_A, 45.0, M_A, MP_A, MPP_A, M_D, xi, (2.5 - M_A) / (M_D - M_A)))
g = law.gates()
check("gates() が要求キーを持つ (端点残差/P jump/min M'/M'' 反転/max|M'''|/J_axis)",
      all(k in g for k in ("endpoint_residual", "P_jumps", "P_jump_max", "min_dMdx",
                           "mpp_sign_changes", "max_abs_Mppp", "J_axis")))

# --- 6. 不正入力 ---
for kw in (dict(xi_P=0.0), dict(xi_P=1.0), dict(eta_P=-0.1), dict(eta_P=1.2),
          dict(L_c=-1.0), dict(M_d=M_A - 0.1), dict(x_A=float("nan"))):
    args = dict(x_A=X_A, L_c=45.0, M_A=M_A, Mp_A=MP_A, Mpp_A=MPP_A, M_d=M_D, xi_P=0.1, eta_P=0.4)
    args.update(kw)
    try:
        OnePointC4AxisLaw(**args)
        check(f"不正入力 {kw} 拒否", False)
    except ValueError:
        check(f"不正入力 {kw} 拒否", True)

# --- 7. 無次元 ↔ 物理 ---
law = OnePointC4AxisLaw(X_A, 50.0, M_A, MP_A, MPP_A, M_D, 0.12, 0.45)
check("x_P = x_A + ξ L_c", abs(law.x_P - (X_A + 0.12 * 50.0)) < 1e-12)
check("M_P = M_A + η (M_d − M_A)", abs(law.M_P - (M_A + 0.45 * (M_D - M_A))) < 1e-12)
check("ell_P = ξ L_c、L_1 + L_2 = L_c",
      abs(law.ell_P - 0.12 * 50.0) < 1e-12 and abs(law.L_1 + law.L_2 - 50.0) < 1e-12)

# --- 9. CPG / semi-perfect アンカー ---
Y = {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.022, "N2": 0.7245}
for tag, gam in (("semiperfect", GasSemiPerfect(Y, Tt=1550.0).gamma_star), ("CPG", GasCPG(1.3, 1150.0).gamma_ref)):
    ht = HallThroat(R=3.0, gamma=gam)
    xs_c, *_ = ht.throat_characteristic(n=41)
    xA = float(xs_c[0]); mA, mpA, mppA = ht.axis_anchor(xA)
    lw = OnePointC4AxisLaw(xA, 45.0, mA, mpA, mppA, 6.0, 0.083, 0.38)
    check(f"{tag} アンカーで構成可・端点残差 {lw.gates()['endpoint_residual']:.1e}",
          lw.gates()["endpoint_residual"] < 1e-9)

# --- 10. 有限差分検算 ---
law = OnePointC4AxisLaw(X_A, 45.0, M_A, MP_A, MPP_A, M_D, 0.083, 0.38)
worst = 0.0
for x in (X_A + 0.3, law.x_P - 0.2, law.x_P + 0.2, law.x_E - 2.0):
    h = 1e-4
    d1 = (float(law(x + h)) - float(law(x - h))) / (2 * h)
    d2 = (float(law(x + h)) - 2 * float(law(x)) + float(law(x - h))) / h ** 2
    worst = max(worst, abs(d1 - float(law.deriv(x, 1))), abs(d2 - float(law.deriv(x, 2))) / 10)
check(f"解析微分 vs 有限差分 (最大差 {worst:.2e})", worst < 1e-5)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
