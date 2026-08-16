#!/usr/bin/env python3
"""axis-Mach チェーン A3: Hall + Hermite + 逆 MOC + 壁組立の統合 test。

検証対象 (plans/active/tooling-nozzle-axismach-chain.md §5.3, §7):
  1. Hall アンカー → 許容窓内 L_c → 逆 MOC が壁を返す (E→F 閉包込み)
  2. 壁 QA: 単調・出口壁角 |θ|<0.2°・r_F が 1D 理論 ±3%・壁 M ≈ M_d
  3. 質量流量診断 (mdot_exit/mdot_start ≈ 1)
  4. AxisMachCFDWall: 接合 C1/C2 (直管/U→T/放物線/設計壁)・単調・最小半径 1
  5. x_F − x_E と terminal Mach line 近似の比較 (情報値: 同オーダー)

実行時間 ~25 s (逆 MOC 三角充填が支配)。
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.geometry.axis_law import QuinticHermiteAxisLaw  # noqa: E402
from forge_design.geometry.moc_inverse import inverse_design  # noqa: E402
from forge_design.geometry.transonic import HallThroat  # noqa: E402
from forge_design.geometry.wall_axismach import (  # noqa: E402
    AxisMachCFDWall, area_ratio_isentropic, wall_qa)

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


g, Md, R = 1.4, 4.0, 2.0
ht = HallThroat(R=R, gamma=g)
x0 = ht.x_axis_of_mach(1.05)
M0, Mp0, Mpp0 = ht.axis_anchor(x0)
lo, hi = QuinticHermiteAxisLaw.admissible_Lc_range(x0, M0, Mp0, Mpp0, Md)
check(f"Hall アンカー窓 ({lo:.3f}, {hi:.3f}) 非退化", hi > 2.0)
Lc = hi * 0.98
law = QuinticHermiteAxisLaw(x0, Lc, M0, Mp0, Mpp0, Md)
check("ゲート合格", not law.gates()["violations"])

rF_pred = float(np.sqrt(area_ratio_isentropic(Md, g)))
x_end = law.x_E + 2.3 * rF_pred * np.sqrt(Md * Md - 1.0)
res = inverse_design(ht, lambda x: float(law(x)), x_axis_end=float(x_end),
                     n_axis=500, n_start=41, gamma=g,
                     th_wall0=float(np.arctan(x0 / R)), M_start=1.05)
wall = res["wall"]
check(f"壁テーブル {len(wall)} 点 (>500)", len(wall) > 500)

qa = wall_qa(wall, Md, law.x_E, g)
print("info QA:", {k: (round(v, 4) if isinstance(v, float) else v)
                   for k, v in qa.items()})
check("壁 QA violations 空", not qa["violations"])
check(f"出口壁角 {qa['theta_exit_deg']:.3f}° < 0.2°", abs(qa["theta_exit_deg"]) < 0.2)
check(f"r_F 誤差 {qa['r_F_err_rel']*100:.2f}% < 3%", abs(qa["r_F_err_rel"]) < 0.03)
check(f"壁出口 M = {qa['M_wall_exit']:.4f} (±0.02)", abs(qa["M_wall_exit"] - Md) < 0.02)
ratio = res["mdot_exit"] / res["mdot_start"]
check(f"質量流量診断 {ratio:.4f} (±2%)", abs(ratio - 1.0) < 0.02)
dEF = qa["xF_minus_xE"] / qa["xF_minus_xE_pred"]
check(f"x_F−x_E / terminal-line 予測 = {dEF:.2f} (同オーダー 0.5–1.5)",
      0.5 < dEF < 1.5)

# --- AxisMachCFDWall 組立 -------------------------------------------------------
w = AxisMachCFDWall(wall[:, :2], R=R, r_U=2.5, L_U=3.5, L_pipe=0.5)
msgs = w.validate()
check(f"AxisMachCFDWall.validate: {msgs or 'クリーン'}", not msgs)
xs = np.linspace(w.x_in, w.x_e, 2000)
rv = w.r(xs)
check("全域 r > 0・有限", bool(np.all(np.isfinite(rv)) and np.all(rv > 0.0)))
check(f"x_in={w.x_in:.2f} < 0 < x0={w.x0:.3f} < x_e={w.x_e:.2f}",
      w.x_in < 0.0 < w.x0 < w.x_e)
# スロートで最小・r=1
i_min = int(np.argmin(rv))
check(f"最小半径 {float(rv[i_min]):.5f} ≈ 1 at x={float(xs[i_min]):.3f} ≈ 0",
      abs(float(rv[i_min]) - 1.0) < 5e-3 and abs(float(xs[i_min])) < 0.2)



# --- スロート整合ゲート + spline リンギング検査 (2026-08-16) ----------------------
from forge_design.geometry.axis_law import QuinticHermiteAxisLaw
from forge_design.geometry.moc_inverse import inverse_design
from forge_design.geometry.transonic import HallThroat
from forge_design.geometry.wall_axismach import throat_curvature_fit

_ht = HallThroat(R=2.0, gamma=1.4)
_xA = float(_ht.throat_characteristic(n=41)[0][0])
_MA, _Mp, _Mpp = _ht.axis_anchor(_xA)


def _design(L_c):
    law = QuinticHermiteAxisLaw(_xA, L_c, _MA, _Mp, _Mpp, 4.0)
    x_end = law.x_E + 2.3 * 3.274 * np.sqrt(15.0)
    return inverse_design(_ht, lambda x: float(law(max(x, _xA))), x_axis_end=x_end,
                          n_axis=500, n_start=41, gamma=1.4, M_start=1.05,
                          exit_mode="characteristic", x_E=law.x_E, M_d=4.0,
                          start_line="throat_char", wall_mode="cplus"), law


_res_ok, _law_ok = _design(9.715)
_qa_ok = wall_qa(_res_ok["wall"], 4.0, _law_ok.x_E, 1.4, R=2.0)
check(f"曲率整合: 生産 L_c で κ₀R = {_qa_ok['kappa0_R_ratio']:.3f} ∈ [0.85,1.15]",
      abs(_qa_ok["kappa0_R_ratio"] - 1.0) <= 0.15 and not _qa_ok["violations"])
_w_ok = AxisMachCFDWall(_res_ok["wall"], R=2.0)
check("リンギング: 生産 L_c の壁 spline はテーブル θ と 0.2° 以内",
      not [m for m in _w_ok.validate() if "リンギング" in m])

_res_ng, _law_ng = _design(4.0)
_qa_ng = wall_qa(_res_ng["wall"], 4.0, _law_ng.x_E, 1.4, R=2.0)
check(f"曲率整合: L_c=4 (設計不整合) を REJECT (κ₀R = {_qa_ng['kappa0_R_ratio']:.2f})",
      any("スロート曲率不整合" in v for v in _qa_ng["violations"]))
_w_ng = AxisMachCFDWall(_res_ng["wall"], R=2.0)
check("リンギング: L_c=4 の壁 spline 振動を検出",
      any("リンギング" in m for m in _w_ng.validate()))
_k0, _resid, _n = throat_curvature_fit(_res_ok["wall"])
check(f"throat_curvature_fit: κ₀={_k0:.3f} が Hall の 1/R に近い (fit 残差 {_resid:.1e})",
      0.4 < _k0 < 0.55 and _resid < 5e-4)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
