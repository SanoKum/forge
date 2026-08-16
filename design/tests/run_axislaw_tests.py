#!/usr/bin/env python3
"""axis-Mach チェーン A1: QuinticHermiteAxisLaw の unit test (素の assert)。

検証対象 (plans/active/tooling-nozzle-axismach-chain.md §5.1, §7):
  1. 両端 6 条件 (M, M', M'' × 両端) の機械精度充足 (乱数 200 試行)
  2. x > x_E で M = M_d・M' = M'' = 0 (一様領域への C² 接続)
  3. 設計ゲート: 許容窓内で単調合格、極端に短い L_c で M'<0 検出
  4. admissible_Lc_range: 窓の内側で合格・両外側で不合格 (上限の存在も確認 —
     M''_A>0 では a2 ∝ L_c² により長すぎる L_c も単調性を破る)
  5. B10 実測アンカー (1.946, 0.290, 0.190) の窓が B8 相当長 12.3 を含まないこと
     (実装時に発見した設計制約のリグレッション)
  6. M'' 符号反転回数: 素直な条件で ≤ 1
  7. 不正入力 (L_c ≤ 0, M_d ≤ M_A, 非有限) の ValueError
  8. KnotQuinticAxisLaw (A6): 始端 3 条件・knot での C² 連続・終端 M'=M''=0・
     区分 2 の構成的単調/凹・M'' 反転 ≤ 1・L_c 窓が単一 quintic の上限を超えて開くこと
     (M6, R=3 相当アンカーで quintic 上限 17.9 に対し knot は L_c=60 でも合格)・
     _axis_grid (スロート側等比細分) の総長/初項/単調
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.geometry.axis_law import QuinticHermiteAxisLaw, KnotQuinticAxisLaw  # noqa: E402
from forge_design.geometry.moc_inverse import _axis_grid  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


rng = np.random.default_rng(41)

# --- 1. 端点条件 (乱数 200 試行) ---------------------------------------------
err = 0.0
for _ in range(200):
    x_A = rng.uniform(-1.0, 2.0)
    L_c = rng.uniform(2.0, 30.0)
    M_A = rng.uniform(1.02, 2.0)
    Mp_A = rng.uniform(0.05, 0.8)
    Mpp_A = rng.uniform(-0.3, 0.5)
    M_d = M_A + rng.uniform(0.5, 4.0)
    law = QuinticHermiteAxisLaw(x_A, L_c, M_A, Mp_A, Mpp_A, M_d)
    x_E = x_A + L_c
    err = max(err,
              abs(float(law(x_A)) - M_A),
              abs(float(law.deriv(x_A, 1)) - Mp_A),
              abs(float(law.deriv(x_A, 2)) - Mpp_A),
              abs(float(law(x_E)) - M_d),
              abs(float(law.deriv(x_E, 1))),
              abs(float(law.deriv(x_E, 2))))
check(f"端点 6 条件 機械精度 (max err = {err:.2e})", err < 1e-9)

# --- 2. x > x_E の一様領域 ------------------------------------------------------
law = QuinticHermiteAxisLaw(0.3, 12.0, 1.05, 0.3, 0.15, 4.0)
xs = np.linspace(law.x_E + 1e-9, law.x_E + 20.0, 50)
check("x>x_E: M=M_d", float(np.max(np.abs(law(xs) - 4.0))) < 1e-12)
check("x>x_E: M'=0", float(np.max(np.abs(law.deriv(xs, 1)))) < 1e-12)
check("x>x_E: M''=0", float(np.max(np.abs(law.deriv(xs, 2)))) < 1e-12)

# --- 3-4. 許容窓 ----------------------------------------------------------------
anch = (0.3, 1.05, 0.3, 0.15, 4.0)     # (x_A, M_A, M'_A, M''_A, M_d)
lo, hi = QuinticHermiteAxisLaw.admissible_Lc_range(*anch)
check(f"許容窓 ({lo:.4f}, {hi:.4f}) が非退化", hi > lo > 0.0)
mid = float(np.sqrt(max(lo, 0.5) * hi))
check("窓内 (幾何平均) で合格",
      QuinticHermiteAxisLaw(anch[0], mid, *anch[1:]).gates()["monotone_ok"])
check("窓上限の直上 3% で不合格 (M''_A>0 は長い L_c が単調性を破る)",
      hi < 199.0
      and not QuinticHermiteAxisLaw(anch[0], hi * 1.03, *anch[1:]).gates()["monotone_ok"])
# 短い側: smoothstep 退化で単調は保たれるが max dM/dx が発散的に増える
g_s1 = QuinticHermiteAxisLaw(0.3, 0.5, 1.05, 0.3, 0.15, 4.0).gates()
g_s2 = QuinticHermiteAxisLaw(0.3, 5.0, 1.05, 0.3, 0.15, 4.0).gates()
check("短い L_c: 単調のまま (smoothstep 退化)", g_s1["monotone_ok"])
check(f"短い L_c: max dM/dx 増大 ({g_s1['max_dMdx']:.2f} > {g_s2['max_dMdx']:.2f})",
      g_s1["max_dMdx"] > 5.0 * g_s2["max_dMdx"])

# --- 5. B10 実測アンカーのリグレッション -----------------------------------------
lo10, hi10 = QuinticHermiteAxisLaw.admissible_Lc_range(1.7042, 1.94576, 0.28960,
                                                       0.19000, 4.0)
check(f"B10 アンカー窓 = ({lo10:.3f}, {hi10:.3f})", hi10 > lo10 > 0.0)
check("B10 アンカー: B8 相当長 L_c=12.3 は窓外 (単一 quintic の設計制約)",
      not (lo10 <= 12.3 <= hi10))

# --- 6. M'' 符号反転 ------------------------------------------------------------
g = QuinticHermiteAxisLaw(0.3, mid, 1.05, 0.3, 0.15, 4.0).gates()
check(f"窓内条件で M'' 反転 {g['mpp_sign_changes']} ≤ 1", g["mpp_quality_ok"])

# --- 7. 不正入力 ---------------------------------------------------------------
for bad in ((0.3, -1.0, 1.05, 0.3, 0.1, 4.0),
            (0.3, 10.0, 4.2, 0.3, 0.1, 4.0),
            (0.3, float("nan"), 1.05, 0.3, 0.1, 4.0)):
    try:
        QuinticHermiteAxisLaw(*bad)
        check(f"不正入力 {bad[1]:.3g}/{bad[2]:.3g} 拒否", False)
    except ValueError:
        check(f"不正入力 {bad[1]:.3g}/{bad[2]:.3g} 拒否", True)

# --- 8. KnotQuinticAxisLaw (A6) --------------------------------------------------
# M6 R=3 相当の Hall アンカー (x_A 0.425, M_A 1.154, M'_A 0.550, M''_A 0.0565)
xA, MA, MpA, MppA, Md = 0.425, 1.154, 0.550, 0.0565, 6.0
err = 0.0
for _ in range(100):
    Lc = rng.uniform(20.0, 80.0)
    MK = rng.uniform(1.6, 3.5)
    k = KnotQuinticAxisLaw(xA, Lc, MA, MpA, MppA, Md, MK)
    e = 1e-9
    err = max(err, abs(k(xA) - MA), abs(k.deriv(xA, 1) - MpA), abs(k.deriv(xA, 2) - MppA),
              abs(k(k.x_E) - Md), abs(k.deriv(k.x_E, 1)), abs(k.deriv(k.x_E, 2)),
              abs(k(k.x_K) - MK),
              abs(k.deriv(k.x_K - e, 1) - k.deriv(k.x_K + e, 1)),
              abs(k.deriv(k.x_K - e, 2) - k.deriv(k.x_K + e, 2)),
              abs(k.deriv(k.x_K + e, 2)),                       # knot で M''=0
              abs(k.L_1 - 2.0 * (MK - MA) / (MpA + k.s_K)))     # 台形則
check(f"knot: 端点/knot C²/M''_K=0/台形則 最大誤差 {err:.2e}", err < 1e-6)
k = KnotQuinticAxisLaw(xA, 50.0, MA, MpA, MppA, Md, 2.5)
xs = np.linspace(k.x_K, k.x_E, 2001)
check("knot: 区分 2 は構成的に M'≥0・M''≤0",
      bool(np.all(k.deriv(xs, 1) >= -1e-12)) and bool(np.all(k.deriv(xs, 2) <= 1e-12)))
check("knot: x > x_E で M=M_d, M'=M''=0",
      abs(k(k.x_E + 5.0) - Md) < 1e-12 and k.deriv(k.x_E + 5.0, 1) == 0.0)
g = k.gates()
check(f"knot: ゲート合格・M'' 反転 {g['mpp_sign_changes']} ≤ 1", not g["violations"] and g["mpp_quality_ok"])
lo_q, hi_q = QuinticHermiteAxisLaw.admissible_Lc_range(xA, MA, MpA, MppA, Md)
lo_k, hi_k = KnotQuinticAxisLaw.admissible_Lc_range(xA, MA, MpA, MppA, Md, 2.5)
check(f"knot: 窓 ({lo_k:.2f}, {hi_k:.0f}) が quintic 上限 {hi_q:.1f} を超えて開く", hi_k > 60.0 > hi_q)
for bad in (dict(M_K=1.0), dict(M_K=6.5), dict(L_c=-1.0)):
    kw = dict(x_A=xA, L_c=50.0, M_A=MA, Mp_A=MpA, Mpp_A=MppA, M_d=Md, M_K=2.5); kw.update(bad)
    try:
        KnotQuinticAxisLaw(**kw)
        check(f"knot 不正入力 {bad} 拒否", False)
    except ValueError:
        check(f"knot 不正入力 {bad} 拒否", True)
# _axis_grid: 総長・初項・単調・等間隔退化
gx = _axis_grid(0.4, 177.0, 1200, 0.03)
d = np.diff(gx)
check(f"_axis_grid: 初項 {d[0]:.4f} / 総長 {gx[-1]:.3f} / 単調 / dx_max {d.max():.3f}",
      abs(d[0] - 0.03) < 1e-12 and abs(gx[-1] - 177.0) < 1e-9 and bool(np.all(d > 0)) and 0.14 < d.max() < 0.16
      and len(gx) == 1200)
check("_axis_grid: dx0=None は等間隔", np.allclose(_axis_grid(0.0, 10.0, 11, None), np.linspace(0, 10, 11)))
check("_axis_grid: dx0 が等間隔より粗いときは等間隔に退化",
      np.allclose(_axis_grid(0.0, 1.0, 11, 0.5), np.linspace(0, 1, 11)))

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
