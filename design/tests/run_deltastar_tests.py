#!/usr/bin/env python3
"""排除厚さ更新計画 (plans/active/tooling-nozzle-deltastar-core-matched-euler.md) の単体テスト
(pytest 非依存・素の assert)。

対象:
  1. `core_matched_deficit` (固定 Euler 基準・コア整合質量欠損 → 半径方向等価排除厚)
     - q_NS = c·q_E なら δ_r = 0
     - 既知の壁側欠損から円環厚さを復元 (一様コア / 非一様コア)
     - Euler 壁と NS 壁が異なっても復元 (Euler 壁の外は壁値外挿)
     - 点数と壁近傍伸長率への収束
     - コア振幅差は α で除去、コア形状差はゲート (core_rms) で検出
     - コア範囲 25/30/35 % 感度
  2. 積分法初期推定器 (`feedback/deltastar_integral.py`) — 実装後に追加
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.metrics.deltastar import core_matched_deficit  # noqa: E402

FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def wall_clustered_grid(rw, n=97, first_frac=4.5e-5, stretch=None):
    """壁側に幾何級数で集中した半径格子 (mesh2d._radial_fracs と同趣旨)。"""
    # 壁からの距離 y を等比級数で
    if stretch is None:
        # first_frac·rw が最初の区間になる伸長率を二分法で
        lo, hi = 1.0001, 1.5
        for _ in range(100):
            g = 0.5 * (lo + hi)
            tot = first_frac * (g ** (n - 1) - 1) / (g - 1)
            lo, hi = (g, hi) if tot < 1.0 else (lo, g)
        stretch = g
    dy = first_frac * stretch ** np.arange(n - 1)
    y = np.concatenate([[0.0], np.cumsum(dy)])
    y = y / y[-1]                       # 0 (壁) → 1 (軸)
    return (rw * (1.0 - y))[::-1]       # 軸 → 壁


def profile_with_annulus_deficit(r, q_core_fn, rw_ns, delta_r):
    """壁側の円環 [rw_ns-δ_r, rw_ns] を完全欠損 (q=0)、内側は q_core_fn(r)。
    → 等価排除厚は厳密に δ_r。"""
    q = q_core_fn(r).copy()
    q[r > rw_ns - delta_r] = 0.0
    return q


def profile_smooth_bl(r, q_core_fn, rw_ns, delta99, n_pow=7.0):
    """べき乗則境界層 (q/q_e = (y/δ)^{1/n})。等価排除厚は数値積分で別途求める。"""
    y = rw_ns - r
    f = np.where(y < delta99, (np.maximum(y, 0) / delta99) ** (1.0 / n_pow), 1.0)
    return q_core_fn(r) * f


def equiv_delta_exact(r_fine, q_ref, q_ns, rw):
    """参照解: 十分細かい格子で D と r_eff を直接計算。"""
    D = 2 * np.pi * np.trapezoid((q_ref - q_ns) * r_fine, r_fine)
    seg = 0.5 * (q_ref[1:] * r_fine[1:] + q_ref[:-1] * r_fine[:-1]) * np.diff(r_fine)
    F = 2 * np.pi * np.concatenate([[0.0], np.cumsum(seg[::-1])])
    return rw - float(np.interp(D, F, r_fine[::-1]))


# --- 1. q_NS = c·q_E → δ_r = 0 ------------------------------------------------------
rw = 5.0
r = wall_clustered_grid(rw)
qE = lambda r_: 1.0 + 0.2 * (r_ / rw) ** 2          # 非一様コア (壁側が高い)
res = core_matched_deficit(r, 1.7 * qE(r), qE(r), rw)
check("比例プロファイル: δ_r = 0", abs(res["delta_r"]) < 1e-12, f"{res['delta_r']:.2e}")
check("比例プロファイル: α = 1.7", abs(res["alpha"] - 1.7) < 1e-12)
check("比例プロファイル: コア RMS = 0", res["core_rms_noaxis"] < 1e-12)

# --- 2. 既知の円環欠損 (一様コア) ---------------------------------------------------
qU = lambda r_: np.ones_like(r_)
for d_true in (0.01, 0.05, 0.2):
    qN = profile_with_annulus_deficit(r, qU, rw, d_true)
    res = core_matched_deficit(r, qN, qU(r), rw)
    check(f"一様コア 円環欠損 δ_r={d_true}", abs(res["delta_r"] - d_true) < 2e-3 * rw,
          f"got {res['delta_r']:.5f}")

# --- 3. 非一様コア + 円環欠損 (α 一致・形状一致) --------------------------------------
for d_true in (0.02, 0.1):
    qN = profile_with_annulus_deficit(r, lambda r_: 1.3 * qE(r_), rw, d_true)
    res = core_matched_deficit(r, qN, qE(r), rw)
    check(f"非一様コア (×1.3) 円環欠損 δ_r={d_true}", abs(res["delta_r"] - d_true) < 2e-3 * rw,
          f"got {res['delta_r']:.5f}, α {res['alpha']:.4f}")
    check(f"  α が振幅差 1.3 を吸収", abs(res["alpha"] - 1.3) < 1e-9)
    check(f"  欠損は壁側に集中 (outer_share>0.8)", res["outer_share"] > 0.8, f"{res['outer_share']:.3f}")

# --- 4. Euler 壁 ≠ NS 壁 (NS が δ_in だけ太い; Euler 壁の外は壁値外挿) -----------------
rw_e = 5.0
for d_in in (0.05, 0.3):
    rw_n = rw_e + d_in
    rn = wall_clustered_grid(rw_n)
    qE_map = qE(np.minimum(rn, rw_e))                # Euler 壁の外は壁値
    d_true = 0.12
    qN = profile_with_annulus_deficit(rn, qE, rw_n, d_true)
    res = core_matched_deficit(rn, qN, qE_map, rw_e)
    check(f"壁半径差 δ_in={d_in}: 円環欠損 {d_true} を復元", abs(res["delta_r"] - d_true) < 3e-3 * rw_n,
          f"got {res['delta_r']:.5f}")
    check(f"  delta_in = {d_in}", abs(res["delta_in"] - d_in) < 1e-12)

# --- 5. 滑らかな境界層: 点数・伸長率への収束 ----------------------------------------
rw = 5.0; d99 = 0.4
r_fine = np.linspace(0, rw, 200001)
q_ref_f = qE(r_fine); q_ns_f = profile_smooth_bl(r_fine, qE, rw, d99)
d_exact = equiv_delta_exact(r_fine, q_ref_f, q_ns_f, rw)
errs = []
for n in (65, 97, 193, 385):
    rn = wall_clustered_grid(rw, n=n)
    res = core_matched_deficit(rn, profile_smooth_bl(rn, qE, rw, d99), qE, rw)
    errs.append(abs(res["delta_r"] - d_exact) / d_exact)
check("べき乗則 BL: n=97 で参照解 ±3 %", errs[1] < 0.03, f"exact {d_exact:.5f}, rel errs {['%.2e' % e for e in errs]}")
check("べき乗則 BL: 点数で単調収束", errs[0] > errs[1] > errs[2] > errs[3] * 0.999,
      f"{['%.2e' % e for e in errs]}")
errs_s = []
for ff in (2e-4, 4.5e-5, 1e-5):
    rn = wall_clustered_grid(rw, n=97, first_frac=ff)
    res = core_matched_deficit(rn, profile_smooth_bl(rn, qE, rw, d99), qE, rw)
    errs_s.append(abs(res["delta_r"] - d_exact) / d_exact)
check("べき乗則 BL: 壁近傍伸長率を変えても ±3 %", max(errs_s) < 0.03, f"{['%.2e' % e for e in errs_s]}")

# --- 6. コア形状差はゲートで検出 (振幅差は検出しない) --------------------------------
r = wall_clustered_grid(5.0)
qN_amp = 1.02 * qE(r)                                          # 振幅のみ 2 %
res_amp = core_matched_deficit(r, qN_amp, qE(r), 5.0)
qN_shape = qE(r) * (1.0 + 0.03 * (r / 5.0))                    # 形状 (勾配) 差 3 %
res_shape = core_matched_deficit(r, qN_shape, qE(r), 5.0)
check("振幅差 2 %: コア RMS ≈ 0 (α が吸収)", res_amp["core_rms_noaxis"] < 1e-12)
check("形状差 3 %: コア RMS > 0.1 % で検出", res_shape["core_rms_noaxis"] > 1e-3,
      f"{res_shape['core_rms_noaxis']:.4f}")

# --- 7. コア範囲 25/30/35 % 感度 (円環欠損なら不変) --------------------------------------
qN = profile_with_annulus_deficit(r, lambda r_: 1.1 * qE(r_), 5.0, 0.08)
ds = [core_matched_deficit(r, qN, qE(r), 5.0, core_frac=fc)["delta_r"] for fc in (0.25, 0.30, 0.35)]
check("コア範囲 25/30/35 %: δ_r 不変 (<0.5 %)", max(abs(d - ds[1]) for d in ds) / ds[1] < 5e-3,
      f"{['%.5f' % d for d in ds]}")

# --- 8. 負の欠損 (NS がコアより多く流す) は negative_deficit を返す -------------------
res_neg = core_matched_deficit(r, qE(r) * (1.0 + 0.05 * (r / 5.0) ** 4), qE(r), 5.0)
check("負の欠損: reason=negative_deficit, δ_r<0", "negative_deficit" in res_neg["reason"] and res_neg["delta_r"] < 0,
      f"{res_neg['delta_r']:.4f}")

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
