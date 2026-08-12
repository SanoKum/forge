#!/usr/bin/env python3
"""forge_design.opt 単体テスト (design/.venv-opt/bin/python で実行)。

山場は EHVI 解析式の検証: (1) 空前線の解析解、(2) σ→0 の決定的極限 = HVI、
(3) Monte-Carlo 照合 (植え込み乱数)、(4) ZDT1 ミニ BO ループの hypervolume 進展。
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.opt.doe import lhs  # noqa: E402
from forge_design.opt.ehvi import ehvi2d, hypervolume2d, nondominated_mask  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


# --- hypervolume / 非劣解 -----------------------------------------------------
F = np.array([[1.0, 3.0], [2.0, 2.0], [3.0, 1.0], [2.5, 2.5]])
check("nondominated: 劣位点を除外", list(nondominated_mask(F)) == [True, True, True, False])
check("hypervolume2d 既知値 (=6)", abs(hypervolume2d(F, (4.0, 4.0)) - 6.0) < 1e-12)
check("hypervolume2d ref 外は寄与 0", hypervolume2d([[5.0, 5.0]], (4.0, 4.0)) == 0.0)

# --- EHVI: 空前線の解析解 -------------------------------------------------------
from scipy.stats import norm  # noqa: E402

e_ana = (norm.pdf(1.0) + 1.0 * norm.cdf(1.0)) ** 2  # E[(1-Z)+]², Z~N(0,1)
check("EHVI 空前線 = E[(r1-Z)+]E[(r2-Z)+]",
      abs(ehvi2d((0, 0), (1, 1), np.empty((0, 2)), (1.0, 1.0)) - e_ana) < 1e-10)

# --- EHVI: σ→0 極限 = 決定的 HVI ------------------------------------------------
front = np.array([[1.0, 3.0], [2.0, 2.0], [3.0, 1.0]])
ref = (4.0, 4.0)
for z in ([1.5, 1.5], [0.5, 3.5], [3.5, 0.5], [2.5, 2.5], [5.0, 5.0]):
    hvi = hypervolume2d(np.vstack([front, z]), ref) - hypervolume2d(front, ref)
    e = ehvi2d(z, (1e-9, 1e-9), front, ref)
    check(f"EHVI σ→0 = HVI @z={z} ({e:.4f} vs {hvi:.4f})", abs(e - hvi) < 1e-6)

# --- EHVI: Monte-Carlo 照合 -----------------------------------------------------
rng = np.random.default_rng(20260813)
N = 400_000
for mu, sg in ([(1.5, 1.5), (0.6, 0.5)], [(3.0, 3.0), (1.0, 1.0)],
               [(0.5, 0.8), (0.3, 0.9)], [(4.5, 4.5), (1.2, 1.2)]):
    Z = rng.normal(mu, sg, size=(N, 2))
    # ベクトル化 HVI: 階段 (f1 昇順) に対する z の改善面積
    P = front[np.argsort(front[:, 0])]
    x_edges = np.concatenate([[-np.inf], P[:, 0], [ref[0]]])
    h = np.concatenate([[np.inf], P[:, 1]])  # 各ストリップの階段高
    hvi_mc = np.zeros(N)
    for i in range(len(x_edges) - 1):
        wlo = np.clip(Z[:, 0], x_edges[i], None)
        w = np.clip(x_edges[i + 1], None, ref[0]) - np.clip(wlo, None, ref[0])
        hgt = np.clip(np.minimum(h[i], ref[1]) - Z[:, 1], 0.0, None)
        hvi_mc += np.clip(w, 0.0, None) * hgt
    e_mc = float(np.mean(hvi_mc))
    e = ehvi2d(mu, sg, front, ref)
    tol = max(0.02 * e_mc, 3e-3)
    check(f"EHVI MC 照合 μ={mu} σ={sg} ({e:.4f} vs {e_mc:.4f})", abs(e - e_mc) < tol)

# --- LHS 決定性 -----------------------------------------------------------------
b = np.array([[0.0, 1.0], [10.0, 20.0], [-1.0, 1.0]])
x1, x2 = lhs(b, 16, seed=42), lhs(b, 16, seed=42)
check("LHS 決定性 (同一 seed → bit 同一)", np.array_equal(x1, x2))
check("LHS bounds 内", bool(np.all(x1 >= b[:, 0]) and np.all(x1 <= b[:, 1])))
check("LHS 層化 (第 1 軸: 各 1/16 区間に 1 点)",
      len(np.unique((x1[:, 0] * 16).astype(int))) == 16)

# --- KRG サロゲート (要 smt) -----------------------------------------------------
from forge_design.opt.surrogate import KrigingSet  # noqa: E402

bounds = np.array([[0.0, 1.0]] * 3)
Xtr = lhs(bounds, 24, seed=1)
Ytr = np.column_stack([np.sin(3 * Xtr[:, 0]) + Xtr[:, 1] ** 2,
                       (Xtr[:, 0] - 0.5) ** 2 + 0.5 * Xtr[:, 2]])
ks = KrigingSet(bounds).fit(Xtr, Ytr)
mu_tr, sg_tr = ks.predict(Xtr)
check("KRG 学習点を補間 (max err < 1e-5)", float(np.max(np.abs(mu_tr - Ytr))) < 1e-5)
check("KRG 学習点で σ≈0", float(np.max(sg_tr)) < 1e-4)

# --- infill 提案: 決定性・境界 ----------------------------------------------------
from forge_design.opt.moo import propose_infill  # noqa: E402

ref2 = (3.0, 3.0)
p1 = propose_infill(ks, Xtr, Ytr, ref2, bounds, n_infill=2, seed=5,
                    pop_size=32, n_gen=20, n_explore=256)
p2 = propose_infill(ks, Xtr, Ytr, ref2, bounds, n_infill=2, seed=5,
                    pop_size=32, n_gen=20, n_explore=256)
check("infill 決定性 (同一 seed → 同一提案)", np.array_equal(p1, p2))
check("infill 提案数 = 2, bounds 内",
      p1.shape == (2, 3) and bool(np.all(p1 >= bounds[:, 0]) and np.all(p1 <= bounds[:, 1])))

# --- ZDT1 ミニ BO ループ (照合ベンチ) ---------------------------------------------
def zdt1(X):
    X = np.atleast_2d(X)
    f1 = X[:, 0]
    g = 1.0 + 9.0 * np.mean(X[:, 1:], axis=1)
    f2 = g * (1.0 - np.sqrt(f1 / g))
    return np.column_stack([f1, f2])

D = 4
bz = np.array([[0.0, 1.0]] * D)
refz = (11.0, 11.0)
hv_true = 11.0 * 11.0 - (11.0 * 1.0 - (10.0 + 2.0 / 3.0))  # = 120.6667 (真の前線)
Xe = lhs(bz, 24, seed=3)
Fe = zdt1(Xe)
hv0 = hypervolume2d(Fe, refz)
for it in range(14):
    kz = KrigingSet(bz).fit(Xe, Fe)
    xs = propose_infill(kz, Xe, Fe, refz, bz, n_infill=2, seed=100 + it,
                        pop_size=48, n_gen=40, n_explore=512)
    Xe = np.vstack([Xe, xs])
    Fe = np.vstack([Fe, zdt1(xs)])
hv1 = hypervolume2d(Fe, refz)
print(f"     ZDT1: HV {hv0:.3f} -> {hv1:.3f} (真値 {hv_true:.3f}), 評価 {len(Fe)} 点")
check("ZDT1 BO: HV が有意に進展 (>hv0+1)", hv1 > hv0 + 1.0)
check("ZDT1 BO: 真値の 99% 到達", hv1 > 0.99 * hv_true)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
