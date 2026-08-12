"""2 目的 EHVI (期待ハイパーボリューム改善) の解析式と付随ユーティリティ。

規約: 全て**最小化** (η は -η で渡す)。予測分布は独立 Gaussian
N(μ1,σ1²)×N(μ2,σ2²) (Kriging の周辺分布)。

EHVI 閉形式は 2-D 区分和 (Hupkens–Deutz–Yang–Emmerich 2015 / Yang et al. 2017):
Pareto 前線を f1 降順に並べ、番兵 y0=(r1,−∞)・y_{n+1}=(−∞,r2) を足した
n+1 個の縦ストリップで

  EHVI = Σ_i [ (y_{i-1}^1 − y_i^1)·Φ((y_i^1−μ1)/σ1)·ψ(y_i^2,y_i^2)
             + (ψ(y_{i-1}^1,y_{i-1}^1) − ψ(y_{i-1}^1,y_i^1))·ψ(y_i^2,y_i^2) ]

  ψ(a,b;μ,σ) = σ·φ((b−μ)/σ) + (a−μ)·Φ((b−μ)/σ) = E[(a−Z)·1{Z<b}]

空前線では EHVI = E[(r1−Z1)+]·E[(r2−Z2)+] に退化する (単体テストの解析照合)。
数値検証は design/tests/run_opt_tests.py の Monte-Carlo 照合 (植え込み乱数)。
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm


def nondominated_mask(F) -> np.ndarray:
    """F (n,m) の非劣解 mask (最小化)。重複点は全て残す。"""
    F = np.asarray(F, dtype=float)
    n = F.shape[0]
    mask = np.ones(n, dtype=bool)
    for i in range(n):
        dom = np.all(F <= F[i], axis=1) & np.any(F < F[i], axis=1)
        mask[i] = not bool(np.any(dom))
    return mask


def hypervolume2d(F, ref) -> float:
    """2 目的 hypervolume (最小化, 参照点 ref)。ref 外の点は寄与しない。"""
    P = np.asarray(F, dtype=float).reshape(-1, 2)
    ref = np.asarray(ref, dtype=float)
    if P.size == 0:
        return 0.0
    P = P[nondominated_mask(P)]
    P = P[np.all(P < ref, axis=1)]
    if P.shape[0] == 0:
        return 0.0
    P = P[np.argsort(P[:, 0])]  # f1 昇順 → f2 降順 (階段)
    hv, prev2 = 0.0, ref[1]
    for f1, f2 in P:
        if f2 < prev2:  # 重複 f1 の劣位側はスキップ
            hv += (ref[0] - f1) * (prev2 - f2)
            prev2 = f2
    return float(hv)


def _psi(a, b, mu, s):
    """E[(a − Z)·1{Z < b}], Z ~ N(mu, s²)。"""
    t = (b - mu) / s
    return s * norm.pdf(t) + (a - mu) * norm.cdf(t)


def ehvi2d(mu, sigma, F, ref) -> float:
    """候補の予測 N(mu, sigma²) (独立 2 目的) の EHVI (最小化)。

    F: 評価済み点の目的値 (非劣解でなくてよい — 内部で前線化)。空も可。
    """
    mu = np.asarray(mu, dtype=float).ravel()
    s = np.maximum(np.asarray(sigma, dtype=float).ravel(), 1e-12)
    ref = np.asarray(ref, dtype=float)
    P = np.asarray(F, dtype=float).reshape(-1, 2)
    if P.shape[0]:
        P = P[nondominated_mask(P)]
        P = P[np.all(P < ref, axis=1)]
    if P.shape[0]:
        P = P[np.argsort(-P[:, 0])]  # f1 降順 (= f2 昇順)
        y1 = np.concatenate([[ref[0]], P[:, 0], [-np.inf]])
        y2 = np.concatenate([[-np.inf], P[:, 1], [ref[1]]])
    else:
        y1 = np.array([ref[0], -np.inf])
        y2 = np.array([-np.inf, ref[1]])

    total = 0.0
    for i in range(1, len(y1)):
        psi2 = _psi(y2[i], y2[i], mu[1], s[1])
        if np.isfinite(y1[i]):
            a_term = (y1[i - 1] - y1[i]) * norm.cdf((y1[i] - mu[0]) / s[0]) * psi2
        else:
            a_term = 0.0  # Φ(−∞)=0 が線形発散に勝つ (極限 0)
        b_term = (_psi(y1[i - 1], y1[i - 1], mu[0], s[0])
                  - _psi(y1[i - 1], y1[i], mu[0], s[0])) * psi2
        total += a_term + b_term
    return float(max(total, 0.0))
