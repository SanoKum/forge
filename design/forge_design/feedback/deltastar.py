"""排除厚さ δ* の経験式 (v1 用の決定的初期補正 — 親計画 §4.7 v1)。

乱流平板相関 + Eckert 参照温度 (断熱壁) の局所適用:
    δ*(s) = 0.0463 · s / Re*_s^{1/5} · H_c(Me)
    Re*_s = ρ* u_e s / μ*  (ρ*, μ* は参照温度 T* で評価)
    T*/Te = 1 + 0.115 (γ-1) Me²  (Eckert, 断熱, 回復係数 r≈0.89)
    H_c(Me) = 1 + 0.4 · (γ-1)/2 · Me²  (圧縮性の排除厚さ増大, 粗い係数)

v1 の役割は「帰還の初期壁の質」だけであり、相関の精度は v2 の収束速度にのみ
影響する (親計画 B3 の論法)。s はスロートからの壁弧長。
"""
from __future__ import annotations

import numpy as np


def _sutherland(T):
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


def dstar_flatplate(s_m, M, Pt: float, Tt: float,
                    gamma: float = 1.4, cp: float = 1004.5) -> np.ndarray:
    """乱流平板 + Eckert 参照温度の δ* [m]。s_m: 弧長 [m], M: 縁 Mach (配列)。

    `deltastar_offset` のコアを関数化したもの (数値は同一)。A13 では s を
    **入口からの弧長**で渡す (境界層の発達履歴 — スロートで δ*>0 になる)。"""
    g = gamma
    R_gas = cp * (g - 1.0) / g
    M = np.maximum(np.asarray(M, dtype=float), 1e-3)
    s_m = np.maximum(np.asarray(s_m, dtype=float), 1e-9)
    Te = Tt / (1.0 + 0.5 * (g - 1.0) * M ** 2)
    pe = Pt * (Te / Tt) ** (g / (g - 1.0))
    ue = M * np.sqrt(g * R_gas * Te)
    Tstar = Te * (1.0 + 0.115 * (g - 1.0) * M ** 2)
    rho_star = pe / (R_gas * Tstar)
    mu_star = _sutherland(Tstar)
    Re_s = np.maximum(rho_star * ue * s_m / mu_star, 1e3)
    Hc = 1.0 + 0.4 * 0.5 * (g - 1.0) * M ** 2
    return 0.0463 * s_m * Re_s ** (-0.2) * Hc


def deltastar_offset(wall, rt_m: float, Pt: float, Tt: float,
                     gamma: float = 1.4, cp: float = 1004.5) -> np.ndarray:
    """逆設計壁 wall=(n,4)[x,r,θ,M] (r*=1 無次元) に δ* 法線オフセットを加える。

    戻り値: (n,2) [x, r] — 物理壁 (無次元)。"""
    g = gamma
    R_gas = cp * (g - 1.0) / g
    x, r, th, M = wall[:, 0], wall[:, 1], wall[:, 2], wall[:, 3]
    # 壁弧長 (スロート近似: テーブル先頭から)
    s = np.concatenate([[1e-6], np.cumsum(np.hypot(np.diff(x), np.diff(r)))]) * rt_m
    Te = Tt / (1.0 + 0.5 * (g - 1.0) * M ** 2)
    pe = Pt * (Te / Tt) ** (g / (g - 1.0))
    ue = M * np.sqrt(g * R_gas * Te)
    Tstar = Te * (1.0 + 0.115 * (g - 1.0) * M ** 2)
    rho_star = pe / (R_gas * Tstar)
    mu_star = _sutherland(Tstar)
    Re_s = np.maximum(rho_star * ue * s / mu_star, 1e3)
    Hc = 1.0 + 0.4 * 0.5 * (g - 1.0) * M ** 2
    dstar = 0.0463 * s * Re_s ** (-0.2) * Hc / rt_m  # 無次元へ戻す
    return np.c_[x - dstar * np.sin(th), r + dstar * np.cos(th)]
