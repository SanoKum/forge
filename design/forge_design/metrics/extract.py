"""res_*.h5 からの目的関数抽出 (③ベル最小セット)。

規約 (サーベイ B4.5 / case/29 postprocess_thrust.py):
- 出口面の状態は境界ダンプでなく res h5 を BCONDS/<outlet>/iCells で引く。
- 軸対称の PLANES/surfArea は 2D 辺長 (dr) — dA = 2πr·dr を手で構成する。
- P は roe から再構成 (CPG)。
- 全指標は最小化形でなく生値も併記 (η は最大化量なので optimizer へは -η)。
"""
from __future__ import annotations

import h5py
import numpy as np

from ..evaluate.ic import invert_area_ratio


def thrust_metrics(mesh_h5, res_h5, gamma, Pt, p_ambient, r_throat,
                   outlet_physid=2) -> dict:
    """推力・C_F・η (= C_F / C_F,ideal 同面積比)。

    出口面の運動量法 F = ∫(ρu² + (P - Pa)) dA。壁摩擦・壁圧の影響は出口状態に
    既に織り込まれているため**別加算しない** (加算すると二重計上 — サーベイ
    B4.5 の「+壁摩擦」はこの抽出法では不採用。子 plan §9 参照)。
    """
    g = gamma
    with h5py.File(mesh_h5, "r") as nz:
        ip = nz[f"/BCONDS/{outlet_physid}/iPlanes"][:]
        ic = nz[f"/BCONDS/{outlet_physid}/iCells"][:]
        pc = nz["/PLANES/centCoords"][:].reshape(-1, 3)
        rf = pc[ip, 1]
        dr = nz["/PLANES/surfArea"][:][ip]  # 軸対称: 2D 辺長 (dr)
        dA = 2.0 * np.pi * rf * dr
    with h5py.File(res_h5, "r") as f:
        ro = f["/VALUE/ro"][:][ic]
        u = f["/VALUE/roUx"][:][ic] / ro
        v = f["/VALUE/roUy"][:][ic] / ro
        roe = f["/VALUE/roe"][:][ic]
        P = (g - 1.0) * (roe - 0.5 * ro * (u * u + v * v))
    F = float(np.sum((ro * u * u + (P - p_ambient)) * dA))
    At = np.pi * r_throat ** 2
    Ae = float(np.sum(dA))
    CF = F / (Pt * At)
    CF_id = cf_ideal(Ae / At, g, p_ambient / Pt)
    return {
        "thrust_N": F,
        "CF": CF,
        "CF_ideal": CF_id,
        "eta_cf": CF / CF_id,
        "eps_measured": Ae / At,
        "mdot_kg_s": float(np.sum(ro * u * dA)),
    }


def cf_ideal(eps, gamma, pa_over_pt) -> float:
    """1D 等エントロピーの理想推力係数 (同面積比・同背圧)。"""
    g = gamma
    Me = float(invert_area_ratio(np.array([eps]), np.array([True]), g)[0])
    pe_pt = (1.0 + 0.5 * (g - 1.0) * Me * Me) ** (-g / (g - 1.0))
    gterm = np.sqrt(
        (2.0 * g * g / (g - 1.0))
        * (2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))
        * (1.0 - pe_pt ** ((g - 1.0) / g))
    )
    return float(gterm + (pe_pt - pa_over_pt) * eps)


def _annular_areas(r_sorted: np.ndarray) -> np.ndarray:
    """半径昇順のサンプル列に対する環状面積 π(r_out²−r_in²)。

    半径方向が非一様なので単純な 2πr·dr では重みが崩れる (セル幅が違う)。
    境界は隣接サンプルの中点、両端は片側幅を鏡像で外挿する。"""
    r = np.asarray(r_sorted, dtype=float)
    mid = 0.5 * (r[:-1] + r[1:])
    lo = np.concatenate([[max(2.0 * r[0] - mid[0], 0.0)], mid])
    hi = np.concatenate([mid, [2.0 * r[-1] - mid[-1]]])
    return np.pi * (np.maximum(hi, 0.0) ** 2 - np.maximum(lo, 0.0) ** 2)


def test_core_radius(x_plane: float, x_d: float, M_design: float,
                     r_wall: float, gamma: float = 1.4) -> float:
    r"""**有効菱形 (テストコア) の幾何定義** — サーベイ B4.5。

    一様出口設計では、軸が $M_d$ に達する点 $(x_d,0)$ から出る C⁺ より下が
    一様域。一様域では傾きが $\tan\mu_d$ 一定の直線なので

        r_core(x) = (x − x_d)·tan μ_d,   μ_d = arcsin(1/M_d)

    (壁半径で上限クリップ)。設計上リップ平面でちょうど壁に達する。
    x, x_d, r_wall は同一単位 (無次元 r* でも実寸でも可)。"""
    mu = np.arcsin(1.0 / max(float(M_design), 1.0 + 1e-12))
    return float(min(max((float(x_plane) - float(x_d)) * np.tan(mu), 0.0), float(r_wall)))


def exit_uniformity(mesh_h5, res_h5, M_design, x_d=None, x_plane=None,
                    core_radius=None, gamma=1.4) -> dict:
    r"""出口一様性 $\varepsilon_M$ / $\varepsilon_\theta$ (サーベイ B4.5 定義)。

    - $\varepsilon_M$ = テストコア上の**質量流束重み RMS 偏差**
      $\sqrt{\langle (M-M_d)^2\rangle_w}/M_d$。重みは $w_i=\rho_i u_{x,i} A_i$ で
      $A_i$ は**環状面積** $\pi(r_{out}^2-r_{in}^2)$ (半径方向が非一様なので
      $\rho u$ だけではメッシュ密度に依存する)。RMS が主指標、max は副指標。
    - $\varepsilon_\theta$ = コア上の $\max|\theta|$ [deg]。$\varepsilon_M$ とは
      **独立の目的** (マッハだけ見ると軸ズレ流れを見逃す)。

    測定面: `x_plane` 未指定なら**最下流のセル列** (x.max() の列) をそのまま使う。
    テストコア: `core_radius` 未指定かつ `x_d` 指定時は `test_core_radius()` の
    幾何定義。どちらも無ければ壁半径 (=コア制限なし) とし、その旨を返り値に残す。
    """
    with h5py.File(mesh_h5, "r") as nz:
        cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res_h5, "r") as f:
        Ux, Uy = f["/VALUE/Ux"][:], f["/VALUE/Uy"][:]
        son, ro = f["/VALUE/sonic"][:], f["/VALUE/ro"][:]
    x, r = cc[:, 0], cc[:, 1]
    M = np.hypot(Ux, Uy) / np.maximum(son, 1e-9)
    th_deg = np.degrees(np.arctan2(Uy, Ux))
    # 測定面 = 最下流のセル列 (構造格子なので x がほぼ揃う)
    xp = float(x.max()) if x_plane is None else float(x_plane)
    xs = np.unique(np.round(x, 12))
    tol = 0.25 * float(np.min(np.diff(xs))) if len(xs) > 1 else 1e-9
    m = np.abs(x - xp) <= tol
    if m.sum() < 8:
        raise ValueError(f"測定面 x={xp:.6g} のセルが不足 ({int(m.sum())})")
    o = np.argsort(r[m])
    rr, MM, tt = r[m][o], M[m][o], th_deg[m][o]
    w_rho_u, A = ro[m][o] * Ux[m][o], _annular_areas(r[m][o])
    r_wall = float(rr.max())
    if core_radius is not None:
        rc, how = float(core_radius), "explicit"
    elif x_d is not None:
        rc, how = test_core_radius(xp, x_d, M_design, r_wall, gamma), "effective_rhombus"
    else:
        rc, how = r_wall, "no_core_limit"
    core = rr <= rc
    if core.sum() < 4:
        raise ValueError(f"コア内セルが不足 ({int(core.sum())}, r_core={rc:.4g})")
    w = w_rho_u[core] * A[core]
    Md = float(M_design)
    eM = float(np.sqrt(np.average((MM[core] - Md) ** 2, weights=w)) / Md)
    return {"x_plane": xp, "r_core": rc, "core_def": how, "r_wall": r_wall,
            "eps_M_rms": eM,
            "eps_M_max": float(np.max(np.abs(MM[core] - Md)) / Md),
            "eps_theta_max_deg": float(np.max(np.abs(tt[core]))),
            "M_core_massflux_avg": float(np.average(MM[core], weights=w)),
            "n_core_cells": int(core.sum())}


def axis_mach(mesh_h5, res_h5, axis_band) -> tuple:
    """軸近傍セル帯 (cy < axis_band [m]) の (x, M) を x 昇順で返す。"""
    with h5py.File(mesh_h5, "r") as nz:
        cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res_h5, "r") as f:
        Ux = f["/VALUE/Ux"][:]
        Uy = f["/VALUE/Uy"][:]
        son = f["/VALUE/sonic"][:]
    m = cc[:, 1] < axis_band
    x = cc[m, 0]
    M = np.hypot(Ux[m], Uy[m]) / np.maximum(son[m], 1e-9)
    o = np.argsort(x)
    return x[o], M[o]
