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


def exit_uniformity(mesh_h5, res_h5, M_design, core_frac: float = 0.85,
                    x_plane=None, band=None) -> dict:
    """出口一様性 $\\varepsilon_M$ / $\\varepsilon_\\theta$ (サーベイ B4.5 定義)。

    - $\\varepsilon_M$ = テストコア上の**質量流束重み RMS 偏差**
      $\\sqrt{\\langle (M-M_d)^2\\rangle_{\\rho u}}/M_d$ (RMS を主指標、max は副指標)
    - $\\varepsilon_\\theta$ = コア上の $\\max|\\theta|$ [deg] (軸からの流れ角。
      $\\varepsilon_M$ とは**独立の目的** — マッハだけ見ると軸ズレ流れを見逃す)
    - 重みは質量流束 $\\rho u$ (模型が実際に「見る」流れ)。幾何・積分系 (推力) の
      面積重みとは区別する規約。

    測定面は既定でリップ直前 (x_plane=None → x_max−0.3·r*)。テストコアは
    有効菱形の簡易代理として半径の core_frac 倍以内 (正式な菱形幾何は将来課題)。
    """
    with h5py.File(mesh_h5, "r") as nz:
        cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res_h5, "r") as f:
        Ux, Uy = f["/VALUE/Ux"][:], f["/VALUE/Uy"][:]
        son, ro = f["/VALUE/sonic"][:], f["/VALUE/ro"][:]
    x, r = cc[:, 0], cc[:, 1]
    M = np.hypot(Ux, Uy) / np.maximum(son, 1e-9)
    th_deg = np.degrees(np.arctan2(Uy, Ux))
    r_ref = float(r.max())
    xp = float(x.max() - 0.03 * r_ref) if x_plane is None else float(x_plane)
    bw = 0.02 * r_ref if band is None else float(band)
    m = np.abs(x - xp) < bw
    if m.sum() < 8:
        raise ValueError(f"測定面 x={xp:.4g} にセルが不足 ({m.sum()})")
    rr, MM, tt = r[m], M[m], th_deg[m]
    w = ro[m] * Ux[m]
    core = rr < core_frac * rr.max()
    Mw = float(np.average(MM[core], weights=w[core]))
    eM = float(np.sqrt(np.average((MM[core] - M_design) ** 2, weights=w[core])) / M_design)
    return {"x_plane_m": xp, "core_frac": core_frac,
            "eps_M_rms": eM,
            "eps_M_max": float(np.max(np.abs(MM[core] - M_design)) / M_design),
            "eps_theta_max_deg": float(np.max(np.abs(tt[core]))),
            "M_core_massflux_avg": Mw,
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
