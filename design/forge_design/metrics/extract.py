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
                   outlet_physid=2, wall_physid=3) -> dict:
    """推力・C_F・η (= C_F / C_F,ideal 同面積比)。壁摩擦は twall があれば加算。"""
    g = gamma
    with h5py.File(mesh_h5, "r") as nz:
        ip = nz[f"/BCONDS/{outlet_physid}/iPlanes"][:]
        ic = nz[f"/BCONDS/{outlet_physid}/iCells"][:]
        pc = nz["/PLANES/centCoords"][:].reshape(-1, 3)
        rf = pc[ip, 1]
        dr = nz["/PLANES/surfArea"][:][ip]
        dA = 2.0 * np.pi * rf * dr
        # 壁摩擦 (軸方向) — BCONDS 出力があれば
        fric = 0.0
        wgrp = f"/BCONDS/{wall_physid}"
        if wgrp + "/VALUE/twall_x" in nz:
            twx = nz[wgrp + "/VALUE/twall_x"][:]
            wip = nz[wgrp + "/iPlanes"][:]
            wdA = 2.0 * np.pi * pc[wip, 1] * nz["/PLANES/surfArea"][:][wip]
            fric = float(np.sum(twx * wdA))
    with h5py.File(res_h5, "r") as f:
        ro = f["/VALUE/ro"][:][ic]
        u = f["/VALUE/roUx"][:][ic] / ro
        v = f["/VALUE/roUy"][:][ic] / ro
        roe = f["/VALUE/roe"][:][ic]
        P = (g - 1.0) * (roe - 0.5 * ro * (u * u + v * v))
    F = float(np.sum((ro * u * u + (P - p_ambient)) * dA)) + fric
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
