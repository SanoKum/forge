"""準 1 次元等エントロピー初期値の貼り付け (uniform IC は圧力比で発散するため必須)。

case/29.bell_vs_conical/mesh/set_isentropic_ic.py と同方式だが、contour CSV
でなく NozzleWall から直接 r_w(x) を取る。x<0 亜音速根 / x>=0 超音速根。
"""
from __future__ import annotations

import h5py
import numpy as np


def area_ratio(M, g):
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** (
        (g + 1.0) / (2.0 * (g - 1.0))
    )


def invert_area_ratio(AR, supersonic, g):
    AR = np.maximum(AR, 1.0 + 1e-12)
    lo = np.where(supersonic, 1.0, 1e-4) * np.ones_like(AR)
    hi = np.where(supersonic, 50.0, 1.0) * np.ones_like(AR)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        f = area_ratio(mid, g) - AR
        go_hi = np.where(supersonic, f < 0.0, f > 0.0)
        lo = np.where(go_hi, mid, lo)
        hi = np.where(go_hi, hi, mid)
    return 0.5 * (lo + hi)


def paste_isentropic_ic(h5path, wall, scale, Pt, Tt, gamma, cp,
                        k_init=1.0, omega_init=18000.0, gas=None) -> dict:
    """wall: NozzleWall (無次元), scale: r* [m]。出口諸元 dict を返す。

    SST の roK/roOmega も貼る (既定 0 のままだと omega=0 で序盤 NaN →
    EOS 床洗浄後プラトーの指紋 — run_0001 で確認)。

    `gas` (semi-perfect, 2026-08-17): 面積比→M・T・P・ρ・roe を NASA-9 テーブルで
    作る。**roe は TP の内部エネルギー e(T)=h(T)−RT** で構成しないと forge (TP) が
    step 0 で温度を誤再構成する ([[wys-tp-divergence-is-cold-not-multispecies]] の
    IC 不整合と同じ罠)。エンタルピー基準は forge 内蔵 DB (絶対基準, thermoHrefTemp=0)
    と同一の NASA-9 なので整合する。
    """
    g = gamma
    R_gas = cp * (g - 1.0) / g
    with h5py.File(h5path, "r+") as f:
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
        xn = cc[:, 0] / scale  # 無次元軸位置
        # 物理壁 (A13) はスロートが (x_throat, r_throat) ≠ (0, 1) に動く
        x_thr = float(getattr(wall, "x_throat", 0.0))
        r_thr = float(getattr(wall, "r_throat", 1.0))
        AR = np.maximum(wall.r(xn) / r_thr, 1.0) ** 2
        if gas is not None and getattr(gas, "kind", "cpg") == "semiperfect":
            M = _invert_area_ratio_gas(AR, xn >= x_thr, gas)
            T = np.where(M >= 1.0, gas.T_of_M(np.maximum(M, 1.0)),
                         np.interp(np.minimum(M, 1.0), gas._Mu, gas._Tu))
            # 等エントロピー P/Pt = exp(∫ cp/(R T) dT) (Tt→T)
            P = Pt * _isentropic_pressure_ratio_gas(T, gas)
            R_gas = gas.R
            ro = P / (R_gas * T)
            gam_loc = gas.gamma(T)
            u = M * np.sqrt(gam_loc * R_gas * T)
            e_int = gas.h_mass(T) - R_gas * T                 # TP 内部エネルギー
            roe = ro * e_int + 0.5 * ro * u * u
        else:
            M = invert_area_ratio(AR, xn >= x_thr, g)
            fac = 1.0 + 0.5 * (g - 1.0) * M * M
            T = Tt / fac
            P = Pt / fac ** (g / (g - 1.0))
            ro = P / (R_gas * T)
            u = M * np.sqrt(g * R_gas * T)
            roe = P / (g - 1.0) + 0.5 * ro * u * u
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * u).astype(np.float32)
        f["/VALUE/roUy"][:] = np.zeros_like(ro, dtype=np.float32)
        f["/VALUE/roUz"][:] = np.zeros_like(ro, dtype=np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for name, val in (("roK", ro * k_init), ("roOmega", ro * omega_init)):
            v32 = val.astype(np.float32)
            if f"/VALUE/{name}" in f:
                f[f"/VALUE/{name}"][:] = v32
            else:
                f.create_dataset(f"/VALUE/{name}", data=v32)
    ie = int(np.argmax(xn))
    return {"M_exit_1d": float(M[ie]), "P_exit_1d": float(P[ie]), "T_exit_1d": float(T[ie])}

def _invert_area_ratio_gas(AR, supersonic, gas):
    """A/A* → M (ガスモデルのテーブル、超音速/亜音速枝)。"""
    AR = np.asarray(AR, dtype=float)
    M = np.empty_like(AR)
    sup = np.asarray(supersonic, bool)
    # 超音速枝: gas._AR は M 増加で単調増
    M[sup] = np.interp(AR[sup], gas._AR, gas._M)
    # 亜音速枝: gas._ARu は M 増加 (0→1) で単調減 → 反転して補間
    o = np.argsort(gas._ARu)
    M[~sup] = np.interp(AR[~sup], gas._ARu[o], gas._Mu[o])
    return M


def _isentropic_pressure_ratio_gas(T, gas):
    """P/Pt = exp( ∫_{Tt}^{T} cp/(R T') dT' ) を NASA-9 で数値積分 (T ごと)。"""
    T = np.asarray(T, dtype=float)
    Tg = np.linspace(gas.Tt, float(np.min(T)) * 0.98, 4000)
    cp = gas.cp_mass(Tg)
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (cp[1:] / Tg[1:] + cp[:-1] / Tg[:-1])
                                             * np.diff(Tg) / gas.R)])
    return np.exp(np.interp(T, Tg[::-1], integ[::-1]))
