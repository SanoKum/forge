#!/usr/bin/env python3
r"""準 1 次元等エントロピー流れを forge の H5 に初期値として貼り付ける.

uniform IC (P=101325 vs 入口 Pt=4 MPa) では 40 倍の圧力比で起動時に強い過渡
衝撃が立ち発散する。代わりに各セルの軸位置 x における壁半径 R(x) から局所
面積比 A/A* を求め、面積-Mach 関係を反転して局所 Mach を得、Pt,Tt からの
等エントロピー関係で (ro, P, u) を貼る。これで初期場が解の近傍になり起動が安定。

  - 収束部 (x<0): 亜音速根、発散部 (x>=0): 超音速根。スロート x=0 で M=1。
  - 速度は軸方向 (+x) のみ。Uy=Uz=0。

使い方 (run ディレクトリから):
    python ../mesh/set_isentropic_ic.py nozzle.h5 ../mesh/contour_bell.csv \
        --pt 4.0e6 --tt 1500 --rt 10.0
"""
from __future__ import annotations

import argparse

import h5py
import numpy as np


def area_ratio(M, g):
    """A/A* (M)。"""
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** (
        (g + 1.0) / (2.0 * (g - 1.0))
    )


def invert_area_ratio(AR, supersonic, g):
    """A/A* から Mach をベクトル化二分法で反転 (branch を mask で選択)。"""
    AR = np.maximum(AR, 1.0 + 1e-12)
    lo = np.where(supersonic, 1.0, 1e-4) * np.ones_like(AR)
    hi = np.where(supersonic, 50.0, 1.0) * np.ones_like(AR)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        f = area_ratio(mid, g) - AR
        # 亜音速枝では area_ratio は M で単調減少、超音速枝では単調増加
        if_inc = supersonic
        go_hi = np.where(if_inc, f < 0.0, f > 0.0)
        lo = np.where(go_hi, mid, lo)
        hi = np.where(go_hi, hi, mid)
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("contour", help="contour_*.csv (x_mm,y_mm[,zone]) で R(x) を与える")
    ap.add_argument("--pt", type=float, default=4.0e6, help="全圧 [Pa]")
    ap.add_argument("--tt", type=float, default=1500.0, help="全温 [K]")
    ap.add_argument("--rt", type=float, default=10.0, help="スロート半径 [mm]")
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--cp", type=float, default=1004.5)
    args = ap.parse_args()

    g, cp = args.gamma, args.cp
    R_gas = cp * (g - 1.0) / g

    cont = np.loadtxt(args.contour, delimiter=",", skiprows=1, usecols=(0, 1))
    cx, cy = cont[:, 0], cont[:, 1]  # mm
    order = np.argsort(cx)
    cx, cy = cx[order], cy[order]

    with h5py.File(args.h5, "r+") as f:
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
        x_m = cc[:, 0]            # m
        x_mm = x_m * 1000.0
        Rwall_mm = np.interp(x_mm, cx, cy)     # 局所壁半径 [mm]
        AR = (Rwall_mm / args.rt) ** 2         # A/A*
        supersonic = x_mm >= 0.0
        M = invert_area_ratio(AR, supersonic, g)

        fac = 1.0 + 0.5 * (g - 1.0) * M * M
        T = args.tt / fac
        P = args.pt / fac ** (g / (g - 1.0))
        ro = P / (R_gas * T)
        a = np.sqrt(g * R_gas * T)
        u = M * a

        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * u).astype(np.float32)
        f["/VALUE/roUy"][:] = np.zeros_like(ro, dtype=np.float32)
        f["/VALUE/roUz"][:] = np.zeros_like(ro, dtype=np.float32)
        f["/VALUE/roe"][:] = (P / (g - 1.0) + 0.5 * ro * u * u).astype(np.float32)

    # 出口 (最大 x) の等エントロピー諸元を報告 (outlet 背圧設定の目安)
    ie = int(np.argmax(x_mm))
    print(f"isentropic IC pasted: {len(M)} cells")
    print(f"  throat  M~1: P={args.pt/fac[np.argmin(np.abs(x_mm))]**(g/(g-1)):.0f} (near)")
    print(f"  exit  x={x_mm[ie]:.2f}mm: M={M[ie]:.3f}, P={P[ie]:.1f} Pa, T={T[ie]:.1f} K, u={u[ie]:.1f} m/s")
    print(f"  -> outlet 背圧 Ps は exit 静圧 {P[ie]:.0f} Pa より低く設定すること")


if __name__ == "__main__":
    main()
