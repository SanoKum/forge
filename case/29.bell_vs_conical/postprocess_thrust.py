#!/usr/bin/env python3
r"""出口面で推力を積分し Rao ベルと C1 コニカルを比較する (軸対称).

**面積・値の取り方 (重要)**: forge の軸対称は 2D 平面方程式 + 軸対称ソース項で解くため、
`/PLANES/surfArea` は回転面積ではなく **2D エッジ長 (dr)** を格納する。物理的な出口
環状面積は手動で回転させる:  dA = 2*pi*r_face*dr。

出口面の状態は境界ダンプ res_outlet_*.h5 ではなく **隣接内部セル値** (res_<step>.h5 の
/VALUE) を用いる (ダンプの COORD/CONNE には巻き込み face と 1-offset があり信頼できない;
セル値 + ソルバ face 幾何なら mdot が両ノズルで一致 = チョーク整合)。

  mdot   = ∫ rho u_x dA
  F_mom  = ∫ rho u_x^2 dA              (軸方向運動量フラックス; 発散損失込み)
  F_pres = ∫ (P - Pa) dA
  F      = F_mom + F_pres
  lambda = F_mom / ∫ rho |U| u_x dA    (発散効率; 完全軸流で 1)
"""
from __future__ import annotations

import sys
import h5py
import numpy as np

GAMMA = 1.4
R_GAS = 1004.5 * (GAMMA - 1.0) / GAMMA  # 287.06
PA = 20000.0  # 周囲圧 [Pa] (= outlet 背圧)
OUTLET_PHYS = 2


def thrust(run, step):
    nz = h5py.File(f"{run}/nozzle.h5", "r")
    ip = nz[f"/BCONDS/{OUTLET_PHYS}/iPlanes"][:]
    ic = nz[f"/BCONDS/{OUTLET_PHYS}/iCells"][:]
    rf = nz["/PLANES/centCoords"][:].reshape(-1, 3)[ip, 1]
    dr = nz["/PLANES/surfArea"][:][ip]          # 2D エッジ長
    dA = 2.0 * np.pi * rf * dr                  # 回転環状面積

    f = h5py.File(f"{run}/res_{step}.h5", "r")
    ro = f["/VALUE/ro"][:][ic]
    u = f["/VALUE/roUx"][:][ic] / ro
    v = f["/VALUE/roUy"][:][ic] / ro
    roe = f["/VALUE/roe"][:][ic]
    P = (GAMMA - 1.0) * (roe - 0.5 * ro * (u * u + v * v))
    T = P / (ro * R_GAS)
    M = np.hypot(u, v) / np.sqrt(GAMMA * R_GAS * T)

    mdot = np.sum(ro * u * dA)
    F_mom = np.sum(ro * u * u * dA)
    F_pres = np.sum((P - PA) * dA)
    lam = F_mom / np.sum(ro * np.hypot(u, v) * u * dA)
    w = ro * u * dA
    return dict(mdot=mdot, F_mom=F_mom, F_pres=F_pres, F=F_mom + F_pres,
                Ae=dA.sum(), lam=lam,
                M_avg=np.sum(M * w) / np.sum(w),
                ang_avg=np.sum(np.degrees(np.arctan2(v, u)) * w) / np.sum(w),
                Pe_avg=np.sum(P * w) / np.sum(w))


def main():
    # 既定: clustered メッシュの inviscid/viscous 比較
    sets = sys.argv[1:] or [
        "INVISCID:run_0006_bell_inv_cl:run_0007_conical_inv_cl:8000",
        "VISCOUS(lam):run_0004_bell_visc:run_0005_conical_visc:12000",
    ]
    print(f"exit-plane thrust  (Pa={PA:.0f} Pa, ambient).  Ae=pi*Re^2=3367 mm^2")
    for spec in sets:
        tag, bell, cone, step = spec.split(":")
        b, c = thrust(bell, step), thrust(cone, step)
        print(f"\n=== {tag}  (step {step}) ===")
        print(f"{'metric':<20}{'Rao bell':>12}{'conical':>12}")
        rows = [("mdot [kg/s]", "mdot", "{:.4f}", 1),
                ("F total [N]", "F", "{:.2f}", 1),
                ("  F_mom [N]", "F_mom", "{:.2f}", 1),
                ("  F_pres [N]", "F_pres", "{:.2f}", 1),
                ("exit area [mm^2]", "Ae", "{:.1f}", 1e6),
                ("div.eff. lambda", "lam", "{:.4f}", 1),
                ("exit M (mavg)", "M_avg", "{:.3f}", 1),
                ("exit ang [deg]", "ang_avg", "{:.2f}", 1)]
        for lab, k, fmt, sc in rows:
            print(f"{lab:<20}{fmt.format(b[k]*sc):>12}{fmt.format(c[k]*sc):>12}")
        print(f"  -> thrust gain (bell vs cone): {(b['F']-c['F'])/c['F']*100:+.2f} %"
              f"   mdot match Δ={(b['mdot']-c['mdot'])/c['mdot']*100:+.2f}%")


if __name__ == "__main__":
    main()
