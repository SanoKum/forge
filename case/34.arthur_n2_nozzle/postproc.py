#!/usr/bin/env python3
"""Arthur (1952) source-flow ノズル dry 膨張の後処理。

run_dir の最終 res_*.h5 から
  - 中心線 (|y|<eps) の Mach / 静圧比 P/P0 を x に対して抽出
  - 壁近傍の静圧比 P/P0 を x に対して抽出 (Arthur が測ったのは壁面静圧)
し、面積比 A/A* = y_wall(x)/y_t (source-flow 解析値) に基づく 1 次元
等エントロピー理論と比較する。Arthur の Fig.4 (M vs A/A*) / Fig.11
(centerline vs wall) に対応。

usage: python3 postproc.py <run_dir>
"""
import sys, glob, os
import numpy as np
import h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

GAM = 1.4
P0  = 844037.0      # 貯気全圧 [Pa]
T0  = 290.0         # 貯気全温 [K]
YT  = 0.000127      # スロート半高 [m]
TANA = 0.0962992    # tan(5.5 deg)

def latest_res(run_dir):
    fs = [f for f in glob.glob(os.path.join(run_dir, "res_*.h5")) if "nan" not in f]
    return max(fs, key=lambda f: int(f.split("res_")[1].split(".")[0]))

def M_from_AR(ar):
    """A/A* -> 超音速枝 Mach (等エントロピー)。"""
    if ar <= 1.0000001: return 1.0
    f = lambda M: (1.0/M)*((2.0/(GAM+1))*(1+0.5*(GAM-1)*M*M))**((GAM+1)/(2*(GAM-1))) - ar
    return brentq(f, 1.0+1e-9, 60.0)

def main(run_dir):
    res = latest_res(run_dir)
    mesh = os.path.join(run_dir, "arthur_nozzle.h5")
    with h5py.File(mesh, "r") as m:
        cc = m["/CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res, "r") as r:
        P  = r["/VALUE/P"][:];   T = r["/VALUE/T"][:]
        Ux = r["/VALUE/Ux"][:];  Uy = r["/VALUE/Uy"][:];  Uz = r["/VALUE/Uz"][:]
        son= r["/VALUE/sonic"][:]
    n = cc.shape[0]
    P, T, Ux, Uy, Uz, son = (a[:n] for a in (P, T, Ux, Uy, Uz, son))
    x, y = cc[:,0], cc[:,1]
    mach = np.sqrt(Ux**2+Uy**2+Uz**2)/son
    print(f"res file: {os.path.basename(res)}  cells={n}")

    # --- 各 x ビンで最下段(中心線)・最上段(壁)セルを抽出 ---
    nb = 200
    xed = np.linspace(0, x.max(), nb+1); xcb = 0.5*(xed[:-1]+xed[1:])
    Pcl = np.full(nb, np.nan); Mcl = np.full(nb, np.nan)
    Pw  = np.full(nb, np.nan); Mw  = np.full(nb, np.nan)
    for i in range(nb):
        sel = (x>=xed[i])&(x<xed[i+1])
        if sel.sum()>0:
            idx = np.where(sel)[0]
            jcl = idx[np.argmin(y[sel])]   # 中心線側
            jw  = idx[np.argmax(y[sel])]   # 壁側
            Pcl[i]=P[jcl]/P0; Mcl[i]=mach[jcl]
            Pw[i] =P[jw]/P0;  Mw[i] =mach[jw]
    okc = np.isfinite(Pcl); okw = np.isfinite(Pw)
    xcl = xcb

    # --- 解析 (source-flow) 面積比と等エントロピー理論 ---
    xth = np.linspace(1e-5, x.max(), 400)
    arw = np.sqrt(YT**2 + (xth*TANA)**2)/YT
    Mthy = np.array([M_from_AR(a) for a in arw])
    Pthy = (1+0.5*(GAM-1)*Mthy**2)**(-GAM/(GAM-1))

    # 出口値
    ar_exit = np.sqrt(YT**2+(x.max()*TANA)**2)/YT
    M_exit_iso = M_from_AR(ar_exit)
    # 中心線出口 (最下流 5% 平均)
    tail = okc & (xcl > 0.95*x.max())
    print(f"area ratio at exit A/A* = {ar_exit:.1f}")
    print(f"isentropic exit Mach    = {M_exit_iso:.3f}")
    print(f"forge centerline exit M = {np.nanmean(Mcl[tail]):.3f}  (P/P0={np.nanmean(Pcl[tail]):.2e})")
    print(f"isentropic exit P/P0    = {(1+0.5*(GAM-1)*M_exit_iso**2)**(-GAM/(GAM-1)):.2e}")

    # --- 図 ---
    fig, ax = plt.subplots(1, 2, figsize=(12,4.5))
    ax[0].plot(arw, Mthy, 'k-', lw=2, label='isentropic (1D)')
    ax[0].plot(np.sqrt(YT**2+(xcl[okc]*TANA)**2)/YT, Mcl[okc], 'C0.', ms=3, label='forge centerline')
    ax[0].plot(np.sqrt(YT**2+(xcb[okw]*TANA)**2)/YT, Mw[okw], 'C3x', ms=4, label='forge wall')
    ax[0].set_xlabel('area ratio  A/A*'); ax[0].set_ylabel('Mach'); ax[0].grid(alpha=.3)
    ax[0].legend(); ax[0].set_title('Mach vs A/A*  (cf. Arthur Fig.4)')

    ax[1].semilogy(xth*1e3, Pthy, 'k-', lw=2, label='isentropic (1D)')
    ax[1].semilogy(xcl[okc]*1e3, Pcl[okc], 'C0.', ms=3, label='forge centerline')
    ax[1].semilogy(xcb[okw]*1e3, Pw[okw], 'C3x', ms=4, label='forge wall')
    ax[1].set_xlabel('x from throat [mm]'); ax[1].set_ylabel('P / P0'); ax[1].grid(alpha=.3, which='both')
    ax[1].legend(); ax[1].set_title('static pressure ratio  (cf. Arthur Fig.11)')
    fig.tight_layout()
    out = os.path.join(run_dir, "postproc_arthur.png")
    fig.savefig(out, dpi=130); print("wrote", out)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv)>1 else ".")
