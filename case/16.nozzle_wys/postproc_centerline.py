#!/usr/bin/env python3
"""Wyslouzil ノズル centerline 後処理。

run_dir の最終 res_*.h5 から中心線 (|y|<eps) の Mach / 静圧比を取り出し、
メッシュ局所高さから求めた面積比 A/A* に基づく等エントロピー理論と比較する。
2D/3D 両対応 (3D は |y|,|z-zc|<eps の中心軸近傍を抽出)。
"""
import sys, glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

GAM = 1.4

def latest_res(run_dir):
    fs = glob.glob(os.path.join(run_dir, "res_*.h5"))
    fs = [f for f in fs if "nan" not in f]
    return max(fs, key=lambda f: int(f.split("res_")[1].split(".")[0]))

def mesh_h5(run_dir):
    cand = [p for p in glob.glob(os.path.join(run_dir, "nozzle_*.h5"))
            if "res_" not in os.path.basename(p)]
    if not cand:
        raise FileNotFoundError("mesh h5 not found")
    return cand[0]

def area_ratio_M(ar, super=True):
    """A/A* -> Mach (等エントロピー)。super=True で超音速枝。"""
    def f(M):
        return (1.0/M)*((2.0/(GAM+1))*(1+0.5*(GAM-1)*M*M))**((GAM+1)/(2*(GAM-1))) - ar
    if ar <= 1.0000001:
        return 1.0
    lo, hi = (1.0, 50.0) if super else (1e-4, 1.0)
    try:
        return brentq(f, lo+1e-6 if super else lo, hi)
    except Exception:
        return np.nan

def main(run_dir, label):
    res = latest_res(run_dir)
    mh = mesh_h5(run_dir)
    with h5py.File(mh, "r") as m:
        cc = m["/CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res, "r") as r:
        P = r["/VALUE/P"][:]; T = r["/VALUE/T"][:]
        Ux = r["/VALUE/Ux"][:]; Uy = r["/VALUE/Uy"][:]; Uz = r["/VALUE/Uz"][:]
        son = r["/VALUE/sonic"][:]
    x, y, z = cc[:,0], cc[:,1], cc[:,2]
    speed = np.sqrt(Ux**2 + Uy**2 + Uz**2)
    mach = speed / son
    Pt_in = 101325.0

    # --- 局所チャンネル高さ -> 面積比 (x ビンごとに max-min y) ---
    nb = 240
    xe = np.linspace(x.min(), x.max(), nb+1)
    xc_b = 0.5*(xe[:-1]+xe[1:])
    h_b = np.full(nb, np.nan)
    for i in range(nb):
        sel = (x >= xe[i]) & (x < xe[i+1])
        if sel.sum() > 3:
            h_b[i] = y[sel].max() - y[sel].min()
    ok = np.isfinite(h_b)
    xc_b, h_b = xc_b[ok], h_b[ok]
    hstar = np.nanmin(h_b)
    ar_b = h_b / hstar
    x_throat = xc_b[np.argmin(h_b)]
    M_iso = np.array([area_ratio_M(a, super=(xx > x_throat)) for a, xx in zip(ar_b, xc_b)])
    P_iso = Pt_in * (1 + 0.5*(GAM-1)*M_iso**2)**(-GAM/(GAM-1))

    # --- 中心線抽出 ---
    zc = 0.5*(z.min()+z.max())
    dz = z.max()-z.min()
    epsz = max(dz*0.15, 1e-6)
    epsy = (y.max()-y.min())*0.02
    sel = (np.abs(y) < epsy) & (np.abs(z - zc) < epsz)
    xs = x[sel]; order = np.argsort(xs)
    xs = xs[order]; Ms = mach[sel][order]; Ps = P[sel][order]; Ts = T[sel][order]

    # --- plot ---
    fig, ax = plt.subplots(1, 2, figsize=(13,4.6))
    ax[0].plot(xs*1e3, Ms, '.', ms=3, color='tab:blue', label=f'forge {label} (centerline)')
    ax[0].plot(xc_b*1e3, M_iso, '-', color='k', lw=1.5, label='isentropic (area ratio)')
    ax[0].axvline(x_throat*1e3, ls=':', color='gray', label='throat')
    ax[0].set_xlabel('x [mm]'); ax[0].set_ylabel('Mach'); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)
    ax[0].set_title(f'Centerline Mach — {label}')

    ax[1].plot(xs*1e3, Ps/1e3, '.', ms=3, color='tab:red', label=f'forge {label}')
    ax[1].plot(xc_b*1e3, P_iso/1e3, '-', color='k', lw=1.5, label='isentropic')
    ax[1].axvline(x_throat*1e3, ls=':', color='gray')
    ax[1].set_xlabel('x [mm]'); ax[1].set_ylabel('static P [kPa]'); ax[1].grid(alpha=0.3); ax[1].legend(fontsize=8)
    ax[1].set_title(f'Centerline static pressure — {label}')
    fig.tight_layout()
    out = os.path.join(run_dir, "centerline_compare.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)
    print(f"[{label}] res={os.path.basename(res)} nCellsCL={sel.sum()}")
    print(f"  throat x={x_throat*1e3:.2f}mm  exit area ratio={ar_b[-1]:.3f}")
    print(f"  exit Mach  forge={Ms[-5:].mean():.3f}   isentropic={M_iso[-1]:.3f}")
    print(f"  exit P[kPa] forge={Ps[-5:].mean()/1e3:.2f}  isentropic={P_iso[-1]/1e3:.2f}")
    print(f"  min static P = {Ps.min()/1e3:.2f} kPa   max Mach = {Ms.max():.3f}   min T = {Ts.min():.1f} K")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
