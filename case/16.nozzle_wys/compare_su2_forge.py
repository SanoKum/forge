#!/usr/bin/env python3
"""forge ↔ SU2 中心線 静圧比 (p/p0) 比較 — Wyslouzil 2D ノズル.

3 物理モデル (inviscid Euler / laminar viscous / turbulent SST) について、
forge (SLAU, セル中心) と SU2 (ROE, 節点) の中心線 (y≈0) 静圧比を重ねる。
p0 = Pt_in = 101325 Pa。面積比等エントロピー線も参考に重ねる。
"""
import sys, glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

GAM = 1.4
PT  = 101325.0
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------- forge 中心線抽出 ----------
def forge_centerline(run_dir):
    fs = [f for f in glob.glob(os.path.join(run_dir, "res_*.h5")) if "nan" not in f]
    res = max(fs, key=lambda f: int(f.split("res_")[1].split(".")[0]))
    mh = [p for p in glob.glob(os.path.join(run_dir, "nozzle_*.h5")) if "res_" not in os.path.basename(p)][0]
    with h5py.File(mh, "r") as m:
        cc = m["/CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res, "r") as r:
        P = r["/VALUE/P"][:]
    x, y, z = cc[:, 0], cc[:, 1], cc[:, 2]
    epsy = (y.max() - y.min()) * 0.02
    zc = 0.5 * (z.min() + z.max()); dz = z.max() - z.min()
    sel = (np.abs(y) < epsy) & (np.abs(z - zc) < max(dz * 0.6, 1e-6))
    xs = x[sel]; ps = P[sel]
    o = np.argsort(xs)
    return xs[o], ps[o], os.path.basename(res)

# ---------- SU2 中心線抽出 (restart_flow.csv, 節点) ----------
def su2_centerline(run_dir):
    f = os.path.join(run_dir, "restart_flow.csv")
    with open(f) as fh:
        hdr = [c.strip().strip('"').strip() for c in fh.readline().split(",")]
    data = np.loadtxt(f, delimiter=",", skiprows=1)
    def ci(*names):
        for n in names:
            if n in hdr: return hdr.index(n)
        raise KeyError(names)
    ix, iy = ci("x"), ci("y")
    x, y = data[:, ix], data[:, iy]
    if "Pressure" in hdr:
        P = data[:, ci("Pressure")]
    else:
        # restart は保存量のみ → 理想気体で p=(γ-1)(ρE-0.5|ρu|^2/ρ) を計算
        rho = data[:, ci("Density")]
        mx, my = data[:, ci("Momentum_x")], data[:, ci("Momentum_y")]
        rhoE = data[:, ci("Energy")]
        P = (GAM - 1.0) * (rhoE - 0.5 * (mx * mx + my * my) / rho)
    epsy = (y.max() - y.min()) * 0.02
    sel = np.abs(y) < epsy
    xs = x[sel]; ps = P[sel]
    o = np.argsort(xs)
    return xs[o], ps[o]

# ---------- 面積比等エントロピー (forge メッシュ高さから) ----------
def isentropic_curve(run_dir):
    mh = [p for p in glob.glob(os.path.join(run_dir, "nozzle_*.h5")) if "res_" not in os.path.basename(p)][0]
    with h5py.File(mh, "r") as m:
        cc = m["/CELLS/centCoords"][:].reshape(-1, 3)
    x, y = cc[:, 0], cc[:, 1]
    nb = 240
    xe = np.linspace(x.min(), x.max(), nb + 1); xc = 0.5 * (xe[:-1] + xe[1:])
    h = np.full(nb, np.nan)
    for i in range(nb):
        s = (x >= xe[i]) & (x < xe[i + 1])
        if s.sum() > 3: h[i] = y[s].max() - y[s].min()
    ok = np.isfinite(h); xc, h = xc[ok], h[ok]
    hstar = np.nanmin(h); ar = h / hstar; xt = xc[np.argmin(h)]
    def ar2M(a, sup):
        if a <= 1.0000001: return 1.0
        f = lambda M: (1.0/M)*((2.0/(GAM+1))*(1+0.5*(GAM-1)*M*M))**((GAM+1)/(2*(GAM-1))) - a
        try: return brentq(f, 1.0+1e-6, 50.0) if sup else brentq(f, 1e-4, 1.0)
        except Exception: return np.nan
    M = np.array([ar2M(a, xx > xt) for a, xx in zip(ar, xc)])
    p = PT * (1 + 0.5*(GAM-1)*M**2) ** (-GAM/(GAM-1))
    return xc, p / PT

# ---------- main ----------
CASES = [
    ("Euler (inviscid)",  "run_0001_slau_2d_imp",  "run_0020_su2_euler", "tab:blue"),
    ("Laminar (viscous)", "run_0003_slau_2d_visc", "run_0023_su2_lam_conv", "tab:green"),
    ("Turbulent (SST)",   "run_0004_slau_2d_sst",  "run_0024_su2_sst_conv", "tab:red"),
]

def main():
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), sharey=True)
    xc_iso, pp_iso = isentropic_curve(os.path.join(HERE, "run_0001_slau_2d_imp"))
    print(f"{'model':22s} {'x[cm]':>6s} {'forge p/p0':>11s} {'SU2 p/p0':>10s} {'diff%':>7s}")
    for ax, (lbl, fdir, sdir, col) in zip(axes, CASES):
        fx, fp, fres = forge_centerline(os.path.join(HERE, fdir))
        sx, sp = su2_centerline(os.path.join(HERE, sdir))
        ax.plot(xc_iso*1e2, pp_iso, '-', color='k', lw=1.0, alpha=0.5, label='isentropic (area)')
        ax.plot(fx*1e2, fp/PT, '.', ms=3, color=col, label=f'forge')
        ax.plot(sx*1e2, sp/PT, '-', color='magenta', lw=1.2, label='SU2')
        ax.set_title(lbl); ax.set_xlabel('x [cm]'); ax.grid(alpha=0.3); ax.legend(fontsize=9)
        ax.set_ylim(0, 1.05)
        # 比較サンプル点で数値出力
        for xq_cm in [-2.0, 0.0, 2.0, 4.0, 6.0, 9.0]:
            xq = xq_cm * 1e-2
            def interp(xx, pp):
                xx = np.asarray(xx); pp = np.asarray(pp)
                if xq < xx.min() or xq > xx.max(): return np.nan
                # 近傍平均で滑らかに
                m = np.abs(xx - xq) < 2e-3
                return pp[m].mean()/PT if m.any() else np.interp(xq, xx, pp)/PT
            fv = interp(fx, fp); sv = interp(sx, sp)
            d = (fv - sv)/sv*100 if (sv==sv and sv!=0) else float('nan')
            print(f"{lbl:22s} {xq_cm:6.1f} {fv:11.4f} {sv:10.4f} {d:7.2f}")
    axes[0].set_ylabel('static pressure ratio  p / p0')
    fig.suptitle('Wyslouzil 2D nozzle — centerline static pressure ratio: forge vs SU2', y=1.02)
    fig.tight_layout()
    out = os.path.join(HERE, "compare_su2_forge.png")
    fig.savefig(out, dpi=130, bbox_inches='tight')
    print("wrote", out)

if __name__ == "__main__":
    main()
