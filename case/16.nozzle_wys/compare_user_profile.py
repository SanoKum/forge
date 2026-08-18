#!/usr/bin/env python3
"""ユーザ指定形状 (2026-08-19) の run 群の軸中心線 Mach / 静圧 / 液滴質量分率 g を、凝縮 off/on × Euler/NS で比較。
usage: python3 compare_user_profile.py --out PNG --csv CSV  RUN_LABEL=RUN_DIR ...
   例: python3 compare_user_profile.py --out compare_user_profile.png "Euler dry=run_0189_user_cell_euler_dry" ...
Mach = |u|/sonic、p/p0 は Pt=59070 Pa。実験 (Fig.3 isentrope / 1 kPa) と 1D 等エントロピー (γ=1.4, ユーザ形状) も重ねる。"""
import sys, os, glob, argparse
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "mesh"))
from nozzle_user_profile import y_user_mm
P0 = 59070.0

def latest(run):
    rs = sorted(glob.glob(os.path.join(run, "res_[0-9]*.h5")), key=lambda f: int(os.path.basename(f)[4:-3]))
    return rs[-1]

def centerline(run, band=0.04, dx=0.5):
    """|y| < band*y_wall(x) のセル/ノードを dx[mm] ビンで平均。戻り x[mm], M, P, T, g, name"""
    r = latest(run)
    m = [p for p in glob.glob(os.path.join(run, "nozzle_*.h5")) if "res_" not in os.path.basename(p)][0]
    cc = h5py.File(m)["/CELLS/centCoords"][:].reshape(-1, 3)
    v = h5py.File(r)["/VALUE"]
    x = cc[:, 0] * 1e3; y = cc[:, 1] * 1e3
    yw = np.array([y_user_mm(t) for t in x])
    sel = np.abs(y) < band * yw
    P = np.array(v["P"])[sel]; T = np.array(v["T"])[sel]; s = np.array(v["sonic"])[sel]
    U = np.sqrt(np.array(v["Ux"])[sel]**2 + np.array(v["Uy"])[sel]**2)
    g = np.array(v["g_0"])[sel] if "g_0" in v else np.zeros_like(P)
    xs = x[sel]
    xb = np.arange(-60, 95 + dx, dx)
    idx = np.digitize(xs, xb)
    out = []
    for i in range(1, len(xb)):
        mm = idx == i
        if mm.any():
            out.append([xs[mm].mean(), (U[mm]/s[mm]).mean(), P[mm].mean(), T[mm].mean(), g[mm].mean()])
    a = np.array(out)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], os.path.basename(r)

def load_exp():
    r = []
    for ln in open(os.path.join(HERE, "wyslouzil_fig3_pp0.csv")):
        p = ln.strip().split(",")
        if p and p[0] and p[0][0].isdigit(): r.append([float(p[0]), float(p[1]), float(p[2])])
    a = np.array(r); return a[:, 0]*10, a[:, 1], a[:, 2]

def isen_1d(gam=1.4):
    from scipy.optimize import brentq
    xs = np.linspace(-59.9, 95, 400); M = []
    for xx in xs:
        ar = y_user_mm(xx)/2.5
        f = lambda m: (1/m)*((2/(gam+1))*(1+(gam-1)/2*m*m))**((gam+1)/(2*(gam-1))) - ar
        if ar <= 1+1e-9: M.append(1.0)
        elif xx > 0: M.append(brentq(f, 1.0, 10))
        else: M.append(brentq(f, 1e-6, 1.0))
    M = np.array(M); return xs, M, (1+(gam-1)/2*M**2)**(-gam/(gam-1))

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+")
ap.add_argument("--out", default="compare_user_profile.png")
ap.add_argument("--csv", default="compare_user_profile_centerline.csv")
ap.add_argument("--title", default="Wyslouzil nozzle, user profile 2026-08-19 (TP split MIXDRY+H2O, Kw+HK)")
a = ap.parse_args()
styles = {"Euler dry": ("C0", "--"), "Euler cond": ("C0", "-"), "NS dry": ("C3", "--"), "NS cond": ("C3", "-")}
fig, ax = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
xi, Mi, pi = isen_1d()
ax[0].plot(xi, Mi, color="gray", lw=1, ls=":", label="1D isentropic (γ=1.4)")
ax[1].plot(xi, pi*P0/1e3, color="gray", lw=1, ls=":", label="1D isentropic (γ=1.4)")
ex, eiso, econd = load_exp()
ax[1].plot(ex, eiso*P0/1e3, "s", color="gray", ms=6, mfc="none", label="exp isentrope (Fig.3)")
ax[1].plot(ex, econd*P0/1e3, "o", color="k", ms=6, mfc="none", label="exp cond 1 kPa (Fig.3)")
rows = {}
for spec in a.runs:
    lab, run = spec.split("=", 1)
    x, M, P, T, g, name = centerline(run)
    c, ls = styles.get(lab, ("k", "-"))
    ax[0].plot(x, M, color=c, ls=ls, lw=1.6, label=f"{lab} ({run}, {name})")
    ax[1].plot(x, P/1e3, color=c, ls=ls, lw=1.6, label=lab)
    ax[2].plot(x, g, color=c, ls=ls, lw=1.6, label=lab)
    rows[lab] = (x, M, P, T, g)
    i = np.argmin(np.abs(x - 94.0)); j = np.argmax(g > 1e-3) if (g > 1e-3).any() else None
    print(f"{lab:11s} {run}: exit(x≈94mm) M={M[i]:.3f} P={P[i]:.0f} Pa (p/p0 {P[i]/P0:.4f}) T={T[i]:.1f} K g={g[i]:.4f} g_max={g.max():.4f}"
          + (f" onset(g>1e-3) x={x[j]:.1f} mm" if j is not None else ""))
ax[0].set_ylabel("Mach (centerline)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3); ax[0].set_title(a.title)
ax[1].set_ylabel("static pressure [kPa]"); ax[1].legend(fontsize=7); ax[1].grid(alpha=.3); ax[1].set_yscale("log")
ax[2].set_ylabel("liquid mass fraction g"); ax[2].set_xlabel("x from throat [mm]"); ax[2].legend(fontsize=7); ax[2].grid(alpha=.3)
for A in ax: A.axvline(0, color="gray", lw=.5)
plt.tight_layout(); plt.savefig(a.out, dpi=130); print("saved", a.out)
# CSV (共通 x グリッドへ補間)
xg = np.arange(-59.5, 95.0, 0.5)
hdr = ["x_mm"]; cols = [xg]
for lab, (x, M, P, T, g) in rows.items():
    for nm, arr in (("M", M), ("P_Pa", P), ("T_K", T), ("g", g)):
        hdr.append(f"{lab.replace(' ', '_')}_{nm}"); cols.append(np.interp(xg, x, arr, left=np.nan, right=np.nan))
np.savetxt(a.csv, np.column_stack(cols), delimiter=",", header=",".join(hdr), comments="", fmt="%.6g")
print("saved", a.csv)
