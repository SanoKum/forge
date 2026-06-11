#!/usr/bin/env python3
"""case 29: residual / Mach 場 / 出口プロファイル / 近壁 vis_turb の図を生成。

正準 run (theta_n=30, theta_e=8, clustered mesh):
  inviscid : run_0006_bell_inv_cl / run_0007_conical_inv_cl   (step 8000)
  laminar  : run_0004_bell_visc   / run_0005_conical_visc     (step 12000)
  RANS SST : run_0012_bell_rans_long / run_0013_conical_rans_long (step 40000)
"""
from __future__ import annotations
import csv
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

GAMMA = 1.4
R_GAS = 1004.5 * (GAMMA - 1.0) / GAMMA

# (bell_run, conical_run, bell_step, conical_step)
PAIRS = {
    "inviscid": ("run_0006_bell_inv_cl", "run_0007_conical_inv_cl", "8000", "8000"),
    "laminar": ("run_0004_bell_visc", "run_0005_conical_visc", "12000", "12000"),
    "RANS SST": ("run_0014_bell_rans_dev", "run_0015_conical_rans_dev", "20000", "30000"),
}


def field(run, step):
    cc = h5py.File(f"{run}/nozzle.h5", "r")["/CELLS/centCoords"][:].reshape(-1, 3)
    f = h5py.File(f"{run}/res_{step}.h5", "r")
    ro = f["/VALUE/ro"][:]; u = f["/VALUE/roUx"][:] / ro; v = f["/VALUE/roUy"][:] / ro
    P = (GAMMA - 1.0) * (f["/VALUE/roe"][:] - 0.5 * ro * (u * u + v * v))
    T = P / (ro * R_GAS)
    M = np.hypot(u, v) / np.sqrt(GAMMA * R_GAS * T)
    return cc[:, 0] * 1000, cc[:, 1] * 1000, M


# ---- residual plots (all pairs) ----
for tag, (br, cr, stb, stc) in PAIRS.items():
    for name, run in [("bell", br), ("conical", cr)]:
        rows = [r for r in csv.DictReader(open(f"{run}/residual_history.csv"))
                if r["phase"] == "outer_end"]
        step = np.array([int(r["step"]) for r in rows])
        fig, ax = plt.subplots(figsize=(7, 4))
        for key, lab in [("rms_roUx", "roUx"), ("rms_roUy", "roUy"), ("rms_roe", "roe")]:
            y = np.array([float(r[key]) if r[key] not in ("nan", "inf", "-inf") else np.nan
                          for r in rows])
            ax.semilogy(step, np.abs(y), label=lab, lw=1)
        ax.set_xlabel("step"); ax.set_ylabel("rms residual"); ax.grid(alpha=0.3)
        ax.set_title(f"case 29 {name} {tag} residual"); ax.legend()
        fig.tight_layout(); fig.savefig(f"{run}/residual_history.png", dpi=120)
        plt.close(fig)
print("wrote residual_history.png for all runs")

# ---- Mach field (inviscid clustered: cleanest core) ----
br, cr, st, _ = PAIRS["inviscid"]
xb, yb, Mb = field(br, st); xc, yc, Mc = field(cr, st)
fig, ax = plt.subplots(figsize=(12, 5))
lv = np.linspace(0, 4.5, 46)
cf = ax.tricontourf(mtri.Triangulation(xb, yb), Mb, levels=lv, cmap="turbo", extend="max")
ax.tricontourf(mtri.Triangulation(xc, -yc), Mc, levels=lv, cmap="turbo", extend="max")
ax.axhline(0, color="k", lw=0.5)
ax.text(-40, 18, "Rao bell", fontsize=12, weight="bold")
ax.text(-40, -20, "conical 15deg", fontsize=12, weight="bold")
ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("r [mm]")
ax.set_title("case 29: Mach field — Rao bell (top) vs 15deg conical (bottom), same throat/length/exit")
fig.colorbar(cf, ax=ax, label="Mach", shrink=0.8)
fig.tight_layout(); fig.savefig("mach_comparison.png", dpi=120); plt.close(fig)
print("wrote mach_comparison.png")

# ---- exit profiles (inviscid clustered) using solver face geom + adjacent cells ----
fig, axx = plt.subplots(1, 2, figsize=(11, 4.5))
for name, run, col in [("Rao bell", br, "tab:blue"), ("conical 15deg", cr, "tab:orange")]:
    nz = h5py.File(f"{run}/nozzle.h5", "r")
    ip = nz["/BCONDS/2/iPlanes"][:]; ic = nz["/BCONDS/2/iCells"][:]
    r = nz["/PLANES/centCoords"][:].reshape(-1, 3)[ip, 1] * 1000
    f = h5py.File(f"{run}/res_{st}.h5", "r")
    ro = f["/VALUE/ro"][:][ic]; u = f["/VALUE/roUx"][:][ic] / ro; v = f["/VALUE/roUy"][:][ic] / ro
    P = (GAMMA - 1) * (f["/VALUE/roe"][:][ic] - 0.5 * ro * (u * u + v * v)); T = P / (ro * R_GAS)
    M = np.hypot(u, v) / np.sqrt(GAMMA * R_GAS * T); ang = np.degrees(np.arctan2(v, u))
    o = np.argsort(r)
    axx[0].plot(M[o], r[o], col, label=name); axx[1].plot(ang[o], r[o], col, label=name)
axx[0].set_xlabel("Mach"); axx[0].set_ylabel("r [mm]"); axx[0].set_title("exit Mach profile")
axx[1].set_xlabel("flow angle [deg]"); axx[1].set_ylabel("r [mm]"); axx[1].set_title("exit flow angle")
for a in axx: a.grid(alpha=0.3); a.legend()
fig.tight_layout(); fig.savefig("exit_profiles.png", dpi=120); plt.close(fig)
print("wrote exit_profiles.png")

# ---- near-wall vis_turb/vis_lam vs wall_dist (RANS, several x stations) ----
br, cr, stb, stc = PAIRS["RANS SST"]
cc = h5py.File(f"{br}/nozzle.h5", "r")["/CELLS/centCoords"][:].reshape(-1, 3) * 1000
f = h5py.File(f"{br}/res_{stb}.h5", "r")
vt = f["/VALUE/vis_turb"][:]; vl = f["/VALUE/vis_lam"][:]; wd = f["/VALUE/wall_dist"][:] * 1000
fig, ax = plt.subplots(figsize=(7, 5))
for xs in [10, 25, 40, 60, 80]:
    m = np.abs(cc[:, 0] - xs) < 1.5
    o = np.argsort(wd[m])
    ax.plot((vt[m] / np.maximum(vl[m], 1e-30))[o], wd[m][o], label=f"x={xs}mm")
ax.set_xlabel("vis_turb / vis_lam"); ax.set_ylabel("wall distance [mm]")
ax.set_ylim(0, 3); ax.set_title("case 29 RANS bell: near-wall turbulent viscosity (converged)")
ax.grid(alpha=0.3); ax.legend()
fig.tight_layout(); fig.savefig("rans_visturb_profile.png", dpi=120); plt.close(fig)
print("wrote rans_visturb_profile.png")
