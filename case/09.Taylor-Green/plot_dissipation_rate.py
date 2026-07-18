#!/usr/bin/env python3
"""L3: 粘性 TGV の散逸率カーブ ε*(t*) = −d(K/K0)/dt* を res 時系列から計算・比較。

t* = t·V0/L (V0=M0=0.4, L=1)。Re=1600 の DNS 参照 (Brachet 1983 / DeBonis 2013 AIAA,
512³ spectral) ではピークが t*≈9。32³ LES では鈍る/やや前後するが、
「ピーク位置がおよそ t* 8-10・カーブが滑らか・keepDiss on/off の差が小さい」を合否とする。
usage: plot_dissipation_rate.py RUN_DIR [RUN_DIR ...] [--dt 0.007] [--v0 0.4] [--out PNG]
"""
import argparse, glob, os, re
import numpy as np, h5py

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+")
ap.add_argument("--dt", type=float, default=0.007)
ap.add_argument("--v0", type=float, default=0.4)
ap.add_argument("--out", default="dissipation_rate_L3.png")
a = ap.parse_args()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

for run in a.runs:
    fs = sorted(glob.glob(os.path.join(run, "res_*.h5")),
                key=lambda p: int(re.search(r"res_(\d+)\.h5", p).group(1)))
    fs = [p for p in fs if "nan" not in p]
    steps, K = [], []
    for p in fs:
        with h5py.File(p, "r") as f:
            ro = np.array(f["VALUE/ro"], dtype=np.float64)
            u = [np.array(f[f"VALUE/U{c}"], dtype=np.float64) for c in "xyz"]
        steps.append(int(re.search(r"res_(\d+)\.h5", p).group(1)))
        K.append(float(np.sum(0.5*ro*(u[0]**2+u[1]**2+u[2]**2))))
    steps = np.array(steps); K = np.array(K)
    tstar = steps*a.dt*a.v0
    Kn = K/K[0]
    # 中央差分で ε* = −dKn/dt*
    eps = -np.gradient(Kn, tstar)
    lab = os.path.basename(run.rstrip("/"))
    ax1.plot(tstar, Kn, "o-", ms=3, label=lab)
    ax2.plot(tstar, eps, "o-", ms=3, label=lab)
    ipk = int(np.argmax(eps))
    print(f"{lab}: K/K0(end)={Kn[-1]:.4f}  ε*ピーク={eps[ipk]:.4f} @ t*={tstar[ipk]:.2f}")

ax2.axvline(9.0, color="gray", ls=":", lw=1, label="DNS peak t*≈9")
ax1.set_xlabel(r"$t^*$"); ax1.set_ylabel(r"$K/K_0$"); ax1.grid(alpha=0.3); ax1.legend(fontsize=7)
ax2.set_xlabel(r"$t^*$"); ax2.set_ylabel(r"$\varepsilon^*=-d(K/K_0)/dt^*$")
ax2.grid(alpha=0.3); ax2.legend(fontsize=7)
fig.suptitle("viscous TGV Re=1600 + WALE (L3)")
fig.tight_layout(); fig.savefig(a.out, dpi=130)
print(f"saved {a.out}")
