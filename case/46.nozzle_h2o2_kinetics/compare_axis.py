#!/usr/bin/env python3
"""forge 軸線 (y=0 ノード) の T, M, Y_OH, Y_NO を Q1D 参照 (finite-rate / frozen) と重ねる。
  .venv-chem/bin/python compare_axis.py run_dir[,run_dir2...] res_step out.png"""
import sys, h5py, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
runs = sys.argv[1].split(","); step = sys.argv[2]; out = sys.argv[3]
names = ["H2","O2","H","O","OH","H2O","HO2","H2O2","N2","N","NO","NO2","HNO","AR","CO2"]
import os
qf = os.environ.get("Q1D_FR", "q1d_finite_rate.csv"); qz = os.environ.get("Q1D_FZ", "q1d_frozen.csv")
q1 = np.genfromtxt(qf, delimiter=",", names=True); q0 = np.genfromtxt(qz, delimiter=",", names=True)
fig, ax = plt.subplots(2, 2, figsize=(12, 8)); ax = ax.ravel()
for lab, q, ls in (("Q1D finite-rate", q1, "k-"), ("Q1D frozen", q0, "k--")):
    ax[0].plot(q["x"], q["T"], ls, label=lab); ax[1].plot(q["x"], q["M"], ls, label=lab)
    ax[2].plot(q["x"], q["Y_OH"], ls, label=lab); ax[3].plot(q["x"], q["Y_NO"], ls, label=lab)
res = {}
for r in runs:
    with h5py.File(f"{r}/nozzle.h5") as m: xyz = m["MESH/COORD"][:].reshape(-1, 3)
    ax_nodes = np.where(np.abs(xyz[:, 1]) < 1e-9)[0]; o = np.argsort(xyz[ax_nodes, 0]); ax_nodes = ax_nodes[o]; x = xyz[ax_nodes, 0]
    with h5py.File(f"{r}/res_{step}.h5") as f:
        v = f["VALUE"]; T = v["T"][:][ax_nodes]; P = v["P"][:][ax_nodes]; U = v["Ux"][:][ax_nodes]
        Y = {n: v[f"Y{i}"][:][ax_nodes] for i, n in enumerate(names)}
        M = (U/np.maximum(v["sonic"][:][ax_nodes], 1e-30)) if "sonic" in v else (U/np.sqrt(np.maximum(v["gamma"][:][ax_nodes]*v["Rmix"][:][ax_nodes]*T, 1e-30)) if "Rmix" in v else None)
    ax[0].plot(x, T, ".", ms=3, label=f"forge {r}"); ax[2].plot(x, Y["OH"], ".", ms=3, label=f"forge {r}"); ax[3].plot(x, Y["NO"], ".", ms=3, label=f"forge {r}")
    if M is not None: ax[1].plot(x, M, ".", ms=3, label=f"forge {r}")
    i_e = -1; res[r] = (x[i_e], T[i_e], (M[i_e] if M is not None else np.nan), Y["OH"][i_e], Y["NO"][i_e])
for i, (yl, lg) in enumerate((("T [K]", False), ("M", False), ("Y_OH", True), ("Y_NO", True))):
    ax[i].set_xlabel("x [m]"); ax[i].set_ylabel(yl); ax[i].grid(alpha=.3); ax[i].legend(fontsize=8)
    if lg: ax[i].set_yscale("log")
fig.tight_layout(); fig.savefig(out, dpi=120)
print("exit (Q1D fr):", f"T={q1['T'][-1]:.1f} M={q1['M'][-1]:.3f} Y_OH={q1['Y_OH'][-1]:.3e}", "| (Q1D frozen):", f"T={q0['T'][-1]:.1f} M={q0['M'][-1]:.3f} Y_OH={q0['Y_OH'][-1]:.3e}")
for r, v in res.items(): print(f"exit forge {r}: x={v[0]:.3f} T={v[1]:.1f} M={v[2]:.3f} Y_OH={v[3]:.3e} Y_NO={v[4]:.3e}")
