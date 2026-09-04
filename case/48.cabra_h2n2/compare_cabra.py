#!/usr/bin/env python3
"""Cabra H2/N2: forge の中心軸 (T, Y_H2, Y_O2, Y_H2O vs z/d) と半径分布 T (z/d=1,9,11,14,26) を実験 (図読み取り) と重ねる。
  .venv-chem/bin/python compare_cabra.py run_dir step out.png"""
import sys, h5py, numpy as np, pathlib, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
run, step, out = sys.argv[1], sys.argv[2], sys.argv[3]; HERE = pathlib.Path(__file__).parent; d = 4.57e-3
with h5py.File(f"{run}/cabra.h5") as m: xyz = m["MESH/COORD"][:].reshape(-1, 3)
with h5py.File(f"{run}/res_{step}.h5") as f:
    v = f["VALUE"]; N = v["T"].shape[0]; T = v["T"][:]; Y = {k: v[f"Y{i}"][:] for i, k in enumerate(["H2","O2","H","O","OH","H2O","HO2","H2O2","N2"])}
x = xyz[:N, 0]; r = xyz[:N, 1]
ax_ = np.where(np.abs(r) < 1e-9)[0]; ax_ = ax_[np.argsort(x[ax_])]; zd = x[ax_]/d
ex = np.loadtxt(HERE / "exp_centerline.csv", delimiter=",", comments="#"); exr = np.loadtxt(HERE / "exp_radial_T.csv", delimiter=",", comments="#")
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
ax[0].plot(ex[:, 0], ex[:, 1], "ko", mfc="none", label="exp T"); ax[0].plot(zd, T[ax_], "-", label="forge T"); ax[0].set_xlim(0, 40); ax[0].set_xlabel("z/d"); ax[0].set_ylabel("T [K]"); ax[0].legend(); ax[0].grid(alpha=.3)
for j, (sp, c) in enumerate((("H2", "C0"), ("O2", "C1"), ("H2O", "C2"))):
    ax[1].plot(ex[:, 0], ex[:, 2+j], "o", mfc="none", color=c, label=f"exp {sp}"); ax[1].plot(zd, Y[sp][ax_], "-", color=c, label=f"forge {sp}")
ax[1].set_xlim(0, 40); ax[1].set_xlabel("z/d"); ax[1].set_ylabel("Y (centerline)"); ax[1].legend(fontsize=7, ncol=2); ax[1].grid(alpha=.3)
xs_all = np.unique(np.round(x, 7))
for j, zdv in enumerate((1, 9, 11, 14, 26)):
    xc = xs_all[np.argmin(np.abs(xs_all - zdv*d))]; sel = np.abs(x - xc) < 1e-7; o = np.argsort(r[sel])
    ax[2].plot(exr[:, 0], exr[:, 1+j], "o", mfc="none", color=f"C{j}", ms=4); ax[2].plot(r[sel][o]*1000, T[sel][o], "-", color=f"C{j}", label=f"z/d={zdv}")
ax[2].set_xlim(0, 30); ax[2].set_xlabel("r [mm]"); ax[2].set_ylabel("T [K]"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
fig.suptitle(f"{run} step {step}"); fig.tight_layout(); fig.savefig(out, dpi=120)
Tl = T[ax_]; i = np.argmax(Tl > 600) if (Tl > 600).any() else -1
print(f"wrote {out}; centerline T max={Tl.max():.0f} K, T>600K first at z/d={zd[i]:.1f} (exp: z/d≈14; lift-off H/d≈10 by OH)" if i >= 0 else f"wrote {out}; no ignition on axis (Tmax {Tl.max():.0f})")
oh = Y["OH"]; m_ = oh > 1e-4; print(f"OH>1e-4 first at z/d={(x[m_].min()/d if m_.any() else float('nan')):.1f} (exp lift-off H/d≈10)")
