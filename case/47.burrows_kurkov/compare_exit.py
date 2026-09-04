#!/usr/bin/env python3
"""x=35.6 cm 出口の forge プロファイル (モル分率, 全温) を Burrows–Kurkov の実験 (図読み取り) と重ねる。
  .venv-chem/bin/python compare_exit.py run_dir step out.png"""
import sys, h5py, numpy as np, yaml, pathlib, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
run, step, out = sys.argv[1], sys.argv[2], sys.argv[3]; HERE = pathlib.Path(__file__).parent
names = ["H2", "O2", "H", "O", "OH", "H2O", "HO2", "H2O2", "N2"]
db = yaml.safe_load(open(HERE / "species_db.yaml")); W = np.array([db[n]["MW"] for n in names])
with h5py.File(f"{run}/bk.h5") as m: xyz = m["MESH/COORD"][:].reshape(-1, 3)
with h5py.File(f"{run}/res_{step}.h5") as f:
    v = f["VALUE"]; N = v["T"].shape[0]
    T = v["T"][:]; P = v["P"][:]; U = np.sqrt(v["Ux"][:]**2 + v["Uy"][:]**2); a = v["sonic"][:]
    Y = np.array([v[f"Y{i}"][:] for i in range(len(names))]).T
    cp = v["cp"][:] if "cp" in v else None
x = xyz[:N, 0]; y = xyz[:N, 1]
xs = 0.356; sel = np.where(np.abs(x - xs) < 2e-4)[0]
if len(sel) < 10: xs = x.max(); sel = np.where(np.abs(x - xs) < 2e-4)[0]
o = np.argsort(y[sel]); sel = sel[o]; yc = y[sel]*100
Xm = Y[sel]/W; Xm /= Xm.sum(1, keepdims=True)
M = U[sel]/a[sel]
Ttot = T[sel] + (0.5*U[sel]**2/cp[sel] if cp is not None else 0.0)   # 近似 (cp 一定): 図の T_ref=2380 K で正規化
ex = np.loadtxt(HERE / "exp_exit_profiles.csv", delimiter=","); et = np.loadtxt(HERE / "exp_exit_ttot.csv", delimiter=",")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
for j, (sp, mk) in enumerate((("H2", "o"), ("O2", "s"), ("N2", "D"), ("H2O", "^"))):
    ax[0].plot(ex[:, 0], ex[:, 1+j], mk, mfc="none", label=f"exp {sp}"); ax[0].plot(yc, Xm[:, names.index(sp)], "-", label=f"forge {sp}")
ax[0].set_xlim(0, 4); ax[0].set_ylim(0, 1); ax[0].set_xlabel("y [cm]"); ax[0].set_ylabel("mole fraction"); ax[0].legend(fontsize=7, ncol=2); ax[0].grid(alpha=.3)
ax[1].plot(et[:, 0], et[:, 1], "ko", mfc="none", label="exp (Fig.13)"); ax[1].plot(yc, Ttot/2380.0, "-", label="forge T_tot/2380K"); ax[1].set_xlim(0, 4); ax[1].set_xlabel("y [cm]"); ax[1].set_ylabel("T_tot/T_ref"); ax[1].legend(); ax[1].grid(alpha=.3)
ax[2].plot(yc, M, "-", label="forge M"); ax[2].set_xlim(0, 4); ax[2].set_xlabel("y [cm]"); ax[2].set_ylabel("Mach"); ax[2].legend(); ax[2].grid(alpha=.3)
fig.suptitle(f"{run} step {step}: exit x=35.6 cm"); fig.tight_layout(); fig.savefig(out, dpi=120)
i = np.argmax(Xm[:, names.index("H2O")]); print(f"wrote {out}; H2O peak X={Xm[i, names.index('H2O')]:.3f} at y={yc[i]:.2f} cm; Ttot max={Ttot.max():.0f} K at y={yc[np.argmax(Ttot)]:.2f} cm; OH max={Y[sel, names.index('OH')].max():.2e}")
