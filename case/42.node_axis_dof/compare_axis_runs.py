"""軸 M(x) の run 間比較: node 軸直読 / 偶関数外挿 / cell 参照 (偶関数外挿) / 設計目標。"""
import json, sys, numpy as np, h5py
from pathlib import Path
sys.path.insert(0, "/home/sano/work/forge-axisnode/design")
from forge_design.evaluate.axis_extract import _columns, even_fit_axis
from scipy.interpolate import make_smoothing_spline

SC = 0.01
REF = Path("run_0001_ref41_0057_dirichlet")
target = np.loadtxt(REF / "target_axis_M.csv", delimiter=",", skiprows=1)
xt, Mt = target[:, 0] / SC, target[:, 1]

def profile(rd, kind):
    rd = Path(rd)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split('_')[1]))[-3:]
    with h5py.File(rd / "nozzle.h5") as nz:
        pos = nz["/MESH/COORD"][:].reshape(-1, 3) if kind == "node" else nz["/CELLS/centCoords"][:].reshape(-1, 3)
    Ms = []
    for rf in res:
        with h5py.File(rf) as f:
            Ms.append(np.hypot(f["/VALUE/Ux"][:], f["/VALUE/Uy"][:]) / f["/VALUE/sonic"][:])
    M = np.mean(Ms, axis=0).astype(float)
    ux, cols = _columns(pos)
    x, fit, ax = [], [], []
    for xv, idx in zip(ux, cols):
        if len(idx) < 6: continue
        r = pos[idx, 1]
        if kind == "node":
            if r[0] > 1e-12: continue
            x.append(xv / SC); ax.append(M[idx[0]]); fit.append(even_fit_axis(r, M[idx], 4, 2, True))
        else:
            x.append(xv / SC); ax.append(np.nan); fit.append(even_fit_axis(r, M[idx], 4, 2, False))
    o = np.argsort(x)
    return np.asarray(x)[o], np.asarray(fit)[o], np.asarray(ax)[o]

def dM_metric(x, M):
    """設計目標との比較 (runner と同じ平滑化スプライン, 目標 x 範囲)."""
    spl = make_smoothing_spline(x, M, lam=1e-5)
    m = (xt >= x.min()) & (xt <= x.max())
    d = spl(xt[m]) - Mt[m]
    return float(np.abs(d).max()), float(np.sqrt(np.mean(d**2))), float(np.abs(d).max() / 4 * 100)

runs = [("run_0001_ref41_0057_dirichlet", "node", "Dirichlet (production, =case/41 run_0057)"),
        ("run_0002_axisdof_euler", "node", "axis DOF, centroid値位置 (nodeAxisDirichlet 0)"),
        ("run_0003_axisdof_van_euler", "node", "axis DOF + nodeValueAtNode 1 (GG)"),
        ("run_0004_axisdof_van_lsq", "node", "axis DOF + nodeValueAtNode 1 + gradLSQ 2"),
        ("run_0005_cell_ref", "cell", "cell 参照 (同メッシュ)")]
P = {}
for rd, kind, label in runs:
    if not Path(rd).exists(): continue
    P[rd] = profile(rd, kind) + (label,)
xc, fc, _, _ = P["run_0005_cell_ref"]
print(f"{'run':34s} {'‖ΔM‖∞ vs target (axis-read)':>28s} {'(evenfit)':>10s} | max|axis-evenfit| | max|evenfit - cell evenfit| (x<22.5)")
for rd, (x, fit, ax, label) in P.items():
    dm_fit = dM_metric(x, fit)
    dm_ax = dM_metric(x, ax) if not np.isnan(ax).any() else (np.nan,)*3
    fci = np.interp(x, xc, fc); m = x < 22.5
    dcell = np.abs(fit - fci)[m].max()
    dax = np.abs(ax - fit)[m].max() if not np.isnan(ax).any() else np.nan
    print(f"{rd:34s} {dm_ax[2]:>14.3f}%Md {dm_fit[2]:>14.3f}%Md | {dax:>14.4f} | {dcell:.4f}   [{label}]")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
for rd, (x, fit, ax, label) in P.items():
    axs[0].plot(x, fit - np.interp(x, xc, fc), label=label + " [evenfit − cell evenfit]")
    if not np.isnan(ax).any(): axs[1].plot(x, ax - fit, label=label + " [axis node − evenfit]")
axs[0].axhline(0, c='k', lw=.5); axs[1].axhline(0, c='k', lw=.5)
axs[0].set_ylabel("ΔM (evenfit vs cell ref)"); axs[1].set_ylabel("ΔM (axis node − own evenfit)"); axs[1].set_xlabel("x / r_t")
axs[0].set_ylim(-0.02, 0.02); axs[1].set_ylim(-0.04, 0.02)
axs[0].legend(fontsize=7); axs[1].legend(fontsize=7); plt.tight_layout(); plt.savefig("axis_runs_compare.png", dpi=120)
