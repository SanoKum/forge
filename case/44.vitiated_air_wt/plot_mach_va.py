"""全 run の軸 M (目標 vs CFD) と Mach コンタを figs/ にまとめる。
実行: design/.venv-opt/bin/python case/44.vitiated_air_wt/plot_mach_va.py"""
import json, sys
from pathlib import Path
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.metrics.extract import axis_mach
CASE = Path("/home/sano/work/forge/case/44.vitiated_air_wt"); FIG = CASE / "figs"; FIG.mkdir(exist_ok=True)
S = 0.210505; MD = 4.19
runs = sorted([r for r in json.loads((CASE / "study_cfd_va.json").read_text()) if r.get("ok")],
              key=lambda r: (r["R"], r["L_U"], r["L_c"]))

def load(run):
    rd = CASE / run
    info = json.loads((rd / "prepare_info.json").read_text())
    with h5py.File(rd / "nozzle.h5") as n: xy = n["MESH/COORD"][:].reshape(-1, 3)[:, :2]
    with h5py.File(rd / "res_12000.h5") as f:
        M = np.hypot(f["VALUE/Ux"][:], f["VALUE/Uy"][:]) / f["VALUE/sonic"][:]
    ni, nj = info["mesh"]["ni"], info["mesh"]["nj"]
    # COORD は i-major (i*nj + j; j=0 が軸, j=nj-1 が壁) → (nj, ni) に転置して X[-1] を壁行にする
    X = xy[:, 0].reshape(ni, nj).T; Y = xy[:, 1].reshape(ni, nj).T; Mg = M.reshape(ni, nj).T
    assert np.all(np.diff(X, axis=1) > 0) and np.all(np.diff(Y, axis=0) > 0)
    tgt = np.loadtxt(rd / "target_axis_M.csv", delimiter=",", skiprows=1)
    xa, Ma = axis_mach(rd / "nozzle.h5", rd / "res_12000.h5", axis_band=0.09 * S)
    return info, X, Y, Mg, tgt, xa, Ma

# --- 1. 軸 M: 目標 vs CFD (15 パネル) ---
fig, axs = plt.subplots(3, 5, figsize=(20, 10), sharey=True)
for ax, r in zip(axs.flat, runs):
    info, X, Y, Mg, tgt, xa, Ma = load(r["run"])
    ok = np.isfinite(tgt[:, 1])
    ax.plot(xa / S, Ma, "r", lw=1.2, label="CFD (node Euler TP)")
    ax.plot(tgt[ok, 0] / S, tgt[ok, 1], "k--", lw=1.2, label="target (quintic law)")
    ax.axvline(info["x_A"], c="gray", ls=":", lw=0.8); ax.axvline(info["x_E"], c="gray", ls=":", lw=0.8)
    ax.axhline(MD, c="b", ls=":", lw=0.8)
    ax2 = ax.twinx(); d = np.loadtxt(CASE / r["run"] / "achieved_vs_target.csv", delimiter=",", skiprows=1)
    ax2.plot(d[:, 0] / S, 100 * (d[:, 2] - d[:, 1]) / MD, "g", lw=0.8); ax2.set_ylim(-0.4, 0.4)
    ax2.axhline(0, c="g", lw=0.3); ax2.tick_params(axis="y", colors="g", labelsize=7)
    ax.set_title(f"{r['tag']}  ‖ΔM‖∞={100*r['dM']:.3f}%  M_exit={r['M_exit_axis']:.4f}", fontsize=9)
    ax.set_xlim(0, xa.max() / S); ax.grid(alpha=0.3); ax.set_xlabel("x / r_t", fontsize=8)
axs[0, 0].set_ylabel("axis M"); axs[0, 0].legend(fontsize=7, loc="lower right")
fig.suptitle("case/44 M4.19 vitiated air — axis Mach: target (black dashed) vs CFD (red); green = (M_cfd−M_tgt)/M_d [%] (right axis, ±0.4%); dotted x_A / x_E", fontsize=11)
fig.tight_layout(); fig.savefig(FIG / "axis_mach_all.png", dpi=110); plt.close(fig)

# --- 2. Mach コンタ: 推奨 run (大) ---
def contour(ax, X, Y, Mg, title, wall=None, levels=np.arange(0, 4.41, 0.1), cbar=True):
    cs = ax.contourf(X / S, Y / S, Mg, levels=levels, cmap="turbo", extend="max")
    ax.contour(X / S, Y / S, Mg, levels=np.arange(0.5, 4.5, 0.5), colors="k", linewidths=0.4)
    ax.contour(X / S, Y / S, Mg, levels=[0.995 * MD], colors="w", linewidths=1.0)
    ax.plot(X[-1] / S, Y[-1] / S, "k", lw=1); ax.set_aspect("equal"); ax.set_title(title, fontsize=9)
    ax.set_xlabel("x / r_t"); ax.set_ylabel("r / r_t")
    return cs
for tag in ("run_0005_va_R2_LU6_Lc8", "run_0013_va_R2_LU6_Lc9"):
    info, X, Y, Mg, tgt, xa, Ma = load(tag)
    fig, ax = plt.subplots(2, 1, figsize=(16, 8), gridspec_kw={"height_ratios": [1.6, 1]})
    cs = contour(ax[0], X, Y, Mg, f"{tag}: Mach (node Euler TP, res_12000; white line = 0.995 M_d (4.169), black = 0.5 step)")
    fig.colorbar(cs, ax=ax[0], fraction=0.02, pad=0.01, label="M")
    ax[0].axvline(info["x_E"], c="w", ls=":", lw=0.8); ax[0].text(info["x_E"], 3.9, " x_E", color="w", fontsize=8)
    ok = np.isfinite(tgt[:, 1])
    ax[1].plot(xa / S, Ma, "r", lw=1.2, label="CFD axis M"); ax[1].plot(tgt[ok, 0] / S, tgt[ok, 1], "k--", label="target")
    ax[1].axhline(MD, c="b", ls=":", lw=0.8); ax[1].set_xlim(ax[0].get_xlim()); ax[1].set_ylabel("axis M"); ax[1].set_xlabel("x / r_t")
    ax[1].grid(alpha=0.3); ax[1].legend(loc="lower right")
    ax1b = ax[1].twinx(); d = np.loadtxt(CASE / tag / "achieved_vs_target.csv", delimiter=",", skiprows=1)
    ax1b.plot(d[:, 0] / S, 100 * (d[:, 2] - d[:, 1]) / MD, "g", lw=0.9); ax1b.set_ylim(-0.3, 0.3); ax1b.set_ylabel("(M_cfd−M_tgt)/M_d [%]", color="g")
    fig.tight_layout(); fig.savefig(FIG / f"mach_contour_{tag}.png", dpi=120); plt.close(fig)
    # 拡張部だけ拡大 (x 0..x_F)
    fig, ax = plt.subplots(figsize=(16, 4.5))
    m = X[0] >= -0.5
    cs = contour(ax, X[:, m], Y[:, m], Mg[:, m], f"{tag}: Mach in the expansion (x ≥ 0)")
    fig.colorbar(cs, ax=ax, fraction=0.02, pad=0.01, label="M"); fig.tight_layout()
    fig.savefig(FIG / f"mach_contour_{tag}_expansion.png", dpi=120); plt.close(fig)

# --- 3. Mach コンタ 全 run 一覧 ---
fig, axs = plt.subplots(15, 1, figsize=(14, 30), constrained_layout=True)
for k, (ax, r) in enumerate(zip(axs, runs)):
    info, X, Y, Mg, tgt, xa, Ma = load(r["run"])
    m = X[0] >= -0.5
    cs = contour(ax, X[:, m], Y[:, m], Mg[:, m], f"{r['run']}   ‖ΔM‖∞={100*r['dM']:.3f}%  ε_M={100*r['eps_M']:.3f}%  ε_θ={r['eps_th']:.4f}°  L={r['L_total_m']:.2f} m")
    ax.set_xlim(-0.5, 27.5)
    if k < len(runs) - 1: ax.set_xlabel("")
fig.colorbar(cs, ax=axs, fraction=0.015, pad=0.01, label="M", shrink=0.4)
fig.suptitle("case/44 — Mach contours of all 15 runs (res_12000; white = 0.995 M_d, black = 0.5 step)")
fig.savefig(FIG / "mach_contour_all.png", dpi=90); plt.close(fig)

# --- 4. 全 run の CFD 軸 M を 1 枚に重ね (目標は run ごとに違うので CFD のみ + 推奨の目標) ---
fig, ax = plt.subplots(figsize=(12, 5))
for r in runs:
    info, X, Y, Mg, tgt, xa, Ma = load(r["run"])
    ax.plot(xa / S, Ma, lw=0.9, label=r["tag"], alpha=0.9)
info, X, Y, Mg, tgt, xa, Ma = load("run_0005_va_R2_LU6_Lc8"); ok = np.isfinite(tgt[:, 1])
ax.plot(tgt[ok, 0] / S, tgt[ok, 1], "k--", lw=1.5, label="target R2_LU6_Lc8")
ax.axhline(MD, c="b", ls=":"); ax.set_xlabel("x / r_t"); ax.set_ylabel("axis M"); ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=3)
ax.set_title("case/44 — CFD axis Mach of all 15 runs")
fig.tight_layout(); fig.savefig(FIG / "axis_mach_overlay.png", dpi=120)
print("wrote", sorted(p.name for p in FIG.iterdir()))
