#!/usr/bin/env python3
"""SST enhanced wall treatment (wallTreatmentSST) の検証まとめ図を生成する。

2 枚出力:
  ewt_comparison.png  : 壁処理の進化 (mode0 → ω-blend → full WF) の Cf 比較
                        左=Cf-Re_x (modeled τ_w=ρu_τ²), 右=Cf/Schlichting @ x=0.6 棒グラフ
  ewt_uplus_yplus.png : full WF 3 帯 (y⁺₁≈0.35/10/30) の u⁺-y⁺ collapse

Cf は壁関数整合の modeled τ_w=ρu_τ² で評価する (tools/wall_law_modeled.py を流用)。
使い方: python3 tools/plot_ewt_summary.py   (case/26 直下で実行)
"""
import importlib.util, os
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m", os.path.join(HERE, "wall_law_modeled.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
MU, RHO, UINF, KAP = m.MU, m.RHO_INF, m.U_INF, m.KAPPA
ROOT = os.path.dirname(HERE)

def load(run):
    h5 = h5py.File(m.load_latest(os.path.join(ROOT, run)), "r")
    C = m.cell_centroids(h5); x, y = C[:, 0], C[:, 1]
    return x, y, h5["VALUE/Ux"][:], h5["VALUE/ro"][:], h5["VALUE/vis_lam"][:]

def column(run, xs):
    x, y, Ux, ro, vis = load(run)
    xc = np.unique(np.round(x[x > 1e-3], 6)); xc = xc[np.argmin(np.abs(xc - xs))]
    c = np.where(np.abs(x - xc) < 1e-4)[0]; c = c[np.argsort(y[c])]
    return y[c], Ux[c], ro[c], vis[c]

def cfcorr(run, xs=0.6):
    y, Ux, ro, vis = column(run, xs); nu = vis[0] / ro[0]
    ut = m.solve_utau(Ux[0], y[0], nu); Rex = RHO * UINF * xs / MU
    return (ro[0] * ut**2 / (0.5 * RHO * UINF**2)) / (0.0592 * Rex**-0.2)

# ---------- Figure 1: progression ----------
runs = [
    ("run_0007_slau_muscl_innersweep", "y+~0.35 wall-resolved (ref)", "#000000", "-"),
    ("run_0010_ewt_yp30_mode0", "y+~30 mode0 (low-Re)",          "#d62728", ":"),
    ("run_0011_ewt_yp30_mode1", "y+~30 mode1 (ω-blend, k=0)",    "#d62728", "--"),
    ("run_0017_ewt_yp30_wf",    "y+~30 full WF",                 "#d62728", "-"),
    ("run_0012_ewt_yp10_mode0", "y+~10 mode0 (low-Re)",          "#ff7f0e", ":"),
    ("run_0013_ewt_yp10_mode1", "y+~10 mode1 (ω-blend, k=0)",    "#ff7f0e", "--"),
    ("run_0018_ewt_yp10_wf",    "y+~10 full WF",                 "#ff7f0e", "-"),
]
fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 5.8))
for run, lab, col, ls in runs:
    x, y, Ux, ro, vis = load(run)
    xcols = np.unique(np.round(x[x > 1e-3], 6)); rex = []; cf = []
    for xc in xcols:
        c = np.where(np.abs(x - xc) < 1e-4)[0]; c = c[np.argsort(y[c])]
        if c.size < 5: continue
        nu = vis[c][0] / ro[c][0]; ut = m.solve_utau(Ux[c][0], y[c][0], nu)
        rex.append(RHO * UINF * xc / MU); cf.append(ro[c][0] * ut**2 / (0.5 * RHO * UINF**2))
    rex = np.array(rex); cf = np.array(cf); o = np.argsort(rex)
    lw = 2.2 if "full WF" in lab else (1.6 if "ref" in lab else 1.2)
    axA.plot(rex[o], cf[o], ls, color=col, lw=lw, label=lab)
re = np.logspace(6, 6.62, 80)
axA.plot(re, 0.0592 * re**-0.2, color="gray", lw=2.5, alpha=0.5, label="Schlichting")
axA.set_xscale("log"); axA.set_yscale("log"); axA.set_xlim(1e6, 4.2e6); axA.set_ylim(2e-4, 7e-3)
axA.set_xlabel("$Re_x$"); axA.set_ylabel(r"$C_f$ (modeled $\tau_w=\rho u_\tau^2$)")
axA.set_title("Skin friction: low-Re collapses; ω-blend overshoots y+10; full WF on correlation")
axA.legend(fontsize=7.2, loc="lower left"); axA.grid(True, which="both", alpha=0.3)

bands = [("fine", ["run_0007_slau_muscl_innersweep", "run_0016_ewt_fine_wf"], ["wall-resolved", "full WF"]),
         ("buffer", ["run_0012_ewt_yp10_mode0", "run_0013_ewt_yp10_mode1", "run_0018_ewt_yp10_wf"], ["mode0", "ω-blend k=0", "full WF"]),
         ("log", ["run_0010_ewt_yp30_mode0", "run_0011_ewt_yp30_mode1", "run_0017_ewt_yp30_wf"], ["mode0", "ω-blend k=0", "full WF"])]
colmap = {"mode0": "#d62728", "ω-blend k=0": "#ff7f0e", "full WF": "#2ca02c", "wall-resolved": "#1f77b4"}
xpos = 0; xt = []; xl = []
for band, rs, labs in bands:
    for run, lb in zip(rs, labs):
        v = cfcorr(run); axB.bar(xpos, v, color=colmap[lb], width=0.8)
        axB.text(xpos, v + 0.02, f"{v:.2f}", ha="center", fontsize=7.5); xt.append(xpos); xl.append(lb); xpos += 1
    xpos += 0.8
axB.axhline(1.0, color="gray", lw=2, ls="--", label="Schlichting (=1.0)")
axB.set_xticks(xt); axB.set_xticklabels(xl, rotation=40, ha="right", fontsize=7)
axB.set_ylabel(r"$C_f / C_{f,\mathrm{Schlichting}}$ @ x=0.6"); axB.set_ylim(0, 2.0)
axB.set_title("Cf/correlation: full WF (green) lands on 1.0 at every band")
axB.legend(fontsize=8, loc="upper right"); axB.grid(True, axis="y", alpha=0.3)
axB.text(0.5, 1.9, "fine y+0.35", ha="center", fontsize=8, color="#444")
axB.text(3.8, 1.9, "buffer y+10", ha="center", fontsize=8, color="#444")
axB.text(7.8, 1.9, "log y+30", ha="center", fontsize=8, color="#444")
fig.suptitle("SST wall treatment progression (case/26 flat plate): mode0 → ω-blend → full wall-function", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(ROOT, "ewt_comparison.png"), dpi=130)
print("wrote", os.path.join(ROOT, "ewt_comparison.png"))

# ---------- Figure 2: u+-y+ collapse ----------
fig2, ax = plt.subplots(figsize=(8.5, 6.2))
yl = np.logspace(-1, 3.3, 300)
ax.plot(yl, yl, "k:", lw=1.2, label=r"$u^+=y^+$ (viscous)")
ax.plot(yl, np.log(yl) / KAP + 5.0, "k-", lw=1.2, label=r"$\frac{1}{0.41}\ln y^+ + 5.0$ (log law)")
ax.plot(yl, m.reichardt_uplus(yl), "-", color="gray", lw=3, alpha=0.45, label="Reichardt (blended)")
wf = [("run_0016_ewt_fine_wf", "y+₁≈0.35 (fine)", "#1f77b4", "o"),
      ("run_0018_ewt_yp10_wf", "y+₁≈10 (buffer)", "#ff7f0e", "s"),
      ("run_0017_ewt_yp30_wf", "y+₁≈30 (log)",   "#2ca02c", "^")]
for run, lab, c, mk in wf:
    y, Ux, ro, vis = column(run, 0.6); nu = vis[0] / ro[0]; ut = m.solve_utau(Ux[0], y[0], nu)
    yp = ro * ut * y / vis; up = Ux / ut; sel = yp < 3000
    ax.plot(yp[sel], up[sel], mk, color=c, ms=5, mfc="none", mew=1.3, label=f"{lab}, $u_\\tau$={ut:.2f}")
ax.set_xscale("log"); ax.set_xlim(0.1, 3000); ax.set_ylim(0, 28)
ax.set_xlabel(r"$y^+$", fontsize=12); ax.set_ylabel(r"$u^+$", fontsize=12)
ax.set_title("SST full wall-function (wallTreatmentSST=1): $u^+$–$y^+$ @ x=0.6 m\n"
             "first-cell $y^+$≈0.35/10/30 all collapse onto the universal profile", fontsize=11)
ax.legend(fontsize=9, loc="upper left"); ax.grid(True, which="both", alpha=0.3)
fig2.tight_layout(); fig2.savefig(os.path.join(ROOT, "ewt_uplus_yplus.png"), dpi=140)
print("wrote", os.path.join(ROOT, "ewt_uplus_yplus.png"))
