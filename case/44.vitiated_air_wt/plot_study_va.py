"""study_cfd_va.json + study_design_va.json → 表 (markdown) と図 study_cfd_va.png。
実行: design/.venv-opt/bin/python case/44.vitiated_air_wt/plot_study_va.py"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
CASE = Path("/home/sano/work/forge/case/44.vitiated_air_wt")
cfd = [r for r in json.loads((CASE / "study_cfd_va.json").read_text()) if r.get("ok")]
des = json.loads((CASE / "study_design_va.json").read_text())
print("| run | R | L_U | L_c | θ_max° | κ₀R | x_F [r_t] | 全長 [m] | ‖ΔM‖∞/M_d | ε_M | ε_θ° | over | M_exit |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in sorted(cfd, key=lambda r: (r["R"], r["L_U"], r["L_c"])):
    print(f"| `{r['run']}` | {r['R']:g} | {r['L_U']:g} | {r['L_c']:g} | {r['theta_max']:.2f} | {r['kappa0R']:.3f} | {r['x_F']:.2f} | {r['L_total_m']:.3f} | "
          f"**{100*r['dM']:.3f} %** | {100*r['eps_M']:.3f} % | {r['eps_th']:.4f} | {100*r['over']:+.3f} % | {r['M_exit_axis']:.4f} |")
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for R, c in ((1.5, "C0"), (2.0, "C1"), (3.0, "C2")):
    s = sorted([r for r in cfd if r["R"] == R and r["L_U"] == 6.0], key=lambda r: r["L_c"])
    if not s: continue
    ax[0].plot([r["L_c"] for r in s], [100*r["dM"] for r in s], "o-", c=c, label=f"R={R:g}")
    ax[1].plot([r["L_c"] for r in s], [100*r["eps_M"] for r in s], "o-", c=c, label=f"R={R:g}")
    ax[2].plot([r["L_c"] for r in s], [r["eps_th"] for r in s], "o-", c=c, label=f"R={R:g}")
lu = sorted([r for r in cfd if r["R"] == 2.0 and r["L_c"] == 8.0], key=lambda r: r["L_U"])
for r in lu:
    if r["L_U"] != 6.0:
        ax[0].plot(r["L_c"], 100*r["dM"], "s", c="C1", mfc="none"); ax[0].annotate(f"L_U{r['L_U']:g}", (r["L_c"], 100*r["dM"]), fontsize=7)
        ax[1].plot(r["L_c"], 100*r["eps_M"], "s", c="C1", mfc="none"); ax[1].annotate(f"L_U{r['L_U']:g}", (r["L_c"], 100*r["eps_M"]), fontsize=7)
        ax[2].plot(r["L_c"], r["eps_th"], "s", c="C1", mfc="none"); ax[2].annotate(f"L_U{r['L_U']:g}", (r["L_c"], r["eps_th"]), fontsize=7)
ax[0].axhline(0.5, ls="--", c="k", lw=0.8); ax[0].set_ylabel("‖ΔM‖∞ / M_d [%]"); ax[0].set_ylim(0, None)
ax[1].set_ylabel("exit ε_M rms [%]"); ax[2].set_ylabel("exit ε_θ max [deg]")
for a_ in ax: a_.set_xlabel("L_c [r_t]"); a_.grid(alpha=0.3); a_.legend(fontsize=8)
fig.suptitle("case/44 M4.19 vitiated air — node Euler TP, 12000 step (L_U=6; squares = L_U 4/9 at R2/L_c8)")
fig.tight_layout(); fig.savefig(CASE / "study_cfd_va.png", dpi=130)
# 軸 M の達成 vs 目標
fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
for r in sorted(cfd, key=lambda r: (r["R"], r["L_U"], r["L_c"])):
    d = np.loadtxt(CASE / r["run"] / "achieved_vs_target.csv", delimiter=",", skiprows=1)
    ax[0].plot(d[:, 0]/0.210505, 100*(d[:, 2]-d[:, 1])/4.19, lw=0.9, label=r["tag"])
    if r["tag"] == "R2_LU6_Lc8":
        ax[1].plot(d[:, 0]/0.210505, d[:, 1], "k--", label="target"); ax[1].plot(d[:, 0]/0.210505, d[:, 2], "r", lw=0.8, label="CFD "+r["tag"])
ax[0].set_xlabel("x / r_t"); ax[0].set_ylabel("(M_cfd − M_target)/M_d [%]"); ax[0].legend(fontsize=6, ncol=3); ax[0].grid(alpha=0.3)
ax[1].set_xlabel("x / r_t"); ax[1].set_ylabel("axis M"); ax[1].legend(); ax[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(CASE / "study_axisM_va.png", dpi=130)
print("wrote study_cfd_va.png study_axisM_va.png")
