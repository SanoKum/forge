#!/usr/bin/env python3
"""TP split (MIXDRY+H2O, run_0170) vs 既存 CPG [N2,H2O] (run_0050) の中心線 p/p0・g 比較
(計画 tooling-nozzle-tp-split-h2o-condensation)。実験 (Fig.3 1 kPa) も重ねる。
usage: python3 compare_tpsplit_wys.py [RUN_TP] [RUN_CPG] [--out PNG]"""
import sys, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/sano/work/forge/case/16.nozzle_wys")
from compare_cond_models import centerline, binavg, load_exp
run_tp = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "run_0170_fig3_2d_sst_kwhk_tpsplit"
run_cpg = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else "run_0050_fig3_2d_sst_kwhk"
out = sys.argv[sys.argv.index("--out")+1] if "--out" in sys.argv else "compare_tpsplit_wys.png"
xe = np.linspace(-2, 8, 101)
res = {}
import glob, h5py, os
def centerline_planar(run):
    from compare_cond_models import latest
    r = latest(run); m = [p for p in glob.glob(f"{run}/nozzle_*.h5") if "res_" not in os.path.basename(p)][0]
    cc = h5py.File(m)["/CELLS/centCoords"][:].reshape(-1, 3); v = h5py.File(r)["/VALUE"]
    P = np.array(v["P"]); g = np.array(v["g_0"]) if "g_0" in v else np.zeros_like(P)
    sel = np.abs(cc[:, 1]) < 0.0025 * 0.06
    x = cc[sel, 0] * 100; o = np.argsort(x); return x[o], P[sel][o] / 59070.0, g[sel][o], os.path.basename(r)
runs = [("cell TP split SST [MIXDRY,H2O]", run_tp, centerline, "C3"), ("cell CPG SST [N2,H2O] (run_0050)", run_cpg, centerline, "C0")]
extra = [("cell CPG inviscid (run_0184)", "run_0184_wysinv_cell_cpg_cond", centerline_planar, "C7"),
         ("node CPG inviscid (run_0185)", "run_0185_wysinv_node_cpg_cond", centerline_planar, "C2"),
         ("node TP split inviscid (run_0186)", "run_0186_wysinv_node_tpsplit_cond", centerline_planar, "C1")]
for tag, run, fn, c in runs + [e for e in extra if os.path.isdir(e[1])]:
    x, pp0, g, name = fn(run)
    xc, pb = binavg(x, pp0, xe); _, gb = binavg(x, g, xe)
    res[tag] = (xc, pb, gb, name, c)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
try:
    ex = load_exp()
    for k, v in ex.items():
        ax[0].plot(v[0], v[1], "o" if "1" in k else "s", ms=3, mfc="none", color="k", label=f"exp {k}")
except Exception:
    pass
for tag, (xc, pb, gb, name, c) in res.items():
    ls = "--" if "inviscid" in tag else "-"
    ax[0].plot(xc, pb, color=c, lw=1.4, ls=ls, label=f"{tag}"); ax[1].plot(xc, gb, color=c, lw=1.4, ls=ls, label=tag)
ax[0].set_xlabel("x from throat [cm]"); ax[0].set_ylabel("p/p0 (centerline)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3); ax[0].set_ylim(0, 0.6)
ax[1].set_xlabel("x from throat [cm]"); ax[1].set_ylabel("g (liquid mass fraction)"); ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(out, dpi=130)
xc, pT, gT, _, _ = res["cell TP split SST [MIXDRY,H2O]"]; _, pC, gC, _, _ = res["cell CPG SST [N2,H2O] (run_0050)"]
m = np.isfinite(pT) & np.isfinite(pC) & (xc > 0)
print(f"p/p0 差 (x>0): max|Δ|={np.nanmax(np.abs(pT[m]-pC[m])):.4f}  rel max={np.nanmax(np.abs(pT[m]-pC[m])/pC[m])*100:.2f}%")
print(f"g 差: max|Δ|={np.nanmax(np.abs(gT[m]-gC[m])):.4f}  g_max TP={np.nanmax(gT[m]):.4f} CPG={np.nanmax(gC[m]):.4f}")
i=np.argmax(np.nan_to_num(gT)>0.001); j=np.argmax(np.nan_to_num(gC)>0.001)
print(f"onset (g>1e-3): TP x={xc[i]:.2f} cm, CPG x={xc[j]:.2f} cm")
print("saved", out)
