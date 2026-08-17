#!/usr/bin/env python3
"""M4.2 ベスト (R2/L_U6/L_c8) H₂O 5 %: pseudo (run_0063) vs split dry (run_0064) の軸 M 一致と、
split+凝縮 (run_0065) の軸 M/T/g/p 分布 (計画 tooling-nozzle-tp-split-h2o-condensation)。
usage: design/.venv-opt/bin/python case/42.isobutane_wt/compare_h2o5_split.py"""
import sys, json, numpy as np, h5py, matplotlib; matplotlib.use("Agg")
from matplotlib import font_manager as _fm
for _f in ("/home/sano/.fonts/NotoSansCJKjp-Regular.otf",):
    try: _fm.fontManager.addfont(_f)
    except Exception: pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]; matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.metrics.extract import axis_mach
from pathlib import Path
CASE = Path("/home/sano/work/forge/case/42.isobutane_wt")
RUNS = {"pseudo (MIX)": "run_0063_ib_R2_LU6_Lc8_h2o5_pseudo", "split dry [MIXDRY,H2O]": "run_0064_ib_R2_LU6_Lc8_h2o5_split",
        "split + 凝縮 Kw+HK": "run_0065_ib_R2_LU6_Lc8_h2o5_split_cond"}
def axis_fields(run, keys):
    d = CASE / run; info = json.load(open(d / "prepare_info.json")); s = info["scale_m"]
    with h5py.File(d / "nozzle.h5") as m: cc = m["/CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(d / "res_12000.h5") as f: vals = {k: f["/VALUE/" + k][:] for k in keys if k in f["/VALUE"]}
    ax = cc[:, 1] < 0.09 * s * 0.1 * 10   # 軸近傍
    ax = cc[:, 1] < 0.9 * s * 0.1
    x = cc[ax, 0] / s; o = np.argsort(x)
    return x[o], {k: v[ax][o] for k, v in vals.items()}, s
out = {}
xa, fa, s = axis_fields(RUNS["pseudo (MIX)"], ["P", "T"]); xb, fb, _ = axis_fields(RUNS["split dry [MIXDRY,H2O]"], ["P", "T"])
xA, MA = axis_mach(CASE / RUNS["pseudo (MIX)"] / "nozzle.h5", CASE / RUNS["pseudo (MIX)"] / "res_12000.h5", axis_band=0.09 * s)
xB, MB = axis_mach(CASE / RUNS["split dry [MIXDRY,H2O]"] / "nozzle.h5", CASE / RUNS["split dry [MIXDRY,H2O]"] / "res_12000.h5", axis_band=0.09 * s)
xC, MC = axis_mach(CASE / RUNS["split + 凝縮 Kw+HK"] / "nozzle.h5", CASE / RUNS["split + 凝縮 Kw+HK"] / "res_12000.h5", axis_band=0.09 * s)
dM = np.abs(np.interp(xA, xB, MB) - MA); m = xA > 0
out["pseudo_vs_split_maxdM"] = float(dM[m].max())
print(f"pseudo vs split (dry): 軸 M 差 max {dM[m].max():.2e} (x>0), 平均 {dM[m].mean():.2e}")
xc, fc, _ = axis_fields(RUNS["split + 凝縮 Kw+HK"], ["P", "T", "g_0", "Y1"])
gmax = float(np.nanmax(fc["g_0"])); i_on = int(np.argmax(fc["g_0"] > 1e-4)); x_on = xc[i_on] if fc["g_0"].max() > 1e-4 else float("nan")
print(f"split+cond: g_max(軸)={gmax:.4f} (Y_H2O 0.0501 の {100*gmax/0.05011:.0f} %), onset(g>1e-4) x={x_on:.2f} r_t, 出口軸 M {MC[np.argmax(xC)]:.4f} vs dry {MB[np.argmax(xB)]:.4f}")
Td = np.interp(xc, xb, fb["T"]); Pd = np.interp(xc, xb, fb["P"])
print(f"  出口付近: T dry {Td[-1]:.1f} K → cond {fc['T'][-1]:.1f} K, P dry {Pd[-1]/1e3:.2f} kPa → cond {fc['P'][-1]/1e3:.2f} kPa")
out.update({"g_max_axis": gmax, "x_onset_rt": float(x_on), "M_exit_dry": float(MB[np.argmax(xB)]), "M_exit_cond": float(MC[np.argmax(xC)]),
            "T_exit_dry": float(Td[-1]), "T_exit_cond": float(fc["T"][-1]), "P_exit_dry": float(Pd[-1]), "P_exit_cond": float(fc["P"][-1])})
json.dump(out, open(CASE / "compare_h2o5_split.json", "w"), indent=1)
fig, ax = plt.subplots(2, 2, figsize=(12, 7.5))
ax[0, 0].plot(xA / s, MA, "C0", label="pseudo"); ax[0, 0].plot(xB / s, MB, "C1--", label="split dry"); ax[0, 0].plot(xC / s, MC, "C3", label="split+凝縮"); ax[0, 0].set_ylabel("軸 M"); ax[0, 0].legend(fontsize=8)
ax[0, 1].plot(xA / s, np.interp(xA, xB, MB) - MA, "C1"); ax[0, 1].set_ylabel("M(split) − M(pseudo)"); ax[0, 1].set_title(f"dry 等価性: max|ΔM|={dM[m].max():.1e}")
ax[1, 0].plot(xb, fb["T"], "C1--", label="dry"); ax[1, 0].plot(xc, fc["T"], "C3", label="凝縮"); ax[1, 0].set_ylabel("軸 T [K]"); ax[1, 0].legend(fontsize=8)
ax2 = ax[1, 0].twinx(); ax2.plot(xc, fc["g_0"], "C2", lw=1); ax2.set_ylabel("g (液相質量分率)", color="C2")
ax[1, 1].plot(xb, fb["P"] / 1e3, "C1--", label="dry"); ax[1, 1].plot(xc, fc["P"] / 1e3, "C3", label="凝縮"); ax[1, 1].set_ylabel("軸 P [kPa]"); ax[1, 1].set_yscale("log"); ax[1, 1].legend(fontsize=8)
for a in ax.flat: a.grid(alpha=.3); a.set_xlabel("x / r_t")
plt.suptitle("イソブタン M4.2 (R2/L_U6/L_c8), H₂O 5 %: pseudo / split / split+凝縮 (node, TP)")
plt.tight_layout(); plt.savefig(CASE / "compare_h2o5_split.png", dpi=130); print("saved", CASE / "compare_h2o5_split.png")
