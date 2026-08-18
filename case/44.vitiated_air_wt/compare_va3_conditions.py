"""va3 L_c8 形状 (1 つの実機ノズル) を **新条件 (Pt 1.139 MPa / Tt 1161 K / 新組成)** と
**一個前の条件 (va2: Pt 1.137 MPa / Tt 1058 K / 旧組成)** で回した 6 run の軸中心分布を比較する。

出力:
  - `axis_compare_va3_vs_va3old.csv` … 6 run の軸分布を 1 ファイルに (x/r_t 昇順、run ごとに列群)
  - `study_va3_conditions.json`       … 主要ステーションの値
  - `figs/va3_cond_axis_compare.png`  … 軸 M / T / P / S / g / (T−Tsat) の 6 段比較
  - `figs/va3_cond_exit_profile.png`  … 物理出口断面 (2 r_t 上流) の M(r), T(r), g(r)
実行: design/.venv-opt/bin/python case/44.vitiated_air_wt/compare_va3_conditions.py
"""
import json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

try:
    fm.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
    matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
except Exception:
    pass

CASE = Path("/home/sano/work/forge/case/44.vitiated_air_wt")
FIG = CASE / "figs"
tz = np.trapezoid

# (ラベル, run, 条件タグ, H2O 質量分率, 色, 線種)
RUNS = [
    ("新条件 dry",   "run_0091_va3_M4.19_Lc8_dry",   "new", 0.0376954, "C0", "-"),
    ("新条件 noneq", "run_0092_va3_M4.19_Lc8_noneq", "new", 0.0376954, "C1", "--"),
    ("新条件 eq",    "run_0094_va3_M4.19_Lc8_eq",    "new", 0.0376954, "C2", "-."),
    ("旧条件 dry",   "run_0095_va3old_Lc8_dry",      "old", 0.0286647, "C3", "-"),
    ("旧条件 noneq", "run_0096_va3old_Lc8_noneq",    "old", 0.0286647, "C4", "--"),
    ("旧条件 eq",    "run_0097_va3old_Lc8_eq",       "old", 0.0286647, "C5", "-."),
]


def axis_cols(run):
    rd = CASE / run
    txt = (rd / "axis_values.csv").read_text().splitlines()
    hdr = txt[1].split(",")
    d = np.loadtxt(txt[2:], delimiter=",")
    q = {k: d[:, j] for j, k in enumerate(hdr)}
    q["M"] = np.hypot(q["Ux"], q["Uy"]) / q["sonic"]
    if "g_0" not in q:
        q["g_0"] = np.zeros_like(q["M"])
    # 凝縮 off の run は forge の condS/condTsat を持たないので後処理値を使う
    q["S"] = q.get("condS_0", q["S_post"])
    q["Tsat"] = q.get("condTsat_0", q["Tsat_post"])
    return q, json.loads((rd / "prepare_info.json").read_text())


def exit_plane(run, dx=2.0):
    """物理出口 (出口 BC artefact を避けて dx r_t 上流) の断面。"""
    rd = CASE / run
    info = json.loads((rd / "prepare_info.json").read_text())
    S = info["scale_m"]
    nc = h5py.File(rd / "nozzle.h5")["MESH/COORD"][:].reshape(-1, 3)
    f = h5py.File(rd / "res_12000.h5")
    cols = np.unique(nc[:, 0])
    xc = cols[np.argmin(abs(cols - (nc[:, 0].max() - dx * S)))]
    idx = np.where(np.abs(nc[:, 0] - xc) < 1e-9)[0]
    idx = idx[np.argsort(nc[idx, 1])]
    g = f["VALUE/g_0"][:][idx] if "VALUE/g_0" in f else np.zeros(len(idx))
    out = dict(x_rt=xc / S, r_rt=nc[idx, 1] / S,
               M=np.hypot(f["VALUE/Ux"][:][idx], f["VALUE/Uy"][:][idx]) / f["VALUE/sonic"][:][idx],
               T=f["VALUE/T"][:][idx], P=f["VALUE/P"][:][idx], ro=f["VALUE/ro"][:][idx],
               Ux=f["VALUE/Ux"][:][idx], g=g)
    w = out["ro"] * out["Ux"] * out["r_rt"]
    for k in ("M", "T", "P", "g"):
        out[k + "_avg"] = float(tz(w * out[k], out["r_rt"]) / tz(w, out["r_rt"]))
    core = out["r_rt"] < 0.85 * out["r_rt"].max()
    out["core_spread"] = float(100 * (out["M"][core].max() - out["M"][core].min()) / out["M_avg"])
    return out


rows = []
axis = {}
for lab, run, tag, Yw, c, ls in RUNS:
    q, info = axis_cols(run)
    e = exit_plane(run)
    axis[run] = (lab, q, info, c, ls)
    x = q["x_over_rt"]
    sup = x > 0
    on = np.where(q["g_0"] > 1e-4)[0]
    rows.append(dict(label=lab, run=run, cond=tag, Yw=Yw,
                     x_E=info["x_E"], M_xE=float(np.interp(info["x_E"], x, q["M"])),
                     T_min=float(q["T"][sup].min()), S_max=float(q["S"][sup].max()),
                     subcool_max=float((q["Tsat"] - q["T"])[sup].max()),
                     g_exit=float(q["g_0"][-1]), g_frac=float(q["g_0"][-1] / Yw),
                     onset_x=float(x[on[0]]) if len(on) else float("nan"),
                     onset_T=float(q["T"][on[0]]) if len(on) else float("nan"),
                     x_exit=e["x_rt"], M_exit_axis=float(e["M"][0]), T_exit_axis=float(e["T"][0]),
                     P_exit_axis=float(e["P"][0]),
                     M_exit_avg=e["M_avg"], T_exit_avg=e["T_avg"], P_exit_avg=e["P_avg"],
                     g_exit_avg=e["g_avg"] / Yw, core_spread=e["core_spread"]))

print("| 条件 | 種別 | run | x_E の軸 M | 評価断面 x/r_t | 軸 M | 軸 T [K] | 軸 P [Pa] | **質量平均 M** | 平均 T [K] | 平均 P [Pa] | コア M 幅 | 軸 T min | S max | 過冷却 max [K] | g/Y_H2O (平均) | onset x (T) |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for w in rows:
    on = "—" if np.isnan(w["onset_x"]) else f"{w['onset_x']:.1f} ({w['onset_T']:.0f} K)"
    print(f"| {'新' if w['cond']=='new' else '旧'} | {w['label'].split()[-1]} | `{w['run']}` | {w['M_xE']:.4f} | {w['x_exit']:.1f} | "
          f"{w['M_exit_axis']:.3f} | {w['T_exit_axis']:.1f} | {w['P_exit_axis']:.0f} | **{w['M_exit_avg']:.3f}** | {w['T_exit_avg']:.1f} | {w['P_exit_avg']:.0f} | "
          f"{w['core_spread']:.1f} % | {w['T_min']:.1f} | {w['S_max']:.2f} | {w['subcool_max']:.1f} | {100*w['g_exit_avg']:.1f} % | {on} |")
json.dump([{k: (float(v) if isinstance(v, (np.floating, np.integer)) else v) for k, v in w.items()} for w in rows],
          open(CASE / "study_va3_conditions.json", "w"), indent=1)

# --- 6 run の軸分布を 1 CSV に ---
ref = axis[RUNS[0][1]][1]["x_over_rt"]
cols = {"x_over_rt": ref, "x_m": axis[RUNS[0][1]][1]["x_m"]}
for lab, run, tag, Yw, c, ls in RUNS:
    q = axis[run][1]
    tagn = run.split("_", 2)[2]
    for k, kk in (("M", "M"), ("T", "T"), ("P", "P"), ("ro", "rho"), ("S", "S"), ("Tsat", "Tsat"), ("g_0", "g")):
        cols[f"{kk}__{tagn}"] = np.interp(ref, q["x_over_rt"], q[k])
hdr = ",".join(cols)
np.savetxt(CASE / "axis_compare_va3_vs_va3old.csv", np.column_stack(list(cols.values())),
           delimiter=",", header="va3 L_c8 形状の軸中心分布: 新条件 (Pt 1.139 MPa/Tt 1161 K/新組成) vs 旧条件 (Pt 1.137 MPa/Tt 1058 K/旧組成) x {dry,noneq,eq}\n" + hdr,
           comments="# ", fmt="%.9g")

# --- 図 1: 軸分布 (左=全域, 右=下流ズーム) ---
fig, axs = plt.subplots(4, 2, figsize=(15, 15))
XZ = (5, 25)


def draw(ax, key, zoom=False):
    for lab, run, tag, Yw, c, ls in RUNS:
        _, q, info, cc, lls = axis[run]
        x = q["x_over_rt"]
        y = q[key]
        ax.plot(x, y, lls, color=cc, lw=1.3, label=lab)
        if key == "T" and zoom:      # 飽和温度 (液基準) を細点線で重ねる
            ax.plot(x, q["Tsat"], ":", color=cc, lw=0.9)
    ax.grid(alpha=.3)
    if zoom:
        ax.set_xlim(*XZ)
        m = (axis[RUNS[0][1]][1]["x_over_rt"] > XZ[0])
        vals = np.concatenate([axis[r[1]][1][key][m] for r in RUNS])
        lo, hi = float(vals.min()), float(vals.max())
        pad = 0.08 * (hi - lo) + 1e-12
        ax.set_ylim(lo - pad, hi + pad)
    else:
        ax.set_xlim(-1, 25)
    for xv in (axis[RUNS[0][1]][2]["x_A"], axis[RUNS[0][1]][2]["x_E"]):
        ax.axvline(xv, c="gray", ls=":", lw=.7)


for row, (key, ylab) in enumerate([("M", "軸 Mach"), ("T", "T [K]"), ("P", "P [Pa]")]):
    for col in (0, 1):
        draw(axs[row, col], key, zoom=bool(col))
        axs[row, col].set_ylabel(ylab + ("  (下流ズーム)" if col else ""))
    if key == "P":
        axs[row, 0].set_yscale("log")
axs[0, 0].axhline(4.19, c="k", lw=.6, ls=":")
axs[0, 1].axhline(4.19, c="k", lw=.6, ls=":")
axs[0, 0].legend(fontsize=9, ncol=2, loc="lower right")
axs[1, 1].set_ylabel("T [K] (下流ズーム; 細点線 = T_sat)")
# 4 段目: S (対数) と g
for col, (key, ylab, log) in enumerate([("S", "過飽和度 S (液基準)", True), ("g_0", "凝縮質量分率 g", False)]):
    ax = axs[3, col]
    for lab, run, tag, Yw, c, ls in RUNS:
        _, q, info, cc, lls = axis[run]
        m = q["x_over_rt"] > 0
        ax.plot(q["x_over_rt"][m], (q[key] / (Yw if key == "g_0" else 1.0))[m], lls, color=cc, lw=1.3, label=lab)
    ax.set_xlim(0, 25)
    ax.grid(alpha=.3)
    ax.set_ylabel(ylab + (" / Y_H2O (凝縮率)" if key == "g_0" else ""))
    if log:
        ax.set_yscale("log")
        ax.axhline(1, c="k", lw=.6, ls=":")
    for xv in (axis[RUNS[0][1]][2]["x_A"], axis[RUNS[0][1]][2]["x_E"]):
        ax.axvline(xv, c="gray", ls=":", lw=.7)
for ax in axs[-1]:
    ax.set_xlabel("x / r_t   (スロート x=0、灰点線 = x_A / x_E)")
fig.suptitle("va3 L_c8 形状 (同一ノズル) の軸中心分布 — 入口条件 新 (Pt 1.139 MPa, Tt 1161 K, 新組成) vs 旧 (Pt 1.137 MPa, Tt 1058 K, 旧組成)", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.985))
fig.savefig(FIG / "va3_cond_axis_compare.png", dpi=130)
plt.close(fig)

# --- 図 2: 出口断面プロファイル ---
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
for lab, run, tag, Yw, c, ls in RUNS:
    e = exit_plane(run)
    axs[0].plot(e["M"], e["r_rt"], ls, color=c, lw=1.3, label=f"{lab} (avg {e['M_avg']:.3f})")
    axs[1].plot(e["T"], e["r_rt"], ls, color=c, lw=1.3)
    axs[2].plot(e["g"] / Yw * 100, e["r_rt"], ls, color=c, lw=1.3)
for ax, xl in zip(axs, ("Mach", "T [K]", "凝縮率 g / Y_H2O [%]")):
    ax.set_xlabel(xl)
    ax.set_ylabel("r / r_t")
    ax.grid(alpha=.3)
axs[0].legend(fontsize=8)
fig.suptitle("物理出口断面 (出口 BC の 2 r_t 上流) の半径方向プロファイル", fontsize=11)
fig.tight_layout()
fig.savefig(FIG / "va3_cond_exit_profile.png", dpi=130)
print("wrote figs/va3_cond_axis_compare.png, figs/va3_cond_exit_profile.png, axis_compare_va3_vs_va3old.csv, study_va3_conditions.json")
