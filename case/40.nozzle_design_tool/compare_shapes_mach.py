"""Rao点・同L最適・前線各Lの壁形状重ね書き + マッハコンタ多パネル比較。"""
import sys
sys.path.insert(0, "design")
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
matplotlib.rcParams["font.family"] = "Noto Sans CJK JP"
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

from forge_design.geometry.wall import NozzleWall, WallParams

BASE = "case/40.nozzle_design_tool/"
RT_MM = 10.0

# (label, (thn,the,L), run_dir, res_file, color)
CASES = [
    ("Rao 80% チャート点\n(θn33.8° θe12.1° L5.97)",
     (33.77, 12.13, 5.9713), BASE + "run_0075_rao_check/rao80_p1", "res_12000.h5", "#eda100"),
    ("同L最適 (θn25.8° θe10.0° L5.97)",
     (25.75, 10.00, 5.9713), BASE + "run_0076_fixedL597/r2_exploit_p1", "res_12000.h5", "#008300"),
    ("前線 L6.49 (θn26.0° θe4.1°)",
     (26.02, 4.15, 6.4871), BASE + "run_0074_moo_top_r1/inf_03_0_p1", "res_12000.h5", "#2a78d6"),
    ("前線 L7.56 (θn27.3° θe6.5°)",
     (27.28, 6.47, 7.5644), BASE + "run_0074_moo_top_r1/doe_009_p1", "res_12000.h5", "#e87ba4"),
    ("前線 L8.64 (θn23.9° θe4.0°)",
     (23.87, 4.00, 8.6371), BASE + "run_0074_moo_top_r1/inf_02_0_p1", "res_12000.h5", "#eb6834"),
    ("前線 L10.05 (θn21.2° θe4.0°)",
     (21.23, 4.01, 10.0527), BASE + "run_0074_moo_top_r1/inf_12_0_p1", "res_12000.h5", "#1baf7a"),
]

# --- 1) 壁形状 重ね書き -------------------------------------------------------
fig, (ax, axz) = plt.subplots(1, 2, figsize=(13.5, 5.6),
                              gridspec_kw={"width_ratios": [2.6, 1.0]})
for label, (thn, the, L), rd, rf, c in CASES:
    wp = WallParams(theta_a=np.deg2rad(thn), theta_e=np.deg2rad(the),
                    L_div=L, bell_type="top")
    w = NozzleWall(wp)
    xs = np.linspace(w.x_in, wp.L_div, 900)
    for a in (ax, axz):
        a.plot(xs * RT_MM, w.r(xs) * RT_MM, color=c, lw=1.8)
    ax.plot([], [], color=c, lw=1.8, label=label)
    ax.plot([w.x_a * RT_MM], [w.r_a * RT_MM], "o", color=c, ms=4.5,
            mec="white", mew=0.6, zorder=5)

for a in (ax, axz):
    a.axhline(0.0, color="#888", lw=0.8, ls="-.", zorder=0)
    a.set_aspect("equal")
    a.grid(True, color="#e6e6e3", lw=0.6, zorder=0)
    a.set_axisbelow(True)
    a.tick_params(colors="#555", labelsize=8.5)
    for s in a.spines.values():
        s.set_color("#cccccc")
ax.set_xlabel("x [mm] (スロート=0)", fontsize=9.5, color="#333")
ax.set_ylabel("r [mm]", fontsize=9.5, color="#333")
ax.set_title("Rao点・同L最適・前線各L の壁形状比較 (rt=10mm, ε=9)", fontsize=11, color="#222")
ax.set_xlim(-45, 108); ax.set_ylim(-2.5, 34)
ax.legend(loc="upper left", fontsize=8, frameon=False)
axz.set_title("スロート〜Pa 近傍拡大 (● = 円弧終端 Pa)", fontsize=9.5, color="#333")
axz.set_xlim(-4, 12); axz.set_ylim(9.8, 15.5)
axz.set_xlabel("x [mm]", fontsize=9.5, color="#333")
fig.tight_layout()
out1 = BASE + "run_0074_moo_top_r1/shape_compare.png"
fig.savefig(out1, dpi=150, facecolor="white")
print(out1)

# --- 2) マッハコンタ 多パネル -------------------------------------------------
fig2, axes = plt.subplots(3, 2, figsize=(13.5, 12.5), sharex=False)
axes = axes.ravel()
vmax = 0.0
fields = []
for label, (thn, the, L), rd, rf, c in CASES:
    with h5py.File(f"{rd}/nozzle.h5", "r") as nz:
        cc = nz["CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(f"{rd}/{rf}", "r") as f:
        Ux, Uy, son = f["VALUE/Ux"][:], f["VALUE/Uy"][:], f["VALUE/sonic"][:]
    M = np.hypot(Ux, Uy) / np.maximum(son, 1e-9)
    vmax = max(vmax, float(np.percentile(M, 99.5)))
    fields.append((cc, M))

for ax_i, ((label, (thn, the, L), rd, rf, c), (cc, M)) in zip(axes, zip(CASES, fields)):
    x, r = cc[:, 0] * 1e3, cc[:, 1] * 1e3  # m -> mm (無次元でなく実寸格納)
    wp = WallParams(theta_a=np.deg2rad(thn), theta_e=np.deg2rad(the),
                    L_div=L, bell_type="top")
    w = NozzleWall(wp)
    # 構造格子への線形補間 (散布データの Delaunay 縁アーティファクト回避:
    # 凸包三角形分割は開いた出口境界付近で細長三角形を作り縞状の描画ノイズを生む)
    xg = np.linspace(w.x_in * RT_MM, w.x_e * RT_MM, 320)
    rg = np.linspace(0.0, float(r.max()), 200)
    XG, RG = np.meshgrid(xg, rg)
    MG = griddata((x, r), M, (XG, RG), method="linear")
    r_wall_g = w.r(np.clip(XG / RT_MM, w.x_in, w.x_e)) * RT_MM
    MG = np.ma.masked_where(RG > r_wall_g, MG)
    cf = ax_i.contourf(XG, RG, MG, levels=np.linspace(0, vmax, 41), cmap="turbo",
                       extend="max")
    ax_i.contour(XG, RG, MG, levels=[1.0], colors="white", linewidths=1.3)
    xs = np.linspace(w.x_in, wp.L_div, 400)
    ax_i.plot(xs * RT_MM, w.r(xs) * RT_MM, color="k", lw=1.0)
    ax_i.set_title(label, fontsize=9, color=c)
    ax_i.set_xlim(-40, 105); ax_i.set_ylim(0, 33)
    ax_i.set_aspect("equal")
    ax_i.tick_params(colors="#555", labelsize=7.5)
    for s in ax_i.spines.values():
        s.set_color("#cccccc")
fig2.suptitle("マッハ数コンタ比較 (M=1 白線, 同一カラースケール 0–%.1f)" % vmax,
              fontsize=11.5, color="#222", y=0.995)
fig2.subplots_adjust(right=0.90, hspace=0.32, wspace=0.12, top=0.95)
cax = fig2.add_axes([0.92, 0.15, 0.018, 0.7])
fig2.colorbar(cf, cax=cax, label="Mach")
out2 = BASE + "run_0074_moo_top_r1/mach_contour_compare.png"
fig2.savefig(out2, dpi=150, facecolor="white")
print(out2)
