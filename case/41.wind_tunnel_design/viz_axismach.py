#!/usr/bin/env python3
r"""axis-Mach チェーンの形状・格子・軸マッハ分布の可視化 (run_0047 = pass3)。

出力 (case/41.wind_tunnel_design/ 直下):
  axismach_geometry_mesh.png  — 形状 + 構造化格子 (全体 / スロート近傍拡大)
  axismach_field.png          — Mach 等値線 + 壁 + 終端特性線
  axismach_axis_mach.png      — target vs 実測の軸 Mach と誤差
  axismach_chars.png          — MOC 特性線ネット + 終端 C⁺ (npz があるときのみ)
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).parent
RUN = HERE / "run_0048_axismach_charexit"   # 終端特性線で出口を決めた最終形
RT = 0.01           # r_throat [m] — 図は r_t 単位で描く
NI, NJ = 357, 65

C_WALL, C_MESH, C_TGT, C_ACH = "#1b3a6b", "#9fb3c8", "#c2410c", "#0f766e"


def load_grid():
    with h5py.File(RUN / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3) / RT
    return nc[:, 0].reshape(NI, NJ), nc[:, 1].reshape(NI, NJ)


def load_mach(res="res_12000.h5"):
    with h5py.File(RUN / res) as f:
        M = np.hypot(f["/VALUE/Ux"][:], f["/VALUE/Uy"][:]) / np.maximum(f["/VALUE/sonic"][:], 1e-9)
    return M.reshape(NI, NJ)


def fig_geometry_mesh(X, Y):
    fig, axes = plt.subplots(2, 1, figsize=(13, 8.2),
                             gridspec_kw={"height_ratios": [1.35, 1]})
    for ax, (xlo, xhi), title, (si, sj) in (
            (axes[0], (X.min(), X.max()),
             "全体: 入口直管 → U→T 5次 Hermite → 骨接放物線 → 逆 MOC 設計壁", (4, 2)),
            (axes[1], (-4.2, 3.0),
             "スロート近傍拡大 (T=0 で曲率 1/R=0.5 を共有・下流に円弧なし)", (1, 1))):
        m = (X[:, 0] >= xlo - 0.5) & (X[:, 0] <= xhi + 0.5)
        Xs, Ys = X[m][::si], Y[m][::si]
        for i in range(Xs.shape[0]):
            ax.plot(Xs[i], Ys[i], lw=0.3, color=C_MESH, zorder=1)
        for j in range(0, NJ, sj if sj > 1 else 2):
            ax.plot(X[m, j], Y[m, j], lw=0.3, color=C_MESH, zorder=1)
        ax.plot(X[m, -1], Y[m, -1], lw=2.0, color=C_WALL, zorder=3, label="壁")
        ax.plot(X[m, 0], Y[m, 0], lw=1.2, color="#334155", ls="--", zorder=3, label="軸")
        ax.axvline(0.0, color="#b91c1c", lw=1.0, ls=":", zorder=4)
        ax.text(0.05, 0.93, "T (スロート)", transform=ax.transAxes, fontsize=9, color="#b91c1c")
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(0, None)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("$x/r_t$")
        ax.set_ylabel("$r/r_t$")
        ax.grid(alpha=0.15)
    axes[0].legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "axismach_geometry_mesh.png", dpi=125)
    print("saved axismach_geometry_mesh.png")


def fig_field(X, Y, M, npz=None):
    fig, ax = plt.subplots(figsize=(13, 4.4))
    cf = ax.contourf(X, Y, M, levels=np.linspace(0.2, 4.05, 40), cmap="turbo", extend="both")
    ax.contour(X, Y, M, levels=np.arange(0.5, 4.01, 0.5), colors="k", linewidths=0.35, alpha=0.5)
    ax.plot(X[:, -1], Y[:, -1], lw=2.0, color="w")
    if npz is not None:
        tp = npz["term_path"]
        ax.plot(tp[:, 0], tp[:, 1], lw=2.0, color="w", ls="--",
                label=f"終端 C⁺ (E→F): $x_E$={float(npz['x_E']):.2f} → $x_F$={float(npz['x_F_char']):.2f}")
        ax.plot([float(npz["x_E"])], [0.0], "wo", ms=6)
        ax.legend(loc="upper left", fontsize=9, framealpha=0.85)
    ax.set_aspect("equal")
    ax.set_xlabel("$x/r_t$")
    ax.set_ylabel("$r/r_t$")
    ax.set_title("run_0048 の Mach 場 (出口 = 終端特性線) — node Euler, 12000 step", fontsize=11)
    fig.colorbar(cf, ax=ax, label="M", pad=0.01, fraction=0.03)
    fig.tight_layout()
    fig.savefig(HERE / "axismach_field.png", dpi=125)
    print("saved axismach_field.png")


def fig_axis_mach(X, Y, M, npz=None):
    x_ax, M_ax = X[:, 0], M[:, 0]
    tgt = np.loadtxt(RUN / "target_axis_M.csv", delimiter=",", skiprows=1)
    tx, tM = tgt[:, 0] / RT, tgt[:, 1]
    ok = np.isfinite(tM)
    import json
    info = json.loads((RUN / "prepare_info.json").read_text())
    xA, xE, Md = info["x_A"], info["x_E"], info["Md"]

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                             gridspec_kw={"height_ratios": [2.1, 1]})
    ax = axes[0]
    ax.plot(x_ax, M_ax, lw=1.8, color=C_ACH, label="実測 (CFD, 軸ノード)")
    ax.plot(tx[ok], tM[ok], lw=2.0, color=C_TGT, ls="--",
            label="target (5次 Hermite 軸 Mach law)")
    ax.axhline(Md, color="#64748b", lw=0.8, ls=":")
    for xv, lab, c in ((xA, "$x_A$ (アンカー = $x_{reach,CFD}$)", "#7c3aed"),
                       (xE, "$x_E$ ($M=M_d$ 到達)", "#b91c1c")):
        ax.axvline(xv, color=c, lw=1.0, ls="-.")
        ax.text(xv + 0.15, 1.35, lab, fontsize=9, color=c, rotation=90, va="bottom")
    if npz is not None:
        xF = float(npz["x_F_char"])
        ax.axvline(xF, color="#0369a1", lw=1.0, ls="-.")
        ax.text(xF + 0.15, 1.35, "$x_F$ (終端特性線)", fontsize=9, color="#0369a1",
                rotation=90, va="bottom")
    ax.set_ylabel("軸上 M")
    ax.set_ylim(0.8, 4.3)
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(alpha=0.2)
    ax.set_title("軸中心 Mach 分布: 指定 (target) vs 実現 (CFD) — run_0048", fontsize=12)

    ax2 = axes[1]
    m = (x_ax >= tx[ok].min()) & (x_ax <= tx[ok].max())
    dM = M_ax[m] - np.interp(x_ax[m], tx[ok], tM[ok])
    ax2.plot(x_ax[m], dM, lw=1.6, color="#334155")
    ax2.axhline(0, color="#94a3b8", lw=0.8)
    for s in (+1, -1):
        ax2.axhline(s * 0.005 * Md, color="#c2410c", lw=0.9, ls="--")
    ax2.text(x_ax[m][-1], 0.005 * Md, " ゲート ±0.5% $M_d$", fontsize=9,
             color="#c2410c", va="bottom", ha="right")
    ax2.set_xlabel("$x/r_t$")
    ax2.set_ylabel("$\\Delta M$ = 実測 − target")
    ax2.grid(alpha=0.2)
    ax2.set_ylim(-0.03, 0.03)
    fig.tight_layout()
    fig.savefig(HERE / "axismach_axis_mach.png", dpi=125)
    print(f"saved axismach_axis_mach.png  (max|dM| = {np.abs(dM).max():.4f})")


def fig_chars(npz, X, Y):
    pts, wall = npz["pts"], npz["wall"]
    i_lip = int(npz["i_lip"])
    xF_c = float(npz["x_F_char"])
    xE = float(npz["x_E"])
    fig, ax = plt.subplots(figsize=(13, 5.4))
    m = pts[:, 0] <= 27.0
    sub = pts[m][::14]
    sc = ax.scatter(sub[:, 0], sub[:, 1], c=sub[:, 3], s=0.5, cmap="turbo",
                    vmin=1.0, vmax=4.0, linewidths=0, alpha=0.55)
    mw = wall[:, 0] <= 27.0
    ax.plot(wall[mw, 0], wall[mw, 1], lw=1.2, color="#64748b", ls=":",
            label="壁流線 (未カット — F より下流は円筒)")
    ax.plot([wall[i_lip, 0]], [wall[i_lip, 1]], "s", color=C_WALL, ms=9, zorder=6,
            label=f"[旧] 壁角しきい値 0.05° で切断  $x_F$={wall[i_lip,0]:.2f}")
    ax.plot([wall[i_lip, 0], wall[i_lip, 0]], [0, wall[i_lip, 1]], lw=1.0,
            color=C_WALL, ls="--", alpha=0.7)
    mc = wall[:, 0] <= xF_c
    ax.plot(wall[mc, 0], wall[mc, 1], lw=2.6, color="#b91c1c",
            label=f"[新] 終端特性線まで  $x_F$={xF_c:.2f}  ($\\theta_w$=8e-11°)")
    tp = npz["term_path"]
    ax.plot(tp[:, 0], tp[:, 1], lw=2.2, color="#b91c1c", ls="--",
            label="終端 C+ 特性線 (E→F) = 一様域の上流境界")
    ax.plot([xE], [0.0], "o", color="#b91c1c", ms=8, zorder=5)
    ax.annotate("E ($M=M_d$ 到達)", (xE, 0.0), textcoords="offset points",
                xytext=(-18, 14), color="#b91c1c", fontsize=10, fontweight="bold")
    ax.plot([xF_c], [float(np.interp(xF_c, wall[:, 0], wall[:, 1]))], "o",
            color="#b91c1c", ms=8, zorder=5)
    ax.annotate("F (物理出口)", (xF_c, float(np.interp(xF_c, wall[:, 0], wall[:, 1]))),
                textcoords="offset points", xytext=(-30, 12), color="#b91c1c",
                fontsize=10, fontweight="bold")
    ax.annotate("", xy=(xF_c, 0.18), xytext=(wall[i_lip, 0], 0.18),
                arrowprops=dict(arrowstyle="<->", color="#0f172a", lw=1.3))
    ax.text(0.5 * (xF_c + wall[i_lip, 0]), 0.30,
            f"しきい値切断で失われていた {xF_c - wall[i_lip,0]:.2f} $r_t$",
            ha="center", fontsize=9.5, color="#0f172a")
    ax.set_xlim(-0.5, 26.5)
    ax.set_ylim(0, 3.9)
    ax.set_xlabel("$x/r_t$")
    ax.set_ylabel("$r/r_t$")
    ax.set_title("ノズル長さの決め方: 終端特性線 E→F が出口を定める (色 = 逆 MOC の M)",
                 fontsize=11.5)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.92)
    fig.colorbar(sc, ax=ax, label="M", pad=0.01, fraction=0.025)
    fig.tight_layout()
    fig.savefig(HERE / "axismach_chars.png", dpi=125)
    print("saved axismach_chars.png")


def main() -> int:
    import matplotlib.font_manager as fm
    for f in ("/home/sano/.fonts/NotoSansCJKjp-Regular.otf",
              "/home/sano/.fonts/NotoSansCJKjp-Bold.otf"):
        if Path(f).exists():
            fm.fontManager.addfont(f)
    plt.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    X, Y = load_grid()
    M = load_mach()
    npz = None
    f = HERE / "axismach_field.npz"
    if f.exists():
        npz = np.load(f)
    fig_geometry_mesh(X, Y)
    fig_field(X, Y, M, npz)
    fig_axis_mach(X, Y, M, npz)
    if npz is not None:
        fig_chars(npz, X, Y)
    return 0


if __name__ == "__main__":
    sys.exit(main())
