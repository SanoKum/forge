#!/usr/bin/env python3
"""wall-driven チェーン Phase W2 可視化: r / θ / κ / dκ/ds と境界・違反ケース。

使い方: python plot_walldriven_demo.py <outdir>
出力: walldriven_nominal.png (代表形状 R_t=2 vs 5), walldriven_limits.png
(単調条件 μ / 変曲条件 L_D の境界挙動と κ' 跳び診断)。
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.geometry.wall_walldriven import (  # noqa: E402
    ThroatExpansionPoly, UpstreamThroatPoly, WallDrivenThroatRegion)

for f in Path.home().glob(".fonts/NotoSansCJKjp-*.otf"):
    font_manager.fontManager.addfont(str(f))
plt.rcParams.update({
    "font.family": "Noto Sans CJK JP", "font.size": 10,
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": "#c9c8c4", "axes.labelcolor": "#0b0b0b",
    "text.color": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e8e7e4", "grid.linewidth": 0.8,
    "axes.axisbelow": True, "lines.linewidth": 2.0,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False,
})
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"  # categorical slots 1-3

outdir = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
outdir.mkdir(parents=True, exist_ok=True)

# --- 図 1: 代表形状 (r_U=2.5, θ_D=12°, L_D=0.9 上限) の R_t=2 vs 5 ---------------
r_U, th_D = 2.5, np.deg2rad(12.0)
walls = {}
for Rt, col in ((2.0, C1), (5.0, C2)):
    L_D = 0.9 * 3.0 * Rt * np.tan(th_D)
    walls[Rt] = (WallDrivenThroatRegion(r_U=r_U, R_t=Rt, L_U=3.5,
                                        L_D=L_D, theta_D=th_D), col)

fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)
(ax_r, ax_t), (ax_k, ax_dk) = axes
for Rt, (w, col) in walls.items():
    s = w.sample(1201)
    lab = f"$R_t={Rt:g}$"
    ax_r.plot(s["x"], s["r"], color=col, label=lab)
    ax_r.plot([w.dn.L_D], [w.dn.r_D], "o", ms=6, color=col)
    ax_t.plot(s["x"], np.rad2deg(s["theta"]), color=col, label=lab)
    ax_k.plot(s["x"], s["kappa"], color=col, label=lab)
    ax_dk.plot(s["x"], s["dkappa_ds"], color=col, label=lab)
ax_r.set_ylabel("$r/r_t$")
ax_r.set_title("壁形状 (● = D 点、r_D は従属量)", fontsize=10, loc="left")
ax_t.set_ylabel(r"$\theta$ [deg]")
ax_t.set_title(r"壁面角 (D で最大 $\theta_D=12°$, $\theta'=0$)", fontsize=10, loc="left")
ax_k.set_ylabel(r"$\kappa \, r_t$")
ax_k.set_title(r"曲率 (T で $1/R_t$ に連続、D で 0)", fontsize=10, loc="left")
ax_dk.set_ylabel(r"$d\kappa/ds \cdot r_t^2$")
ax_dk.set_title(r"曲率変化率 (T で有限の跳び = 診断値)", fontsize=10, loc="left")
for ax in axes.flat:
    ax.axvline(0.0, color="#c9c8c4", lw=1.0, ls=":")
    ax.set_xlabel("$x/r_t$ (T = 0)")
    ax.legend(loc="best")
fig.suptitle("wall-driven スロート域多項式壁 — 代表形状 "
             "($r_U=2.5$, $L_U=3.5$, $\\theta_D=12°$, $L_D=0.9\\,L_{D,max}$)",
             fontsize=11)
fig.savefig(outdir / "walldriven_nominal.png", dpi=150)
plt.close(fig)

# --- 図 2: 設計条件の境界挙動 ----------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)
(ax_a, ax_b), (ax_c, ax_d) = axes

# (a)(b) U→T 単調条件 μ ≤ 20
dr = r_U - 1.0
for mu, col in ((4.08, C1), (20.0, C3), (22.0, C2)):
    L_U = np.sqrt(mu * 2.0 * dr)
    u = UpstreamThroatPoly(r_U=r_U, R_t=2.0, L_U=L_U)
    xi = np.linspace(0, 1, 2001)
    lab = f"$\\mu={mu:g}$" + (" (違反)" if mu > 20 else " (境界)" if mu == 20 else "")
    for ax in (ax_a, ax_b):
        ax.plot(xi, u.r(-L_U + xi * L_U, 1), color=col, label=lab)
ax_a.axhline(0, color="#c9c8c4", lw=1.0)
ax_a.set_title("U→T 勾配 $dr/dx$ — 単調収縮条件 $\\mu\\leq20$", fontsize=10, loc="left")
ax_a.set_ylabel("$dr/dx$")
ax_b.set_xlim(0, 0.25)
ax_b.set_ylim(-0.02, 0.02)
ax_b.axhline(0, color="#c9c8c4", lw=1.0)
ax_b.set_title("同 U 端拡大 — $\\mu>20$ で $dr/dx>0$ (半径の逆行)", fontsize=10, loc="left")
for ax in (ax_a, ax_b):
    ax.set_xlabel(r"$\xi=(x+L_U)/L_U$")
    ax.legend(loc="best")

# (c) T→D 変曲条件 L_D ≤ 3 R_t tanθ_D
Rt = 2.0
Lmax = 3.0 * Rt * np.tan(th_D)
ax_cz = ax_c.inset_axes([0.44, 0.46, 0.52, 0.46])
for fac, col in ((0.9, C1), (1.0, C3), (1.1, C2)):
    d = ThroatExpansionPoly(R_t=Rt, L_D=fac * Lmax, theta_D=th_D)
    xi = np.linspace(0, 1, 1201)
    lab = f"$L_D={fac:g}\\,L_{{D,max}}$" + (" (違反)" if fac > 1 else "")
    ax_c.plot(xi, d.r(xi * d.L_D, 2), color=col, label=lab)
    ax_cz.plot(xi, d.r(xi * d.L_D, 2), color=col)
ax_c.axhline(0, color="#c9c8c4", lw=1.0)
ax_cz.set_xlim(0.6, 1.0)
ax_cz.set_ylim(-0.012, 0.05)
ax_cz.axhline(0, color="#c9c8c4", lw=1.0)
ax_cz.set_title("末端拡大: 違反で $r''<0$", fontsize=8)
ax_cz.tick_params(labelsize=8)
ax_c.set_title("T→D $r''$ — $L_D>3R_t\\tan\\theta_D$ で D 手前に余分な変曲点",
               fontsize=10, loc="left")
ax_c.set_xlabel(r"$\xi=x/L_D$")
ax_c.set_ylabel("$r''$")
ax_c.legend(loc="lower left", fontsize=9)

# (d) T 近傍の dκ/ds — κ' 跳び診断
for Rt2, (w, col) in walls.items():
    x = np.linspace(-0.8, 0.8, 1601)
    ax_d.plot(x, w.dkappa_ds(x), color=col,
              label=f"$R_t={Rt2:g}$ (跳び {w.kappa_prime_jump_at_throat:+.3f})")
ax_d.axvline(0.0, color="#c9c8c4", lw=1.0, ls=":")
ax_d.set_title(r"T 近傍の $d\kappa/ds$ — κ は連続、κ′ は有限跳び (常時診断)",
               fontsize=10, loc="left")
ax_d.set_xlabel("$x/r_t$")
ax_d.set_ylabel(r"$d\kappa/ds \cdot r_t^2$")
ax_d.legend(loc="best")

fig.suptitle("設計条件の境界挙動 — 検証済み必要十分条件と違反時の実挙動", fontsize=11)
fig.savefig(outdir / "walldriven_limits.png", dpi=150)
plt.close(fig)
print("saved:", outdir / "walldriven_nominal.png", outdir / "walldriven_limits.png")
