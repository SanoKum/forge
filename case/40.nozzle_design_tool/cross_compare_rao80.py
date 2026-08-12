"""rao80 内部衝撃波の forge vs SU2 / RANS vs Euler 2x2 比較 + 軸プロファイル + 筋トレース。"""
import sys, struct, json
sys.path.insert(0, "/home/sano/work/forge/design")
from pathlib import Path
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

BASE = Path("/home/sano/work/forge/case/40.nozzle_design_tool")
W = NozzleWall(WallParams(theta_a=np.deg2rad(33.77), theta_e=np.deg2rad(12.13),
                          L_div=5.9713, bell_type="top"))
G, RGAS = 1.4, 287.058


def load_forge(run_dir):
    rd = BASE / run_dir
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    with h5py.File(rd / "nozzle.h5") as nz:
        cc = nz["CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res) as f:
        Ux, Uy, son = f["VALUE/Ux"][:], f["VALUE/Uy"][:], f["VALUE/sonic"][:]
    return cc[:, 0] * 1e3, cc[:, 1] * 1e3, np.hypot(Ux, Uy) / np.maximum(son, 1e-9)


def load_su2(run_dir):
    raw = open(BASE / run_dir / "restart_flow.dat", "rb").read()
    magic, nf, npts, e1, e2 = struct.unpack("iiiii", raw[:20])
    off = 20 + 33 * nf
    d = np.frombuffer(raw[off:off + 8 * nf * npts], np.float64).reshape(npts, nf)
    x, y, rho, mx, my, E = (d[:, i] for i in range(6))
    u, v = mx / rho, my / rho
    p = (G - 1) * (E - 0.5 * rho * (u * u + v * v))
    T = p / (rho * RGAS)
    M = np.hypot(u, v) / np.sqrt(G * RGAS * np.maximum(T, 1.0))
    return x * 1e3, y * 1e3, M


def axis_profile(x, r, M):
    m = r < 0.9
    o = np.argsort(x[m])
    return x[m][o], M[m][o]


def jump_trace(x, r, M):
    out = []
    for xt in np.linspace(8, 59, 35):
        mm = np.abs(x - xt) < 0.7
        rr, mv = r[mm], M[mm]
        rwall = float(W.r(np.clip(xt / 10, W.x_in, W.x_e)) * 10)
        keep = (rr < 0.92 * rwall) & (rr > 2.0)
        rr, mv = rr[keep], mv[keep]
        if len(rr) < 6:
            continue
        o = np.argsort(rr)
        rr, mv = rr[o], mv[o]
        dM = np.diff(mv)
        i = int(np.argmin(dM))
        if dM[i] < -0.04:
            out.append((xt, rr[i], dM[i]))
    return np.array(out)


CASES = [
    ("forge RANS-SST (node)", "run_0075_rao_check/rao80_p1", load_forge, "#2a78d6"),
    ("SU2 RANS-SST (同一メッシュ)", "run_0078_su2_rao80", load_su2, "#eb6834"),
    ("forge Euler (cell, slip)", "run_0085_forge_euler_cell", load_forge, "#1baf7a"),
    ("SU2 Euler (同一メッシュ)", "run_0086_su2_euler", load_su2, "#eda100"),
]

fields = [(lab, fn(rd), c) for lab, rd, fn, c in CASES]

# --- 2x2 コンタ ---------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.0))
vmax = 5.1
for ax, (lab, (x, r, M), c) in zip(axes.ravel(), fields):
    xg = np.linspace(W.x_in * 10, W.x_e * 10, 340)
    rg = np.linspace(0, 30, 210)
    XG, RG = np.meshgrid(xg, rg)
    MG = griddata((x, r), M, (XG, RG), method="linear")
    MG = np.ma.masked_where(RG > W.r(np.clip(XG / 10, W.x_in, W.x_e)) * 10, MG)
    cf = ax.contourf(XG, RG, MG, levels=np.linspace(0, vmax, 41), cmap="turbo",
                     extend="max")
    ax.contour(XG, RG, MG, levels=[1.0], colors="white", linewidths=1.2)
    tr = jump_trace(x, r, M)
    if len(tr):
        ax.plot(tr[:, 0], tr[:, 1], ".", color="white", ms=3.5, alpha=0.85)
    xs = np.linspace(W.x_in, W.x_e, 300)
    ax.plot(xs * 10, W.r(xs) * 10, "k-", lw=1.0)
    ax.set_title(lab, fontsize=9.5, color=c)
    ax.set_xlim(-42, 62)
    ax.set_ylim(0, 31)
    ax.set_aspect("equal")
    ax.tick_params(colors="#555", labelsize=8)
fig.suptitle("rao80 (θn33.8° θe12.1° L5.97, ε=9) — forge vs SU2 / RANS vs Euler の内部衝撃波 2×2\n"
             "(白点 = 各 x 断面の最大 Mach 降下位置 = 内部衝撃波トレース。4 者とも同一軌跡)",
             fontsize=10.5, color="#222")
fig.subplots_adjust(right=0.90, hspace=0.28, wspace=0.10, top=0.88)
cax = fig.add_axes([0.92, 0.15, 0.016, 0.65])
fig.colorbar(cf, cax=cax, label="Mach")
out = BASE / "run_0078_su2_rao80/cross_compare_2x2.png"
fig.savefig(out, dpi=150, facecolor="white")
print(out)

# --- 軸プロファイル + 出口プロファイル ------------------------------------------
fig2, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.6))
for lab, (x, r, M), c in fields:
    xs, Ms = axis_profile(x, r, M)
    a1.plot(xs, Ms, color=c, lw=1.7, label=lab)
    me = x > (W.x_e * 10 - 1.2)
    o = np.argsort(r[me])
    a2.plot(r[me][o], M[me][o], color=c, lw=1.7)
a1.axhline(3.806, color="#888", lw=0.9, ls="--")
a1.annotate("1D 等エントロピー Me=3.81", xy=(-38, 3.86), fontsize=8.5, color="#666")
a1.set_xlabel("x [mm]", fontsize=9.5)
a1.set_ylabel("軸上 Mach", fontsize=9.5)
a1.set_title("軸プロファイル (kernel 過膨張のコード間一致)", fontsize=10)
a1.legend(fontsize=8, frameon=False, loc="upper left")
a2.axvline(0, color="#999", lw=0.5)
a2.set_xlabel("r [mm] (出口面)", fontsize=9.5)
a2.set_ylabel("Mach", fontsize=9.5)
a2.set_title("出口半径プロファイル (二層構造 = 内部衝撃波の内外)", fontsize=10)
for a in (a1, a2):
    a.grid(True, color="#e6e6e3", lw=0.6)
    a.set_axisbelow(True)
    a.tick_params(colors="#555", labelsize=8.5)
    for s in a.spines.values():
        s.set_color("#cccccc")
fig2.tight_layout()
out2 = BASE / "run_0078_su2_rao80/cross_compare_profiles.png"
fig2.savefig(out2, dpi=150, facecolor="white")
print(out2)

# --- サマリ表 -----------------------------------------------------------------
print("\n=== サマリ ===")
for lab, (x, r, M), c in fields:
    xs, Ms = axis_profile(x, r, M)
    tr = jump_trace(x, r, M)
    print(f"{lab:34s} axisM_exit={Ms[-1]:.3f}  jump@exit_dM={tr[-1,2] if len(tr) else 0:+.3f} "
          f"r_jump@exit={tr[-1,1] if len(tr) else 0:.1f}mm")
