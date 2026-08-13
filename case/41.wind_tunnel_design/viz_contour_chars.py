"""① 初期 vs 最終ノズルコンタ比較、② 逆設計場の特性線網の可視化。"""
import sys
sys.path.insert(0, "design")
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
matplotlib.rcParams["font.family"] = "Noto Sans CJK JP"
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator

from forge_design.geometry.transonic import SauerThroat
from forge_design.geometry.bezier import MachBezier
from forge_design.geometry.moc_inverse import inverse_design
from forge_design.geometry.moc_kernel import pm_mach

BASE = "case/41.wind_tunnel_design/"
RT_MM = 10.0

# --- ① 初期 (v1 設計壁) vs 最終 (v2 pass13 収束壁) ------------------------------
w0 = np.loadtxt(BASE + "run_0001_v1_euler/wall_design.csv", delimiter=",", skiprows=1)
x0w, r0w = w0[:, 0] * 1e3, w0[:, 1] * 1e3  # m→mm
wf = np.loadtxt(BASE + "run_0003_v2_tr/pass_12/wall_next.csv", delimiter=",", skiprows=1)
xfw, rfw = wf[:, 0] * RT_MM, wf[:, 1] * RT_MM  # rt単位→mm

fig, (a1, a2) = plt.subplots(2, 1, figsize=(11.5, 7.2), sharex=True,
                             gridspec_kw={"height_ratios": [2.0, 1.0]})
a1.plot(x0w, r0w, color="#eda100", lw=1.8, ls="--", label="初期コンタ (v1 逆 MOC 設計)")
a1.plot(xfw, rfw, color="#2a78d6", lw=1.8, label="最終コンタ (v2 帰還 13 パス収束)")
a1.axhline(0, color="#888", lw=0.7, ls="-.")
a1.set_ylabel("r [mm]", fontsize=10)
a1.set_title("初期 vs 最終ノズルコンタ (Md=4 風洞, rt=10mm) — 差は下段", fontsize=11)
a1.legend(fontsize=9, frameon=False, loc="lower right")
a1.set_aspect("equal")
dr = np.interp(x0w, xfw, rfw) - r0w
a2.plot(x0w, dr * 1e3, color="#1baf7a", lw=1.6)
a2.axhline(0, color="#999", lw=0.6)
a2.set_xlabel("x [mm] (スロート=0)", fontsize=10)
a2.set_ylabel("Δr = 最終 − 初期 [μm]", fontsize=10)
for a in (a1, a2):
    a.grid(True, color="#e6e6e3", lw=0.6)
    a.set_axisbelow(True)
    a.tick_params(colors="#555", labelsize=8.5)
    for s in a.spines.values():
        s.set_color("#cccccc")
fig.tight_layout()
fig.savefig(BASE + "contour_initial_vs_final.png", dpi=150, facecolor="white")
print(BASE + "contour_initial_vs_final.png")

# --- ② 特性線網の可視化 (M4 設計の逆マーチ場をトレース) --------------------------
G, Md, R = 1.4, 4.0, 2.0
st = SauerThroat(R=R, gamma=G)
x0, rr, MM, tt = st.starting_line(M_start=1.05, n=41)
h = 1e-4
M0 = float(st.mach(x0, 0.0))
dM0 = (st.mach(x0 + h, 0.0) - st.mach(x0 - h, 0.0)) / (2 * h)
xd = 14.0
bz = MachBezier.from_constraints(x0, xd, start=(M0, dM0), free_cp=[2.0, 2.9, 3.6],
                                 end=(Md, 0.0, 0.0))
target = lambda x: Md if x >= xd else float(np.atleast_1d(bz(float(x)))[0])  # noqa: E731
eps = (1 / Md) * ((2 + (G - 1) * Md * Md) / (G + 1)) ** ((G + 1) / (2 * (G - 1)))
x_end = xd + 2.3 * np.sqrt(eps) / np.tan(np.arcsin(1 / Md))
res = inverse_design(st, target, x_axis_end=float(x_end), n_axis=500, n_start=41,
                     gamma=G, th_wall0=float(np.arcsin(min(x0 / R, 1.0))))
pts, wall = res["pts"], res["wall"]
xy = np.array([[p.x, p.r] for p in pts])
thv = np.array([p.th for p in pts])
nuv = np.array([p.nu for p in pts])
itp_th = LinearNDInterpolator(xy, thv)
itp_nu = LinearNDInterpolator(xy, nuv)


def _M(x, r):
    v = itp_nu(x, r)
    return pm_mach(max(float(v), 1e-9), G) if np.isfinite(v) else np.nan


def trace(x, r, fam, ds=0.04, nmax=2000):
    """fam=+1: C+ (軸→壁方向), fam=-1: C- (壁→軸方向) を dx 前進で追う。"""
    out = [(x, r)]
    for _ in range(nmax):
        t = itp_th(x, r)
        if not np.isfinite(t):
            break
        M = _M(x, r)
        if not np.isfinite(M):
            break
        mu = np.arcsin(1.0 / max(M, 1.0 + 1e-9))
        slope = np.tan(float(t) + fam * mu)
        xm, rm = x + 0.5 * ds, r + 0.5 * ds * slope
        t2, M2 = itp_th(xm, rm), _M(xm, rm)
        if not (np.isfinite(t2) and np.isfinite(M2)):
            break
        mu2 = np.arcsin(1.0 / max(M2, 1.0 + 1e-9))
        x, r = x + ds, r + ds * np.tan(float(t2) + fam * mu2)
        if r < 0:
            out.append((x, 0.0))
            break
        out.append((x, r))
        # 壁の外に出たら終了
        if r > np.interp(x, wall[:, 0], wall[:, 1]) + 1e-6:
            break
    return np.array(out)


fig2, ax = plt.subplots(figsize=(13.0, 5.0))
# C+ (軸の各点から壁へ) — 青
for xa in np.linspace(x0 + 0.3, 20.0, 16):
    tr = trace(float(xa), 1e-4, +1)
    ax.plot(tr[:, 0], tr[:, 1], color="#2a78d6", lw=0.8, alpha=0.85)
# C- (壁の各点から軸へ) — 橙
for i in np.linspace(5, len(wall) - 10, 16).astype(int):
    tr = trace(float(wall[i, 0]), float(wall[i, 1]) - 1e-4, -1)
    ax.plot(tr[:, 0], tr[:, 1], color="#eb6834", lw=0.8, alpha=0.85)
ax.plot(wall[:, 0], wall[:, 1], "k-", lw=1.8, label="逆設計壁 (流線抽出)")
ax.plot([x0, x0], [0, float(rr[-1])], color="#008300", lw=2.2, label="starting line (Sauer 遷音速解)")
ax.axhline(0, color="#888", lw=0.8, ls="-.")
ax.plot([], [], color="#2a78d6", lw=1.2, label="C+ 特性線 (軸→壁)")
ax.plot([], [], color="#eb6834", lw=1.2, label="C- 特性線 (壁→軸)")
ax.set_xlabel("x / rt (スロート=0)", fontsize=10)
ax.set_ylabel("r / rt", fontsize=10)
ax.set_title("逆 MOC 設計の特性線網 (Md=4 モード F) — 壁の波が C- に乗って軸へ届く様子",
             fontsize=11)
ax.set_xlim(-0.5, 26)
ax.set_ylim(0, 3.6)
ax.set_aspect("equal")
ax.legend(fontsize=8.5, frameon=False, loc="upper left")
ax.grid(True, color="#e6e6e3", lw=0.6)
ax.set_axisbelow(True)
for s in ax.spines.values():
    s.set_color("#cccccc")
fig2.tight_layout()
fig2.savefig(BASE + "characteristics_net.png", dpi=150, facecolor="white")
print(BASE + "characteristics_net.png")
