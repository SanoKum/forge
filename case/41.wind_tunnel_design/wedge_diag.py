"""楔 (starting line 直後・壁近傍) の被覆診断 + 解像度収束。"""
import sys
sys.path.insert(0, "design")
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
matplotlib.rcParams["font.family"] = "Noto Sans CJK JP"
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

from forge_design.geometry.transonic import SauerThroat
from forge_design.geometry.bezier import MachBezier
from forge_design.geometry.moc_inverse import InverseMOC, inverse_design
from forge_design.geometry.moc_kernel import _Pt, pm_nu

G, Md, R = 1.4, 4.0, 2.0
xd = 14.0
eps = (1/Md)*((2+(G-1)*Md*Md)/(G+1))**((G+1)/(2*(G-1)))
x_end = xd + 2.3*np.sqrt(eps)/np.tan(np.arcsin(1/Md))


def design(n_start=41, n_axis=500, dx_wall=0.02):
    st = SauerThroat(R=R, gamma=G)
    x0, rr, MM, tt = st.starting_line(M_start=1.05, n=n_start)
    h = 1e-4
    M0 = float(st.mach(x0, 0.0))
    dM0 = (st.mach(x0+h, 0.) - st.mach(x0-h, 0.)) / (2*h)
    bz = MachBezier.from_constraints(x0, xd, start=(M0, dM0),
                                     free_cp=[2.0, 2.9, 3.6], end=(Md, 0., 0.))
    tgt = lambda x: Md if x >= xd else float(np.atleast_1d(bz(float(x)))[0])  # noqa
    res = inverse_design(st, tgt, x_axis_end=float(x_end), n_axis=n_axis,
                         n_start=n_start, gamma=G, dx_wall=dx_wall,
                         th_wall0=float(np.arcsin(x0/R)))
    return st, x0, res


# ---- 解像度収束 ----------------------------------------------------------------
print("解像度収束 (壁足まわりの量):")
print(f"{'n_start':>8}{'n_axis':>8}{'dx_wall':>9} | {'tanθ[1]':>9}{'tanθ[2]':>9} | "
      f"{'r(x0+0.2)':>10}{'r(x0+1)':>9}{'lip r':>8}")
base = None
for ns, na, dxw in [(41, 400, 0.02), (81, 400, 0.02),
                    (41, 800, 0.02), (41, 400, 0.01), (81, 800, 0.01)]:
    st, x0, res = design(ns, na, dxw)
    w = res["wall"]
    t1, t2 = np.tan(w[1, 2]), np.tan(w[2, 2])
    r02 = float(np.interp(x0+0.2, w[:, 0], w[:, 1]))
    r10 = float(np.interp(x0+1.0, w[:, 0], w[:, 1]))
    print(f"{ns:8d}{na:8d}{dxw:9.3f} | {t1:9.4f}{t2:9.4f} | {r02:10.5f}{r10:9.5f}{w[-1,1]:8.4f}")
    if base is None:
        base = (ns, na, dxw, res)
tan_arc = np.tan(np.arcsin(x0/R))
print(f"{'':>25} | 円弧 tanθ = {tan_arc:.4f}")

# ---- 被覆図 --------------------------------------------------------------------
st, x0, res = design(41, 500, 0.02)
pts, wall = res["pts"], res["wall"]
px = np.array([p.x for p in pts]); pr = np.array([p.r for p in pts])
tri = Delaunay(np.c_[px, pr])
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
for ax, (xl, yl, ttl) in zip(axes, [
        ((x0-0.05, x0+1.2), (0.95, 1.35), "楔 (starting line 直後・壁近傍) の被覆"),
        ((x0-0.02, x0+0.25), (1.00, 1.10), "さらに拡大 (壁足の 3 点)")]):
    ax.triplot(px, pr, tri.simplices, color="#cfcfe0", lw=0.4, zorder=1)
    ax.plot(px, pr, ".", color="#333", ms=2.5, zorder=3, label="MOC 計算点")
    ax.plot(wall[:, 0], wall[:, 1], "-", color="#eb6834", lw=2.0, zorder=5,
            label="抽出された壁流線")
    ax.plot(wall[:6, 0], wall[:6, 1], "o", color="#eb6834", ms=5, zorder=6)
    xs_arc = np.linspace(x0-0.3, x0, 60)
    ax.plot(xs_arc, 1+R-np.sqrt(R**2-xs_arc**2), "-", color="#008300", lw=2.2,
            zorder=5, label="円弧 (上流)")
    x0r = st.starting_line(M_start=1.05, n=41)[1]
    ax.plot([x0, x0], [0, x0r[-1]], "-", color="#2a78d6", lw=2.2, zorder=4,
            label="starting line")
    ax.set_xlim(*xl); ax.set_ylim(*yl)
    ax.set_xlabel("x / rt", fontsize=9.5); ax.set_ylabel("r / rt", fontsize=9.5)
    ax.set_title(ttl, fontsize=10)
    ax.grid(True, color="#eeeeec", lw=0.5); ax.set_axisbelow(True)
    ax.tick_params(colors="#555", labelsize=8.5)
    for s in ax.spines.values():
        s.set_color("#cccccc")
axes[0].legend(fontsize=8.5, frameon=False, loc="upper left")
fig.tight_layout()
out = "case/41.wind_tunnel_design/wedge_coverage.png"
fig.savefig(out, dpi=150, facecolor="white")
print("\n" + out)
