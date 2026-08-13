"""「上端の縁 その 2 (C− 側)」がどこから来るのかの説明図。

左列: 軸に注文を書く下流端を 2.5 / 4.0 / 6.0 と変えると、縁③がそのまま右へ動く。
右: 縁③の外側の点を計算するには「書いていない軸データ」が要る、という仕組み。
"""
import sys
sys.path.insert(0, "design")
import numpy as np
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
matplotlib.rcParams["font.family"] = "Noto Sans CJK JP"
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

from forge_design.geometry.transonic import SauerThroat
from forge_design.geometry.bezier import MachBezier
from forge_design.geometry.moc_inverse import InverseMOC
from forge_design.geometry.moc_kernel import _Pt, pm_nu, pm_mach

BASE = "case/41.wind_tunnel_design/"
G, Md, R = 1.4, 4.0, 2.0
NS, DX = 12, 0.2296
C_CP, C_CM, C_AX, C_SL, C_W = "#2a78d6", "#eb6834", "#7a3fbf", "#008300", "#c0392b"
C_EP, C_EM = "#0b4f9e", "#a1400f"

st = SauerThroat(R=R, gamma=G)
x0, rr, MM, tt = st.starting_line(M_start=1.05, n=NS)
h = 1e-4
M0 = float(st.mach(x0, 0.0))
dM0 = (st.mach(x0 + h, 0.0) - st.mach(x0 - h, 0.0)) / (2 * h)
bz = MachBezier.from_constraints(x0, 14.0, start=(M0, dM0), free_cp=[2.0, 2.9, 3.6],
                                 end=(Md, 0.0, 0.0))
target = lambda x: Md if x >= 14.0 else float(np.atleast_1d(bz(float(x)))[0])  # noqa: E731
r_foot = float(rr[-1])
inv = InverseMOC(gamma=G, delta=1.0)


def build(x_end):
    ax_x = np.arange(x0 + DX, x_end + 1e-9, DX)[::-1]
    init = [_Pt(float(x), 0.0, 0.0, float(pm_nu(float(target(x)), G)), G) for x in ax_x]
    init += [_Pt(float(x0), float(rr[i]), float(tt[i]), float(pm_nu(float(MM[i]), G)), G)
             for i in range(NS)]
    init[-1].th = float(np.arcsin(x0 / R))
    pts, pa, pb = list(init), [None] * len(init), [None] * len(init)
    front = list(range(len(init)))
    while len(front) >= 2:
        new = []
        for i in range(len(front) - 1):
            ib, ia = front[i], front[i + 1]
            try:
                P = inv._interior(pts[ia], pts[ib])
            except RuntimeError:
                continue
            if not (np.isfinite(P.x) and np.isfinite(P.r)) or P.r < -1e-12:
                continue
            pts.append(P)
            pa.append(ia)
            pb.append(ib)
            new.append(len(pts) - 1)
        front = new
    return init, pts, pa, pb, float(ax_x[0])


def field(pts):
    xy = np.column_stack([[p.x for p in pts], [p.r for p in pts]])
    th, nu = np.array([p.th for p in pts]), np.array([p.nu for p in pts])
    lin_t, lin_n = LinearNDInterpolator(xy, th), LinearNDInterpolator(xy, nu)
    nea_t, nea_n = NearestNDInterpolator(xy, th), NearestNDInterpolator(xy, nu)

    def st_(x, r):
        t, v = float(lin_t(x, r)), float(lin_n(x, r))
        if not np.isfinite(t):
            t = float(nea_t(x, r))
        if not np.isfinite(v):
            v = float(nea_n(x, r))
        return t, pm_mach(max(v, 1e-9), G)
    return st_


def trace(st_, x, r, fam, ds, xlo, xhi, nmax=1500):
    out = [(x, r)]
    for _ in range(nmax):
        t, M = st_(x, r)
        mu = np.arcsin(1.0 / max(M, 1.0 + 1e-9))
        x, r = x + ds, r + ds * np.tan(t + fam * mu)
        if r < 0 or r > 2.8 or not (xlo <= x <= xhi):
            out.append((x, max(r, 0.0)))
            break
        out.append((x, r))
    return np.array(out)


CASES = [2.55, 4.16, 6.00]
fig = plt.figure(figsize=(14.4, 8.6))
gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.12], hspace=0.42, wspace=0.13)
BOX = dict(fc="white", ec="#c8c8c8", lw=0.7, pad=2.6, alpha=0.97)


def note(a, txt, xy, xyt, fs=9.2, ha="left"):
    a.annotate(txt, xy=xy, xytext=xyt, fontsize=fs, ha=ha, bbox=BOX, zorder=14,
               arrowprops=dict(arrowstyle="->", color="#666", lw=1.05, shrinkA=3, shrinkB=4))


# ② の全長 (最長ケースの場で 1 回だけ引く) — 3 枚に共通で重ねる
_i, _p, _pa, _pb, _xl = build(CASES[-1])
E2_FULL = trace(field(_p), x0, r_foot, +1, 0.02, x0 - 0.02, CASES[-1] + 0.02)

keep = None
for k, xe in enumerate(CASES):
    a = fig.add_subplot(gs[k, 0])
    init, pts, pa, pb, x_last = build(xe)
    X = np.array([p.x for p in pts])
    Y = np.array([p.r for p in pts])
    st_ = field(pts)
    e_cp = trace(st_, x0, r_foot, +1, 0.02, x0 - 0.02, xe + 0.02)
    e_cm = trace(st_, x_last, 0.0, -1, -0.02, x0 - 0.02, xe + 0.02)[::-1]
    d = np.interp(e_cp[:, 0], e_cm[:, 0], e_cm[:, 1], left=np.inf, right=-np.inf) - e_cp[:, 1]
    ix = int(np.argmax(d < 0))
    e_cp, e_cm = e_cp[:ix + 1], e_cm[e_cm[:, 0] >= e_cp[ix, 0]]
    wall = inv.wall_streamline(pts, x0, r_foot, dx=0.05)
    for i in range(len(init), len(pts)):
        a.plot([X[pa[i]], X[i]], [Y[pa[i]], Y[i]], color=C_CP, lw=0.4, alpha=0.40, zorder=2)
        a.plot([X[pb[i]], X[i]], [Y[pb[i]], Y[i]], color=C_CM, lw=0.4, alpha=0.40, zorder=2)
    a.plot(X, Y, ".", color="#333", ms=1.3, zorder=3)
    xa = np.linspace(-R * np.sin(np.deg2rad(20)), x0, 40)
    a.plot(xa, 1 + R - np.sqrt(R ** 2 - xa ** 2), color="#111", lw=2.4, zorder=6)
    a.plot([x0, x0], [0, r_foot], color=C_SL, lw=2.4, zorder=6)
    a.plot([x0, x_last], [0, 0], color=C_AX, lw=3.4, zorder=6)
    a.plot([x_last], [0], "o", color=C_AX, ms=7, zorder=8, mec="white", mew=1.1)
    a.plot(E2_FULL[:, 0], E2_FULL[:, 1], ":", color=C_EP, lw=1.6, alpha=0.75, zorder=8)
    a.plot(e_cp[:, 0], e_cp[:, 1], color=C_EP, lw=2.6, zorder=9)
    a.plot(e_cm[:, 0], e_cm[:, 1], color=C_EM, lw=2.6, zorder=9)
    a.plot([e_cp[-1, 0]], [e_cp[-1, 1]], "o", color="#111", ms=6.5, zorder=11,
           mec="white", mew=1.2)
    a.plot(wall[:, 0], wall[:, 1], "--", color=C_W, lw=2.0, zorder=10)
    a.set_aspect("equal")
    a.set_xlim(-0.85, 6.35)
    a.set_ylim(-0.05, 2.25)
    a.grid(True, color="#efefec", lw=0.5)
    a.set_axisbelow(True)
    a.tick_params(colors="#555", labelsize=8)
    for s in a.spines.values():
        s.set_color("#cccccc")
    a.set_title(f"軸に注文を書いた下流端 = x {x_last:.1f}   →   壁は x {wall[-1, 0]:.1f} まで決まる",
                fontsize=10, loc="left")
    if k == 2:
        a.set_xlabel("x / rt", fontsize=9)
    if k == 0:
        note(a, "紫の太線 = 目標 M(x) を書いた区間", (0.5 * (x0 + x_last), 0.0), (1.15, 1.62),
             fs=8.8)
        note(a, "青の点線 = ②の続き。②は 3 枚とも同じ 1 本で、\n"
                "③にぶつかる所 (●) で切れる — 切れる位置だけが違う",
             (float(E2_FULL[len(E2_FULL) // 2, 0]), float(E2_FULL[len(E2_FULL) // 2, 1])),
             (3.05, 1.20), fs=8.8)
    if k == 1:
        note(a, "●= ②と③の交点。ここで場が尽きる", (float(e_cp[-1, 0]), float(e_cp[-1, 1])),
             (2.05, 1.62), fs=8.8)
    keep = (init, pts, pa, pb, x_last, st_, e_cp, e_cm, wall, X, Y)

# ============ 右: 縁③の正体 ============
init, pts, pa, pb, x_last, st_, e_cp, e_cm, wall, X, Y = keep
a = fig.add_subplot(gs[:, 1])
for i in range(len(init), len(pts)):
    a.plot([X[pa[i]], X[i]], [Y[pa[i]], Y[i]], color=C_CP, lw=0.4, alpha=0.28, zorder=2)
    a.plot([X[pb[i]], X[i]], [Y[pb[i]], Y[i]], color=C_CM, lw=0.4, alpha=0.28, zorder=2)
a.plot(X, Y, ".", color="#333", ms=1.4, zorder=3)
xa = np.linspace(-R * np.sin(np.deg2rad(20)), x0, 40)
a.plot(xa, 1 + R - np.sqrt(R ** 2 - xa ** 2), color="#111", lw=2.6, zorder=6)
a.plot([x0, x0], [0, r_foot], color=C_SL, lw=2.6, zorder=6)
a.plot([x0, x_last], [0, 0], color=C_AX, lw=3.6, zorder=6)
a.plot([x_last], [0], "o", color=C_AX, ms=8, zorder=8, mec="white", mew=1.2)
a.plot([x_last, 8.2], [0, 0], color="#bbbbbb", lw=3.6, ls=(0, (2, 2)), zorder=5)
a.plot(e_cp[:, 0], e_cp[:, 1], color=C_EP, lw=2.8, zorder=9)
a.plot(e_cm[:, 0], e_cm[:, 1], color=C_EM, lw=2.8, zorder=9)

# 縁③のすぐ内側の点 (計算できる) と すぐ外側の点 (計算できない)
for qx, ok in ((3.10, True), (3.55, False)):
    qr = float(np.interp(qx, e_cm[:, 0], e_cm[:, 1])) + (-0.16 if ok else 0.30)
    fp = trace(st_, qx, qr, +1, -0.02, -1.0, 8.5)      # C+ の足 (上流・下)
    fm = trace(st_, qx, qr, -1, +0.02, -1.0, 8.5)      # C− の足 (下流・下)
    col = "#1b7f4b" if ok else "#c0392b"
    a.plot(fp[:, 0], fp[:, 1], color=col, lw=1.9, zorder=11)
    a.plot(fm[:, 0], fm[:, 1], color=col, lw=1.9, zorder=11)
    a.plot([qx], [qr], "o", color=col, ms=9, zorder=12, mec="white", mew=1.3)
    a.plot([fp[-1, 0]], [0], "o", color=col, ms=8, zorder=12, mec="white", mew=1.2)
    a.plot([fm[-1, 0]], [0], "o" if ok else "X", color=col, ms=8 if ok else 11,
           zorder=12, mec="white", mew=1.2)
    print(("in " if ok else "out"), "Q=", round(qx, 2), round(qr, 2),
          " C+足 x=", round(float(fp[-1, 0]), 2), " C−足 x=", round(float(fm[-1, 0]), 2))

note(a, "緑の点 = 縁③の内側。\n2 本の足がどちらも「注文を書いた区間」に着地する → 計算できる",
     (3.10, float(np.interp(3.10, e_cm[:, 0], e_cm[:, 1])) - 0.16), (-0.80, 2.52))
note(a, "赤の点 = 縁③の外側。\nC− の足が x≈7 に着地する。\n"
        "そこには注文を書いていない (灰の破線) → 計算できない",
     (3.55, float(np.interp(3.55, e_cm[:, 0], e_cm[:, 1])) + 0.30), (3.95, 2.22))
note(a, "縁③ = 「C− の足がちょうど最後の軸点に着地する」点をつないだ線。\n"
        "だから③は C− そのものであって、②を伸ばしたものではない",
     (float(e_cm[len(e_cm) // 2, 0]), float(e_cm[len(e_cm) // 2, 1])), (0.35, 1.48))
a.set_aspect("equal")
a.set_xlim(-0.85, 8.30)
a.set_ylim(-0.05, 3.15)
a.grid(True, color="#efefec", lw=0.5)
a.set_axisbelow(True)
a.tick_params(colors="#555", labelsize=8.5)
for s in a.spines.values():
    s.set_color("#cccccc")
a.set_xlabel("x / rt", fontsize=9.5)
a.set_ylabel("r / rt", fontsize=9.5)
a.set_title("縁③の正体 — どの点も「上流側の C+」と「下流側の C−」の 2 本足で立っている",
            fontsize=10.5, loc="left", pad=12)

fig.suptitle("上端の縁③ はどこから来るか — ②を伸ばしたものではなく、"
             "「軸に注文を書いた下流端」が作る別の線", fontsize=12.5, y=1.005)
fig.savefig(BASE + "moc_edge3.png", dpi=150, facecolor="white", bbox_inches="tight")
print(BASE + "moc_edge3.png")
