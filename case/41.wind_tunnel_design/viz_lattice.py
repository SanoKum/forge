"""逆 MOC の「実際の計算格子」可視化: 単位セル (C+片, C-片) を親子リンク付きで描く。"""
import sys
sys.path.insert(0, "design")
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
matplotlib.rcParams["font.family"] = "Noto Sans CJK JP"
import matplotlib.pyplot as plt

from forge_design.geometry.transonic import SauerThroat
from forge_design.geometry.bezier import MachBezier
from forge_design.geometry.moc_inverse import InverseMOC
from forge_design.geometry.moc_kernel import _Pt, pm_nu

G, Md, R = 1.4, 4.0, 2.0
st = SauerThroat(R=R, gamma=G)
x0, rr, MM, tt = st.starting_line(M_start=1.05, n=10)   # 粗くして見えるように
h = 1e-4
M0 = float(st.mach(x0, 0.0))
dM0 = (st.mach(x0 + h, 0.0) - st.mach(x0 - h, 0.0)) / (2 * h)
xd = 14.0
bz = MachBezier.from_constraints(x0, xd, start=(M0, dM0), free_cp=[2.0, 2.9, 3.6],
                                 end=(Md, 0.0, 0.0))
target = lambda x: Md if x >= xd else float(np.atleast_1d(bz(float(x)))[0])  # noqa

inv = InverseMOC(gamma=G, delta=1.0)
# 初期データ: 軸 (x0→6, 22点) + starting line (10点) — 粗解像のデモ
ax_pts = np.linspace(x0, 6.0, 22)[1:][::-1]
init = [_Pt(float(x), 0.0, 0.0, float(pm_nu(float(target(x)), G)), G) for x in ax_pts]
init += [_Pt(float(x0), float(rr[i]), float(tt[i]), float(pm_nu(float(MM[i]), G)), G)
         for i in range(10)]
init[-1].th = float(np.arcsin(x0 / R))

# fill と同じロジック + 親リンク記録
edges_cp, edges_cm = [], []   # (A→P) = C+片, (B→P) = C-片
front = list(init)
pts_all = list(init)
while len(front) >= 2:
    new = []
    for i in range(len(front) - 1):
        P = None
        for A, B in ((front[i + 1], front[i]), (front[i], front[i + 1])):
            try:
                cand = inv._interior(A, B)
            except RuntimeError:
                continue
            if (np.isfinite(cand.x) and np.isfinite(cand.r) and cand.r > -1e-12
                    and cand.x >= min(A.x, B.x) - 1e-9):
                P = cand
                edges_cp.append(((A.x, A.r), (P.x, P.r)))
                edges_cm.append(((B.x, B.r), (P.x, P.r)))
                break
        if P is not None:
            new.append(P)
    front = new
    pts_all.extend(new)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.4),
                             gridspec_kw={"width_ratios": [1.6, 1.0]})
for a in (a1, a2):
    for (p0, p1) in edges_cp:
        a.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#2a78d6", lw=0.6, alpha=0.7)
    for (p0, p1) in edges_cm:
        a.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#eb6834", lw=0.6, alpha=0.7)
    xs = np.array([p.x for p in pts_all])
    rs = np.array([p.r for p in pts_all])
    a.plot(xs, rs, ".", color="#333", ms=2.0, zorder=5)
    # 初期データ線
    a.plot([x0, x0], [0, float(rr[-1])], color="#008300", lw=3.0, zorder=6)
    a.plot([x0, 6.0], [0, 0], color="#7a3fbf", lw=3.0, zorder=6)
    a.axhline(0, color="#888", lw=0.5, ls="-.")
    a.set_aspect("equal")
    a.grid(True, color="#eeeeec", lw=0.5)
    a.set_axisbelow(True)
    a.tick_params(colors="#555", labelsize=8.5)
    for s in a.spines.values():
        s.set_color("#cccccc")
a1.plot([], [], color="#008300", lw=3, label="starting line (Sauer 遷音速解: M,θ 既知)")
a1.plot([], [], color="#7a3fbf", lw=3, label="軸 (r=0: 目標 M(x) と θ=0 を指定)")
a1.plot([], [], color="#2a78d6", lw=1.2, label="C+ 小片 (下側の既知点 → 新点)")
a1.plot([], [], color="#eb6834", lw=1.2, label="C- 小片 (もう一方の既知点 → 新点)")
a1.plot([], [], ".", color="#333", ms=5, label="計算された点 (M,θ が決まった場所)")
a1.legend(fontsize=8, frameon=False, loc="upper left")
a1.set_xlim(-0.1, 6.2); a1.set_ylim(-0.05, 2.6)
a1.set_xlabel("x / rt", fontsize=9.5); a1.set_ylabel("r / rt", fontsize=9.5)
a1.set_title("逆 MOC の実計算格子 (粗解像デモ): 2 本のデータ線から編み上がる", fontsize=10)
a2.set_xlim(0.1, 1.3); a2.set_ylim(-0.02, 1.15)
a2.set_xlabel("x / rt", fontsize=9.5)
a2.set_title("スロート近傍の拡大 (楔領域)", fontsize=10)
fig.tight_layout()
out = "case/41.wind_tunnel_design/moc_lattice_demo.png"
fig.savefig(out, dpi=150, facecolor="white")
print(out)
