"""逆 MOC の用語を図の上で名指しする説明図。

A: 用語の地図 (壁足 / 上端の縁 = C+ と C− / 三角形が閉じる場所 / 壁の流線)
B: 「単位セル 1 回」の拡大
C: 「その C+ 線 1 本の上の点だけ」= 流束閉包が要る情報のすべて

格子は production `InverseMOC.fill` と同じ役割固定ロジックを親リンク付きで
再現した粗解像デモ。"""
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
from forge_design.geometry.moc_inverse import InverseMOC, _mass_flux_density
from forge_design.geometry.moc_kernel import _Pt, pm_nu, pm_mach

BASE = "case/41.wind_tunnel_design/"
G, Md, R = 1.4, 4.0, 2.0
NS, XEND, NAX = 12, 6.0, 26
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

inv = InverseMOC(gamma=G, delta=1.0)
ax_x = np.linspace(x0, XEND, NAX)[1:][::-1]
init = [_Pt(float(x), 0.0, 0.0, float(pm_nu(float(target(x)), G)), G) for x in ax_x]
init += [_Pt(float(x0), float(rr[i]), float(tt[i]), float(pm_nu(float(MM[i]), G)), G)
         for i in range(NS)]
init[-1].th = float(np.arcsin(x0 / R))
r_foot = float(rr[-1])

pts, pa, pb = list(init), [None] * len(init), [None] * len(init)
front = list(range(len(init)))
while len(front) >= 2:
    new = []
    for i in range(len(front) - 1):
        ib, ia = front[i], front[i + 1]      # B=右(下流)=C− 担体, A=左=C+ 担体
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
X = np.array([p.x for p in pts])
Y = np.array([p.r for p in pts])
wall = inv.wall_streamline(pts, x0, r_foot, dx=0.05)
mstar = 0.0
sl = init[len(ax_x):]
for P0, P1 in zip(sl[:-1], sl[1:]):
    thm = 0.5 * (P0.th + P1.th)
    mstar += (float(_mass_flux_density(0.5 * (P0.M + P1.M), G))
              * (np.cos(thm) * (P1.r - P0.r) - np.sin(thm) * (P1.x - P0.x))
              * 2.0 * np.pi * 0.5 * (P0.r + P1.r))

_xy = np.column_stack([X, Y])
_th, _nu = np.array([p.th for p in pts]), np.array([p.nu for p in pts])
itp_th, itp_nu = LinearNDInterpolator(_xy, _th), LinearNDInterpolator(_xy, _nu)
nea_th, nea_nu = NearestNDInterpolator(_xy, _th), NearestNDInterpolator(_xy, _nu)


def state(x, r):
    """線形補間、hull 外/縁で nan なら最近傍で埋める (縁を描くための保険)。"""
    t, v = float(itp_th(x, r)), float(itp_nu(x, r))
    if not np.isfinite(t):
        t = float(nea_th(x, r))
    if not np.isfinite(v):
        v = float(nea_nu(x, r))
    return t, pm_mach(max(v, 1e-9), G)


def trace(x, r, fam, ds, nmax=1200):
    """fam=+1: C+, fam=−1: C−。ds<0 で上流へ。"""
    out = [(x, r)]
    for _ in range(nmax):
        t, M = state(x, r)
        mu = np.arcsin(1.0 / max(M, 1.0 + 1e-9))
        x, r = x + ds, r + ds * np.tan(t + fam * mu)
        if r < 0 or r > 2.6 or not (x0 - 0.02 <= x <= XEND + 0.02):
            break
        out.append((x, r))
    return np.array(out)


e_cp = trace(x0, r_foot, +1, 0.02)                 # 壁足から出る C+ (上端の縁 左)
e_cm = trace(XEND, 0.0, -1, -0.02)[::-1]           # 最下流の軸点へ降りる C− (上端の縁 右)
# 交点 (三角形の頂点) で両方を切る
d = np.interp(e_cp[:, 0], e_cm[:, 0], e_cm[:, 1], left=np.inf, right=-np.inf) - e_cp[:, 1]
i_x = int(np.argmax(d < 0))
apex = (float(e_cp[i_x, 0]), float(e_cp[i_x, 1]))
e_cp = e_cp[:i_x + 1]
e_cm = e_cm[e_cm[:, 0] >= apex[0]]
print("apex", apex)

# C+ 親チェーン: 流束が mdot* に届く 1 本を選ぶ
best = None
for k in range(len(init) - NS):
    ch = [k]
    while True:
        kid = [i for i in range(len(init), len(pts)) if pa[i] == ch[-1]]
        if not kid:
            break
        ch.append(kid[0])
    acc, hit = 0.0, None
    for i0, i1 in zip(ch[:-1], ch[1:]):
        P0, P1 = pts[i0], pts[i1]
        thm = 0.5 * (P0.th + P1.th)
        d = (float(_mass_flux_density(0.5 * (P0.M + P1.M), G))
             * (np.cos(thm) * (P1.r - P0.r) - np.sin(thm) * (P1.x - P0.x))
             * 2.0 * np.pi * 0.5 * (P0.r + P1.r))
        if acc + d >= mstar and hit is None:
            s = (mstar - acc) / max(d, 1e-30)
            hit = (P0.x + s * (P1.x - P0.x), P0.r + s * (P1.r - P0.r))
        acc += d
    if hit is not None:
        best = (ch, hit)
        break
chain, hit = best
print("chain from x=%.2f len=%d closure=(%.2f, %.2f) wall_r=%.2f"
      % (X[chain[0]], len(chain), hit[0], hit[1],
         float(np.interp(hit[0], wall[:, 0], wall[:, 1]))))

fig = plt.figure(figsize=(14.0, 9.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.16, 0.84], hspace=0.30, wspace=0.16)
aA, aB, aC = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])


def net(a, lw=0.5, alpha=0.5, ms=1.8):
    for i in range(len(init), len(pts)):
        a.plot([X[pa[i]], X[i]], [Y[pa[i]], Y[i]], color=C_CP, lw=lw, alpha=alpha, zorder=2)
        a.plot([X[pb[i]], X[i]], [Y[pb[i]], Y[i]], color=C_CM, lw=lw, alpha=alpha, zorder=2)
    a.plot(X, Y, ".", color="#333", ms=ms, zorder=3)


def frame(a):
    xa = np.linspace(-R * np.sin(np.deg2rad(22)), x0, 60)
    a.plot(xa, 1 + R - np.sqrt(R ** 2 - xa ** 2), color="#111", lw=3.2, zorder=6)
    a.plot([x0, x0], [0, r_foot], color=C_SL, lw=3.2, zorder=6)
    a.plot([x0, XEND], [0, 0], color=C_AX, lw=3.2, zorder=6)
    a.plot([x0], [r_foot], "o", color=C_W, ms=9, zorder=9, mec="white", mew=1.3)
    a.set_aspect("equal")
    a.grid(True, color="#efefec", lw=0.5)
    a.set_axisbelow(True)
    a.tick_params(colors="#555", labelsize=8.5)
    for s in a.spines.values():
        s.set_color("#cccccc")


BOX = dict(fc="white", ec="#c8c8c8", lw=0.7, pad=2.8, alpha=0.97)


def note(a, txt, xy, xyt, col="#555", fs=9.2, ha="left"):
    """引出線は必ずグレー (図中の線と混同しないため)。"""
    a.annotate(txt, xy=xy, xytext=xyt, fontsize=fs, ha=ha, bbox=BOX, zorder=14,
               arrowprops=dict(arrowstyle="->", color="#666", lw=1.05,
                               shrinkA=3, shrinkB=4))


# ================= A. 用語の地図 =================
net(aA)
frame(aA)
for _e, _c in ((e_cp, C_EP), (e_cm, C_EM)):
    aA.plot(_e[:, 0], _e[:, 1], color="white", lw=5.2, zorder=9, solid_capstyle="round")
    aA.plot(_e[:, 0], _e[:, 1], color=_c, lw=3.0, zorder=10, solid_capstyle="round")
aA.plot(wall[:, 0], wall[:, 1], "-", color="white", lw=4.6, zorder=10)
aA.plot(wall[:, 0], wall[:, 1], "--", color=C_W, lw=2.4, zorder=11)
note(aA, "① 壁足 — スロート円弧 (太い黒) が終わる 1 点。\n"
         "上流の壁はここまで実在し、設計壁はここから始まる",
     (x0, r_foot), (-1.28, 2.42))
note(aA, "② 上端の縁 その 1 — 壁足から出た C+ 線 (濃青)。\n"
         "この線より上は 2 本のデータ線だけでは計算が届かない",
     (float(e_cp[len(e_cp) // 2, 0]), float(e_cp[len(e_cp) // 2, 1])), (0.35, 3.02),
     col=C_EP)
note(aA, "③ 上端の縁 その 2\n一番下流の軸点 (ここでは x=6) へ\n降りてくる C− 線 (茶)",
     (4.40, float(np.interp(4.40, e_cm[:, 0], e_cm[:, 1]))), (4.15, 2.34),
     col=C_EM)
note(aA, "④ 「三角形が閉じる」のはココ = ②と③の交点。\n"
         "②・③・軸の 3 辺で囲まれた内側が計算できる範囲",
     apex, (1.95, 1.86))
note(aA, "⑤ 壁の流線 (これが答えの壁)\n壁足から dr/dx = tanθ で積分\nするだけ。傾きが θ < θ+μ なので\n"
         "必ず縁のすぐ下を這い、③に\n突き当たった所 (x≈2.4) で止まる",
     (0.95, float(np.interp(0.95, wall[:, 0], wall[:, 1]))), (-1.28, 0.05), col=C_W)
note(aA, "⑥ 「データが痩せた所」= ②のすぐ下のこの辺り。\n"
         "縁の上には点が 1 つも無いので補間が片側だけになる",
     (2.30, float(np.interp(2.30, e_cm[:, 0], e_cm[:, 1])) - 0.09), (3.45, 1.30))
aA.set_xlim(-1.30, 6.45)
aA.set_ylim(-0.06, 3.50)
aA.set_xlabel("x / rt  (スロート = 0)", fontsize=9.5)
aA.set_ylabel("r / rt", fontsize=9.5)
aA.set_title("A. 用語の地図 — 太黒=既に在る壁(円弧) / 緑=starting line / 紫=軸 / "
             "青=C+ 小片 / 橙=C− 小片 / 黒点=計算で M,θ が決まった点", fontsize=10.5)

# ================= B. 単位セル 1 回 =================
ic = next(i for i in range(len(init), len(pts)) if 2.4 < X[i] < 3.2 and 0.30 < Y[i] < 0.75)
A, B, P = pts[pa[ic]], pts[pb[ic]], pts[ic]
net(aB, lw=1.1, alpha=0.28, ms=6.0)
aB.plot([A.x, P.x], [A.r, P.r], color=C_CP, lw=4.2, zorder=9, solid_capstyle="round")
aB.plot([B.x, P.x], [B.r, P.r], color=C_CM, lw=4.2, zorder=9, solid_capstyle="round")
aB.plot([A.x, B.x], [A.r, B.r], "o", color="#111", ms=11, zorder=10, mec="white", mew=1.4)
aB.plot([P.x], [P.r], "o", color=C_W, ms=13, zorder=10, mec="white", mew=1.5)
aB.set_xlim(A.x - 0.34, P.x + 0.60)
aB.set_ylim(min(A.r, B.r) - 0.30, P.r + 0.46)
aB.set_aspect("equal")
aB.grid(True, color="#efefec", lw=0.5)
aB.set_axisbelow(True)
aB.tick_params(colors="#555", labelsize=8.5)
for s in aB.spines.values():
    s.set_color("#cccccc")
note(aB, "既知点 A (M,θ 既知)", (A.x, A.r), (A.x - 0.30, A.r + 0.26))
note(aB, "既知点 B (M,θ 既知)", (B.x, B.r), (B.x - 0.10, B.r - 0.25))
note(aB, "新しい点 P\nA から C+ を、B から C− を伸ばした交点。\n"
         "ここで P の M と θ が 1 組だけ決まる。\n"
         "この 1 回が「単位セル 1 回」で、\n格子はこの反復だけでできている",
     (P.x, P.r), (P.x + 0.05, P.r + 0.10), col=C_W, fs=9.4)
aB.set_xlabel("x / rt", fontsize=9.5)
aB.set_ylabel("r / rt", fontsize=9.5)
aB.set_title("B. 「単位セル 1 回」= 既知 2 点から新しい 1 点を作る操作 (A の拡大)",
             fontsize=10.5)

# ================= C. C+ 線 1 本 =================
net(aC, lw=0.5, alpha=0.14, ms=1.6)
frame(aC)
cx, cy = X[chain], Y[chain]
keep = cy <= hit[1] + 1e-9
aC.plot(np.append(cx[keep], hit[0]), np.append(cy[keep], hit[1]),
        color=C_EP, lw=2.4, zorder=9)
aC.plot(cx[keep], cy[keep], "o", color=C_EP, ms=5.4, zorder=10, mec="white", mew=0.9)
aC.plot([cx[0]], [cy[0]], "o", color=C_AX, ms=10, zorder=11, mec="white", mew=1.3)
aC.plot([hit[0]], [hit[1]], "*", color=C_W, ms=17, zorder=11, mec="white", mew=0.8)
aC.plot(wall[:, 0], wall[:, 1], "--", color=C_W, lw=1.8, zorder=8)
note(aC, "軸の 1 点から出発", (cx[0], cy[0]), (cx[0] - 1.85, 0.22), col=C_AX)
note(aC, "この線に乗っている ● だけで\n「線を横切る流量」を軸から順に足せる。\n"
         "周りの場は 1 点も要らない",
     (cx[len(cx) // 3], cy[len(cx) // 3]), (2.55, 0.30), col=C_EP, fs=9.3)
note(aC, "足し算がスロート流量に達した★で切る → そこが壁点\n"
         "(★が破線=流線の壁に乗っているのが両者一致の確認)",
     hit, (-1.25, 2.05), col=C_W, fs=9.3)
aC.set_xlim(-1.30, 5.20)
aC.set_ylim(-0.06, 2.60)
aC.set_xlabel("x / rt", fontsize=9.5)
aC.set_ylabel("r / rt", fontsize=9.5)
aC.set_title("C. 「その C+ 線 1 本の上の点だけ」= 流束閉包が必要とする情報のすべて",
             fontsize=10.5)

fig.savefig(BASE + "moc_glossary.png", dpi=150, facecolor="white", bbox_inches="tight")
print(BASE + "moc_glossary.png")
