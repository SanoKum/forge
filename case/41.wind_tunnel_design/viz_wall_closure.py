"""逆 MOC で「壁がどう決まるか」の可視化:
左 = 流管の縁 (流線) としての壁、右 = 各 C+ 線の質量流束閉包で決まる壁点。
両者が同一曲線であることを数値で確認する。
"""
import sys
sys.path.insert(0, "design")
import numpy as np
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz   # numpy<2 互換 (この環境の numpy は 1.26)
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
font_manager.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
matplotlib.rcParams["font.family"] = "Noto Sans CJK JP"
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator

from forge_design.geometry.transonic import SauerThroat
from forge_design.geometry.bezier import MachBezier
from forge_design.geometry.moc_inverse import inverse_design, _mass_flux_density
from forge_design.geometry.moc_kernel import pm_mach

BASE = "case/41.wind_tunnel_design/"
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
pts, wall, mdot_star = res["pts"], res["wall"], res["mdot_start"]
xy = np.array([[p.x, p.r] for p in pts])
itp_th = LinearNDInterpolator(xy, np.array([p.th for p in pts]))
itp_nu = LinearNDInterpolator(xy, np.array([p.nu for p in pts]))


def state(x, r):
    t, v = itp_th(x, r), itp_nu(x, r)
    if not (np.isfinite(t) and np.isfinite(v)):
        return np.nan, np.nan
    return float(t), pm_mach(max(float(v), 1e-9), G)


def stream(x, r, ds=0.03, nmax=4000):
    """流線 dr/dx = tanθ (RK2)。"""
    out = [(x, r)]
    for _ in range(nmax):
        t, _M = state(x, r)
        if not np.isfinite(t):
            break
        t2, _ = state(x + 0.5 * ds, r + 0.5 * ds * np.tan(t))
        if not np.isfinite(t2):
            break
        x, r = x + ds, r + ds * np.tan(t2)
        if x > wall[-1, 0]:
            break
        out.append((x, r))
    return np.array(out)


def cplus_flux(xa, ds=0.02, nmax=6000):
    """軸点 (xa,0) から C+ 線を追い、累積質量流束を返す。
    戻り: 線 (n,2), 累積流束 (n,), 閉包点 (x,r) or None。"""
    x, r = float(xa), 1e-5
    line, acc = [(x, r)], [0.0]
    close = None
    for _ in range(nmax):
        t, M = state(x, r)
        if not np.isfinite(t):
            break
        mu = np.arcsin(1.0 / max(M, 1.0 + 1e-9))
        xn, rn = x + ds, r + ds * np.tan(t + mu)
        t2, M2 = state(x + 0.5 * ds, r + 0.5 * ds * np.tan(t + mu))
        if np.isfinite(t2):
            mu2 = np.arcsin(1.0 / max(M2, 1.0 + 1e-9))
            xn, rn = x + ds, r + ds * np.tan(t2 + mu2)
        tn, Mn = state(xn, rn)
        if not np.isfinite(tn):
            break
        thm, Mm, rm = 0.5 * (t + tn), 0.5 * (M + Mn), 0.5 * (r + rn)
        df = (float(_mass_flux_density(Mm, G))
              * (np.cos(thm) * (rn - r) - np.sin(thm) * (xn - x)) * 2.0 * np.pi * rm)
        if acc[-1] + df >= mdot_star and close is None:
            s = (mdot_star - acc[-1]) / max(df, 1e-30)
            close = (x + s * (xn - x), r + s * (rn - r))
            line.append(close)
            acc.append(mdot_star)
            break
        acc.append(acc[-1] + df)
        line.append((xn, rn))
        x, r = xn, rn
    return np.array(line), np.array(acc), close


fig, (a1, a2) = plt.subplots(2, 1, figsize=(13.0, 6.6), sharex=True)

# ---- 上: 流管の縁としての壁 -------------------------------------------------
a1.fill_between(wall[:, 0], 0, wall[:, 1], color="#dbeafe", alpha=0.55, zorder=0)
for r_s in np.linspace(0.12, 0.95, 7) * float(rr[-1]):
    tr = stream(x0, float(r_s))
    if len(tr) > 3:
        a1.plot(tr[:, 0], tr[:, 1], color="#7fa9d8", lw=0.8, alpha=0.9)
a1.plot(wall[:, 0], wall[:, 1], "k-", lw=2.4, zorder=6,
        label="逆設計壁 = starting line 壁足から出た流線 (この流管に mdot* が全部入る)")
a1.plot([x0, x0], [0, float(rr[-1])], color="#008300", lw=2.6, zorder=6,
        label=f"starting line: ここを横切る流量が mdot* = {mdot_star:.4f}")
a1.plot([x0], [float(rr[-1])], "o", color="#c0392b", ms=7, zorder=7,
        label="壁足 (円弧の終端 = 流線の出発点。ここで既存壁と接続)")
a1.axhline(0, color="#7a3fbf", lw=2.2, ls="-")
a1.annotate("軸 r=0: 対称 → 流束ゼロ", xy=(6.0, 0.0), xytext=(6.0, 0.55),
            fontsize=9, color="#7a3fbf",
            arrowprops=dict(arrowstyle="->", color="#7a3fbf", lw=1.0))
iw = int(0.42 * len(wall))
a1.annotate("壁 = 流線 → 流束ゼロ (非粘性すべり壁の条件そのもの)",
            xy=(wall[iw, 0], wall[iw, 1]), xytext=(wall[iw, 0] - 6.0, wall[iw, 1] - 1.35),
            fontsize=9, arrowprops=dict(arrowstyle="->", color="#333", lw=1.0))
a1.set_ylabel("r / rt", fontsize=10)
a1.set_title("① 壁の定義: 「スロートを通った質量 mdot* が全部入る流管の縁」= 流線 "
             "(非粘性では任意の流線を固体壁に置き換えても流れは変わらない)", fontsize=10.5)
a1.legend(fontsize=8.5, frameon=False, loc="upper left")

# ---- 下: C+ 線の質量流束閉包 -------------------------------------------------
xas = np.linspace(x0 + 0.25, 17.0, 13)
err = []
for xa in xas:
    line, acc, close = cplus_flux(float(xa))
    frac = np.clip(acc / mdot_star, 0, 1)
    sc = a2.scatter(line[:, 0], line[:, 1], c=frac, cmap="viridis", s=3.5,
                    vmin=0, vmax=1, zorder=4)
    if close is not None:
        a2.plot([close[0]], [close[1]], "*", color="#c0392b", ms=11, zorder=7)
        rw = float(np.interp(close[0], wall[:, 0], wall[:, 1]))
        err.append(abs(close[1] - rw) / rw * 100)
a2.plot(wall[:, 0], wall[:, 1], "k-", lw=2.0, zorder=6, label="流線抽出の壁 (上段と同じ)")
a2.plot([], [], "*", color="#c0392b", ms=11,
        label=f"流束閉包点: 軸から積算した流量が mdot* に達した所 (壁との差 平均 {np.mean(err):.2f}%)")
a2.plot([x0, x0], [0, float(rr[-1])], color="#008300", lw=2.6, zorder=6)
a2.axhline(0, color="#7a3fbf", lw=2.2)
cb = fig.colorbar(sc, ax=a2, pad=0.01, fraction=0.035)
cb.set_label("軸からの累積質量流量 / mdot*", fontsize=9)
cb.ax.tick_params(labelsize=8)
a2.set_xlabel("x / rt  (スロート = 0)", fontsize=10)
a2.set_ylabel("r / rt", fontsize=10)
a2.set_title("② 壁の求め方 (古典法): 軸の各点から C+ 線を立て、線を横切る流量を軸から積算し、"
             "mdot* に達した点で切る → その点列が壁", fontsize=10.5)
a2.legend(fontsize=8.5, frameon=False, loc="upper left")

for a in (a1, a2):
    a.grid(True, color="#eeeeec", lw=0.5)
    a.set_axisbelow(True)
    a.set_aspect("equal")
    a.set_xlim(-0.3, wall[-1, 0] + 0.4)
    a.set_ylim(-0.08, wall[-1, 1] + 0.45)
    a.tick_params(colors="#555", labelsize=8.5)
    for s in a.spines.values():
        s.set_color("#cccccc")
fig.tight_layout()
out = BASE + "moc_wall_closure.png"
fig.savefig(out, dpi=150, facecolor="white")
print(out, "wall_end=", wall[-1], "mdot*=", mdot_star, "err%=", np.round(err, 3))
