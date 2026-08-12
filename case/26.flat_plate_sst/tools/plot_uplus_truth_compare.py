#!/usr/bin/env python3
"""平板の u+ - y+ / Cf 比較 (cell 壁解像 "真値" vs cell 壁関数 vs node 壁関数 1次/2次)。

「真値」の中身を明示するための図。本 case には実験値がないので、基準は 2 つある:
  (1) forge 自身の cell 壁解像 run (y+~0.35) の分子勾配 tau_w = mu*U1/y1  ... 内部基準
  (2) Schlichting 相関 Cf = 0.0592 Re_x^-0.2                              ... 経験式
どちらも「真値」と呼びうるので両方描く。

u+ - y+ は 2 通りに正規化する:
  A) 各 run 自身の u_tau  -> プロファイル形状の collapse を見る (絶対値の差は消える)
  B) 真値 run の u_tau    -> 絶対量の差 (Cf 欠損) がそのまま見える

使い方: python3 tools/plot_uplus_truth_compare.py [x_station]
"""
import sys, os, glob
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm
for _p in ("/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf",
           os.path.expanduser("~/.fonts/NotoSansCJKjp-Regular.otf")):
    if os.path.exists(_p):
        try:
            _fm.fontManager.addfont(_p)
        except Exception:
            pass
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "Droid Sans Fallback", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MU, GAM, CP = 1.8e-5, 1.4, 1004.5
R = CP * (GAM - 1.0) / GAM
PT, TT, MACH = 100000.0, 288.15, 0.2
PS = PT * (1.0 + 0.5 * (GAM - 1.0) * MACH ** 2) ** (-GAM / (GAM - 1.0))
TS = TT * (1.0 + 0.5 * (GAM - 1.0) * MACH ** 2) ** (-1.0)
RHO_INF = PS / (R * TS)
U_INF = MACH * np.sqrt(GAM * R * TS)
KAPPA, BLOG = 0.41, 5.0

RUNS = [
    ("run_0007_slau_muscl_innersweep", "res_120000.h5", "cell  y+0.35 壁解像 (内部基準)", "k", "-",  "molec"),
    ("run_0017_ewt_yp30_wf",           None,            "cell  y+30 壁関数",              "tab:blue",   "-",  "model"),
    ("run_0038_node_yp30_1st_long",    "res_40000.h5",  "node y+30 壁関数 1次 (押し出し)", "tab:orange", "--", "model"),
    ("run_0040_node_yp30_planar_2nd_long", "res_40000.h5", "node y+30 壁関数 2次 (平面2D)", "tab:red",  "-",  "model"),
]


def reichardt(yp):
    return np.log(1.0 + KAPPA * yp) / KAPPA + 7.8 * (
        1.0 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-yp / 3.0))


def solve_utau(Ut, y, nu):
    """Reichardt 則の逆解き (ソルバが課す modeled tau_w = rho u_tau^2 と整合)。"""
    utau = np.sqrt(nu * Ut / y)
    for _ in range(80):
        utau = max(utau, 1e-12)
        f = Ut / utau - reichardt(utau * y / nu)
        d = 1e-6 * utau
        df = ((Ut / (utau + d) - reichardt((utau + d) * y / nu)) - f) / d
        if abs(df) < 1e-20:
            break
        utau -= f / df
    return max(utau, 0.0)


def load(run, fname):
    d = os.path.join(HERE, run)
    p = os.path.join(d, fname) if fname else max(
        [f for f in glob.glob(os.path.join(d, "res_*.h5"))
         if "wall" not in os.path.basename(f) and os.path.basename(f) != "res_0.h5"],
        key=os.path.getmtime)
    h = h5py.File(p, "r")
    coord = h["MESH/COORD"][:].reshape(-1, 3)
    n = h["VALUE/ro"].shape[0]
    if coord.shape[0] == n:                    # node モード: DOF = ノードそのもの
        C = coord
    else:                                      # cell モード: 要素重心
        conne = h["MESH/CONNE"][:]
        conn = conne.reshape(n, conne.size // n)[:, 1:]
        if conn.min() == 1:
            conn = conn - 1
        C = coord[conn].mean(axis=1)
    return os.path.basename(p), C, h["VALUE/Ux"][:], h["VALUE/ro"][:], h["VALUE/vis_lam"][:]


def column(C, Ux, ro, vis, xs):
    x, y = C[:, 0], C[:, 1]
    xc_all = np.unique(np.round(x[x > 1e-6], 6))
    xc = xc_all[np.argmin(np.abs(xc_all - xs))]
    col = np.where(np.abs(x - xc) < 1e-4)[0]
    col = col[np.argsort(y[col])]
    yy, uu, rr, vv = y[col], Ux[col], ro[col], vis[col]
    keep = yy > 1e-12                          # node は壁ノード (y=0, u=0) を落とす
    return xc, yy[keep], uu[keep], rr[keep], vv[keep]


def main():
    xs = float(sys.argv[1]) if len(sys.argv) > 1 else 0.6
    data = []
    for run, fname, label, color, ls, mode in RUNS:
        base, C, Ux, ro, vis = load(run, fname)
        xc, y, u, r, v = column(C, Ux, ro, vis, xs)
        nu = v[0] / r[0]
        if mode == "molec":
            tau = v[0] * u[0] / y[0]
            utau = np.sqrt(tau / r[0])
        else:
            utau = solve_utau(u[0], y[0], nu)
            tau = r[0] * utau ** 2
        cf = tau / (0.5 * RHO_INF * U_INF ** 2)
        data.append(dict(label=label, color=color, ls=ls, y=y, u=u, nu=nu,
                         utau=utau, cf=cf, yp1=utau * y[0] / nu, xc=xc, file=base))

    truth = data[0]
    print(f"x = {xs} m   (U_inf = {U_INF:.3f} m/s, rho_inf = {RHO_INF:.4f})")
    print(f"Schlichting Cf = 0.0592 Re_x^-0.2 = "
          f"{0.0592*(RHO_INF*U_INF*xs/MU)**-0.2:.4e}")
    print(f"{'run':34s} {'y+_1':>7s} {'u_tau':>8s} {'Cf':>11s} {'/内部基準':>10s} {'/Schl':>8s}")
    cf_schl = 0.0592 * (RHO_INF * U_INF * xs / MU) ** -0.2
    for d in data:
        print(f"{d['label']:34s} {d['yp1']:7.2f} {d['utau']:8.3f} {d['cf']:11.4e} "
              f"{d['cf']/truth['cf']:10.3f} {d['cf']/cf_schl:8.3f}")

    fig, axs = plt.subplots(2, 2, figsize=(13.5, 10.4))
    axes = axs.ravel()
    ypl = np.logspace(-1, 4, 300)
    for ax, key, ttl in [
            (axes[0], "own",   "A) 各 run 自身の $u_\\tau$ で正規化 — 形状の collapse\n"
                               "(ここだけ見ると揃って見えるので絶対量の判定に使わない)"),
            (axes[1], "truth", "B) 内部基準 run_0007 の $u_\\tau$ で正規化 — 絶対量\n"
                               "(★ = 各 run の壁から最初の DOF)")]:
        ax.semilogx(ypl, ypl, "0.6", lw=1, label="$u^+=y^+$")
        ax.semilogx(ypl, np.log(ypl) / KAPPA + BLOG, "0.6", ls="--", lw=1,
                    label=f"$u^+=\\ln y^+/{KAPPA}+{BLOG}$")
        for d in data:
            ut = d["utau"] if key == "own" else truth["utau"]
            yp = (d["utau"] * d["y"] / d["nu"]) if key == "own" else (truth["utau"] * d["y"] / truth["nu"])
            ax.semilogx(yp, d["u"] / ut, color=d["color"], ls=d["ls"], marker="o", ms=2.5,
                        lw=1.4, label=d["label"])
            if key == "truth":
                ax.plot(yp[0], d["u"][0] / ut, marker="*", ms=17, color=d["color"],
                        mec="k", mew=.6, zorder=5)
        ax.set_xlim(0.1, 3e4); ax.set_ylim(0, 32)
        ax.set_xlabel("$y^+$"); ax.set_ylabel("$u^+$")
        ax.set_title(ttl, fontsize=9.5); ax.grid(alpha=.3); ax.legend(fontsize=7.5, loc="upper left")

    # C) 対数則からのずれ (内部基準 u_tau 正規化) — Cf 欠損の出所を直接見る
    ax = axes[2]
    ax.axhline(0, color="0.6", ls="--", lw=1)
    for d in data:
        yp = truth["utau"] * d["y"] / truth["nu"]
        m = (yp > 20) & (yp < 2000)
        ax.semilogx(yp[m], d["u"][m] / truth["utau"] - (np.log(yp[m]) / KAPPA + BLOG),
                    color=d["color"], ls=d["ls"], marker="o", ms=3, lw=1.5, label=d["label"])
        yp1 = yp[0]
        if yp1 > 20:
            ax.plot(yp1, d["u"][0] / truth["utau"] - (np.log(yp1) / KAPPA + BLOG),
                    marker="*", ms=17, color=d["color"], mec="k", mew=.6, zorder=5)
    ax.set_xlabel("$y^+$ (内部基準 $u_\\tau$)"); ax.set_ylabel("$u^+ - u^+_{\\rm log}$")
    ax.set_title("C) 対数則からのずれ — ★ (最初の DOF) の位置が壁関数の入力\n"
                 "node は y+52 で対数則を下回る → 壁関数が低い $u_\\tau$ を推定する",
                 fontsize=9.5)
    ax.grid(alpha=.3); ax.legend(fontsize=7.5)

    ax = axes[3]
    rex = np.logspace(5.5, 6.8, 100)
    ax.loglog(rex, 0.0592 * rex ** -0.2, "0.4", ls="--", lw=1.2,
              label="Schlichting $0.0592\\,Re_x^{-0.2}$ (経験式)")
    for run, fname, label, color, ls, mode in RUNS:
        xsts, cfs = [], []
        base, C, Ux, ro, vis = load(run, fname)
        for st in (0.30, 0.60, 0.89):
            xc, y, u, r, v = column(C, Ux, ro, vis, st)
            nu = v[0] / r[0]
            if mode == "molec":
                tau = v[0] * u[0] / y[0]
            else:
                tau = r[0] * solve_utau(u[0], y[0], nu) ** 2
            xsts.append(RHO_INF * U_INF * xc / MU)
            cfs.append(tau / (0.5 * RHO_INF * U_INF ** 2))
        ax.loglog(xsts, cfs, color=color, ls=ls, marker="s", ms=5, lw=1.5, label=label)
    ax.set_xlabel("$Re_x$"); ax.set_ylabel("$C_f$")
    ax.set_title("D) $C_f$ vs $Re_x$ — 基準は 2 つある", fontsize=10)
    ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7.5)

    fig.suptitle(f"case/26 平板 x={xs} m — 「真値」の中身と node 壁関数の Cf 欠損", fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "uplus_truth_compare.png")
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
