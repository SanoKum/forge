#!/usr/bin/env python3
r"""Rao 推力最適化放物線 (TOP) ベルノズルと等軸長コニカルノズルの輪郭生成.

RRS "Making Correct Parabolic Nozzles" (2023-01-28) の作図法を実装する。
要点は「放物線は *canted* (傾けた / 回転した) でなければならない」こと。
軸並行の放物線は 6 拘束に対し自由度 5 で過拘束になり破綻するが、
**2 次 Bézier 曲線**は始点・終点・両端接線の 4 条件 + 制御点で
ちょうど「傾いた放物線」を表し、過拘束を解消する (RRS 式 7.1-7.2 の
回転放物線と数学的に等価)。本スクリプトは G.V.R. Rao(1960) の
近似輪郭を Bézier で構成する標準手法を用いる。

幾何構成 (上流→下流):
  1. chamber       : 半径 Rc の直管 (長さ Lc)
  2. conv cone     : 収縮半角 beta の直線収縮部
  3. throat up-arc : 半径 1.5*Rt の円弧 (収縮部とスロートに接する)
  4. throat dn-arc : 半径 0.382*Rt の円弧 (スロート→放物線始点 N, 角 theta_n)
  5. parabola      : N→E の 2 次 Bézier (始点接線 theta_n, 終点接線 theta_e)

コニカルは 1-4 を共有し、5 を「N から x=Ln まで直線」に置換 (等軸長)。
既定では出口半径 Re も一致させる (等膨張比) ので、半角は theta_n/theta_e
ではなく N→E を結ぶ直線で決まる (Re を一致させる純粋な形状比較)。

出力 (mm 単位, x: 軸方向, y: 半径方向):
  contour_bell.csv      (x_mm, y_mm, zone)
  contour_conical.csv   (x_mm, y_mm, zone)
  nozzle_params.txt     主要諸元
  nozzle_preview.png    ベル/コニカル重ね描き
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

DEG = math.pi / 180.0


# ---------------------------------------------------------------------------
# Rao TOP チャート (theta_n, theta_e vs 膨張比 eps) — 80% ベル基準で digitize
# 出典: Rao(1960) のチャート (Sutton/RRS 記事に再掲)。Lf=80% の代表値。
# 他の Lf ではおおよそ theta は緩やかに変わるが、ここでは 80% 値を基準とし
# --thetan / --thetae で明示上書きできるようにする。
# ---------------------------------------------------------------------------
_RAO_EPS = np.array([5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0])
_RAO_TN = np.array([32.5, 34.0, 35.0, 35.8, 37.0, 37.8, 38.5])  # theta_n [deg]
_RAO_TE = np.array([14.0, 11.8, 10.5, 9.7, 8.7, 8.0, 7.5])      # theta_e [deg]


def rao_angles(eps: float) -> tuple[float, float]:
    """膨張比 eps から (theta_n, theta_e) [deg] を log(eps) 線形補間で読む。"""
    le = math.log(eps)
    tn = float(np.interp(le, np.log(_RAO_EPS), _RAO_TN))
    te = float(np.interp(le, np.log(_RAO_EPS), _RAO_TE))
    return tn, te


def area_ratio_from_mach(me: float, gamma: float) -> float:
    """設計出口 Mach から等エントロピー面積比 Ae/At を返す。"""
    g = gamma
    return (1.0 / me) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * me * me)) ** (
        (g + 1.0) / (2.0 * (g - 1.0))
    )


def arc(cx, cy, r, a0, a1, n):
    """中心(cx,cy)半径r, 角度a0->a1 [rad] の円弧点列 (始点含み終点含む)。"""
    a = np.linspace(a0, a1, n)
    return cx + r * np.cos(a), cy + r * np.sin(a)


def build(args):
    r"""幾何を exact primitive (line/arc/bezier) のセグメント列として構築する.

    設計フロー (ユーザ確定方針):
      1. C1 コニカルを「良い形状」として先に決める:
         スロート下流円弧 (Rfd=0.382Rt) を半角 alpha (既定 15deg) の直線コーンに
         *接して* つなぐ (= 円弧の終端角を alpha にとる → 1 階微分連続)。
         出口は膨張比 eps (Me) で固定し、C1 を満たす軸長 Ln を解く。
      2. Rao ベルを同じ Rt・同じ Ln・同じ出口半径 Re で作り直す:
         円弧を角 theta_n まで回し、そこから 2 次 Bézier で (Ln, Re) へ。
         始点接線 theta_n・終点接線 theta_e。円弧終端角 theta_n = Bézier 始点接線
         なので C1 連続。

    返り値: bell/conical それぞれの (sampled points, segment list) と params。
    各 segment は dict:
      {'kind':'line', 'p0':(x,y), 'p1':(x,y), 'zone':..}
      {'kind':'arc',  'p0':, 'p1':, 'center':(cx,cy), 'r':R, 'zone':..}
      {'kind':'bezier','p0':N, 'pc':Q, 'p1':E, 'zone':..}
    """
    Rt = args.rt
    gamma = args.gamma
    Rfd = 0.382 * Rt          # スロート下流円弧半径
    Rfu = 1.5 * Rt            # スロート上流円弧半径
    Rc = args.rc
    beta = args.beta * DEG

    # --- 出口 (膨張比で固定) ---
    eps = args.eps if args.eps is not None else area_ratio_from_mach(args.mach, gamma)
    Re = Rt * math.sqrt(eps)

    # --- 収束側 (両ノズル共通): chamber 直管 / 収縮コーン / スロート上流円弧 ---
    # スロート上流円弧: 中心 (0, Rt+Rfu)、スロート T=(0,Rt) で水平接線、
    # 収縮側端 conv_join は接線角 beta。
    conv_join = (-Rfu * math.sin(beta), (Rt + Rfu) - Rfu * math.cos(beta))
    T = (0.0, Rt)
    # 収縮コーン: conv_join から半角 beta で Rc まで
    cone_start_x = conv_join[0] - (Rc - conv_join[1]) / math.tan(beta)
    cone_start = (cone_start_x, Rc)
    cham_start = (cone_start_x - args.chamber, Rc)

    common = [
        dict(kind="line", p0=cham_start, p1=cone_start, zone="chamber"),
        dict(kind="line", p0=cone_start, p1=conv_join, zone="conv_cone"),
        dict(kind="arc", p0=conv_join, p1=T, center=(0.0, Rt + Rfu), r=Rfu, zone="throat_up"),
    ]

    # --- C1 コニカル: 円弧を alpha まで → 接する直線コーン → 出口 Re ---
    alpha = args.cone_half * DEG
    Nc = (Rfd * math.sin(alpha), (Rt + Rfd) * 1.0 - Rfd * math.cos(alpha))
    # コーンは Nc から半角 alpha。出口 Re に届く軸長 Ln を解く (C1 前提)。
    Ln = Nc[0] + (Re - Nc[1]) / math.tan(alpha)
    Ec = (Ln, Re)
    conical_div = [
        dict(kind="arc", p0=T, p1=Nc, center=(0.0, Rt + Rfd), r=Rfd, zone="throat_dn"),
        dict(kind="line", p0=Nc, p1=Ec, zone="cone"),
    ]

    # --- Rao ベル: 同じ Rt/Ln/Re。円弧を theta_n まで → 2 次 Bézier → (Ln,Re) ---
    tn_deg, te_deg = rao_angles(eps)
    if args.thetan is not None:
        tn_deg = args.thetan
    if args.thetae is not None:
        te_deg = args.thetae
    tn, te = tn_deg * DEG, te_deg * DEG
    N = (Rfd * math.sin(tn), (Rt + Rfd) - Rfd * math.cos(tn))
    E = (Ln, Re)
    mn, me_ = math.tan(tn), math.tan(te)
    Qx = (N[1] - E[1] - mn * N[0] + me_ * E[0]) / (me_ - mn)
    Qy = N[1] + mn * (Qx - N[0])
    Q = (Qx, Qy)
    bell_div = [
        dict(kind="arc", p0=T, p1=N, center=(0.0, Rt + Rfd), r=Rfd, zone="throat_dn"),
        dict(kind="bezier", p0=N, pc=Q, p1=E, zone="parabola"),
    ]

    bell_segs = common + bell_div
    conical_segs = common + conical_div

    # --- CSV / 描画用サンプリング ---
    def sample(seg, n):
        if seg["kind"] == "line":
            x = np.linspace(seg["p0"][0], seg["p1"][0], n)
            y = np.linspace(seg["p0"][1], seg["p1"][1], n)
        elif seg["kind"] == "arc":
            cx, cy = seg["center"]
            a0 = math.atan2(seg["p0"][1] - cy, seg["p0"][0] - cx)
            a1 = math.atan2(seg["p1"][1] - cy, seg["p1"][0] - cx)
            a = np.linspace(a0, a1, n)
            x = cx + seg["r"] * np.cos(a)
            y = cy + seg["r"] * np.sin(a)
        else:  # bezier (quadratic)
            t = np.linspace(0, 1, n)
            P0, Pc, P1 = seg["p0"], seg["pc"], seg["p1"]
            x = (1 - t) ** 2 * P0[0] + 2 * (1 - t) * t * Pc[0] + t**2 * P1[0]
            y = (1 - t) ** 2 * P0[1] + 2 * (1 - t) * t * Pc[1] + t**2 * P1[1]
        return x, y

    def stack(segs):
        xs, ys, zs = [], [], []
        for s in segs:
            x, y = sample(s, 60 if s["kind"] != "line" else 2)
            xs.extend(x); ys.extend(y); zs.extend([s["zone"]] * len(x))
        return np.array(xs), np.array(ys), zs

    bx, by, bz = stack(bell_segs)
    cx_, cy_, cz = stack(conical_segs)

    # 15deg 基準コーン長との比 (= ベル % length 診断)
    Nc15 = (Rfd * math.sin(15 * DEG), (Rt + Rfd) - Rfd * math.cos(15 * DEG))
    L15 = Nc15[0] + (Re - Nc15[1]) / math.tan(15 * DEG)

    params = dict(
        Rt=Rt, Re=Re, Rc=Rc, eps=eps, gamma=gamma,
        Me_design=args.mach if args.eps is None else float("nan"),
        theta_n=tn_deg, theta_e=te_deg,
        Ln=Ln, Nx=N[0], Ny=N[1], Qx=Qx, Qy=Qy,
        cone_half_deg=args.cone_half, cone_Ncx=Nc[0], cone_Ncy=Nc[1],
        percent=100.0 * Ln / L15,
        beta=args.beta, chamber=args.chamber,
        x_start=cham_start[0], x_throat=0.0, x_exit=Ln,
    )
    return (bx, by, bz, bell_segs), (cx_, cy_, cz, conical_segs), params


def write_csv(path, x, y, z):
    with open(path, "w") as f:
        f.write("x_mm,y_mm,zone\n")
        for xi, yi, zi in zip(x, y, z):
            f.write(f"{xi:.6f},{yi:.6f},{zi}\n")


# 各 zone の x 分割数 (壁の曲率に応じて配分)
NX_ZONE = dict(chamber=18, conv_cone=45, throat_up=30, throat_dn=28,
               parabola=150, cone=150)


def write_geo(path, which, segs, p, ny, prog_r=1.0, bump=None):
    """exact primitive (Line/Circle/Bezier) による多ブロック構造軸対称メッシュ.

    各セグメント (chamber/conv_cone/throat_up/throat_dn/parabola|cone) を 1 ブロックとし、
    壁の折れ (kink) は必ずブロック境界に置く。スロート T=(0,Rt) はブロック角点になるので
    セルがスロートを跨がず、負体積を生まない。壁は近似 Spline ではなく厳密プリミティブ:
      line -> Line, arc -> Circle(start,center,end), parabola -> Bezier(N,Q,E)。
    Physical: inlet(1) / outlet(2) / wall(3) / axis(4)。ScalingFactor=0.001 (mm->m)。
    """
    L = []
    em = L.append
    em(f"// case 29 {which} nozzle — axisymmetric multi-block mesh (exact primitives)")
    em("Geometry.PointNumbers = 0; Geometry.LineNumbers = 0;")
    em("Mesh.ScalingFactor = 0.001;   // mm -> m")
    em("lc = 1.0;")
    em(f"ny = {ny};")
    em(f"prog_r = {prog_r};   // 半径方向クラスタリング (1.0=一様; <1 で壁に密)")
    if bump is not None:
        em(f"bump_r = {bump};  // 両端クラスタ (壁+軸): Transfinite Bump")
    em("")

    pid = [0]
    pcache = {}

    def P(x, y):
        key = (round(x, 7), round(y, 7))
        if key not in pcache:
            pid[0] += 1
            i = pid[0]
            em(f"Point({i}) = {{{x:.7f}, {y:.7f}, 0.0, lc}};")
            pcache[key] = i
        return pcache[key]

    lid = [0]

    def newline(stmt):
        lid[0] += 1
        em(f"Line({lid[0]}) = {stmt};")
        return lid[0]

    def newcirc(a, c, b):
        lid[0] += 1
        em(f"Circle({lid[0]}) = {{{a}, {c}, {b}}};")
        return lid[0]

    def newbez(a, c, b):
        lid[0] += 1
        em(f"Bezier({lid[0]}) = {{{a}, {c}, {b}}};")
        return lid[0]

    wall_curves, axis_curves, surfs = [], [], []
    sid = [0]

    em("// ---- blocks (left -> right) ----")
    prev_vert = None  # 直前ブロックの右縦線 (= 次ブロックの左縦線として共有)
    inlet_curve = None
    for s in segs:
        x0, y0 = s["p0"]
        x1, y1 = (s["p1"])
        nx = NX_ZONE.get(s["zone"], 40)
        pw0, pw1 = P(x0, y0), P(x1, y1)
        pa0, pa1 = P(x0, 0.0), P(x1, 0.0)
        # top wall curve (exact)
        if s["kind"] == "line":
            top = newline(f"{{{pw0}, {pw1}}}")
        elif s["kind"] == "arc":
            pc = P(*s["center"])
            top = newcirc(pw0, pc, pw1)
        else:  # bezier
            pc = P(*s["pc"])
            top = newbez(pw0, pc, pw1)
        bot = newline(f"{{{pa0}, {pa1}}}")          # axis
        left = prev_vert if prev_vert is not None else newline(f"{{{pa0}, {pw0}}}")
        right = newline(f"{{{pa1}, {pw1}}}")
        if inlet_curve is None:
            inlet_curve = left
        em(f"Transfinite Curve {{{top}, {bot}}} = {nx + 1};")
        radial_spec = "Bump bump_r" if bump is not None else "Progression prog_r"
        em(f"Transfinite Curve {{{right}}} = ny + 1 Using {radial_spec};")
        if prev_vert is None:
            em(f"Transfinite Curve {{{left}}} = ny + 1 Using {radial_spec};")
        sid[0] += 1
        ll = sid[0]
        em(f"Curve Loop({ll}) = {{{bot}, {right}, -{top}, -{left}}};")
        em(f"Plane Surface({ll}) = {{{ll}}};")
        em(f"Transfinite Surface {{{ll}}} = {{{pa0}, {pa1}, {pw1}, {pw0}}};")
        em(f"Recombine Surface {{{ll}}};")
        wall_curves.append(top)
        axis_curves.append(bot)
        surfs.append(ll)
        prev_vert = right
    outlet_curve = prev_vert

    em("")
    em(f'Physical Curve("inlet", 1)  = {{{inlet_curve}}};')
    em(f'Physical Curve("outlet", 2) = {{{outlet_curve}}};')
    em(f'Physical Curve("wall", 3)   = {{{",".join(map(str, wall_curves))}}};')
    em(f'Physical Curve("axis", 4)   = {{{",".join(map(str, axis_curves))}}};')
    em(f'Physical Surface("fluid", 8) = {{{",".join(map(str, surfs))}}};')
    Path(path).write_text("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rt", type=float, default=10.0, help="スロート半径 [mm]")
    ap.add_argument("--mach", type=float, default=4.0, help="設計出口 Mach (eps 未指定時)")
    ap.add_argument("--eps", type=float, default=None, help="膨張比 Ae/At を直接指定")
    ap.add_argument("--gamma", type=float, default=1.4, help="比熱比")
    ap.add_argument("--percent", type=float, default=80.0, help="ベル長 (15deg コニカル比 %%)")
    ap.add_argument("--thetan", type=float, default=None, help="始点角 theta_n [deg] 上書き")
    ap.add_argument("--thetae", type=float, default=None, help="出口角 theta_e [deg] 上書き")
    ap.add_argument("--rc", type=float, default=30.0, help="チャンバー半径 [mm]")
    ap.add_argument("--beta", type=float, default=30.0, help="収縮半角 [deg]")
    ap.add_argument("--chamber", type=float, default=20.0, help="チャンバー直管長 [mm]")
    ap.add_argument("--cone-half", type=float, default=15.0,
                    help="コニカル半角 [deg] (C1: スロート円弧をこの角で接続)。出口εと共に軸長Lnを決める")
    ap.add_argument("--ndiv", type=int, default=120, help="放物線分割数")
    ap.add_argument("--nx-conv", type=int, default=80, help="収束部 x 分割")
    ap.add_argument("--nx-div", type=int, default=140, help="発散部 x 分割")
    ap.add_argument("--ny", type=int, default=80, help="半径方向分割")
    ap.add_argument("--prog-r", type=float, default=1.0,
                    help="半径方向クラスタリング比 (1.0=一様; 粘性は壁密 <1, 例 0.95)")
    ap.add_argument("--bump-r", type=float, default=None,
                    help="両端クラスタ (壁+軸, Transfinite Bump 係数 <1)。指定時 prog-r より優先")
    ap.add_argument("--outdir", type=str, default=".", help="出力先")
    args = ap.parse_args()

    (bx, by, bz, bell_segs), (cx, cy, cz, conical_segs), p = build(args)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "contour_bell.csv", bx, by, bz)
    write_csv(out / "contour_conical.csv", cx, cy, cz)
    write_geo(out / "bell.geo", "bell", bell_segs, p, args.ny, args.prog_r, args.bump_r)
    write_geo(out / "conical.geo", "conical", conical_segs, p, args.ny, args.prog_r, args.bump_r)

    # --- 諸元テキスト ---
    lines = [
        "Rao TOP bell vs C1 conical nozzle (same throat / length / exit area)",
        "=" * 60,
        f"throat radius Rt        : {p['Rt']:.3f} mm",
        f"exit radius   Re        : {p['Re']:.3f} mm  (both nozzles)",
        f"chamber radius Rc       : {p['Rc']:.3f} mm",
        f"expansion ratio eps     : {p['eps']:.3f}  (Ae/At)",
        f"design exit Mach        : {p['Me_design']:.3f}  (gamma={p['gamma']})",
        f"axial length Ln         : {p['Ln']:.3f} mm  (throat x=0 -> exit; both)",
        "-- conical --",
        f"cone half angle         : {p['cone_half_deg']:.2f} deg (C1: arc tangent to cone)",
        f"  arc->cone join Nc     : ({p['cone_Ncx']:.3f}, {p['cone_Ncy']:.3f}) mm",
        "-- Rao bell --",
        f"theta_n (start angle)   : {p['theta_n']:.2f} deg",
        f"theta_e (exit angle)    : {p['theta_e']:.2f} deg",
        f"bell length / 15deg-cone: {p['percent']:.1f}%",
        f"parabola start N        : ({p['Nx']:.3f}, {p['Ny']:.3f}) mm",
        f"Bezier control Q        : ({p['Qx']:.3f}, {p['Qy']:.3f}) mm",
        "-- shared converging --",
        f"converging half angle   : {p['beta']:.1f} deg",
        f"chamber straight length : {p['chamber']:.1f} mm",
        f"axial extent            : x = {p['x_start']:.3f} .. {p['x_exit']:.3f} mm",
    ]
    txt = "\n".join(lines)
    (out / "nozzle_params.txt").write_text(txt + "\n")
    print(txt)

    # --- 描画 ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axx = plt.subplots(2, 1, figsize=(11, 7), height_ratios=[2, 1])
    ax = axx[0]
    for sign in (1, -1):
        ax.plot(cx, sign * cy, color="tab:orange", lw=1.6,
                label=f"conical {p['cone_half_deg']:.0f}deg (C1)" if sign == 1 else None)
        ax.plot(bx, sign * by, color="tab:blue", lw=1.8,
                label="Rao bell" if sign == 1 else None)
    ax.axvline(0.0, color="gray", ls=":", lw=0.8)
    ax.plot([p["Nx"]], [p["Ny"]], "ko", ms=4)
    ax.annotate("N", (p["Nx"], p["Ny"]), textcoords="offset points", xytext=(2, 6))
    ax.plot([p["Qx"]], [p["Qy"]], "k+", ms=8)
    ax.annotate("Q (Bezier ctrl)", (p["Qx"], p["Qy"]),
                textcoords="offset points", xytext=(4, 0), fontsize=8)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title(f"case 29: Rao bell vs conical  (Rt={p['Rt']:.0f}mm, eps={p['eps']:.2f}, "
                 f"Me={p['Me_design']:.1f}, Ln={p['Ln']:.1f}mm)")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)

    # 半径プロファイル (上半分のみ, 拡大)
    ax2 = axx[1]
    ax2.plot(cx, cy, color="tab:orange", lw=1.6, label="conical")
    ax2.plot(bx, by, color="tab:blue", lw=1.8, label="Rao bell")
    ax2.axvline(0.0, color="gray", ls=":", lw=0.8)
    ax2.set_xlabel("x [mm]")
    ax2.set_ylabel("r [mm]")
    ax2.set_title("wall radius profile (upper half)")
    ax2.legend(loc="upper left")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out / "nozzle_preview.png", dpi=130)
    print(f"\nwrote: {out/'contour_bell.csv'}, {out/'contour_conical.csv'}, "
          f"{out/'nozzle_params.txt'}, {out/'nozzle_preview.png'}")


if __name__ == "__main__":
    main()
