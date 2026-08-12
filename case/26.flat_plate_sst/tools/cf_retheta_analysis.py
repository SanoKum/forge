#!/usr/bin/env python3
"""case/26 平板: $C_f(Re_\\theta)$ を外部相関 (Kármán–Schoenherr) と比較する正式ツール。

README に載る数値はすべて本ツールで再生成できること。scratchpad の使い捨てスクリプトで
数値を出して README に書くのは禁止 (2026-08-12 レビュー指摘)。

**$C_f$ の 4 定義を分離して出す** (混ぜないこと):

  (a) resolved : 解像勾配 tau_w = mu*u_1/y_1                  … 壁関数メッシュでは無意味 (参考値)
  (w) wi-force : W-I 双対面で **実際に運動量残差へ加えた** 接線力 (要 --wi-force、
                 FORGE_WI_FORCE_DIAG=1 で走らせた run が必要)
  (b) solver   : ソルバ自身の壁面 Cf 出力                      … SU2 のみ (forge は本体 res に無い)
  (c) momentum : 運動量積分 Cf = 2(dtheta/dx + PG)             … **壁出力に依存しない収支診断**
  (d) walllaw  : 壁法則の逆解きで u_tau を求めた Cf = 2(u_tau/Ue)^2

**(c) は「定義非依存の壁摩擦」ではない**。定常性・境界層近似・Ue/rho_e の取り方・積分上端・
微分 fit に依存する。あくまで (a)(b)(d) と突き合わせる**収支診断**として使う。

**(d) の壁法則はソルバごとに違う**ので `--wall-law` で明示すること:
  reichardt : forge の SST automatic wall treatment (`wallLaw_reichardt_uplus`)
  spalding  : Spalding 則 (SU2 STANDARD_WALL_FUNCTION は Nichols&Nelson の圧縮性版なので近似)
  none      : (d) を計算しない (壁法則が違うソルバに Reichardt を当てない)

運動量積分は圧力勾配込みの完全形:
  Cf/2 = dtheta/dx + theta*(H + 2 - Me^2)/Ue * dUe/dx,   H = delta*/theta

使い方:
  python3 tools/cf_retheta_analysis.py run_0042_node_yp30_planar_2nd
  python3 tools/cf_retheta_analysis.py run_0049_su2_sst_lowre_fine --su2 --wall-law none
  python3 tools/cf_retheta_analysis.py <run> --ytop 0.05 --fit-window 0.12 --fit-order 3
  python3 tools/cf_retheta_analysis.py <run> --edge local|global
"""
import sys, os, re, glob, struct, argparse
import numpy as np
import h5py

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# --- 自由流 (config と一致させる。solverConfig.yaml / SU2 cfg 双方の基準) ---
MU, GAM, CP = 1.8e-5, 1004.5 * 0.4 / 1.4 * 3.5, 1004.5   # CP は下で使う
CP = 1004.5
GAM = 1.4
R = CP * (GAM - 1.0) / GAM
PT, TT, MACH = 100000.0, 288.15, 0.2
PS = PT * (1.0 + 0.5 * (GAM - 1.0) * MACH ** 2) ** (-GAM / (GAM - 1.0))
TS = TT * (1.0 + 0.5 * (GAM - 1.0) * MACH ** 2) ** (-1.0)
RHO_INF = PS / (R * TS)
A_INF = np.sqrt(GAM * R * TS)
U_INF = MACH * A_INF
KAPPA, BLOG = 0.41, 5.0


# ---------------------------------------------------------------- 外部相関
def cf_karman_schoenherr(re_theta):
    """NASA TMR flat-plate validation と同じ Kármán–Schoenherr 相関 (外部相関)。"""
    ret = np.asarray(re_theta, float)
    l = np.log10(np.where(ret > 1.0, ret, np.nan))
    return 1.0 / (17.08 * l * l + 25.11 * l + 6.012)


def cf_schlichting(re_x, coeff=0.0592):
    """補助基準。0.0592 = Schlichting の経験フィット / 0.0576 = 1/7 乗則の自己整合値 (別物)。"""
    return coeff * np.asarray(re_x, float) ** -0.2


# ---------------------------------------------------------------- 壁法則
def uplus_reichardt(yp):
    return np.log(1.0 + KAPPA * yp) / KAPPA + 7.8 * (
        1.0 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-yp / 3.0))


def uplus_spalding(yp, up):
    """Spalding: y+ = u+ + e^{-kB}(e^{ku+} - 1 - ku+ - (ku+)^2/2 - (ku+)^3/6) の残差。"""
    ku = KAPPA * up
    return up + np.exp(-KAPPA * BLOG) * (
        np.exp(ku) - 1.0 - ku - ku * ku / 2.0 - ku ** 3 / 6.0) - yp


def solve_utau(Ut, y, nu, law="reichardt"):
    if law == "none":
        return float("nan")
    utau = np.sqrt(max(nu * Ut / max(y, 1e-30), 1e-30))
    for _ in range(200):
        utau = max(utau, 1e-12)
        if law == "reichardt":
            f = Ut / utau - uplus_reichardt(utau * y / nu)
        else:  # spalding
            f = uplus_spalding(utau * y / nu, Ut / utau)
        d = 1e-6 * utau
        if law == "reichardt":
            f2 = Ut / (utau + d) - uplus_reichardt((utau + d) * y / nu)
        else:
            f2 = uplus_spalding((utau + d) * y / nu, Ut / (utau + d))
        df = (f2 - f) / d
        if abs(df) < 1e-30:
            break
        step = f / df
        utau -= step
        if abs(step) < 1e-14 * max(utau, 1e-12):
            break
    return max(utau, 0.0)


# ---------------------------------------------------------------- 読み込み
def load_forge(run, fname=None):
    d = os.path.join(HERE, run)
    p = os.path.join(d, fname) if fname else max(
        [f for f in glob.glob(os.path.join(d, "res_*.h5"))
         if "wall" not in os.path.basename(f) and os.path.basename(f) != "res_0.h5"],
        key=lambda f: int(re.findall(r"res_(\d+)\.h5", os.path.basename(f))[0]))
    h = h5py.File(p, "r")
    coord = h["MESH/COORD"][:].reshape(-1, 3)
    n = h["VALUE/ro"].shape[0]
    if coord.shape[0] == n:                       # node: DOF = ノード
        C = coord
    else:                                         # cell: 要素重心
        conne = h["MESH/CONNE"][:]
        conn = conne.reshape(n, conne.size // n)[:, 1:]
        if conn.min() == 1:
            conn = conn - 1
        C = coord[conn].mean(axis=1)
    ro = h["VALUE/ro"][:]
    V = dict(ro=ro, u=h["VALUE/roUx"][:] / ro, mu=h["VALUE/vis_lam"][:],
             P=h["VALUE/P"][:], cf_solver=None)
    return C, V, os.path.basename(p)


def wi_force_band(run, fname, lo, hi):
    """W-I 実力: 帯 [lo,hi] の壁ノードで実際に加えた接線力の和 [N/m 幅] と壁ノード x を返す。

    **規約 (曖昧さを残さないため明示)**:
      - `wi_ftan[i]` は壁ノード i に入射する W-I 双対面の |接線 traction| (面積込み [N]) の和。
        壁ノードの双対幅 dx_i = (x_{i+1}-x_{i-1})/2 で割って **局所 tau_w,applied(x_i)** に直す。
        こうしないと「点和 vs 台形積分」で帯の端が半セルずれ、数 % ではなく 10% 級の偏りが出る。
      - 返すのは (x_i, tau_w_applied(x_i))。呼び出し側で Cf_KS と**同じ台形則**で積分して比を取る。
        したがって比は **同一離散の積分同士の比**である。
      - 端点は双対幅が片側になるノード (帯の最初/最後) を除外する。
    """
    d = os.path.join(HERE, run)
    p = os.path.join(d, fname) if fname else max(
        [f for f in glob.glob(os.path.join(d, "res_*.h5"))
         if "wall" not in os.path.basename(f) and os.path.basename(f) != "res_0.h5"],
        key=lambda f: int(re.findall(r"res_(\d+)\.h5", os.path.basename(f))[0]))
    h = h5py.File(p, "r")
    if "wi_ftan" not in h["VALUE"]:
        raise SystemExit(f"{p} に wi_ftan が無い。FORGE_WI_FORCE_DIAG=1 で走らせた run が要る。")
    C = h["MESH/COORD"][:].reshape(-1, 3)
    x, y = C[:, 0], C[:, 1]
    ft = h["VALUE/wi_ftan"][:]
    fna = h["VALUE/wi_fnrm_abs"][:] if "wi_fnrm_abs" in h["VALUE"] else None
    fr = h["VALUE/wi_ftan_res"][:] if "wi_ftan_res" in h["VALUE"] else None
    w = (np.abs(y) < 1e-12) & (x > 1e-6)
    iw = np.where(w)[0]
    iw = iw[np.argsort(x[iw])]
    eps = 1e-6                       # station は round(x,6) なので端点は許容差で拾う
    xall = x[iw]
    dxall = np.empty_like(xall)
    dxall[1:-1] = 0.5 * (xall[2:] - xall[:-2])
    dxall[0] = xall[1] - xall[0]
    dxall[-1] = xall[-1] - xall[-2]
    tau_all = ft[iw] / np.maximum(dxall, 1e-30)          # 局所 tau_w,applied [Pa]
    m = (xall >= lo - eps) & (xall <= hi + eps)
    m[0] = m[-1] = False                                  # 双対幅が片側のノードを除外
    extra = {}
    if fna is not None and ft[iw][m].sum() > 0:
        extra["fn_ratio"] = float(fna[iw][m].sum() / ft[iw][m].sum())
        extra["fn_ratio_max"] = float((fna[iw][m] / np.maximum(ft[iw][m], 1e-30)).max())
    if fr is not None:
        rs = ft[iw][m] / np.maximum(fr[iw][m], 1e-30)
        extra["scale_min"], extra["scale_med"], extra["scale_max"] = \
            float(rs.min()), float(np.median(rs)), float(rs.max())
    return xall[m], tau_all[m], extra


def load_su2(run, fname="flow.vtu"):
    path = os.path.join(HERE, run, fname)
    raw = open(path, "rb").read()
    hdr = raw[:raw.index(b"<AppendedData")].decode("utf8", "ignore")
    off0 = raw.index(b"_", raw.index(b"<AppendedData")) + 1
    arrs = re.findall(
        r'<DataArray type="(\w+)" Name="([^"]*)" NumberOfComponents=\s*"(\d+)" offset="(\d+)"', hdr)
    tmap = {"Float32": "f4", "Float64": "f8", "Int32": "i4", "Int64": "i8", "UInt8": "u1"}
    out = {}
    for t, n, nc, off in arrs:
        o = off0 + int(off)
        nb = struct.unpack("<Q", raw[o:o + 8])[0]
        a = np.frombuffer(raw[o + 8:o + 8 + nb], dtype=np.dtype("<" + tmap[t]))
        out[n if n else "pts"] = a.reshape(-1, int(nc)) if int(nc) > 1 else a
    C = out["pts"]
    V = dict(ro=out["Density"], u=out["Velocity"][:, 0], mu=out["Laminar_Viscosity"],
             P=out["Pressure"],
             cf_solver=np.abs(out["Skin_Friction_Coefficient"][:, 0]))
    return C, V, fname


# ---------------------------------------------------------------- 抽出
def stations(C, xmin, xmax):
    x = C[:, 0]
    xa = np.unique(np.round(x[x > 1e-6], 6))
    return xa[(xa >= xmin) & (xa <= xmax)]


def column(C, xc):
    x, y = C[:, 0], C[:, 1]
    i = np.sort(np.where(np.abs(x - xc) < 1e-4)[0])
    return i[np.argsort(y[i])]


def profile(C, V, xc, ytop, edge):
    """壁点を必ず含むカラムを返す。cell は (y,u)=(0,0) を補い、node は既にある壁点を使う。"""
    col = column(C, xc)
    y = C[col, 1].astype(float)
    u = V["u"][col].astype(float)
    ro = V["ro"][col].astype(float)
    mu = V["mu"][col].astype(float)
    if y[0] > 1e-12:                               # cell: 壁点を補う
        y = np.r_[0.0, y]; u = np.r_[0.0, u]; ro = np.r_[ro[0], ro]; mu = np.r_[mu[0], mu]
    m = np.ones(len(y), bool) if ytop is None else (y <= ytop)
    j = int(np.argmax(u[m]))
    if edge == "global":
        ue, roe, mue = U_INF, RHO_INF, MU
    else:
        ue, roe, mue = u[m][j], ro[m][j], mu[m][j]
    theta = float(np.trapz((ro[m] / roe) * (u[m] / ue) * (1.0 - u[m] / ue), y[m]))
    dstar = float(np.trapz(1.0 - (ro[m] * u[m]) / (roe * ue), y[m]))
    return dict(y=y, u=u, ro=ro, mu=mu, ue=ue, roe=roe, mue=mue,
                theta=theta, dstar=dstar, col=col, y1=y[1], u1=u[1], nu1=mu[1] / ro[1])


# ---------------------------------------------------------------- 本体
def analyse(C, V, args):
    xa = stations(C, args.xmin, args.xmax)
    prof = {x: profile(C, V, x, args.ytop, args.edge) for x in xa}
    th = np.array([prof[x]["theta"] for x in xa])
    ue = np.array([prof[x]["ue"] for x in xa])
    rows = []
    for xt in args.stations:
        i = int(np.argmin(np.abs(xa - xt)))
        xc = xa[i]
        w = np.abs(xa - xc) <= args.fit_window
        if w.sum() < args.fit_order + 2:
            w = np.abs(xa - xc) <= 2 * args.fit_window
        dthdx = np.polyval(np.polyder(np.polyfit(xa[w], th[w], args.fit_order)), xc)
        duedx = np.polyval(np.polyder(np.polyfit(xa[w], ue[w], args.fit_order)), xc)
        p = prof[xc]
        Me = p["ue"] / np.sqrt(GAM * R * TS)
        H = p["dstar"] / p["theta"]
        pg = p["theta"] * (H + 2.0 - Me * Me) / p["ue"] * duedx
        q = 0.5 * p["roe"] * p["ue"] ** 2
        ret = p["roe"] * p["ue"] * p["theta"] / p["mue"]
        rex = RHO_INF * U_INF * xc / MU
        ks = float(cf_karman_schoenherr(ret))
        cf_a = (p["mu"][1] * p["u1"] / p["y1"]) / q
        cf_b = float(V["cf_solver"][p["col"][0]]) if V["cf_solver"] is not None else float("nan")
        cf_c = 2.0 * (dthdx + pg)
        utau = solve_utau(p["u1"], p["y1"], p["nu1"], args.wall_law)
        cf_d = 2.0 * (utau / p["ue"]) ** 2
        rows.append(dict(x=xc, rex=rex, theta=p["theta"], ret=ret, H=H, ks=ks,
                         a=cf_a, b=cf_b, c=cf_c, d=cf_d, pgfrac=pg / dthdx,
                         yp1=utau * p["y1"] / p["nu1"], nfit=int(w.sum()),
                         schl=float(cf_schlichting(rex))))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--file", default=None, help="res_*.h5 / flow.vtu (既定: 最新)")
    ap.add_argument("--su2", action="store_true", help="SU2 の flow.vtu を読む")
    ap.add_argument("--wall-law", default="reichardt", choices=["reichardt", "spalding", "none"],
                    help="(d) の壁法則。SU2 は Nichols&Nelson なので reichardt を当てないこと")
    ap.add_argument("--ytop", type=float, default=None, help="theta 積分の上端 [m] (既定: 全域)")
    ap.add_argument("--edge", default="local", choices=["local", "global"],
                    help="Ue/rho_e を BL 外縁の局所値 (local) か自由流固定 (global) か")
    ap.add_argument("--fit-window", type=float, default=0.08, help="dtheta/dx 局所 fit の半幅 [m]")
    ap.add_argument("--fit-order", type=int, default=2, help="局所 fit の次数")
    ap.add_argument("--xmin", type=float, default=0.10, help="前縁側の除外境界 [m]")
    ap.add_argument("--xmax", type=float, default=0.95, help="出口側の除外境界 [m]")
    ap.add_argument("--stations", type=float, nargs="+", default=[0.30, 0.45, 0.60, 0.75, 0.90])
    ap.add_argument("--wi-force", metavar="RUN", default=None,
                    help="W-I 実力 (wi_ftan) を持つ run (FORGE_WI_FORCE_DIAG=1 で走らせたもの)。"
                         "--band の帯で ΣF_t と ∫Cf_KS·q dx の比を出す")
    ap.add_argument("--wi-file", default=None, help="--wi-force run 内の res_*.h5 (既定: 最新)")
    ap.add_argument("--band", type=float, nargs=2, default=None, metavar=("LO", "HI"),
                    help="帯平均の範囲 [m]。--wi-force と併用")
    ap.add_argument("--retheta-min", type=float, default=None,
                    help="この Re_theta 未満の station を判定から除外する (K-S 主比較域ゲート)")
    args = ap.parse_args()

    if args.su2:
        C, V, used = load_su2(args.run, args.file or "flow.vtu")
    else:
        C, V, used = load_forge(args.run, args.file)

    print(f"# run       : {args.run}")
    print(f"# file      : {used}")
    print(f"# 設定      : wall-law={args.wall_law}  ytop={args.ytop}  edge={args.edge}  "
          f"fit(window={args.fit_window}, order={args.fit_order})  x=[{args.xmin},{args.xmax}]")
    print(f"# 自由流    : U_inf={U_INF:.4f} m/s  rho_inf={RHO_INF:.5f}  Re/m={RHO_INF*U_INF/MU:.4e}")
    print(f"# 外部相関  : Karman-Schoenherr。主比較域 4000<Re_theta<13000")
    print(f"# Cf 定義   : (a)解像勾配 (b)ソルバ出力 (c)運動量積分[収支診断] (d)壁法則逆解き")
    print()
    hdr = (f'{"x":>6s} {"Re_x":>9s} {"theta":>10s} {"Re_th":>7s} {"KS域":>5s} {"H":>5s} '
           f'{"(a)/KS":>7s} {"(b)/KS":>7s} {"(c)/KS":>7s} {"(d)/KS":>7s} '
           f'{"Cf/Schl":>8s} {"PG/dth":>7s} {"y+1":>6s} {"nfit":>4s}')
    print(hdr)
    for r in analyse(C, V, args):
        ink = "○" if 4000 <= r["ret"] <= 13000 else "-"
        f = lambda v: f"{v:7.3f}" if np.isfinite(v) else f'{"—":>7s}'
        print(f'{r["x"]:6.3f} {r["rex"]:9.3e} {r["theta"]:10.4e} {r["ret"]:7.0f} {ink:>5s} '
              f'{r["H"]:5.2f} {f(r["a"]/r["ks"])} {f(r["b"]/r["ks"])} {f(r["c"]/r["ks"])} '
              f'{f(r["d"]/r["ks"])} {r["d"]/r["schl"]:8.3f} {r["pgfrac"]*100:6.1f}% '
              f'{r["yp1"]:6.1f} {r["nfit"]:4d}')
    if args.wi_force and args.band:
        lo, hi = args.band
        # 先に Re_theta ゲートを掛けて **有効な x 範囲** を確定し、分子・分母を同じ範囲で取る
        xa = stations(C, args.xmin, args.xmax)
        prof = {x: profile(C, V, x, args.ytop, args.edge) for x in xa}
        ret = np.array([prof[x]["roe"] * prof[x]["ue"] * prof[x]["theta"] / prof[x]["mue"] for x in xa])
        q = np.array([0.5 * prof[x]["roe"] * prof[x]["ue"] ** 2 for x in xa])
        m = (xa >= lo) & (xa <= hi)
        if args.retheta_min is not None:
            m &= (ret >= args.retheta_min)
        if m.sum() < 2:
            raise SystemExit("帯 + Re_theta ゲートで station が 2 個未満。--band か --retheta-min を見直す。")
        lo_eff, hi_eff = float(xa[m].min()), float(xa[m].max())
        xw, tauw, extra = wi_force_band(args.wi_force, args.wi_file, lo_eff, hi_eff)
        # 同じ台形則で両辺を積分する (x は壁ノード列に統一)
        ksx = np.interp(xw, xa[m], cf_karman_schoenherr(ret[m]) * q[m])
        F   = float(np.trapz(tauw, xw))
        Fks = float(np.trapz(ksx, xw))
        lo, hi = float(xw.min()), float(xw.max())
        print()
        print(f"# W-I 実力 (有効帯 x=[{lo:.3f},{hi:.3f}]" +
              (f", Re_theta>={args.retheta_min:.0f} でゲート" if args.retheta_min else "") + ")")
        print(f"#   規約: wi_ftan を壁ノード双対幅で割って局所 tau_w,applied に直し、")
        print(f"#         Cf_KS·q を同じ壁ノード列へ内挿して **同一の台形則**で両辺を積分する")
        print(f"#         (点和 vs 台形積分の端点ずれを避けるため)")
        print(f"  ∫tau_w,applied dx  = {F:.6e} N   (壁ノード {len(xw)} 個, x=[{xw.min():.3f},{xw.max():.3f}])")
        print(f"  ∫Cf_KS·q dx        = {Fks:.6e} N   (同一 x 列で台形積分)")
        print(f"  **(w) W-I実力 Cf/KS = {F/Fks:.4f}**")
        if "fn_ratio" in extra:
            print(f"  Σ|F_n|/Σ|F_t| = {extra['fn_ratio']:.4e}  (局所最大 {extra['fn_ratio_max']:.4e})")
        if "scale_min" in extra:
            print(f"  AddTauWall 再スケール係数: min={extra['scale_min']:.4f} "
                  f"median={extra['scale_med']:.4f} max={extra['scale_max']:.4f}")
    print()
    print("注: (c) は壁出力に依存しないが「定義非依存の壁摩擦」ではない。定常性・境界層近似・")
    print("    Ue/rho_e の取り方・積分上端・微分 fit に依存する収支診断である。")
    print("    残差は check_convergence.py で全列を、準定常は check_quasisteady.py")
    print("    --quantity theta,cf_retheta で確認すること (rms_ro / RMS_DENSITY だけで判断しない)。")


if __name__ == "__main__":
    main()
