#!/usr/bin/env python3
"""凍結場での壁法則逆解き比較 — 連成 A/B の**前**にモデル形式感度を測る。

[plan §5.2 ①](../../plans/active/turbulence-reichardt-gap-residual.md) の実体化ツール。

## なぜ連成の前にこれをやるか

壁法則を差し替えると `ransWallFunction_d.cu` の**4 系統が同時に動く**
($u^+$ → $u_\\tau$ 逆解き / $g=du^+/dy^+$ → `wf_pk` / $g$ → E3 `wf_sprod` /
$u_\\tau$ → `roK_wf`・Kader 壁熱流束)。いきなり連成 A/B をやると交絡する。
**凍結場**なら「逆解き則を替えると $u_\\tau$ がどれだけ動くか」だけを取り出せる。

## 2 つのモード

- `--mode resolved` : **壁解像 run** の場を使う。真の $u_\\tau$ は壁第一節点の
  **分子粘性勾配**から独立に決める。指定 $y_{\\rm rep}$ へ補間した速度を各法則で逆解きし、
  $u_\\tau(\\text{法則})/u_\\tau(\\text{内部基準})$ を出す。
  → **どの法則が壁解像解と整合するか** (= 端点比分解の「壁法則段」を法則別に測る)。
- `--mode wf` : **壁関数 run** の場を使う。代表点 (第一内層ノード) の速度を各法則で逆解きし、
  現行 Reichardt に対する比を出す。
  → **連成させたときに $u_\\tau$ が動く向きと大きさの凍結場見積り**。
  まず本体 `utau` を Reichardt で再構成できるか検証してから比を出す。

**平板 (case/26) 前提**: 壁法線を $\\hat{\\mathbf n}=(0,1)$ とし接線速度を $|U_x|$ とする。

使い方:
  python3 tools/wall_law_inverse.py run_0064_B_nasa_ref --mode resolved --yrep 1.7681e-4
  python3 tools/wall_law_inverse.py run_0063_B_nasa_wf  --mode wf
"""
import os, re, glob, argparse
import numpy as np
import h5py

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KAPPA = 0.41
BLOG = 5.0
LAWS = ["reichardt", "spalding", "log"]


def u_plus(yp, law):
    """$u^+(y^+)$。定数は Reichardt: kappa=0.41 + 7.8 項 / Spalding,log: kappa=0.41, B=5.0。"""
    yp = max(yp, 1e-12)
    if law == "reichardt":
        return (np.log(1.0 + KAPPA * yp) / KAPPA
                + 7.8 * (1.0 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-yp / 3.0)))
    if law == "log":
        return np.log(yp) / KAPPA + BLOG
    if law == "spalding":
        up = np.log(yp) / KAPPA + BLOG
        for _ in range(300):
            e = KAPPA * up
            f = up + np.exp(-KAPPA * BLOG) * (np.exp(e) - 1 - e - e**2 / 2 - e**3 / 6) - yp
            df = 1 + np.exp(-KAPPA * BLOG) * KAPPA * (np.exp(e) - 1 - e - e**2 / 2)
            s = f / df
            up -= s
            if abs(s) < 1e-13:
                break
        return up
    raise SystemExit(f"未知の壁法則: {law}")


def solve_utau(Ut, y, nu, law):
    """$U_t/u_\\tau = u^+(u_\\tau y/\\nu)$ を Newton で解く (本体と同じ形)。"""
    ut = np.sqrt(max(nu * Ut / max(y, 1e-30), 1e-30))
    for _ in range(300):
        ut = max(ut, 1e-12)
        f = Ut / ut - u_plus(ut * y / nu, law)
        d = 1e-6 * ut
        df = ((Ut / (ut + d) - u_plus((ut + d) * y / nu, law)) - f) / d
        if abs(df) < 1e-30:
            break
        s = f / df
        ut -= s
        if abs(s) < 1e-14 * max(ut, 1e-12):
            break
    return max(ut, 0.0)


def latest(run, pattern):
    d = run if os.path.isabs(run) else os.path.join(HERE, run)
    fs = [f for f in glob.glob(os.path.join(d, pattern))
          if "outlet" not in os.path.basename(f) and os.path.basename(f) != "res_0.h5"]
    if not fs:
        raise SystemExit(f"{pattern} が見つからない: {d}")
    return max(fs, key=lambda f: int(re.findall(r"(\d+)\.h5", os.path.basename(f))[-1]))


def columns(x, y, xmin, xmax, wall_tol):
    """同一 x の**壁に接続した**法線カラムに束ねる (平板・構造格子前提)。

    最下点が壁 (y=0) から `wall_tol` 以内にあるカラムだけ採用する。
    これを見ないと、壁から浮いた点群 (例 x=0.433 で min(y)=0.086) まで
    カラム扱いしてしまい、$u_\tau$ 再構成が壊れる。
    """
    xr = np.round(x, 7)
    out, skipped = [], []
    for xc in np.unique(xr):
        if xc <= 1e-9 or (xmin is not None and xc < xmin) or (xmax is not None and xc > xmax):
            continue
        i = np.where(xr == xc)[0]
        if len(i) < 5:
            continue
        i = i[np.argsort(y[i])]
        if y[i[0]] > wall_tol:
            skipped.append((float(xc), float(y[i[0]])))
            continue
        out.append((xc, i))
    return out, skipped


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--mode", required=True, choices=["resolved", "wf"])
    ap.add_argument("--yrep", type=float, default=None,
                    help="resolved モードで比較する代表点高さ [m] (wf モードでは第一内層を自動使用)")
    ap.add_argument("--xmin", type=float, default=0.05)
    ap.add_argument("--xmax", type=float, default=None)
    ap.add_argument("--laws", default=",".join(LAWS))
    ap.add_argument("--wall-tol", type=float, default=1e-6,
                    help="カラム最下点がこの高さ以内なら壁接続とみなす [m]")
    ap.add_argument("--ypmin", type=float, default=None, help="この y+ 未満のカラムを除外")
    ap.add_argument("--ypmax", type=float, default=None, help="この y+ を超えるカラムを除外")
    a = ap.parse_args()
    laws = a.laws.split(",")

    h = h5py.File(latest(a.run, "res_[0-9]*.h5"), "r")
    C = h["MESH/COORD"][:].reshape(-1, 3)
    ro = h["VALUE/ro"][:]
    ux = h["VALUE/roUx"][:] / ro
    mu = h["VALUE/vis_lam"][:]
    x, y = C[:, 0], C[:, 1]
    cols, skipped = columns(x, y, a.xmin, a.xmax, a.wall_tol)
    if not cols:
        raise SystemExit("カラムが作れない (xmin/xmax/--wall-tol を確認)")

    print(f"# run={a.run}  mode={a.mode}  カラム {len(cols)} 本  x=[{a.xmin},{a.xmax}]")
    print(f"# 定数: Reichardt kappa=0.41 + 7.8 項 / Spalding,log kappa=0.41, B=5.0")
    print(f"# 平板前提: 壁法線 n=(0,1)、接線速度 = |U_x|")
    if skipped:
        print(f"# 壁未接続で除外したカラム {len(skipped)} 本 "
              f"(例 x={skipped[0][0]:.6f} min(y)={skipped[0][1]:.4g})")

    rows = {L: [] for L in laws}
    ref, xs, yps = [], [], []

    if a.mode == "resolved":
        if a.yrep is None:
            raise SystemExit("--yrep が要る (壁関数 run の代表点高さを渡す)")
        print(f"# 内部基準 u_tau = sqrt(nu*dU/dy|_w) (壁第一節点の分子勾配)、y_rep={a.yrep:.4e}")
        for xc, i in cols:
            j = i[1] if y[i[0]] < 1e-12 else i[0]
            nu = mu[j] / ro[j]
            ut_ref = np.sqrt(nu * abs(ux[j]) / y[j])
            Ut = np.interp(a.yrep, y[i], np.abs(ux[i]))          # 構造カラム上の線形補間
            nu_r = np.interp(a.yrep, y[i], mu[i] / ro[i])
            ref.append(ut_ref); xs.append(xc); yps.append(ut_ref * a.yrep / nu_r)
            for L in laws:
                rows[L].append(solve_utau(Ut, a.yrep, nu_r, L))
    else:
        hw = h5py.File(latest(a.run, "res_wall_*.h5"), "r")
        xw = hw["MESH/COORD"][:].reshape(-1, 3)[:, 0]
        uw = hw["VALUE/utau"][:]
        o = np.argsort(xw)
        print(f"# 基準 = 本体 utau (現行 Reichardt)。まず再構成検証を行う")
        rec = []
        for xc, i in cols:
            j = i[1] if y[i[0]] < 1e-12 else i[0]                # 代表点 = 第一内層ノード
            nu = mu[j] / ro[j]
            Ut, yj = abs(ux[j]), y[j]
            ut_soln = np.interp(xc, xw[o], uw[o])
            ref.append(ut_soln); xs.append(xc); yps.append(ut_soln * yj / nu)
            rec.append(solve_utau(Ut, yj, nu, "reichardt"))
            for L in laws:
                rows[L].append(solve_utau(Ut, yj, nu, L))
        rec = np.array(rec)

    ref = np.array(ref); yps = np.array(yps); xs = np.array(xs)
    keep = np.ones(len(ref), bool)
    if a.ypmin is not None: keep &= yps >= a.ypmin
    if a.ypmax is not None: keep &= yps <= a.ypmax
    if keep.sum() < 3:
        raise SystemExit(f"y+ フィルタ後のカラムが {keep.sum()} 本しかない")
    if keep.sum() != len(ref):
        print(f"# y+ フィルタ [{a.ypmin},{a.ypmax}]: {len(ref)} → {keep.sum()} カラム "
              f"(x=[{xs[keep].min():.3f},{xs[keep].max():.3f}])")
    ref, yps, xs = ref[keep], yps[keep], xs[keep]
    rows = {L: np.array(rows[L])[keep] for L in laws}
    print(f"# y+ 範囲: {yps.min():.1f}–{yps.max():.1f}")
    if a.mode == "wf":
        rec = rec[keep]
        g = ref > 0
        rel = np.abs(rec[g] - ref[g]) / ref[g]
        print(f"# ★ 再構成検証 (y+ フィルタ後): 本体 utau との相対差 "
              f"median={np.median(rel):.3e} max={rel.max():.3e}")
        if rel.max() > 0.05:
            print("#   → 5% 超。定義か物性の取り違えがあるので以降の比は信用しないこと。")
    print()
    base = "内部基準 (分子勾配)" if a.mode == "resolved" else "本体 utau (Reichardt)"
    print(f'{"壁法則":12s} {"u_tau/基準":>12s} {"tau_w/基準":>12s}   ({base} 比)')
    print(f'{"":12s} {"中央値":>12s} {"中央値":>12s}   min / max (u_tau 比)')
    for L in laws:
        r = np.array(rows[L]) / ref
        ru = float(np.median(r))
        print(f"{L:12s} {ru:12.4f} {ru**2:12.4f}   {r.min():.4f} / {r.max():.4f}")
    print()
    print(f'{"x":>7s} {"y+":>7s} {"基準u_tau":>10s} ' + " ".join(f"{L[:9]:>10s}" for L in laws))
    step = max(1, len(xs) // 6)
    for k in range(0, len(xs), step):
        print(f"{xs[k]:7.3f} {yps[k]:7.2f} {ref[k]:10.4f} "
              + " ".join(f"{rows[L][k]:10.4f}" for L in laws))


if __name__ == "__main__":
    main()
