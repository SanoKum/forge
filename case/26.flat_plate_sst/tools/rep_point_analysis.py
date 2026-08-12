#!/usr/bin/env python3
"""node 壁関数の代表点 (Normal_Neighbor) 診断の後処理 — 正式ツール。

README/plan に載る代表点まわりの数値はすべて本ツールで再生成できること。
入力は `FORGE_WF_REP_DIAG=1` で走らせた run の `res_*.h5`
(`rep_id`/`rep_y`/`rep_dist`/`rep_cos`/`rep_toff`/`rep_wdratio`/`rep_nx,ny,nz`)。

## できること

  geom    : 選択の幾何 (best_cos / 接線オフセット / wall_dist[irep] 比) の統計
  split   : **端点比の要因分解** (下記の注意を読むこと)
  cf      : 凍結場カウンターファクト — 法線補間速度を壁法則へ入れて u_tau, tau_w を再計算

## ★ split の解釈上の注意 (レビュー指摘 2026-08-13)

出力する 3 つの比の積が元の比に一致するのは、**中間状態を定義した比の積だから代数的に自明**で
あり、**連成解に対する因果寄与の一意な分解ではない**。あくまで「端点比を、選択 / 速度 / 壁法則
の 3 段に分けて眺めるための分解」として読むこと。

**壁法則へ入れる物性の固定 (曖昧さを残さないため明示)**:
  - `u_tau(粗場)` と `u_tau(基準速度から)` は **どちらも粗場側の nu と y_rep** を使う
    (差が **速度だけ** になるようにするため)。
  - `u_tau(壁解像真値)` は基準 run 自身の分子勾配 tau_w = mu*U_1/y_1 から求める
    (基準解自身の壁応力なので、粗場の nu/y は使わない)。
  → したがって「粗メッシュ速度場」の段は **速度のみ** の寄与だが、「壁法則」の段には
    物性・y の定義差と基準解自身の誤差が混ざる。

## 法線補間の規約

壁ノード $x_W$ から診断出力の **外向き法線** `rep_n` の逆向きに `rep_y` だけ進んだ点で、
逆距離重み (k 近傍, 既定 k=4) により速度を補間する。**外挿はしない** (対象点が凸包外なら
最近傍側に寄るだけなので、`--k` を変えて感度を見ること)。

使い方:
  python3 tools/rep_point_analysis.py geom  <diag_run> [--xmin --xmax]
  python3 tools/rep_point_analysis.py split <diag_run> --ref <ref_run> [--stations ...]
  python3 tools/rep_point_analysis.py cf    <diag_run> [--xmin --xmax]
"""
import sys, os, re, glob, argparse
import numpy as np
import h5py

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KAPPA = 0.41


def uplus_reichardt(yp):
    return np.log(1.0 + KAPPA * yp) / KAPPA + 7.8 * (
        1.0 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-yp / 3.0))


def solve_utau(Ut, y, nu):
    """Reichardt 則の逆解き。forge の SST automatic wall treatment と同じ式。"""
    ut = np.sqrt(max(nu * Ut / max(y, 1e-30), 1e-30))
    for _ in range(200):
        ut = max(ut, 1e-12)
        f = Ut / ut - uplus_reichardt(ut * y / nu)
        d = 1e-6 * ut
        df = ((Ut / (ut + d) - uplus_reichardt((ut + d) * y / nu)) - f) / d
        if abs(df) < 1e-30:
            break
        step = f / df
        ut -= step
        if abs(step) < 1e-14 * max(ut, 1e-12):
            break
    return max(ut, 0.0)


def load(run, fname=None):
    d = run if os.path.isabs(run) else os.path.join(HERE, run)
    p = os.path.join(d, fname) if fname else max(
        [f for f in glob.glob(os.path.join(d, "res_*.h5"))
         if "wall" not in os.path.basename(f) and "outlet" not in os.path.basename(f)
         and os.path.basename(f) != "res_0.h5"],
        key=lambda f: int(re.findall(r"res_(\d+)\.h5", os.path.basename(f))[0]))
    h = h5py.File(p, "r")
    C = h["MESH/COORD"][:].reshape(-1, 3)
    n = h["VALUE/ro"].shape[0]
    if C.shape[0] != n:                       # cell モード
        cn = h["MESH/CONNE"][:]
        c = cn.reshape(n, cn.size // n)[:, 1:]
        if c.min() == 1:
            c = c - 1
        C = C[c].mean(axis=1)
    ro = h["VALUE/ro"][:]
    return dict(C=C, ro=ro, ux=h["VALUE/roUx"][:] / ro, uy=h["VALUE/roUy"][:] / ro,
                mu=h["VALUE/vis_lam"][:], h=h, file=os.path.basename(p))


def wall_nodes(d, xmin, xmax):
    if "rep_id" not in d["h"]["VALUE"]:
        raise SystemExit("rep_id が無い。FORGE_WF_REP_DIAG=1 で走らせた run が要る。")
    rid = d["h"]["VALUE/rep_id"][:]
    x = d["C"][:, 0]
    m = np.where(rid >= 0)[0]
    if xmin is not None:
        m = m[x[m] >= xmin]
    if xmax is not None:
        m = m[x[m] <= xmax]
    return m


def interp_speed(tree, ux, uy, P, k):
    dd, ii = tree.query(P, k=k)
    w = 1.0 / np.maximum(dd, 1e-12)
    w /= w.sum()
    return float(np.sum(w * np.hypot(ux[ii], uy[ii]))), float(np.min(dd))


def cmd_geom(a):
    d = load(a.run, a.file)
    m = wall_nodes(d, a.xmin, a.xmax)
    V = {k: d["h"]["VALUE/" + k][:] for k in
         ["rep_y", "rep_dist", "rep_cos", "rep_toff", "rep_wdratio"]}
    print(f"# run={a.run}  file={d['file']}  壁ノード {len(m)} 個  x=[{a.xmin},{a.xmax}]")
    print(f'{"量":26s} {"min":>12s} {"median":>12s} {"max":>12s}')
    for k, nm in [("rep_y", "法線距離 y"), ("rep_dist", "実距離 |d|"), ("rep_cos", "best_cos"),
                  ("rep_toff", "接線オフセット"), ("rep_wdratio", "wall_dist[irep]/y")]:
        v = V[k][m]
        print(f'{nm:26s} {v.min():12.4e} {np.median(v):12.4e} {v.max():12.4e}')
    r = V["rep_toff"][m] / np.maximum(V["rep_y"][m], 1e-30)
    print(f'{"接線オフセット/y":26s} {r.min():12.4e} {np.median(r):12.4e} {r.max():12.4e}')
    print("\n注: wall_dist[irep]/y は壁距離場と法線射影 y の比。wall_dist 自身も離散量なので")
    print("    「真の壁距離」と断定しない。")


def cmd_split(a):
    """端点比の 3 段分解 (§3.3)。物性は上の docstring の規約に従う。"""
    dc = load(a.run, a.file)
    dr = load(a.ref, a.ref_file)
    from scipy.spatial import cKDTree
    tc = cKDTree(dc["C"][:, :2])
    tr = cKDTree(dr["C"][:, :2])
    rid = dc["h"]["VALUE/rep_id"][:].astype(int)
    ry = dc["h"]["VALUE/rep_y"][:]
    nx = dc["h"]["VALUE/rep_nx"][:]
    ny = dc["h"]["VALUE/rep_ny"][:]
    m = wall_nodes(dc, a.xmin, a.xmax)
    xw = dc["C"][:, 0]
    print(f"# 粗場 run={a.run} ({dc['file']})   基準 run={a.ref} ({dr['file']})")
    print(f"# 物性規約: u_tau(粗場) と u_tau(基準速度から) は **粗場の nu と y_rep** を使う")
    print(f"#           u_tau(基準真値) は基準 run の分子勾配 mu*U_1/y_1 から")
    print(f"# 法線補間: 逆距離重み k={a.k}, 外挿なし")
    print(f'{"x":>8s} {"y_rep":>10s} {"U(irep)":>10s} {"U(粗場法線)":>11s} {"U(基準法線)":>11s} '
          f'{"ut粗場":>8s} {"ut基準速度":>10s} {"ut基準真値":>10s} {"選択比":>7s} {"速度比":>7s} {"則比":>7s}')
    for xt in a.stations:
        j = m[np.argmin(np.abs(xw[m] - xt))]
        I = rid[j]
        y = ry[j]
        nin = np.array([-nx[j], -ny[j]])
        P = dc["C"][j, :2] + nin * y
        U1 = float(np.hypot(dc["ux"][I], dc["uy"][I]))
        U2, _ = interp_speed(tc, dc["ux"], dc["uy"], P, a.k)
        U3, _ = interp_speed(tr, dr["ux"], dr["uy"], P, a.k)
        nu = dc["mu"][I] / dc["ro"][I]                       # ★ 粗場側の nu で固定
        ut1 = solve_utau(U1, y, nu)
        ut3 = solve_utau(U3, y, nu)                          # ★ 同じ nu, y で速度だけ差替え
        # 基準の真値: 基準 run の壁法線第一点の分子勾配
        dd, ii = tr.query(dc["C"][j, :2], k=1)
        # 基準側の同 x の壁法線カラムから第一内点を拾う
        xr = dr["C"][:, 0]
        yr = dr["C"][:, 1]
        colm = np.abs(xr - dc["C"][j, 0]) < 2e-4
        idx = np.where(colm)[0]
        idx = idx[np.argsort(yr[idx])]
        i1 = idx[1] if yr[idx[0]] < 1e-12 else idx[0]
        ut_true = np.sqrt(dr["mu"][i1] * np.hypot(dr["ux"][i1], dr["uy"][i1])
                          / max(yr[i1], 1e-30) / dr["ro"][i1])
        print(f'{dc["C"][j,0]:8.4f} {y:10.3e} {U1:10.3f} {U2:11.3f} {U3:11.3f} '
              f'{ut1:8.4f} {ut3:10.4f} {ut_true:10.4f} '
              f'{U1/U2:7.4f} {ut3/ut1:7.4f} {ut_true/ut3:7.4f}')
    print("\n注: 「選択比 x 速度比 x 則比」が元の比に一致するのは **中間状態を定義した比の積**")
    print("    だからで、代数的に自明。連成解に対する因果寄与の一意な分解ではない。")


def cmd_cf(a):
    """凍結場カウンターファクト: 法線補間速度を壁法則へ入れて u_tau, tau_w を再計算。"""
    d = load(a.run, a.file)
    from scipy.spatial import cKDTree
    t = cKDTree(d["C"][:, :2])
    rid = d["h"]["VALUE/rep_id"][:].astype(int)
    ry = d["h"]["VALUE/rep_y"][:]
    nx = d["h"]["VALUE/rep_nx"][:]
    ny = d["h"]["VALUE/rep_ny"][:]
    m = wall_nodes(d, a.xmin, a.xmax)
    xw = d["C"][:, 0]
    print(f"# run={a.run}  file={d['file']}  **凍結場** (場は再計算しない)")
    print(f"# 現行: 代表点 irep の速度を壁法則へ / 反実: 同じ y の法線上へ補間した速度を壁法則へ")
    print(f"# 物性は同じ (粗場の nu, y_rep)。法線補間は逆距離重み k={a.k}、外挿なし")
    rows = []
    for j in m:
        I = rid[j]
        y = ry[j]
        nu = d["mu"][I] / d["ro"][I]
        U1 = float(np.hypot(d["ux"][I], d["uy"][I]))
        P = d["C"][j, :2] + np.array([-nx[j], -ny[j]]) * y
        U2, _ = interp_speed(t, d["ux"], d["uy"], P, a.k)
        ut1 = solve_utau(U1, y, nu)
        ut2 = solve_utau(U2, y, nu)
        rows.append((xw[j], U1, U2, ut1, ut2, d["ro"][I]))
    R = np.array(rows)
    tau1 = R[:, 5] * R[:, 3] ** 2
    tau2 = R[:, 5] * R[:, 4] ** 2
    print(f'{"":24s} {"min":>10s} {"median":>10s} {"max":>10s}')
    print(f'{"U(法線)/U(irep)":24s} {(R[:,2]/R[:,1]).min():10.4f} '
          f'{np.median(R[:,2]/R[:,1]):10.4f} {(R[:,2]/R[:,1]).max():10.4f}')
    print(f'{"u_tau(法線)/u_tau(現行)":24s} {(R[:,4]/R[:,3]).min():10.4f} '
          f'{np.median(R[:,4]/R[:,3]):10.4f} {(R[:,4]/R[:,3]).max():10.4f}')
    print(f'{"tau_w(法線)/tau_w(現行)":24s} {(tau2/tau1).min():10.4f} '
          f'{np.median(tau2/tau1):10.4f} {(tau2/tau1).max():10.4f}')
    print(f'\n面平均 tau_w: 現行 {tau1.mean():.2f} Pa -> 反実 {tau2.mean():.2f} Pa '
          f'(比 {tau2.mean()/tau1.mean():.4f})')
    print("注: これは **凍結場** の即時効果であり、実際に代表点選択を変えれば場が再平衡する。")
    print("    最終解の予測ではない (親 plan の E3 で凍結場見積り 0.484 に対し実測 0.513 だった前例)。")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["geom", "split", "cf"])
    ap.add_argument("run")
    ap.add_argument("--file", default=None)
    ap.add_argument("--ref", default=None, help="split: 壁解像基準 run")
    ap.add_argument("--ref-file", default=None)
    ap.add_argument("--xmin", type=float, default=None)
    ap.add_argument("--xmax", type=float, default=None)
    ap.add_argument("--stations", type=float, nargs="+", default=[0.45, 0.60, 0.70, 0.80])
    ap.add_argument("--k", type=int, default=4, help="逆距離重み補間の近傍点数")
    a = ap.parse_args()
    if a.cmd == "geom":
        cmd_geom(a)
    elif a.cmd == "split":
        if not a.ref:
            raise SystemExit("split には --ref (壁解像基準 run) が要る")
        cmd_split(a)
    else:
        cmd_cf(a)


if __name__ == "__main__":
    main()
