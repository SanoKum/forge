#!/usr/bin/env python3
"""代表点の法線化カウンターファクト — **要素ベース補間**による正式版。

旧実装 (`rep_point_analysis.py cf` の IDW) は近傍数 k に強く依存し、case/40 の
tau_w 比が k=1/4/8/16 で 1.000/0.679/0.735/0.825 と動いた。したがって IDW の値は
**補間の産物**であって物理量ではない。本ツールはそれを次で置き換える:

  1. **メッシュ接続から対象点を含む要素を特定**する (VIZMESH の primal quad を
     2 三角形に分割し重心座標で包含判定)。
  2. **速度ベクトルを補間**してから、本体と同じ定義で接線成分を取る:
       U_t = |U - (U.n) n|           (ransWallFunction_d.cu と同一)
     旧実装は hypot(ux,uy) を補間しており、これは接線成分ではない。
  3. **包含要素が無ければエラー**にする (真に外挿しない)。IDW の「最近傍側に寄る」は
     外挿を防いだことにならない。
  4. **壁面積重み付き平均**を使う (壁ノードの算術平均ではない)。軸対称なら周長 2*pi*r を掛ける。
  5. まず **本体の utau (res_wall_*.h5) を再構成できるか検証** してから反実値を出す。

使い方:
  python3 tools/rep_normal_interp.py <diag_run> --wall-file res_wall_3_100.h5 \\
      [--xmin 0.010 --xmax 0.070] [--axisym]
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


class Mesh2D:
    """VIZMESH の primal 要素を 2 三角形へ分割し、重心座標で包含判定・補間する。"""

    def __init__(self, coord, conne, nDOF):
        self.P = coord[:, :2]
        stride = conne.size // (conne.size // 5) if False else 5
        e = conne.reshape(-1, 5)
        assert (e[:, 0] == 5).all(), "想定外の要素タイプ (quad=5 のみ対応)"
        q = e[:, 1:]
        if q.max() >= nDOF:
            raise SystemExit("CONNE の節点 index が DOF 数を超える (node モードの res か確認)")
        # quad -> 2 triangles
        self.tri = np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
        c = self.P[self.tri].mean(axis=1)
        from scipy.spatial import cKDTree
        self.tree = cKDTree(c)

    def locate(self, p, ncand=32):
        """p を含む三角形を返す。無ければ None (= 外挿になるので呼び出し側でエラー)。"""
        _, cand = self.tree.query(p, k=min(ncand, len(self.tri)))
        for t in np.atleast_1d(cand):
            a, b, c = self.P[self.tri[t]]
            v0, v1, v2 = b - a, c - a, p - a
            den = v0[0] * v1[1] - v1[0] * v0[1]
            if abs(den) < 1e-30:
                continue
            u = (v2[0] * v1[1] - v1[0] * v2[1]) / den
            v = (v0[0] * v2[1] - v2[0] * v0[1]) / den
            if u >= -1e-9 and v >= -1e-9 and u + v <= 1.0 + 1e-9:
                return self.tri[t], np.array([1.0 - u - v, u, v])
        return None, None

    def interp_vec(self, p, ux, uy):
        nodes, w = self.locate(p)
        if nodes is None:
            return None
        return np.array([float(np.dot(w, ux[nodes])), float(np.dot(w, uy[nodes]))])


def latest(run, pattern, explicit=None):
    d = run if os.path.isabs(run) else os.path.join(HERE, run)
    if explicit:
        return os.path.join(d, explicit)
    fs = [f for f in glob.glob(os.path.join(d, pattern))
          if "outlet" not in os.path.basename(f) and os.path.basename(f) != "res_0.h5"]
    if not fs:
        raise SystemExit(f"{pattern} が見つからない: {d}")
    return max(fs, key=lambda f: int(re.findall(r"(\d+)\.h5", os.path.basename(f))[-1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--file", default=None)
    ap.add_argument("--wall-file", default=None, help="res_wall_<physID>_<step>.h5 (utau 検証用)")
    ap.add_argument("--xmin", type=float, default=None)
    ap.add_argument("--xmax", type=float, default=None)
    ap.add_argument("--axisym", action="store_true", help="壁面積重みに周長 2*pi*r を掛ける")
    a = ap.parse_args()

    h = h5py.File(latest(a.run, "res_[0-9]*.h5", a.file), "r")
    C = h["MESH/COORD"][:].reshape(-1, 3)
    ro = h["VALUE/ro"][:]
    ux = h["VALUE/roUx"][:] / ro
    uy = h["VALUE/roUy"][:] / ro
    mu = h["VALUE/vis_lam"][:]
    if "rep_id" not in h["VALUE"]:
        raise SystemExit("rep_id が無い。FORGE_WF_REP_DIAG=1 の run が要る。")
    rid = h["VALUE/rep_id"][:].astype(int)
    ry = h["VALUE/rep_y"][:]
    nx = h["VALUE/rep_nx"][:]
    ny = h["VALUE/rep_ny"][:]
    mesh = Mesh2D(C, h["MESH/CONNE"][:], ro.shape[0])

    m = np.where(h["VALUE/rep_id"][:] >= 0)[0]
    x = C[:, 0]
    if a.xmin is not None:
        m = m[x[m] >= a.xmin]
    if a.xmax is not None:
        m = m[x[m] <= a.xmax]
    m = m[np.argsort(x[m])]

    # --- 検証: 本体の utau を再構成できるか -----------------------------------
    print(f"# run={a.run}  壁ノード {len(m)} 個  x=[{a.xmin},{a.xmax}]  axisym={a.axisym}")
    print(f"# U_t = |U - (U.n)n| (本体と同一定義)。補間は **要素ベース** (quad->2 triangle, 重心座標)")
    print(f"# 包含要素が無ければエラー (真に外挿しない)")
    ut_rec = np.empty(len(m))
    for i, j in enumerate(m):
        I = rid[j]
        n = np.array([nx[j], ny[j]])
        U = np.array([ux[I], uy[I]])
        Ut = np.linalg.norm(U - np.dot(U, n) * n)
        ut_rec[i] = solve_utau(Ut, ry[j], mu[I] / ro[I])
    if a.wall_file:
        hw = h5py.File(os.path.join(HERE, a.run, a.wall_file)
                       if not os.path.isabs(a.run) else os.path.join(a.run, a.wall_file), "r")
        xw = hw["MESH/COORD"][:].reshape(-1, 3)[:, 0]
        uw = hw["VALUE/utau"][:]
        o = np.argsort(xw)
        ref = np.interp(x[m], xw[o], uw[o])
        good = ref > 0
        rel = np.abs(ut_rec[good] - ref[good]) / np.maximum(ref[good], 1e-30)
        print(f"# ★ 再構成検証: 本体 utau との相対差  median={np.median(rel):.3e}  "
              f"max={rel.max():.3e}  ({good.sum()} 点)")
        if rel.max() > 0.05:
            print("#   → 5% を超える。定義か物性の取り違えがあるので反実値は信用しないこと。")
    else:
        print("# (--wall-file 未指定: 本体 utau との再構成検証をスキップ)")

    # --- 反実: 法線上へ要素ベース補間 -----------------------------------------
    ut_cf = np.full(len(m), np.nan)
    miss = 0
    for i, j in enumerate(m):
        I = rid[j]
        n = np.array([nx[j], ny[j]])
        P = C[j, :2] + (-n) * ry[j]
        V = mesh.interp_vec(P, ux, uy)
        if V is None:
            miss += 1
            continue
        Ut = np.linalg.norm(V - np.dot(V, n) * n)
        ut_cf[i] = solve_utau(Ut, ry[j], mu[I] / ro[I])
    if miss:
        print(f"# ! 包含要素が見つからない壁ノード {miss} 個 (外挿になるので除外)")
    ok = np.isfinite(ut_cf)
    if ok.sum() < 3:
        raise SystemExit("包含要素が見つかる点が少なすぎる")

    # --- 壁面積重み付き平均 ---------------------------------------------------
    xs = x[m][ok]
    ys = C[m, 1][ok]
    ds = np.empty(len(xs))
    ds[1:-1] = 0.5 * np.hypot(xs[2:] - xs[:-2], ys[2:] - ys[:-2])
    ds[0] = np.hypot(xs[1] - xs[0], ys[1] - ys[0])
    ds[-1] = np.hypot(xs[-1] - xs[-2], ys[-1] - ys[-2])
    w = ds * (2.0 * np.pi * np.maximum(ys, 1e-30) if a.axisym else 1.0)
    rod = ro[rid[m][ok]]
    t1 = rod * ut_rec[ok] ** 2
    t2 = rod * ut_cf[ok] ** 2
    r = ut_cf[ok] / ut_rec[ok]
    print()
    print(f'{"":26s} {"min":>10s} {"median":>10s} {"max":>10s}')
    print(f'{"u_tau(法線)/u_tau(現行)":26s} {r.min():10.4f} {np.median(r):10.4f} {r.max():10.4f}')
    print(f'{"tau_w 比":26s} {(t2/t1).min():10.4f} {np.median(t2/t1):10.4f} {(t2/t1).max():10.4f}')
    print(f'\n壁面積重み付き平均 tau_w: 現行 {np.sum(t1*w)/np.sum(w):.2f} Pa '
          f'-> 反実 {np.sum(t2*w)/np.sum(w):.2f} Pa  (比 {np.sum(t2*w)/np.sum(t1*w):.4f})')
    print("\n注: **凍結場**の即時効果。実際に選択を変えれば場が再平衡するので最終解の予測ではない。")
    print("    また、既往の 0.945 (y+1 基準比) と比較するには、同じ範囲・同じ重み・同じ tau_w 定義で")
    print("    基準側も評価し直す必要がある (単純な掛け算はできない)。")


if __name__ == "__main__":
    main()
