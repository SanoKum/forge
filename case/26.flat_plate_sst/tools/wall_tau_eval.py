#!/usr/bin/env python3
"""壁せん断応力の **単一評価関数** — 壁関数 run と壁解像 run を同一定義で比べる。

これまで case/40 の「y+1 基準比 0.945」と、代表点カウンターファクトの「0.610」は
**別々の評価関数**で出しており、単純に掛け算できなかった。本ツールは両者を
**同じ範囲・同じ重み・同じ tau_w 定義**に通すためのもの。

## 定義 (両モード共通)

`res_wall_<physID>_<step>.h5` の **`twall_x/y/z`** を使う。これは壁半割面の traction
(単位面積あたり) であり:

  - **壁解像 run** (`wallTreatmentSST=0`): 解像した粘性 traction = その解の真の壁応力。
  - **壁関数 run**: `viscousFlux_d.cu` が向きは解像値・大きさをモデル `tau_w=rho*u_tau^2`
    へ再スケールした値 = **実際に課しているモデル応力**。

したがって `twall` は「そのソルバ構成が壁に与えている応力」であり、両者の比較に使える。
(`utau` は壁解像 run では 0 なので使えない。)

**接線成分**を取る: `tau_w = |t - (t.n) n|`。法線 `n` は壁面メッシュの隣接節点から
局所接線を作って求める (診断出力に依存しないので任意の run に使える)。

## 重み

壁面に沿った弧長 `ds` (隣接節点の中点間) を重みとし、`--axisym` なら周長 `2*pi*r` を掛ける。
壁ノードの単純算術平均は使わない。

使い方:
  python3 tools/wall_tau_eval.py <run> [--wall-file res_wall_3_9000.h5] \\
      [--xmin 0.010 --xmax 0.070] [--axisym]
  python3 tools/wall_tau_eval.py <runA> <runB> ... --xmin .. --xmax .. --axisym
"""
import os, re, glob, argparse
import numpy as np
import h5py

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def latest_wall(run, explicit=None, phys_id=None):
    d = run if os.path.isabs(run) else os.path.join(HERE, run)
    if explicit:
        return os.path.join(d, explicit)
    pat = re.compile(r"^res_wall_(\d+)_(\d+)\.h5$")
    cand = []
    for f in glob.glob(os.path.join(d, "res_wall_*.h5")):
        m = pat.match(os.path.basename(f))
        if m and (phys_id is None or int(m.group(1)) == phys_id):
            cand.append((int(m.group(2)), f))
    if not cand:
        raise SystemExit(f"res_wall_*.h5 が無い: {d}")
    return max(cand)[1]


def eval_wall_tau(path, xmin, xmax, axisym):
    """壁面積重み付き平均の接線 traction [Pa] を返す。"""
    with h5py.File(path, "r") as h:
        C = h["MESH/COORD"][:].reshape(-1, 3)
        tx = h["VALUE/twall_x"][:]
        ty = h["VALUE/twall_y"][:]
        tz = h["VALUE/twall_z"][:] if "VALUE/twall_z" in h else np.zeros_like(tx)
    x, y = C[:, 0], C[:, 1]
    o = np.argsort(x)                       # 壁面に沿って並べる (単調な x を仮定)
    x, y, tx, ty, tz = x[o], y[o], tx[o], ty[o], tz[o]
    m = np.ones(len(x), bool)
    if xmin is not None:
        m &= x >= xmin
    if xmax is not None:
        m &= x <= xmax
    idx = np.where(m)[0]
    if len(idx) < 3:
        raise SystemExit(f"範囲内の壁ノードが {len(idx)} 個しかない")
    # 局所接線 -> 法線 (2D)
    tang = np.empty((len(x), 2))
    tang[1:-1] = np.column_stack([x[2:] - x[:-2], y[2:] - y[:-2]])
    tang[0] = [x[1] - x[0], y[1] - y[0]]
    tang[-1] = [x[-1] - x[-2], y[-1] - y[-2]]
    tang /= np.linalg.norm(tang, axis=1, keepdims=True)
    nrm = np.column_stack([-tang[:, 1], tang[:, 0]])
    t2 = np.column_stack([tx, ty])
    tn = np.sum(t2 * nrm, axis=1)
    tau = np.linalg.norm(t2 - tn[:, None] * nrm, axis=1)      # 接線成分の大きさ
    # 弧長重み
    ds = np.empty(len(x))
    ds[1:-1] = 0.5 * np.hypot(x[2:] - x[:-2], y[2:] - y[:-2])
    ds[0] = np.hypot(x[1] - x[0], y[1] - y[0])
    ds[-1] = np.hypot(x[-1] - x[-2], y[-1] - y[-2])
    w = ds * (2.0 * np.pi * np.maximum(y, 1e-30) if axisym else 1.0)
    return float(np.sum(tau[idx] * w[idx]) / np.sum(w[idx])), len(idx), (x[idx].min(), x[idx].max())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--wall-file", default=None)
    ap.add_argument("--phys-id", type=int, default=None)
    ap.add_argument("--xmin", type=float, default=None)
    ap.add_argument("--xmax", type=float, default=None)
    ap.add_argument("--axisym", action="store_true")
    a = ap.parse_args()
    print(f"# 定義: tau_w = |twall - (twall.n)n| (壁半割面 traction の接線成分)")
    print(f"#       壁解像 run では解像値 = 真の壁応力 / 壁関数 run ではモデル値に再スケール済み")
    print(f"# 重み: 壁面弧長 ds" + (" x 周長 2*pi*r (axisym)" if a.axisym else "") +
          f"   範囲 x=[{a.xmin},{a.xmax}]")
    print()
    print(f'{"run":46s} {"tau_w [Pa]":>12s} {"n":>5s} {"x範囲":>21s}')
    vals = []
    for r in a.runs:
        p = latest_wall(r, a.wall_file, a.phys_id)
        v, n, rng = eval_wall_tau(p, a.xmin, a.xmax, a.axisym)
        vals.append((r, v))
        print(f'{os.path.basename(r):46s} {v:12.2f} {n:5d}  [{rng[0]:.4f},{rng[1]:.4f}]  ({os.path.basename(p)})')
    if len(vals) >= 2:
        print()
        base = vals[-1][1]
        for r, v in vals[:-1]:
            print(f'  {os.path.basename(r)} / {os.path.basename(vals[-1][0])} = {v/base:.4f}')


if __name__ == "__main__":
    main()
