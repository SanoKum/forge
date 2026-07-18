#!/usr/bin/env python3
"""input h5 の VALUE を「一様低マッハ流 + 圧力のみ市松摂動」に書き換える (L1 診断 IC)。

検証ラダー L1 (notes/investigations/convection-central-scheme-oscillation-control.md §10):
  ρ=ρ0, u=(u0,0,0) 一様、P(i,j,k) = P0 [1 + ε (-1)^(i+j+k)]  (2Δ odd-even モード)
純 KEEP の中心平均はこのモードを相殺して見えない (null-mode) ため RHS が機械ゼロになる、
を実証し、keepDissType=1 (scalar ES 散逸) で減衰することを確認する。

- パリティはセル重心から算出 (セル順序に依存しない)。一様直交 box 前提。
- 保存量で書く: ro, roUx..roUz, roe = P/(γ-1) + ½ρ|u|²。roK/roOmega は 0。
usage: make_checkerboard_ic.py INPUT.h5 [--eps 1e-3] [--u0 0.1] [--gamma 1.4]
"""
import argparse
import numpy as np, h5py, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "solver_density_cuda", "tools"))
from res_h5_to_vtu import parse_conne

ap = argparse.ArgumentParser()
ap.add_argument("h5")
ap.add_argument("--eps", type=float, default=1e-3)
ap.add_argument("--u0", type=float, default=0.1)
ap.add_argument("--gamma", type=float, default=1.4)
a = ap.parse_args()

with h5py.File(a.h5, "r+") as f:
    coord = np.array(f["MESH/COORD"]).reshape(-1, 3)
    nc = f["VALUE/ro"].shape[0]
    conn, offs, _ = parse_conne(np.array(f["MESH/CONNE"]), nc)
    cc = np.zeros((nc, 3)); s = 0
    for i, o in enumerate(offs):
        cc[i] = coord[conn[s:o]].mean(axis=0); s = o

    # 一様格子前提で重心→整数インデックス (方向ごとに一意重心座標の rank)
    idx = np.zeros((nc, 3), dtype=np.int64)
    for d in range(3):
        u = np.unique(np.round(cc[:, d], 9))
        idx[:, d] = np.searchsorted(u, np.round(cc[:, d], 9))
    parity = np.where((idx.sum(axis=1) % 2) == 0, 1.0, -1.0)

    P0 = 1.0 / a.gamma; ro0 = 1.0
    P = P0 * (1.0 + a.eps * parity)
    roe = P / (a.gamma - 1.0) + 0.5 * ro0 * a.u0**2

    f["VALUE/ro"][:]   = np.float32(ro0)
    f["VALUE/roUx"][:] = np.float32(ro0 * a.u0)
    f["VALUE/roUy"][:] = np.float32(0.0)
    f["VALUE/roUz"][:] = np.float32(0.0)
    f["VALUE/roe"][:]  = roe.astype(np.float32)
    for k in ("roK", "roOmega"):
        if k in f["VALUE"]:
            f["VALUE"][k][:] = np.float32(0.0)

print(f"[checkerboard IC] nc={nc} P0={P0:.6f} eps={a.eps} u0={a.u0} "
      f"(+cells={int((parity>0).sum())}, -cells={int((parity<0).sum())})")
