#!/usr/bin/env python3
"""段差角1セル下流 (x≈2.2) の縦列で P の odd-even (市松) 振幅を測る。

低マッハ checkerboard は段差リップ直後の再循環域 (y<1) に局在するため、その縦列で
y 方向に並べた P の 2 次差分 (odd-even 成分) の rms を「市松振幅 [Pa]」とする。
notes/investigations/backstep-lowmach-checkerboard.md (2) の指標を再実装したもの。

usage: oddeven_corner.py res_XXXX.h5 [--x 2.2] [--ymax 1.0]
"""
import sys, argparse
import h5py, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("h5")
ap.add_argument("--x", type=float, default=2.2, help="対象縦列の x 目標値")
ap.add_argument("--ymax", type=float, default=1.0, help="再循環域 (リップ下) の y 上限")
a = ap.parse_args()

f = h5py.File(a.h5, "r")
P = f["VALUE/P"][:]
conn = f["MESH/CONNE"][:].reshape(-1, 5)[:, 1:5]
coord = f["MESH/COORD"][:].reshape(-1, 3)
cx = coord[conn, 0].mean(1)
cy = coord[conn, 1].mean(1)

# x が最も近い縦列を選ぶ
xcol = np.unique(np.round(cx, 3))
xsel = xcol[np.argmin(np.abs(xcol - a.x))]
m = (np.abs(cx - xsel) < 1e-3) & (cy < a.ymax)
ys = cy[m]; Ps = P[m]
order = np.argsort(ys)
ys = ys[order]; Ps = Ps[order]

# 2 次差分 = odd-even 成分。両端を除く内部点で評価。
oe = Ps[1:-1] - 0.5 * (Ps[:-2] + Ps[2:])
oe_rms = float(np.sqrt(np.mean(oe**2)))
oe_ptp = float(oe.max() - oe.min())
print(f"file        : {a.h5}")
print(f"column x    : {xsel:.3f}  (target {a.x}),  y<{a.ymax}: {m.sum()} cells")
print(f"P range     : {Ps.min():.2f} .. {Ps.max():.2f} Pa")
print(f"P odd-even  : rms={oe_rms:.2f} Pa,  peak-to-peak={oe_ptp:.2f} Pa")
