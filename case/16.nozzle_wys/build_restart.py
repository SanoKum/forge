#!/usr/bin/env python3
"""収束済み res_*.h5 の場を別メッシュ h5 の /VALUE に焼き込んで引き継ぎ初期場を作る。

usage:
  build_restart.py SRC_RES.h5 SRC_MESH.h5 DST_MESH.h5  [--rok K --roomega W]

SRC_MESH と DST_MESH の /CELLS/centCoords を KD-tree で対応付け (同一形状なら厳密一致)、
ro,roUx,roUy,roUz,roe を写す。--rok/--roomega 指定時は RANS 用に roK,roOmega も一様で書く。
DST_MESH h5 はその場で書き換える (上書き)。
"""
import sys, argparse
import numpy as np
import h5py
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("src_res"); ap.add_argument("src_mesh"); ap.add_argument("dst_mesh")
ap.add_argument("--rok", type=float, default=None)
ap.add_argument("--roomega", type=float, default=None)
a = ap.parse_args()

with h5py.File(a.src_mesh, "r") as f:
    src_cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
with h5py.File(a.src_res, "r") as f:
    vals = {k: f["/VALUE/"+k][:] for k in ("ro", "roUx", "roUy", "roUz", "roe")}

tree = cKDTree(src_cc)
with h5py.File(a.dst_mesh, "r+") as f:
    dst_cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
    d, idx = tree.query(dst_cc, k=1)
    print(f"match: nDst={len(dst_cc)} maxDist={d.max():.3e} meanDist={d.mean():.3e}")
    for k, v in vals.items():
        f["/VALUE/"+k][:] = v[idx]
    if a.rok is not None:
        ro = vals["ro"][idx]
        f["/VALUE/roK"][:] = ro * a.rok          # roK = ro * k   (k = TKE [m^2/s^2])
        f["/VALUE/roOmega"][:] = ro * a.roomega  # roOmega = ro * omega [1/s]
        print(f"set roK=ro*{a.rok}, roOmega=ro*{a.roomega}")
print("wrote restart field into", a.dst_mesh)
