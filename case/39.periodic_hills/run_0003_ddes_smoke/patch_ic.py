#!/usr/bin/env python3
"""run_0002 の RANS 収束場 → DDES 初期場パッチ。

1. run_0002 最終 res の場を入力 h5 (同一メッシュ) へ移植
2. バルク運動量を目標 <ρu_x>_V=44.1248 へリスケール (bodyForceCtrl の deadbeat 初撃回避)
3. 剥離域〜チャネル下半分 (y_w < y < 2.5h) に ±6% の 3D 撹乱 (KH 3 次元化の種)
4. roe を新しい運動エネルギーと整合させる
"""
import glob
import re
import sys

import h5py
import numpy as np

sys.path.insert(0, '../mesh')
from make_hill_mesh import wall_y, H, LX, LY_TOP

TARGET = 44.1248  # <ro*u_x>_V [kg/m^2/s]

fs = [f for f in glob.glob('../run_0002_rans_dev/res_*.h5') if re.search(r'res_\d+\.h5$', f)]
fs.sort(key=lambda s: int(re.findall(r'(\d+)\.h5$', s)[0]))
src_path = fs[-1]
print("IC source:", src_path)

src = h5py.File(src_path, 'r')
dst = h5py.File('hill_des_160x100x60.h5', 'r+')
for k in ['ro', 'roUx', 'roUy', 'roUz', 'roe', 'roK', 'roOmega']:
    dst['VALUE/' + k][:] = src['VALUE/' + k][:]

xyz = dst['CELLS/centCoords'][:].reshape(-1, 3)
x, y = xyz[:, 0], xyz[:, 1]
vol = dst['CELLS/volume'][:]
ro = dst['VALUE/ro'][:]
roUx = dst['VALUE/roUx'][:]
roUy = dst['VALUE/roUy'][:]
roUz = dst['VALUE/roUz'][:]

# 旧運動エネルギー (roe から差し引くため)
ek_old = 0.5 * (roUx**2 + roUy**2 + roUz**2) / ro

# --- バルクリスケール ---
M0 = (roUx * vol).sum() / vol.sum()
scale = TARGET / M0
print(f"M/V: {M0:.3f} -> {TARGET} (scale {scale:.4f})")
roUx *= scale
roUy *= scale
roUz *= scale

# --- 3D 撹乱 (壁で 0 に減衰) ---
rng = np.random.default_rng(390003)
yw = np.interp(x, np.linspace(0, LX, 4001), [wall_y(xi) for xi in np.linspace(0, LX, 4001)])
eta = np.clip((y - yw) / (2.5 * H - yw), 0.0, 1.0)  # y_w..2.5h で 0..1
damp = np.where(y < 2.5 * H, np.sin(np.pi * np.clip((y - yw) / (2.5 * H - yw), 0, 1)), 0.0)
u_local = np.abs(roUx / ro)
amp = 0.06 * u_local * damp
roUx += ro * amp * rng.uniform(-1, 1, len(ro))
roUy += ro * amp * rng.uniform(-1, 1, len(ro))
roUz += ro * amp * rng.uniform(-1, 1, len(ro))

ek_new = 0.5 * (roUx**2 + roUy**2 + roUz**2) / ro
dst['VALUE/roUx'][:] = roUx
dst['VALUE/roUy'][:] = roUy
dst['VALUE/roUz'][:] = roUz
dst['VALUE/roe'][:] = dst['VALUE/roe'][:] - ek_old + ek_new  # 内部エネルギー保存で KE を整合

M1 = (roUx * vol).sum() / vol.sum()
print(f"final M/V = {M1:.3f} (target {TARGET}); perturbed nodes: {(damp > 0).sum()}")
src.close()
dst.close()
