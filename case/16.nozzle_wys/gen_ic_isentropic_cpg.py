#!/usr/bin/env python3
"""
quasi-1D 等エントロピー IC を CPG (定数 cp/γ) で生成する。
面積-マッハ関係 A(x)/A* から各セルの M(x) を解き (喉上流=亜音速, 下流=超音速)、
等エントロピー関係で P,T,ρ,u を与える。ρuA=const を満たすため startup ショックが出にくい
(README run_0097/0098 の非整合 IC への対策)。

使い方:
  python3 gen_ic_isentropic_cpg.py <mesh.h5> <dst_ic.h5> [Pt=101325] [Tt=293.15] [gamma=1.4] [cp=1039] [kIC=1.0] [omegaIC=1000]

mesh.h5 を dst_ic.h5 にコピーし /VALUE/{ro,roUx,roUy,roUz,roe} を等エントロピー場で上書き、
roK/roOmega があれば ρ·kIC / ρ·omegaIC を書く (wall_dist 等は保持)。
"""
import sys, shutil, numpy as np, h5py

mesh = sys.argv[1]
dst  = sys.argv[2]
Pt   = float(sys.argv[3]) if len(sys.argv) > 3 else 101325.0
Tt   = float(sys.argv[4]) if len(sys.argv) > 4 else 293.15
gam  = float(sys.argv[5]) if len(sys.argv) > 5 else 1.4
cp   = float(sys.argv[6]) if len(sys.argv) > 6 else 1039.0
kIC     = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
omegaIC = float(sys.argv[8]) if len(sys.argv) > 8 else 1000.0
R = cp * (gam - 1.0) / gam

with h5py.File(mesh, 'r') as h:
    cc    = np.array(h['/CELLS/centCoords']).reshape(-1, 3)
    coord = np.array(h['/MESH/COORD']).reshape(-1, 3)
nC = cc.shape[0]
xc = cc[:, 0]
xn, yn = coord[:, 0], coord[:, 1]

# --- 面積分布 A(x): ノード y 範囲 (2D 単位奥行 → A ∝ 高さ) ---
nb = 250
xe = np.linspace(xn.min(), xn.max(), nb + 1)
xm = 0.5 * (xe[:-1] + xe[1:])
height = np.full(nb, np.nan)
for i in range(nb):
    m = (xn >= xe[i]) & (xn < xe[i + 1])
    if m.sum() > 1:
        height[i] = yn[m].max() - yn[m].min()
good = ~np.isnan(height)
height = np.interp(xm, xm[good], height[good])
height = np.convolve(height, np.ones(5) / 5, mode='same')
Astar = height.min()
x_th = xm[np.argmin(height)]
A_cell = np.interp(xc, xm, height)
AR = np.maximum(A_cell / Astar, 1.0 + 1e-12)

# --- 面積-マッハ関係を分枝別に二分法で解く ---
expo = (gam + 1.0) / (2.0 * (gam - 1.0))
def AR_of_M(M):
    return (1.0 / M) * ((2.0 / (gam + 1.0)) * (1.0 + 0.5 * (gam - 1.0) * M * M)) ** expo

sup = xc > x_th
M = np.empty(nC)
# 亜音速分枝 (M∈(0,1], AR は M に対し単調減少)
lo = np.full(nC, 1e-5); hi = np.full(nC, 1.0)
for _ in range(100):
    mid = 0.5 * (lo + hi)
    f = AR_of_M(mid) - AR            # 減少関数 → f>0 なら M を増やす
    hi = np.where(f > 0, hi, mid)
    lo = np.where(f > 0, mid, lo)
M_sub = 0.5 * (lo + hi)
# 超音速分枝 (M∈[1,Mmax], AR は M に対し単調増加)
lo = np.full(nC, 1.0); hi = np.full(nC, 60.0)
for _ in range(100):
    mid = 0.5 * (lo + hi)
    f = AR_of_M(mid) - AR            # 増加関数 → f>0 なら M を減らす
    hi = np.where(f > 0, mid, hi)
    lo = np.where(f > 0, lo, mid)
M_sup = 0.5 * (lo + hi)
M = np.where(sup, M_sup, M_sub)

# --- 等エントロピー関係 ---
T  = Tt / (1.0 + 0.5 * (gam - 1.0) * M * M)
P  = Pt * (T / Tt) ** (gam / (gam - 1.0))
ro = P / (R * T)
a  = np.sqrt(gam * R * T)
u  = M * a
roe = P / (gam - 1.0) + 0.5 * ro * u * u

print(f"R={R:.3f}  x_throat={x_th*1000:.2f}mm  A/A*[{AR.min():.3f},{AR.max():.3f}]")
print(f"M[{M.min():.3f},{M.max():.3f}]  T[{T.min():.1f},{T.max():.1f}]K  P[{P.min():.0f},{P.max():.0f}]Pa")
print(f"ro[{ro.min():.4f},{ro.max():.4f}]  u[{u.min():.1f},{u.max():.1f}]m/s")
# 質量流束一様性チェック (ρuA, A∝height) — 喉下流の超音速部
mdot = ro * u * np.interp(xc, xm, height)
ds = mdot[sup]
print(f"rho*u*A spread (supersonic): mean={ds.mean():.4e} std/mean={ds.std()/ds.mean():.2%}")

shutil.copy(mesh, dst)
with h5py.File(dst, 'r+') as h:
    h['/VALUE/ro'][...]   = ro.astype(np.float32)
    h['/VALUE/roUx'][...] = (ro * u).astype(np.float32)
    h['/VALUE/roUy'][...] = np.zeros(nC, np.float32)
    h['/VALUE/roUz'][...] = np.zeros(nC, np.float32)
    h['/VALUE/roe'][...]  = roe.astype(np.float32)
    if '/VALUE/roK' in h:
        h['/VALUE/roK'][...] = (ro * kIC).astype(np.float32)
    if '/VALUE/roOmega' in h:
        h['/VALUE/roOmega'][...] = (ro * omegaIC).astype(np.float32)
print(f"wrote {dst}")
