#!/usr/bin/env python3
# 均一メッシュ (1D_shock_tube.h5, x in [0,1]) 上に N2 衝撃管初期値を生成する。
#   x < xdia: 左 (高 P,T, hot driver)、x >= xdia: 右 (低 P,T, cold driven)、u=0。
#   mode=tp:  roe = ro*(e_NASA(T)+ek)  [NASA-9 絶対基準]
#   mode=cpg: roe = P/(gamma-1)        [gamma=1.4 一定]
#
# 使い方: python3 gen_shocktube_ic.py <mesh.h5> <out.h5> <tp|cpg>
import sys, numpy as np, h5py, shutil

mesh = sys.argv[1]; dst = sys.argv[2]; mode = sys.argv[3]  # 'tp' or 'cpg'

# --- 初期不連続条件 (高温で cp(T) 変化が顕著な N2 衝撃管) ---
xdia = 0.5
PL, TL = 1.0e6, 2000.0   # 左 (hot driver)
PR, TR = 1.0e5,  400.0   # 右 (cold driven)

# --- N2 物性 (NASA-9, CEA 係数) ---
RU = 8.314462618; MW = 0.0280134; R = RU/MW; TMID = 1000.0; GAM = 1.4
LO = np.array([2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,
               1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1])
HI = np.array([5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,
               1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1])

def hmass(T):
    """NASA-9 比エンタルピー h(T) [J/kg] (絶対基準)。"""
    T = np.atleast_1d(T).astype(float)
    a = np.where(T[..., None] < TMID, LO, HI); Ti = 1/T; lnT = np.log(T)
    hRT = (-a[..., 0]*Ti*Ti + a[..., 1]*lnT*Ti + a[..., 2] + a[..., 3]*T/2
           + a[..., 4]*T**2/3 + a[..., 5]*T**3/4 + a[..., 6]*T**4/5 + a[..., 7]*Ti)
    return RU*T*hRT/MW

shutil.copy(mesh, dst)
with h5py.File(dst, 'r+') as f:
    co = f['MESH/COORD'][:].reshape(-1, 3); c = f['MESH/CONNE'][:].reshape(-1, 9)
    cx = co[c[:, 1:9], 0].mean(1); n = len(cx)
    P = np.where(cx < xdia, PL, PR); T = np.where(cx < xdia, TL, TR)
    ro = P/(R*T)
    if mode == 'tp':
        e = hmass(T) - R*T          # e = h - R*T (ideal gas), NASA 絶対基準
        roe = ro*e
    else:
        roe = P/(GAM-1.0)
    g = f['VALUE'] if 'VALUE' in f else f.create_group('VALUE')

    def setv(k, v):
        v = np.asarray(v, np.float32)
        if 'VALUE/'+k in f:
            f['VALUE/'+k][...] = v
        else:
            g.create_dataset(k, data=v)
    setv('ro', ro); setv('roUx', np.zeros(n)); setv('roUy', np.zeros(n)); setv('roUz', np.zeros(n))
    setv('roe', roe); setv('T', T); setv('P', P)
    setv('Ux', np.zeros(n)); setv('Uy', np.zeros(n)); setv('Uz', np.zeros(n))
    if 'VALUE/wall_dist' not in f:
        g.create_dataset('wall_dist', data=np.ones(n, np.float32))
    print(f"[{mode}] n={n} xdia={xdia} "
          f"L(P={PL:.3e},T={TL},ro={ro[cx<xdia][0]:.4f}) "
          f"R(P={PR:.3e},T={TR},ro={ro[cx>=xdia][0]:.4f}) "
          f"roe[{roe.min():.3e},{roe.max():.3e}] -> {dst}")
