#!/usr/bin/env python3
# 2 化学種 (M2) 衝撃管初期値を均一メッシュ上に生成する。
#   x < 0.5: 左 (高 P,T, 化学種0)、x >= 0.5: 右 (低 P,T, 化学種1)、u=0。
#   roe は NASA-9 絶対基準: roe = ro*(e_gasX(T)+ek), ek=0。
#   roY0,roY1 は左右に分離 (左 = 純 species0, 右 = 純 species1) -> contact が組成界面。
#
# 使い方: python3 gen_shocktube_species_ic.py <mesh.h5> <out.h5> <gasL> <gasR>
#   例 (識別子トレーサ, 流れは単成分 N2 と一致):     ... N2 N2
#   例 (実 2 ガス O2 駆動 / N2 被駆動):               ... O2 N2
import sys, numpy as np, h5py, shutil

mesh = sys.argv[1]; dst = sys.argv[2]; gasL = sys.argv[3]; gasR = sys.argv[4]

xdia = 0.5
PL, TL = 1.0e6, 2000.0   # 左 (hot driver)
PR, TR = 1.0e5,  400.0   # 右 (cold driven)
RU = 8.314462618; TMID = 1000.0

# NASA-9 (CEA) DB: MW [kg/mol], low[9], high[9]  (thermo_d.cu と一致)
DB = {
 'N2': (0.0280134,
        np.array([2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1]),
        np.array([5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1])),
 'O2': (0.0319988,
        np.array([-3.425563420e4,4.847000970e2,1.119010961,4.293889240e-3,-6.836300520e-7,-2.023372700e-9,1.039040018e-12,-3.391454870e3,1.849699470e1]),
        np.array([-1.037939022e6,2.344830282e3,1.819732036,1.267847582e-3,-2.188067988e-7,2.053719572e-11,-8.193467050e-16,-1.689010929e4,1.738716506e1])),
}

def hmass(gas, T):
    MW, LO, HI = DB[gas]; a = LO if T < TMID else HI; Ti = 1.0/T; lnT = np.log(T)
    hRT = (-a[0]*Ti*Ti + a[1]*lnT*Ti + a[2] + a[3]*T/2 + a[4]*T**2/3
           + a[5]*T**3/4 + a[6]*T**4/5 + a[7]*Ti)
    return RU*T*hRT/MW

def state(gas, P, T):
    MW = DB[gas][0]; R = RU/MW; ro = P/(R*T); e = hmass(gas, T) - R*T
    return ro, ro*e

roL, roeL = state(gasL, PL, TL)
roR, roeR = state(gasR, PR, TR)

shutil.copy(mesh, dst)
with h5py.File(dst, 'r+') as f:
    co = f['MESH/COORD'][:].reshape(-1, 3); c = f['MESH/CONNE'][:].reshape(-1, 9)
    cx = co[c[:, 1:9], 0].mean(1); n = len(cx)
    left = cx < xdia
    ro  = np.where(left, roL, roR)
    roe = np.where(left, roeL, roeR)
    T   = np.where(left, TL, TR)
    P   = np.where(left, PL, PR)
    roY0 = np.where(left, ro, 0.0)   # 左 = 純 species0
    roY1 = np.where(left, 0.0, ro)   # 右 = 純 species1
    g = f['VALUE'] if 'VALUE' in f else f.create_group('VALUE')
    def setv(k, v):
        v = np.asarray(v, np.float32)
        if 'VALUE/'+k in f: f['VALUE/'+k][...] = v
        else: g.create_dataset(k, data=v)
    setv('ro', ro); setv('roUx', np.zeros(n)); setv('roUy', np.zeros(n)); setv('roUz', np.zeros(n))
    setv('roe', roe); setv('T', T); setv('P', P)
    setv('Ux', np.zeros(n)); setv('Uy', np.zeros(n)); setv('Uz', np.zeros(n))
    setv('roY0', roY0); setv('roY1', roY1)
    setv('Y0', roY0/ro); setv('Y1', roY1/ro)
    if 'VALUE/wall_dist' not in f: g.create_dataset('wall_dist', data=np.ones(n, np.float32))
    print(f"[species] gasL={gasL} gasR={gasR} n={n} "
          f"L(P={PL:.2e},T={TL},ro={roL:.4f}) R(P={PR:.2e},T={TR},ro={roR:.4f}) "
          f"roe[{roe.min():.3e},{roe.max():.3e}] -> {dst}")
