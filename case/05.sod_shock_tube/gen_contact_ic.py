#!/usr/bin/env python3
# 多成分 material-contact 初期値: 圧力 p0・速度 u0 を一様にし、組成 (と ρ,T) だけを界面で変える。
# 理想スキームは p,u を永久に一様保持する; スプリアス圧力振動が PEP 不足を測る。
#   x < xdia: gasL (純 species0), x >= xdia: gasR (純 species1)。smooth>0 で tanh 平滑界面。
# 使い方: python3 gen_contact_ic.py <mesh.h5> <out.h5> <gasL> <gasR> <p0> <u0> <T> [smooth_width]
#   例 (He/N2 sharp stationary): ... He N2 1.0e5 0.0 300 0
import sys, numpy as np, h5py, shutil

mesh, dst, gasL, gasR = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
p0 = float(sys.argv[5]); u0 = float(sys.argv[6]); T0 = float(sys.argv[7])
smooth = float(sys.argv[8]) if len(sys.argv) > 8 else 0.0
xdia = float(sys.argv[9]) if len(sys.argv) > 9 else 0.5
RU = 8.314462618; TMID = 1000.0

# NASA-9 (CEA) DB: MW [kg/mol], low[9], high[9] (thermo_d.cu と一致)
DB = {
 'He': (0.0040026,
        np.array([0,0,2.5,0,0,0,0,-745.375,0.928723974], float),
        np.array([0,0,2.5,0,0,0,0,-745.375,0.928723974], float)),
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

def Rgas(gas): return RU/DB[gas][0]

shutil.copy(mesh, dst)
with h5py.File(dst, 'r+') as f:
    co = f['MESH/COORD'][:].reshape(-1, 3); c = f['MESH/CONNE'][:].reshape(-1, 9)
    cx = co[c[:, 1:9], 0].mean(1); n = len(cx)
    # 組成混合率 wL (左成分の質量分率): sharp = step, smooth = tanh
    if smooth > 0.0:
        wL = 0.5*(1.0 - np.tanh((cx - xdia)/smooth))
    else:
        wL = np.where(cx < xdia, 1.0, 0.0)
    YL = wL                      # species0 (gasL) 質量分率
    YR = 1.0 - wL                # species1 (gasR)
    # 一様 p0,T0 → 局所 R_mix, ρ = p0/(R_mix T0)。e = Σ Y_s e_s(T0)。
    Rmix = YL*Rgas(gasL) + YR*Rgas(gasR)
    eL = hmass(gasL, T0) - Rgas(gasL)*T0
    eR = hmass(gasR, T0) - Rgas(gasR)*T0
    emix = YL*eL + YR*eR
    ro = p0/(Rmix*T0)
    ek = 0.5*u0*u0
    roe = ro*(emix + ek)
    def setv(k, v):
        v = np.asarray(v, np.float32); g = f['VALUE']
        if 'VALUE/'+k in f: f['VALUE/'+k][...] = v
        else: g.create_dataset(k, data=v)
    if 'VALUE' not in f: f.create_group('VALUE')
    setv('ro', ro); setv('roUx', ro*u0); setv('roUy', np.zeros(n)); setv('roUz', np.zeros(n))
    setv('roe', roe); setv('T', np.full(n, T0)); setv('P', np.full(n, p0))
    setv('Ux', np.full(n, u0)); setv('Uy', np.zeros(n)); setv('Uz', np.zeros(n))
    setv('roY0', ro*YL); setv('roY1', ro*YR); setv('Y0', YL); setv('Y1', YR)
    if 'VALUE/wall_dist' not in f: f['VALUE'].create_dataset('wall_dist', data=np.ones(n, np.float32))
    print(f"[contact] {gasL}/{gasR} p0={p0:.3e} u0={u0} T0={T0} smooth={smooth} n={n} "
          f"ro[{ro.min():.4f},{ro.max():.4f}] Rmix[{Rmix.min():.1f},{Rmix.max():.1f}] -> {dst}")
