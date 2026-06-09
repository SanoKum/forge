#!/usr/bin/env python3
# 28.cutler ternary [He,O2,N2] の初期場。ドメインは ambient 空気 (O2/N2 by mass) を
# 緩やかな下流速度で初期化。roe=ρ(e_mix(Y,T)+ek) (NASA 絶対基準), roY{s} を書く。
# 組成順は [He, O2, N2] (index 0,1,2)。使い方: python3 gen_ternary_ic.py <mesh.h5> <input.h5>
import sys, numpy as np, h5py, shutil

src = sys.argv[1] if len(sys.argv) > 1 else "run_0011_he_o2_n2/coaxial.h5"
dst = sys.argv[2] if len(sys.argv) > 2 else "run_0011_he_o2_n2/input.h5"

RU = 8.314462618
TMID = 1000.0
# NASA-9 (下段 200-1000K, 上段 1000-6000K)
DB = {
 'O2': (0.0319988, np.array([-3.425563420e4,4.847000970e2,1.119010961,4.293889240e-3,-6.836300520e-7,-2.023372700e-9,1.039040018e-12,-3.391454870e3,1.849699470e1]),
                   np.array([-1.037939022e6,2.344830282e3,1.819732036,1.267847582e-3,-2.188067988e-7,2.053719572e-11,-8.193467050e-16,-1.689010929e4,1.738716506e1])),
 'N2': (0.0280134, np.array([2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1]),
                   np.array([5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1])),
}
def hmass(sp, T):
    MW, LO, HI = DB[sp]
    a = np.where(T[...,None] < TMID, LO, HI); Ti = 1/T; lnT = np.log(T)
    hRT = (-a[...,0]*Ti*Ti + a[...,1]*lnT*Ti + a[...,2] + a[...,3]*T/2
           + a[...,4]*T**2/3 + a[...,5]*T**3/4 + a[...,6]*T**4/5 + a[...,7]*Ti)
    return RU*T*hRT/MW

# ambient 空気 (mass) と緩やかな下流速度
YO2, YN2 = 0.23293, 0.76707
Rmix = RU*(YO2/DB['O2'][0] + YN2/DB['N2'][0])
P0, T0, u0 = 101325.0, 300.0, 100.0

shutil.copy(src, dst)
with h5py.File(dst, 'r+') as f:
    n = f['VALUE/ro'].shape[0]
    T = np.full(n, T0); P = np.full(n, P0)
    ro = P / (Rmix * T)
    ek = 0.5*u0*u0
    hmix = YO2*hmass('O2', T) + YN2*hmass('N2', T)
    e = hmix - Rmix*T
    roe = ro*(e + ek)
    for k, v in [('ro',ro),('roUx',ro*u0),('roUy',ro*0.0),('roUz',ro*0.0),('roe',roe)]:
        f['VALUE/'+k][...] = v.astype(np.float32)
    for nm, Y in [('roY0', np.zeros(n)), ('roY1', np.full(n,YO2)), ('roY2', np.full(n,YN2))]:
        d = 'VALUE/'+nm
        if d in f: del f[d]
        f.create_dataset(d, data=(Y*ro).astype(np.float32))
    print('ternary IC: air ambient ro=%.4f Rmix=%.1f roe=%.4e  Y[He,O2,N2]=[0,%.3f,%.3f] n=%d -> %s'
          % (ro[0], Rmix, roe[0], YO2, YN2, n, dst))
