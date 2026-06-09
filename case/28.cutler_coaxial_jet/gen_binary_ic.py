#!/usr/bin/env python3
# 28.cutler binary [He, N2] の初期場を生成する。
#   ドメインは静止した N2 ambient (P=101325, T=300, u=0)。
#   TP 整合のため ro=P/(R_N2 T), roe=ro*e_N2(T) (NASA 絶対基準) を再構成し、
#   組成 roY0(He)=0, roY1(N2)=ro を VALUE/ に書く (readValueHDF5 がこれを読む)。
# 使い方: python3 gen_binary_ic.py <converted.h5> <input.h5>
import sys, numpy as np, h5py, shutil

src = sys.argv[1] if len(sys.argv) > 1 else "run_0003_he_n2_inviscid/coaxial.h5"
dst = sys.argv[2] if len(sys.argv) > 2 else "run_0003_he_n2_inviscid/input.h5"

RU = 8.314462618
# --- N2 NASA-9 (gen_tp_ic.py と同一係数) ---
MW_N2 = 0.0280134
R_N2 = RU / MW_N2
TMID = 1000.0
LO = np.array([2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,
               1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1])
HI = np.array([5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,
               1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1])
def hmass_N2(T):
    a = np.where(T[...,None] < TMID, LO, HI); Ti = 1/T; lnT = np.log(T)
    hRT = (-a[...,0]*Ti*Ti + a[...,1]*lnT*Ti + a[...,2] + a[...,3]*T/2
           + a[...,4]*T**2/3 + a[...,5]*T**3/4 + a[...,6]*T**4/5 + a[...,7]*Ti)
    return RU*T*hRT/MW_N2

# --- ambient N2 (緩やかな下流速度で初期化) ---
#   outlet_statPress は Un<0 (backflow) 分岐で全温・全圧から逆算するが、静止場 (u=0) だと
#   出口セルが Un=0 近傍で backflow に落ち、Pt/Tt 未指定なら 0/0 で NaN 化する。
#   下流向きの初期速度を与えて出口を常に outgoing に保ち、startup の発散を防ぐ。
P0, T0, u0 = 101325.0, 300.0, 100.0

shutil.copy(src, dst)
with h5py.File(dst, 'r+') as f:
    n = f['VALUE/ro'].shape[0]
    T = np.full(n, T0); P = np.full(n, P0)
    ro = P / (R_N2 * T)
    ek = 0.5 * u0 * u0
    e = hmass_N2(T) - R_N2 * T          # NASA 絶対内部エネルギー
    roe = ro * (e + ek)
    upd = {'ro': ro, 'roUx': ro*u0, 'roUy': ro*0.0, 'roUz': ro*0.0, 'roe': roe}
    for k, v in upd.items():
        f['VALUE/'+k][...] = v.astype(np.float32)
    # 組成: He=0, N2=ro (ambient は純 N2)
    roY0 = np.zeros(n, dtype=np.float32)   # He
    roY1 = ro.astype(np.float32)           # N2
    for name, arr in [('roY0', roY0), ('roY1', roY1)]:
        if 'VALUE/'+name in f: del f['VALUE/'+name]
        f.create_dataset('VALUE/'+name, data=arr)
    print(f"binary IC: N2 ambient ro={ro[0]:.4f} roe={roe[0]:.4e} "
          f"roY0(He)=0 roY1(N2)={roY1[0]:.4f} n={n} -> {dst}")
