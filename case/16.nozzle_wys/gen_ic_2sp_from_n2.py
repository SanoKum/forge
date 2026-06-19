#!/usr/bin/env python3
"""
TP N2 単独収束解 (run_0056 等) に H2O 組成を均一付加して N2+H2O 2 成分 IC を作る。
- ro, roUx/Uy/Uz, roK, roOmega, wall_dist は収束解から継承
- Y_H2O = Y1_in で均一設定、Y_N2 = 1 - Y1_in
- roe は収束解の T と新組成で NASA-9 TP 基準に再構成 (gen_tp_ic_wys.py と同一ロジック)

使い方:
  python3 gen_ic_2sp_from_n2.py <src_n2_res.h5> <dst_2sp.h5> [Y_H2O=0.01095] [Tlo=200] [Href_Tref=0]

Href_Tref > 0 を指定すると、solver の thermoHrefTemp と同じく各化学種のエンタルピー基準を
h_s(Tref)=0 へオフセットする (sensible-enthalpy datum)。solver 側 thermoHrefTemp と必ず一致させること。
"""
import sys, numpy as np, h5py

src  = sys.argv[1]
dst  = sys.argv[2]
Y1_in = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01095
Y0_in = 1.0 - Y1_in
HREF_TREF = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0  # 0: 絶対基準 (従来)

RU = 8.314462618
MW = {"N2": 0.0280134, "H2O": 0.0180153}
TMID = 1000.0
LO = {
 "N2":  np.array([2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1]),
 "H2O": np.array([-3.947960830e4,5.755731020e2,9.317826530e-1,7.222712860e-3,-7.342557370e-6,4.955043490e-9,-1.336933246e-12,-3.303974310e4,1.724205775e1]),
}
HI = {
 "N2":  np.array([5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1]),
 "H2O": np.array([1.034972096e6,-2.412698562e3,4.646110780,2.291998307e-3,-6.836830480e-7,9.426468930e-11,-4.822380530e-15,-1.384286509e4,-7.978148510e0]),
}
TLO_ARG = float(sys.argv[4]) if len(sys.argv) > 4 else 200.0
TLO, THI = TLO_ARG, 6000.0

def _hRT(sp, T):
    a = np.where(T[..., None] < TMID, LO[sp], HI[sp])
    Ti = 1.0 / T; lnT = np.log(T)
    return (-a[...,0]*Ti*Ti + a[...,1]*lnT*Ti + a[...,2] + a[...,3]*T/2
            + a[...,4]*T**2/3 + a[...,5]*T**3/4 + a[...,6]*T**4/5 + a[...,7]*Ti)

def _cpR(sp, T):
    a = np.where(T[..., None] < TMID, LO[sp], HI[sp])
    Ti = 1.0 / T
    return (a[...,0]*Ti*Ti + a[...,1]*Ti + a[...,2] + a[...,3]*T
            + a[...,4]*T**2 + a[...,5]*T**3 + a[...,6]*T**4)

def hmass(sp, T):
    R = RU / MW[sp]
    Tc = np.clip(T, TLO, THI)
    h = R * Tc * _hRT(sp, Tc)
    lo = T < TLO;  hi = T > THI
    if np.any(lo):
        cp_lo = R * _cpR(sp, np.full_like(T, TLO))
        h = np.where(lo, R*TLO*_hRT(sp, np.full_like(T, TLO)) + cp_lo*(T - TLO), h)
    if np.any(hi):
        cp_hi = R * _cpR(sp, np.full_like(T, THI))
        h = np.where(hi, R*THI*_hRT(sp, np.full_like(T, THI)) + cp_hi*(T - THI), h)
    return h

# TP N2 収束解を読む
with h5py.File(src, 'r') as f:
    V = f['VALUE']
    ro   = V['ro'][:].astype(np.float64)
    roUx = V['roUx'][:].astype(np.float64)
    roUy = V['roUy'][:].astype(np.float64)
    roUz = V['roUz'][:].astype(np.float64)
    roK  = V['roK'][:].astype(np.float64)
    roOm = V['roOmega'][:].astype(np.float64)
    T    = V['T'][:].astype(np.float64)
    wall = V['wall_dist'][:].astype(np.float64)

# ek は質量当たり運動エネルギー 0.5|u|^2 (solver 規約: dependentVariables_d.cu intE=roe/ro-ek)。
# roU=ro*u なので 0.5*(roU^2)/ro^2。以前は /ro で割っており ek が係数 ro だけ過大→
# roe が高速域 (出口) で過小→solver 読み戻しで T が 200K 床に張り付き発散していた。
ek = 0.5 * (roUx**2 + roUy**2 + roUz**2) / (ro * ro)

R_N2  = RU / MW["N2"]
R_H2O = RU / MW["H2O"]

# エンタルピー基準オフセット (solver の thermoHrefTemp と一致させる)。
# h_s -> h_s - h_s(Tref) で各種の生成エンタルピー差を除去。e_mix も同じ基準へ。
href = {"N2": 0.0, "H2O": 0.0}
if HREF_TREF > 0.0:
    Tref_arr = np.array([HREF_TREF])
    href = {sp: float(hmass(sp, Tref_arr)[0]) for sp in ("N2", "H2O")}
    print(f"href offset: Tref={HREF_TREF}K  h_ref(N2)={href['N2']/1e6:.4f} MJ/kg  h_ref(H2O)={href['H2O']/1e6:.4f} MJ/kg")

def hmass_off(sp, Tarr):
    return hmass(sp, Tarr) - href[sp]

e_mix = Y0_in*(hmass_off("N2", T) - R_N2*T) + Y1_in*(hmass_off("H2O", T) - R_H2O*T)

# forge の thermo は TLO=200K でクランプするため、e_internal < e_mix(TLO) だと
# T 回復の Newton 解が収束せず NaN になる。TLO 床でクランプして回避する。
T_floor = np.full_like(T, TLO)
e_mix_floor = Y0_in*(hmass_off("N2", T_floor) - R_N2*TLO) + Y1_in*(hmass_off("H2O", T_floor) - R_H2O*TLO)
n_clamp = int(np.sum(e_mix < e_mix_floor))
e_mix = np.maximum(e_mix, e_mix_floor)

roe_new = ro * (e_mix + ek)

roY0 = ro * Y0_in
roY1 = ro * Y1_in

print(f"src : {src}")
print(f"T   : {T.min():.1f} - {T.max():.1f} K")
print(f"e_mix_floor(TLO=200K) = {e_mix_floor.flat[0]/1000:.2f} kJ/kg  (clamped {n_clamp} cells)")
print(f"roe/ro: {(roe_new/ro).min()/1000:.1f} to {(roe_new/ro).max()/1000:.1f} kJ/kg")
print(f"Y0(N2)={Y0_in:.5f}  Y1(H2O)={Y1_in:.5f}")
print(f"dst : {dst}")

with h5py.File(dst, 'w') as f:
    grp = f.create_group('VALUE')
    def w(name, arr):
        grp.create_dataset(name, data=arr.astype(np.float32), compression='gzip')
    w('ro',       ro)
    w('roUx',     roUx)
    w('roUy',     roUy)
    w('roUz',     roUz)
    w('roe',      roe_new)
    w('roY0',     roY0)
    w('roY1',     roY1)
    w('roK',      roK)
    w('roOmega',  roOm)
    w('wall_dist', wall)

print("done.")
