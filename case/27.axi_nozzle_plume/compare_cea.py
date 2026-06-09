#!/usr/bin/env python3
# ノズル中心線の CEA 照合: 計算 (CPG-N2 / TP-N2) の中心線 (T,P,M) を、
# Pt=4e6, Tt=1500K からの NASA-N2 cp(T) 1次元等エントロピー解 (=非反応フローズン CEA) と比較。
#   TP は cp(T) 等エントロピーに乗り、CPG (定数 cp) は高温側でずれるはず。
import sys, numpy as np, h5py, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

Pt, Tt = 4.0e6, 1500.0
RU = 8.314462618; MW = 0.0280134; R = RU/MW; TMID = 1000.0
LO = np.array([2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1])
HI = np.array([5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1])
def _co(T): return np.where(np.atleast_1d(T)[...,None]<TMID, LO, HI)
def cp_mass(T):
    a=_co(T); T=np.atleast_1d(T).astype(float); Ti=1/T
    return RU*(a[...,0]*Ti*Ti+a[...,1]*Ti+a[...,2]+a[...,3]*T+a[...,4]*T*T+a[...,5]*T**3+a[...,6]*T**4)/MW
def h_mass(T):
    a=_co(T); T=np.atleast_1d(T).astype(float); Ti=1/T; lnT=np.log(T)
    hRT=-a[...,0]*Ti*Ti+a[...,1]*lnT*Ti+a[...,2]+a[...,3]*T/2+a[...,4]*T**2/3+a[...,5]*T**3/4+a[...,6]*T**4/5+a[...,7]*Ti
    return RU*T*hRT/MW
def s0_mass(T):
    a=_co(T); T=np.atleast_1d(T).astype(float); Ti=1/T; lnT=np.log(T)
    sR=-a[...,0]*Ti*Ti/2-a[...,1]*Ti+a[...,2]*lnT+a[...,3]*T+a[...,4]*T*T/2+a[...,5]*T**3/3+a[...,6]*T**4/4+a[...,8]
    return RU*sR/MW

def isen_T_from_P(P):
    """NASA-N2 等エントロピー: s0(T)-R ln(P) = s0(Tt)-R ln(Pt) を解いて T(P)。"""
    target = s0_mass(Tt)[0] + R*np.log(P/Pt)   # s0(T) = s0(Tt)+R ln(P/Pt)
    T = np.full_like(P, Tt, dtype=float)
    for _ in range(60):
        f = s0_mass(T) - target; df = cp_mass(T)/T
        dT = f/np.maximum(df,1e-12); dT = np.clip(dT, -0.4*T, 0.4*T)
        T = np.clip(T - dT, 50.0, Tt*1.01)
    return T
def cpg_T_from_P(P, gam=1.4):
    return Tt*(P/Pt)**((gam-1)/gam)

def centerline(h5, meshh5=None):
    import os
    if meshh5 is None: meshh5=os.path.join(os.path.dirname(h5),'axi_nozzle.h5')
    with h5py.File(meshh5,'r') as f:
        cc=f['CELLS/centCoords'][:].reshape(-1,3)
    with h5py.File(h5,'r') as f:
        d={k:f['VALUE/'+k][:] for k in ('P','T','ro','sonic','Ux','Uy','Uz')}
    n=len(d['P']); cc=cc[:n]
    rr=np.sqrt(cc[:,1]**2)
    # 軸近傍 (最小 r バンド) かつ ノズル内部 (x <= 排気端 0.0675 m) のセルを中心線とする。
    # プルーム (x>exit) は衝撃ダイヤモンドで非等エントロピーなので除外。
    NOZ_EXIT=0.0675
    thr=np.percentile(rr, 4)
    m=(rr<=thr)&(cc[:,0]<=NOZ_EXIT)
    x=cc[m,0]; i=np.argsort(x)
    P=d['P'][m][i]; T=d['T'][m][i]; U=np.sqrt(d['Ux'][m]**2+d['Uy'][m]**2+d['Uz'][m]**2)[i]
    Mach=U/np.maximum(d['sonic'][m][i],1e-9)
    return x[i],P,T,Mach

cpg=centerline(sys.argv[1] if len(sys.argv)>1 else 'run_0008_cpgN2_2nd/res_8000.h5')
tp =centerline(sys.argv[2] if len(sys.argv)>2 else 'run_0009_tpN2_2nd/res_8000.h5')

# isentropic references over the pressure range
Pgrid=np.logspace(np.log10(max(min(cpg[1].min(),tp[1].min()),1e3)), np.log10(Pt), 200)
T_nasa=isen_T_from_P(Pgrid); T_cpg=cpg_T_from_P(Pgrid)

fig,ax=plt.subplots(1,3,figsize=(16,5))
# (1) T vs P  (the CEA comparison)
ax[0].plot(Pgrid/1e3, T_nasa,'k-',lw=2,label='NASA-N2 cp(T) isentropic (CEA frozen)')
ax[0].plot(Pgrid/1e3, T_cpg,'k--',lw=1.2,label='CPG isentropic (const cp, γ=1.4)')
ax[0].plot(cpg[1]/1e3, cpg[2],'C0.',ms=3,label='forge CPG-N2 centerline')
ax[0].plot(tp[1]/1e3, tp[2],'C3.',ms=3,label='forge TP-N2 centerline')
ax[0].set_xscale('log'); ax[0].set_xlabel('P [kPa]'); ax[0].set_ylabel('T [K]'); ax[0].grid(alpha=0.3,which='both'); ax[0].legend(fontsize=8)
ax[0].set_title('centerline T(P) vs isentropic')
# (2) T vs x
ax[1].plot(cpg[0],cpg[2],'C0.',ms=3,label='CPG-N2'); ax[1].plot(tp[0],tp[2],'C3.',ms=3,label='TP-N2')
ax[1].set_xlabel('x [m]'); ax[1].set_ylabel('T [K]'); ax[1].grid(alpha=0.3); ax[1].legend(fontsize=8); ax[1].set_title('centerline T(x)')
# (3) Mach vs x
ax[2].plot(cpg[0],cpg[3],'C0.',ms=3,label='CPG-N2'); ax[2].plot(tp[0],tp[3],'C3.',ms=3,label='TP-N2')
ax[2].set_xlabel('x [m]'); ax[2].set_ylabel('Mach'); ax[2].grid(alpha=0.3); ax[2].legend(fontsize=8); ax[2].set_title('centerline Mach(x)')
fig.suptitle('27.axi_nozzle_plume: N2 nozzle Pt=4MPa Tt=1500K — TP cp(T) vs CPG vs NASA isentropic (CEA)',fontsize=11)
fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig('cea_comparison.png',dpi=130); print('wrote cea_comparison.png')

# quantitative: deviation of each from NASA isentropic at matched P
def devN(cl):
    P,T=cl[1],cl[2]; mm=(P>5e3)&(P<Pt)
    Tn=isen_T_from_P(P[mm])
    return np.sqrt(np.mean(((T[mm]-Tn)/Tn)**2))*100
print(f"centerline T vs NASA-N2 isentropic (RMS rel %):  CPG-N2={devN(cpg):.2f}%   TP-N2={devN(tp):.2f}%")
print(f"  -> TP should be << CPG if cp(T) matters (TP on CEA curve, CPG off)")
