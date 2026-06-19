import sys, shutil, numpy as np, h5py
src,dst=sys.argv[1],sys.argv[2]
Ru=8.314462618; MW=0.0280134; R=Ru/MW; gam=1.4
LO=np.array([2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1])
def e_mass(T):  # NASA-9 N2 internal energy [J/kg], T<1000 (LO)
    a=LO; Ti=1/T; lnT=np.log(T)
    hRT=-a[0]*Ti*Ti + a[1]*lnT*Ti + a[2] + a[3]*T/2 + a[4]*T**2/3 + a[5]*T**3/4 + a[6]*T**4/5 + a[7]*Ti
    h=hRT*Ru*T; e=h-Ru*T; return e/MW
shutil.copy(src,dst)
with h5py.File(dst,'r+') as h:
    ro=np.array(h['/VALUE/ro']); u=np.array(h['/VALUE/roUx'])/ro; roe=np.array(h['/VALUE/roe'])
    P=(gam-1)*(roe-0.5*ro*u*u); Rcpg=1039*(gam-1)/gam; T=P/(ro*Rcpg)
    roe_tp=ro*(e_mass(T)+0.5*u*u)
    h['/VALUE/roe'][...]=roe_tp.astype(np.float32)
    print(f"T[{T.min():.1f},{T.max():.1f}]K  roe_cpg/ro[{(roe/ro).min():.0f},{(roe/ro).max():.0f}]  roe_tp/ro[{(roe_tp/ro).min():.0f},{(roe_tp/ro).max():.0f}]")
print("wrote",dst)
