#!/usr/bin/env python3
"""保存スカラー (Bilger 混合分率, 全 H/O 含有種から元素で計算) の断面フラックス ∫ρuξ 2πr dr が軸方向に保存されているかの検査。
定常なら入口の燃料流量 ṁ_jet と一致し続けるはず。  python3 xi_flux_check.py run_dir step [run_dir step ...]"""
import sys, h5py, numpy as np
from scipy.interpolate import griddata
D=0.00457; WH=1.00794e-3; WO=15.9994e-3; WN=14.0067e-3
names=["H2","O2","H","O","OH","H2O","HO2","H2O2","N2"]; nH={"H2":2,"H":1,"OH":1,"H2O":2,"HO2":1,"H2O2":2}; nO={"O2":2,"O":1,"OH":1,"H2O":1,"HO2":2,"H2O2":2}
W={"H2":2*WH,"O2":2*WO,"H":WH,"O":WO,"OH":WH+WO,"H2O":2*WH+WO,"HO2":WH+2*WO,"H2O2":2*WH+2*WO,"N2":2*WN}
def beta(Y):  # Y: dict name→array
    ZH=sum(Y[s]*nH.get(s,0)*WH/W[s] for s in Y); ZO=sum(Y[s]*nO.get(s,0)*WO/W[s] for s in Y); return ZH/(2*WH)-ZO/WO
mF=0.25*W["H2"]+0.75*W["N2"]; mO=0.15*W["O2"]+0.099*W["H2O"]+0.751*W["N2"]
bF=beta({"H2":np.array(0.25*W["H2"]/mF)}); bO=beta({"O2":np.array(0.15*W["O2"]/mO),"H2O":np.array(0.099*W["H2O"]/mO)})
args=sys.argv[1:]
for run,step in zip(args[::2],args[1::2]):
    with h5py.File(f"{run}/cabra.h5") as m: xyz=m["MESH/COORD"][:].reshape(-1,3)
    with h5py.File(f"{run}/res_{step}.h5") as h:
        V=h["VALUE"]; ro=V["ro"][:]; Ux=V["Ux"][:]; N=len(ro); ns=sum(1 for i in range(16) if f"Y{i}" in V)
        Y={names[i]:V[f"Y{i}"][:] for i in range(ns)} if ns==9 else {"H2":V["Y0"][:]}
    x=xyz[:N,0]; y=xyz[:N,1]; xi=(beta(Y)-bO)/(bF-bO) if ns==9 else Y["H2"]/(0.25*W["H2"]/mF)
    r=np.linspace(0,0.1,1601); out=[]
    for z in (0.5,2,5,10,15,20,25,30,35,40,50):
        pts=np.column_stack([np.full_like(r,z*D),r]); R,U,X=(griddata((x,y),a,pts,method="linear") for a in (ro,Ux,xi)); ok=np.isfinite(R*U*X)
        out.append((z, np.trapezoid((R*U*X*2*np.pi*r)[ok],r[ok]), np.trapezoid((R*U*2*np.pi*r)[ok],r[ok]), X[0]))
    print(f"{run} step {step} (ns={ns}):")
    print("   z/d   " + " ".join(f"{z:7.1f}" for z,_,_,_ in out)); print("   m_xi  " + " ".join(f"{a*1e3:7.3f}" for _,a,_,_ in out) + "  [g/s]  (jet fuel flow ~1.5)")
    print("   m     " + " ".join(f"{b*1e3:7.2f}" for _,_,b,_ in out) + "  [g/s]"); print("   xi_ax " + " ".join(f"{c:7.3f}" for _,_,_,c in out))
