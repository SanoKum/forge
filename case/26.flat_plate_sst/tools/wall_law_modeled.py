#!/usr/bin/env python3
"""壁関数 (enhanced wall treatment) 整合の Cf 評価。

postprocess_wall_law.py は tau_w = mu_lam*du/dy (第一セル分子勾配) で Cf を測るが、
これは壁関数 run では不適切 (粗メッシュの第一セルは対数層にあり、壁せん断の大半を
乱流が担うため分子勾配は真のせん断を大幅に過小評価する)。

本ツールはソルバと同じ Reichardt 普遍速度則の逆解きを収束第一セル速度に適用して
modeled u_tau を求め、Cf_model = 2(u_tau/U_inf)^2 で評価する (ソルバが運動量式に
課している modeled tau_w = rho*u_tau^2 と整合)。比較のため分子勾配 Cf も併記する。
"""
import sys, glob, os
import numpy as np
import h5py

MU = 1.8e-5
RHO_INF = 1.1854
U_INF = 67.782
KAPPA = 0.41

def reichardt_uplus(yp):
    return np.log(1.0+KAPPA*yp)/KAPPA + 7.8*(1.0-np.exp(-yp/11.0)-(yp/11.0)*np.exp(-yp/3.0))

def solve_utau(Ut, y, nu):
    # Newton: f(utau)=Ut/utau - uplus(yp(utau))=0
    utau = np.sqrt(nu*Ut/y)
    for _ in range(50):
        utau = max(utau, 1e-12)
        yp = utau*y/nu
        f = Ut/utau - reichardt_uplus(yp)
        # 数値微分
        d = 1e-6*utau
        yp2=(utau+d)*y/nu
        df = ((Ut/(utau+d)-reichardt_uplus(yp2))-f)/d
        if abs(df)<1e-20: break
        utau -= f/df
    return max(utau,0.0)

def load_latest(run_dir):
    files = glob.glob(os.path.join(run_dir, "res_*.h5"))
    files = [f for f in files if os.path.basename(f)!="res_0.h5"] or files
    return max(files, key=lambda f:int(''.join(c for c in os.path.basename(f) if c.isdigit())))

def cell_centroids(h5):
    coord=h5["MESH/COORD"][:].reshape(-1,3); conne=h5["MESH/CONNE"][:]
    nC=h5["VALUE/ro"].shape[0]; stride=conne.size//nC
    nodes=conne.reshape(nC,stride)[:,1:]
    if nodes.min()==1: nodes=nodes-1
    return coord[nodes].mean(axis=1)

def main():
    run_dir=sys.argv[1]
    stations=[float(s) for s in sys.argv[2:]] or [0.3,0.6,0.9]
    h5=h5py.File(load_latest(run_dir),"r")
    C=cell_centroids(h5); x,y=C[:,0],C[:,1]
    Ux=h5["VALUE/Ux"][:]; ro=h5["VALUE/ro"][:]; vis=h5["VALUE/vis_lam"][:]
    xcol=np.unique(np.round(x[x>1e-6],6))
    print(f"# {os.path.basename(load_latest(run_dir))}  ({run_dir})")
    print(" station |   Re_x    | u_tau_mdl | y+_1  | Cf_model | Cf_molec | Cf_corr  | mdl/corr")
    for xs in stations:
        xc=xcol[np.argmin(np.abs(xcol-xs))]
        col=np.where(np.abs(x-xc)<1e-4)[0]; col=col[np.argsort(y[col])]
        yc,uc,rc,vc=y[col],Ux[col],ro[col],vis[col]
        nu=vc[0]/rc[0]
        ut_m=solve_utau(uc[0],yc[0],nu)
        yp1=ut_m*yc[0]/nu
        Cf_model=rc[0]*ut_m**2/(0.5*RHO_INF*U_INF**2)
        tau_mol=vc[0]*uc[0]/yc[0]; Cf_mol=tau_mol/(0.5*RHO_INF*U_INF**2)
        Rex=RHO_INF*U_INF*xc/MU; Cfc=0.0592*Rex**-0.2
        print(f"  {xc:.2f}  | {Rex:.3e} |  {ut_m:6.3f}  | {yp1:5.2f} | {Cf_model:.3e}| {Cf_mol:.3e}| {Cfc:.3e}|  {Cf_model/Cfc:.3f}")

if __name__=="__main__":
    main()
