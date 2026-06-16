#!/usr/bin/env python3
"""メッシュの実面積比から等エントロピー p/p0 を算出し、forge dry と実験に重ねる。

切り分け:
  forge centerline ≈ geometry-isentrope < exp-isentrope  → 形状(面積分布)が実験と違う
  forge centerline < geometry-isentrope                  → forge の過膨張(数値/BL)
usage: check_geometry_isentrope.py <run_dir> <label>
"""
import sys, glob, os
import numpy as np, h5py
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

GAM = 1.4
P0 = 59070.0

def latest(run):
    fs = [f for f in glob.glob(f"{run}/res_*.h5") if "nan" not in f]
    return max(fs, key=lambda f: int(f.split("res_")[1].split(".")[0]))
def mesh(run):
    return [p for p in glob.glob(f"{run}/nozzle_*.h5") if "res_" not in os.path.basename(p)][0]

def ar_to_M(ar, sup):
    if ar <= 1.0 + 1e-7: return 1.0
    f = lambda M: (1.0/M)*((2.0/(GAM+1))*(1+0.5*(GAM-1)*M*M))**((GAM+1)/(2*(GAM-1))) - ar
    try: return brentq(f, 1.0+1e-6, 50.0) if sup else brentq(f, 1e-4, 1.0)
    except Exception: return np.nan

def load_exp():
    rows=[]
    for ln in open("wyslouzil_fig3_pp0.csv"):
        p=ln.strip().split(",")
        if p and p[0] and p[0][0].isdigit(): rows.append([float(p[0]),float(p[1]),float(p[2])])
    a=np.array(rows); return a[:,0],a[:,1],a[:,2]

def main(run, label):
    cc = h5py.File(mesh(run),"r")["/CELLS/centCoords"][:].reshape(-1,3)
    x,y,z = cc[:,0],cc[:,1],cc[:,2]
    P = np.array(h5py.File(latest(run),"r")["/VALUE"]["P"])
    # local full height -> area ratio  (full duct height = 2*ymax in plane)
    nb=240; xe=np.linspace(x.min(),x.max(),nb+1); xcb=0.5*(xe[:-1]+xe[1:]); h=np.full(nb,np.nan)
    for i in range(nb):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>3: h[i]=y[s].max()-y[s].min()
    ok=np.isfinite(h); xcb,h=xcb[ok],h[ok]
    hstar=h.min(); xth=xcb[np.argmin(h)]; ar=h/hstar
    M_iso=np.array([ar_to_M(a, xx>xth) for a,xx in zip(ar,xcb)])
    pp0_geom = (1+0.5*(GAM-1)*M_iso**2)**(-GAM/(GAM-1))
    xcm_geom = (xcb-xth)*100.0
    # forge centerline
    zc=0.5*(z.min()+z.max()); dz=z.max()-z.min()
    sel=(np.abs(y)<(y.max()-y.min())*0.02)&(np.abs(z-zc)<max(dz*0.12,1e-6))
    xs=x[sel]; o=np.argsort(xs); xcm_f=(xs[o]-xth)*100.0; pp0_f=P[sel][o]/P0
    # exp
    ex,eiso,econd=load_exp()

    plt.figure(figsize=(8.5,6))
    plt.plot(xcm_f, pp0_f, '.', ms=3, color='tab:blue', label=f'forge {label} centerline')
    plt.plot(xcm_geom, pp0_geom, '-', color='green', lw=2, label='isentropic from MESH area ratio')
    plt.plot(ex,eiso,'s',color='gray',ms=7,label='Wyslouzil isentrope (exp)')
    plt.plot(ex,econd,'o',color='black',ms=6,mfc='none',label='Wyslouzil cond 1kPa (exp)')
    plt.xlim(-1,9); plt.ylim(0.1,0.55); plt.xlabel('distance from throat [cm]'); plt.ylabel('p/p0')
    plt.grid(alpha=0.3); plt.legend(fontsize=9)
    plt.title('Geometry area-ratio isentrope vs forge vs experiment')
    plt.tight_layout(); plt.savefig("geometry_isentrope_check.png",dpi=130)
    print("saved geometry_isentrope_check.png")
    print(f"throat x={xth*100:.2f}cm  exit area ratio (mesh)={ar[-1]:.3f}  -> isentropic exit M={M_iso[-1]:.3f}")

    print("\n x_cm | exp_iso | geom_iso | forge | (geom-exp)% | (forge-exp)%")
    for xq in [0.1,1.1,2.1,3.1,4.2,5.2,6.0,7.0,8.0]:
        ei=np.interp(xq,ex,eiso)
        mg=np.isfinite(pp0_geom); gg=np.interp(xq,xcm_geom[mg],pp0_geom[mg])
        ff=np.interp(xq,xcm_f,pp0_f)
        print(f"{xq:5.1f} | {ei:.3f}  | {gg:.3f}   | {ff:.3f} |  {100*(gg-ei)/ei:+5.1f}    |  {100*(ff-ei)/ei:+5.1f}")

if __name__=="__main__":
    main(sys.argv[1], sys.argv[2])
