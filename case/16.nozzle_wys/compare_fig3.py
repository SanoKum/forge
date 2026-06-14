import sys, glob, os, numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P0=59070.0
def latest(run):
    fs=[f for f in glob.glob(f"{run}/res_*.h5") if "nan" not in f]
    return max(fs,key=lambda f:int(f.split("res_")[1].split(".")[0]))
def mesh(run):
    return [p for p in glob.glob(f"{run}/nozzle_*.h5") if "res_" not in os.path.basename(p)][0]
def centerline(run):
    cc=h5py.File(mesh(run),"r")["/CELLS/centCoords"][:].reshape(-1,3)
    x,y=cc[:,0],cc[:,1]
    res=latest(run); v=h5py.File(res,"r")["/VALUE"]
    P=np.array(v["P"]); g=np.array(v["g_0"]) if "g_0" in v else np.zeros_like(P)
    T=np.array(v["T"])
    # throat = x-bin of min channel height
    nb=200; xe=np.linspace(x.min(),x.max(),nb+1); xcb=0.5*(xe[:-1]+xe[1:]); h=np.full(nb,np.nan)
    for i in range(nb):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>3: h[i]=y[s].max()-y[s].min()
    ok=np.isfinite(h); xthroat=xcb[ok][np.argmin(h[ok])]
    epsy=(y.max()-y.min())*0.025; sel=np.abs(y)<epsy
    xs=x[sel]; o=np.argsort(xs)
    return (xs[o]-xthroat)*100.0, P[sel][o]/P0, g[sel][o], T[sel][o], os.path.basename(res)

xc_d,pp_d,_,T_d,rd = centerline("run_0007_h2o_cond")
xc_c,pp_c,g_c,T_c,rc = centerline("run_0008_h2o_cond_on")
fig=np.array([[float(t) for t in l.split(",")] for l in open("wyslouzil_fig3_1kPa.csv") if l.strip() and not l.startswith("#") and l[0].isdigit()])

# bin to smooth (centerline has scatter from unstructured cells)
def binavg(x,y,xe):
    xc=0.5*(xe[:-1]+xe[1:]); yb=np.full(len(xc),np.nan)
    for i in range(len(xc)):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>0: yb[i]=np.median(y[s])
    return xc,yb
xe=np.linspace(-1,9,60)
xb_d,pb_d=binavg(xc_d,pp_d,xe); xb_c,pb_c=binavg(xc_c,pp_c,xe)
_,gb_c=binavg(xc_c,g_c,xe); _,Tb_d=binavg(xc_d,T_d,xe); _,Tb_c=binavg(xc_c,T_c,xe)

fig_,(ax,ax2)=plt.subplots(1,2,figsize=(13,5))
ax.plot(xb_d,pb_d,'-',color='tab:blue',label=f'forge dry ({rd})')
ax.plot(xb_c,pb_c,'-',color='tab:red',label=f'forge cond ({rc})')
ax.plot(fig[:,0],fig[:,1],'s--',color='gray',ms=5,label='Wyslouzil isentrope')
ax.plot(fig[:,0],fig[:,2],'o-',color='black',ms=5,label='Wyslouzil cond 1kPa')
ax.set_xlim(-1,8.5); ax.set_ylim(0.15,0.55); ax.set_xlabel('distance from throat [cm]'); ax.set_ylabel('p/p0')
ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.set_title('Wyslouzil Fig.3 comparison (p/p0)')
ax2.plot(xb_c,gb_c*1e3,'-',color='tab:green',label='g (liquid mass frac) ×1e3')
ax2.plot(xb_d,Tb_d,'--',color='tab:blue',label='T dry [K]')
ax2.plot(xb_c,Tb_c,'-',color='tab:red',label='T cond [K]')
ax2.set_xlim(-1,8.5); ax2.set_xlabel('distance from throat [cm]'); ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
ax2.set_title('condensate g and T (latent heating)')
plt.tight_layout(); plt.savefig("fig3_compare.png",dpi=110)
# print numeric comparison at Fig.3 stations
print("x_cm  dry  cond  exp_iso exp_cond   (forge)")
for xq in [2.5,3,4,5,6,8]:
    pd=np.interp(xq,xb_d[np.isfinite(pb_d)],pb_d[np.isfinite(pb_d)])
    pcv=np.interp(xq,xb_c[np.isfinite(pb_c)],pb_c[np.isfinite(pb_c)])
    ei=np.interp(xq,fig[:,0],fig[:,1]); ec=np.interp(xq,fig[:,0],fig[:,2])
    print("%4.1f  %.3f %.3f   %.3f  %.3f   ratio cond/dry=%.3f (exp %.3f)"%(xq,pd,pcv,ei,ec,pcv/pd,ec/ei))
print("saved fig3_compare.png")
