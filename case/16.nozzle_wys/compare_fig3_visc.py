import glob, os, numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P0=59070.0
def latest(run):
    fs=[f for f in glob.glob(f"{run}/res_*.h5") if "nan" not in f]
    return max(fs,key=lambda f:int(f.split("res_")[1].split(".")[0])) if fs else None
def mesh(run):
    return [p for p in glob.glob(f"{run}/nozzle_*.h5") if "res_" not in os.path.basename(p)][0]
def centerline(run):
    res=latest(run)
    if res is None: return None
    cc=h5py.File(mesh(run),"r")["/CELLS/centCoords"][:].reshape(-1,3)
    x,y=cc[:,0],cc[:,1]
    v=h5py.File(res,"r")["/VALUE"]
    P=np.array(v["P"]); g=np.array(v["g_0"]) if "g_0" in v else np.zeros_like(P); T=np.array(v["T"])
    nb=200; xe=np.linspace(x.min(),x.max(),nb+1); xcb=0.5*(xe[:-1]+xe[1:]); h=np.full(nb,np.nan)
    for i in range(nb):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>3: h[i]=y[s].max()-y[s].min()
    ok=np.isfinite(h); xth=xcb[ok][np.argmin(h[ok])]
    epsy=(y.max()-y.min())*0.02; sel=np.abs(y)<epsy
    xs=x[sel]; o=np.argsort(xs)
    return (xs[o]-xth)*100.0, P[sel][o]/P0, g[sel][o], T[sel][o], os.path.basename(res)
def binavg(x,y,xe):
    xc=0.5*(xe[:-1]+xe[1:]); yb=np.full(len(xc),np.nan)
    for i in range(len(xc)):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>0: yb[i]=np.median(y[s])
    return xc,yb
xe=np.linspace(-1,9,55)
runs={'inv_dry':'run_0007_h2o_cond','inv_cond':'run_0008_h2o_cond_on',
      'sst_dry':'run_0009_h2o_sst_dry','sst_cond':'run_0010_h2o_sst_cond'}
D={}
for k,r in runs.items():
    c=centerline(r)
    if c is None: print('MISSING',r); continue
    xb,pb=binavg(c[0],c[1],xe); _,gb=binavg(c[0],c[2],xe); _,Tb=binavg(c[0],c[3],xe)
    D[k]=(xb,pb,gb,Tb,c[4])
fig=np.array([[float(t) for t in l.split(",")] for l in open("wyslouzil_fig3_1kPa.csv") if l.strip() and not l.startswith("#") and l[0].isdigit()])

plt.figure(figsize=(8,6))
sty={'inv_dry':('--','tab:blue','forge inviscid dry'),'inv_cond':('--','tab:red','forge inviscid cond'),
     'sst_dry':('-','navy','forge SST(visc+turb) dry'),'sst_cond':('-','crimson','forge SST(visc+turb) cond')}
for k,(ls,col,lab) in sty.items():
    if k in D: plt.plot(D[k][0],D[k][1],ls,color=col,lw=2,label=f'{lab}')
plt.plot(fig[:,0],fig[:,1],'s',color='gray',ms=6,label='Wyslouzil isentrope (exp)')
plt.plot(fig[:,0],fig[:,2],'o',color='black',ms=6,label='Wyslouzil cond 1kPa (exp)')
plt.xlim(-1,8.5); plt.ylim(0.15,0.55); plt.xlabel('distance from throat [cm]'); plt.ylabel('p/p0')
plt.legend(fontsize=8); plt.grid(alpha=0.3); plt.title('Wyslouzil Fig.3: condensation effect (inviscid vs viscous+SST)')
plt.tight_layout(); plt.savefig("fig3_compare_visc.png",dpi=120)

print("=== p/p0 at Fig.3 stations ===")
print("x_cm | inv_dry inv_cond | sst_dry sst_cond | exp_iso exp_cond | sst_ratio exp_ratio")
for xq in [2.5,3,4,5,6,8]:
    def g(k):
        if k not in D: return np.nan
        xb,pb=D[k][0],D[k][1]; m=np.isfinite(pb); return np.interp(xq,xb[m],pb[m])
    ei=np.interp(xq,fig[:,0],fig[:,1]); ec=np.interp(xq,fig[:,0],fig[:,2])
    sd,sc=g('sst_dry'),g('sst_cond')
    print("%4.1f |  %.3f   %.3f  |  %.3f   %.3f  |  %.3f  %.3f  |  %.3f    %.3f"%(
        xq,g('inv_dry'),g('inv_cond'),sd,sc,ei,ec, (sc/sd if sd==sd else np.nan), ec/ei))
print("saved fig3_compare_visc.png")
