import glob, os, numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P0=59070.0
exp=np.array([[float(t) for t in l.split(",")[:3]] for l in open("wyslouzil_fig3_pp0.csv")
              if l.strip() and l.split(",")[0].replace(".","",1).isdigit()])
def latest(run):
    fs=[f for f in glob.glob(f"{run}/res_*.h5") if "nan" not in f]
    return max(fs,key=lambda f:int(f.split("res_")[1].split(".")[0])) if fs else None
def mesh(run): return [p for p in glob.glob(f"{run}/nozzle_*.h5") if "res_" not in os.path.basename(p)][0]
def binavg(x,y,xe):
    xc=0.5*(xe[:-1]+xe[1:]); yb=np.full(len(xc),np.nan)
    for i in range(len(xc)):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum(): yb[i]=np.median(y[s])
    return xc,yb
def line(run,field="P"):
    res=latest(run)
    if res is None: return None
    cc=h5py.File(mesh(run),"r")["/CELLS/centCoords"][:].reshape(-1,3); x,y=cc[:,0],cc[:,1]
    F=np.array(h5py.File(res,"r")[f"/VALUE/{field}"])
    nb=200; xe=np.linspace(x.min(),x.max(),nb+1); xcb=0.5*(xe[:-1]+xe[1:]); h=np.full(nb,np.nan)
    for i in range(nb):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>3: h[i]=y[s].max()-y[s].min()
    ok=np.isfinite(h); xth=xcb[ok][np.argmin(h[ok])]
    sel=np.abs(y)<(y.max()-y.min())*0.02; xs=x[sel]; o=np.argsort(xs)
    sc=1.0/P0 if field=="P" else 1.0
    return binavg((xs[o]-xth)*100.0,F[sel][o]*sc,np.linspace(-1,9,60))
PAIRS=[("run_0010_h2o_sst_cond","isoCNT+HK, 1-temp","tab:blue","-"),
       ("run_0018_h2o_sst_nokw_hk_2T","isoCNT+HK, 2-temp (Td)","tab:blue","--"),
       ("run_0011_h2o_sst_kw_hk","Kantrowitz+HK, 1-temp","tab:red","-"),
       ("run_0019_h2o_sst_kw_hk_2T","Kantrowitz+HK, 2-temp (Td)","tab:red","--")]
fig,(ax,ax2)=plt.subplots(1,2,figsize=(15,6))
d=line("run_0009_h2o_sst_dry","P"); ax.plot(d[0],d[1],'-',color='gray',lw=1.2,label='dry')
for run,lab,col,ls in PAIRS:
    p=line(run,"P")
    if p: ax.plot(p[0],p[1],ls,color=col,lw=2,label=lab)
    g=line(run,"g_0")
    if g: ax2.plot(g[0],g[1]*1e3,ls,color=col,lw=2,label=lab)
ax.plot(exp[:,0],exp[:,1],'s',color='k',ms=5,mfc='none',label='exp isentrope')
ax.plot(exp[:,0],exp[:,2],'o',color='k',ms=5,label='exp 1.00kPa')
ax.set_xlim(-1,8.5);ax.set_ylim(0.15,0.55);ax.set_xlabel('x from throat [cm]');ax.set_ylabel('p/p0')
ax.legend(fontsize=8);ax.grid(alpha=0.3);ax.set_title('Droplet temperature (1-temp vs 2-temp), Hertz-Knudsen, SST')
ax2.set_xlim(-1,8.5);ax2.set_xlabel('x from throat [cm]');ax2.set_ylabel('condensed mass frac g x1e3')
ax2.legend(fontsize=8);ax2.grid(alpha=0.3);ax2.set_title('condensate g (1-temp solid vs 2-temp dashed)')
plt.tight_layout();plt.savefig('compare_2temp.png',dpi=110)
print("=== p/p0 at exp stations ===")
print(" | ".join(f"{h:>10}" for h in ["x","iso 1T","iso 2T","Kant 1T","Kant 2T","exp_cond"]))
cache={r:line(r,"P") for r,_,_,_ in PAIRS}
for xq,ei,ec in zip(exp[:,0],exp[:,1],exp[:,2]):
    vals=[xq]
    for run,_,_,_ in PAIRS:
        dd=cache[run]; m=np.isfinite(dd[1]); vals.append(np.interp(xq,dd[0][m],dd[1][m]))
    vals.append(ec)
    print(" | ".join(f"{v:10.3f}" for v in vals))
print("saved compare_2temp.png")
