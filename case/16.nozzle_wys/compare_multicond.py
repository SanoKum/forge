import glob, os, numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P0=59070.0
exp=np.array([[float(t) for t in l.split(",")] for l in open("wyslouzil_fig3_multiP.csv")
              if l.strip() and l[0].isdigit()])
ex=exp[:,0]
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
def pline(run):
    res=latest(run)
    if res is None: return None
    cc=h5py.File(mesh(run),"r")["/CELLS/centCoords"][:].reshape(-1,3); x,y=cc[:,0],cc[:,1]
    P=np.array(h5py.File(res,"r")["/VALUE/P"])
    nb=200; xe=np.linspace(x.min(),x.max(),nb+1); xcb=0.5*(xe[:-1]+xe[1:]); h=np.full(nb,np.nan)
    for i in range(nb):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>3: h[i]=y[s].max()-y[s].min()
    ok=np.isfinite(h); xth=xcb[ok][np.argmin(h[ok])]
    sel=np.abs(y)<(y.max()-y.min())*0.02; xs=x[sel]; o=np.argsort(xs)
    return binavg((xs[o]-xth)*100.0,P[sel][o]/P0,np.linspace(-1,9,60))
RUNS=[("run_0013_h2o_sst_kw_gyar","1.00 kPa","tab:red",2),
      ("run_0020_h2o_sst_kg_0p50kPa","0.50 kPa","tab:green",3),
      ("run_0021_h2o_sst_kg_0p26kPa","0.26 kPa","tab:blue",4)]
plt.figure(figsize=(8.5,6.5))
d=pline("run_0009_h2o_sst_dry"); plt.plot(d[0],d[1],'-',color='gray',lw=1.5,label='forge dry (~isentrope)')
plt.plot(ex,exp[:,1],'k-',lw=1,alpha=0.5,label='exp isentrope')
for run,lab,col,ci in RUNS:
    p=pline(run)
    if p: plt.plot(p[0],p[1],'-',color=col,lw=2,label=f'forge {lab}')
    plt.plot(ex,exp[:,ci],'o',color=col,ms=6,label=f'exp {lab}')
plt.xlim(-1,8.5);plt.ylim(0.15,0.55);plt.xlabel('x from throat [cm]');plt.ylabel('p/p0')
plt.legend(fontsize=8,ncol=2);plt.grid(alpha=0.3)
plt.title('Multi-condition generalization: Kantrowitz+Gyarmathy vs Wyslouzil (1.0/0.5/0.26 kPa)')
plt.tight_layout();plt.savefig('compare_multicond.png',dpi=120)
print("=== p/p0 forge vs exp at stations ===")
print(" x_cm | 1.0:forge exp | 0.5:forge exp | 0.26:forge exp")
cache={r:pline(r) for r,_,_,_ in RUNS}
for k,xq in enumerate(ex):
    row=[f"{xq:5.1f}"]
    for (run,_,_,ci) in RUNS:
        d=cache[run]
        fv=np.interp(xq,d[0][np.isfinite(d[1])],d[1][np.isfinite(d[1])]) if d else float('nan')
        row.append(f"{fv:.3f} {exp[k,ci]:.3f}")
    print(" | ".join(row))
print("saved compare_multicond.png")
