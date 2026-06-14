import glob, os, numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P0=59070.0
EXP="wyslouzil_fig3_pp0.csv"
# 新 exp データ (x, isentrope, 1.00kPa)。ヘッダ行 + 3列の数値
exp=np.array([[float(t) for t in l.split(",")[:3]]
              for l in open(EXP) if l.strip() and l.split(",")[0].replace(".","",1).isdigit()])
ex, e_iso, e_cond = exp[:,0], exp[:,1], exp[:,2]

def latest(run):
    fs=[f for f in glob.glob(f"{run}/res_*.h5") if "nan" not in f]
    return max(fs,key=lambda f:int(f.split("res_")[1].split(".")[0])) if fs else None
def mesh(run):
    return [p for p in glob.glob(f"{run}/nozzle_*.h5") if "res_" not in os.path.basename(p)][0]
def binavg(x,y,xe):
    xc=0.5*(xe[:-1]+xe[1:]); yb=np.full(len(xc),np.nan)
    for i in range(len(xc)):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>0: yb[i]=np.median(y[s])
    return xc,yb
def centerline(run):
    res=latest(run)
    if res is None: return None
    cc=h5py.File(mesh(run),"r")["/CELLS/centCoords"][:].reshape(-1,3)
    x,y=cc[:,0],cc[:,1]
    P=np.array(h5py.File(res,"r")["/VALUE/P"])
    nb=200; xe=np.linspace(x.min(),x.max(),nb+1); xcb=0.5*(xe[:-1]+xe[1:]); h=np.full(nb,np.nan)
    for i in range(nb):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>3: h[i]=y[s].max()-y[s].min()
    ok=np.isfinite(h); xth=xcb[ok][np.argmin(h[ok])]
    epsy=(y.max()-y.min())*0.02; sel=np.abs(y)<epsy
    xs=x[sel]; o=np.argsort(xs)
    return binavg((xs[o]-xth)*100.0, P[sel][o]/P0, np.linspace(-1,9,55))

runs={'inv_dry':'run_0007_h2o_cond','inv_cond':'run_0008_h2o_cond_on',
      'sst_dry':'run_0009_h2o_sst_dry','sst_cond':'run_0010_h2o_sst_cond'}
D={k:centerline(r) for k,r in runs.items()}

plt.figure(figsize=(8,6))
sty={'inv_dry':('--','tab:blue','forge inviscid dry'),'inv_cond':('--','tab:red','forge inviscid cond'),
     'sst_dry':('-','navy','forge SST(visc+turb) dry'),'sst_cond':('-','crimson','forge SST(visc+turb) cond')}
for k,(ls,col,lab) in sty.items():
    if D[k]: plt.plot(D[k][0],D[k][1],ls,color=col,lw=2,label=lab)
plt.plot(ex,e_iso,'s',color='gray',ms=6,label='Wyslouzil isentrope (exp, new)')
plt.plot(ex,e_cond,'o',color='black',ms=6,label='Wyslouzil 1.00kPa (exp, new)')
plt.xlim(-1,8.5); plt.ylim(0.15,0.55); plt.xlabel('distance from throat [cm]'); plt.ylabel('p/p0')
plt.legend(fontsize=8); plt.grid(alpha=0.3); plt.title('Wyslouzil Fig.3 (re-digitized): forge vs experiment')
plt.tight_layout(); plt.savefig("fig3_compare_new.png",dpi=120)

print("=== p/p0 at experiment x stations ===")
print("x_cm | inv_dry inv_cond | sst_dry sst_cond | exp_iso exp_cond")
def gv(k,xq):
    if not D[k]: return np.nan
    xb,pb=D[k]; m=np.isfinite(pb); return np.interp(xq,xb[m],pb[m])
for xq,ei,ec in zip(ex,e_iso,e_cond):
    print("%4.2f |  %.3f   %.3f  |  %.3f   %.3f  |  %.3f  %.3f"%(
        xq,gv('inv_dry',xq),gv('inv_cond',xq),gv('sst_dry',xq),gv('sst_cond',xq),ei,ec))
print("saved fig3_compare_new.png")
