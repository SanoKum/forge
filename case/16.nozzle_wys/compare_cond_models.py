#!/usr/bin/env python3
"""凝縮コア4モデルの中心線 p/p0 を Wyslouzil Fig.3 (isentrope/1.00kPa) と比較する。

dry SST (凝縮なし) を基準に、4 モデル (HK / Kw+HK / Gyar / Kw+Gyar) を重ね、
実験の isentrope(■) と cond 1kPa(○) に対する一致を見る。凝縮潜熱で p/p0 が
dry より上振れし 1kPa 実験点に近づく挙動と、onset 位置を確認する。
"""
import glob, os
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

P0 = 59070.0
RUNS = [
    ("run_0048_fig3_2d_sst_dry",   "SST dry",      'k',          '--'),
    ("run_0049_fig3_2d_sst_hk",    "HK",           'tab:blue',   '-'),
    ("run_0050_fig3_2d_sst_kwhk",  "Kw+HK",        'tab:cyan',   '-'),
    ("run_0051_fig3_2d_sst_gyar",  "Gyar",         'tab:orange', '-'),
    ("run_0052_fig3_2d_sst_kwgyar","Kw+Gyar",      'tab:red',    '-'),
]

def latest(run):
    fs = [f for f in glob.glob(f"{run}/res_*.h5") if "nan" not in f]
    return max(fs, key=lambda f: int(f.split("res_")[1].split(".")[0])) if fs else None
def mesh(run):
    return [p for p in glob.glob(f"{run}/nozzle_*.h5") if "res_" not in os.path.basename(p)][0]

def centerline(run):
    res = latest(run)
    if res is None: return None
    cc = h5py.File(mesh(run), "r")["/CELLS/centCoords"][:].reshape(-1, 3)
    x, y, z = cc[:,0], cc[:,1], cc[:,2]
    v = h5py.File(res, "r")["/VALUE"]
    P = np.array(v["P"])
    g = np.array(v["g_0"]) if "g_0" in v else np.zeros_like(P)   # 液相質量分率 (条件で名称差異あり)
    nb=200; xe=np.linspace(x.min(),x.max(),nb+1); xcb=0.5*(xe[:-1]+xe[1:]); h=np.full(nb,np.nan)
    for i in range(nb):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>3: h[i]=y[s].max()-y[s].min()
    ok=np.isfinite(h); xth=xcb[ok][np.argmin(h[ok])]
    zc=0.5*(z.min()+z.max()); dz=z.max()-z.min()
    sel=(np.abs(y)<(y.max()-y.min())*0.02)&(np.abs(z-zc)<max(dz*0.12,1e-6))
    xs=x[sel]; o=np.argsort(xs)
    return (xs[o]-xth)*100.0, P[sel][o]/P0, g[sel][o], os.path.basename(res)

def binavg(x,y,xe):
    xc=0.5*(xe[:-1]+xe[1:]); yb=np.full(len(xc),np.nan)
    for i in range(len(xc)):
        s=(x>=xe[i])&(x<xe[i+1])
        if s.sum()>0: yb[i]=np.median(y[s])
    return xc,yb

def load_exp():
    r=[]
    for ln in open("wyslouzil_fig3_pp0.csv"):
        p=ln.strip().split(",")
        if p and p[0] and p[0][0].isdigit(): r.append([float(p[0]),float(p[1]),float(p[2])])
    a=np.array(r); return a[:,0],a[:,1],a[:,2]

def main():
    ex,eiso,econd=load_exp()
    xe=np.linspace(-1,9,55)
    D={}
    fig,ax=plt.subplots(1,2,figsize=(14,5.2))
    for run,lab,col,ls in RUNS:
        c=centerline(run)
        if c is None: print("MISSING",run); continue
        xb,pb=binavg(c[0],c[1],xe); _,gb=binavg(c[0],c[2],xe)
        D[lab]=(xb,pb,gb,c[3])
        ax[0].plot(xb,pb,ls,color=col,lw=2,label=lab)
        if lab!="SST dry": ax[1].plot(c[0],c[2],ls,color=col,lw=1.5,label=lab)
    ax[0].plot(ex,eiso,'s',color='gray',ms=7,label='exp isentrope')
    ax[0].plot(ex,econd,'o',color='black',ms=6,mfc='none',label='exp cond 1kPa')
    ax[0].set_xlim(-1,8.5); ax[0].set_ylim(0.15,0.55); ax[0].set_xlabel('distance from throat [cm]')
    ax[0].set_ylabel('p/p0'); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8); ax[0].set_title('Centerline p/p0: condensation models vs Fig.3')
    ax[1].set_xlim(-1,8.5); ax[1].set_xlabel('distance from throat [cm]'); ax[1].set_ylabel('condensate g')
    ax[1].grid(alpha=0.3); ax[1].legend(fontsize=8); ax[1].set_title('Condensate mass fraction (onset/growth)')
    fig.tight_layout(); fig.savefig("cond_models_compare.png",dpi=130)
    print("saved cond_models_compare.png")

    print("\n=== p/p0 vs exp cond(1kPa) at stations ===")
    hdr="x_cm| exp_cond |"+ "".join(f" {l:>8}|" for _,l,_,_ in RUNS)
    print(hdr)
    for xq in [2.1,3.1,4.2,5.2,6.2,7.2]:
        ec=np.interp(xq,ex,econd)
        row=f"{xq:4.1f}|  {ec:.3f}  |"
        for _,lab,_,_ in RUNS:
            if lab in D:
                xb,pb=D[lab][0],D[lab][1]; m=np.isfinite(pb); row+=f" {np.interp(xq,xb[m],pb[m]):8.3f}|"
            else: row+="    -   |"
        print(row)

if __name__=="__main__":
    main()
