"""物理出口 (出口 BC artefact を避け --dx r_t 上流の断面) の軸値・質量平均を表にする。
usage: design/.venv-opt/bin/python exit_table_va.py <glob> [--dx 2.0] [--out study_xxx_exit.json]"""
import argparse, json, h5py, numpy as np
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument("glob"); ap.add_argument("--dx",type=float,default=2.0); ap.add_argument("--out",default=None); a=ap.parse_args()
tz=np.trapezoid; rows=[]; Yw=0.0286647
for r in sorted(p for p in Path(".").glob(a.glob) if p.is_dir()):
    Md=float(r.name.split("_M")[1].split("_")[0]); kind=r.name.split("_")[-1]
    info=json.load(open(r/"prepare_info.json")); S=info["scale_m"]
    nc=h5py.File(r/"nozzle.h5")["MESH/COORD"][:].reshape(-1,3); f=h5py.File(r/"res_12000.h5")
    xe=nc[:,0].max()-a.dx*S; cols=np.unique(nc[:,0]); xc=cols[np.argmin(abs(cols-xe))]
    sel=np.abs(nc[:,0]-xc)<1e-9; o=np.argsort(nc[sel,1]); idx=np.where(sel)[0][o]
    rr=nc[idx,1]/S; ro=f["VALUE/ro"][:][idx]; Ux=f["VALUE/Ux"][:][idx]; Uy=f["VALUE/Uy"][:][idx]; so=f["VALUE/sonic"][:][idx]; T=f["VALUE/T"][:][idx]; P=f["VALUE/P"][:][idx]
    M=np.hypot(Ux,Uy)/so; g=f["VALUE/g_0"][:][idx] if "VALUE/g_0" in f else np.zeros_like(M); Sx=f["VALUE/condS_0"][:][idx] if "VALUE/condS_0" in f else np.zeros_like(M)
    w=ro*Ux*rr; Ma=tz(w*M,rr)/tz(w,rr); Ta=tz(w*T,rr)/tz(w,rr); ga=tz(w*g,rr)/tz(w,rr)
    core=rr<0.85*rr.max(); eps=100*(M[core].max()-M[core].min())/Ma
    # 軸履歴: onset (g>1e-4) と S max
    ax=nc[:,1]<1e-9; ia=np.where(ax)[0][np.argsort(nc[ax,0])]; xa=nc[ia,0]/S
    ga_=f["VALUE/g_0"][:][ia] if "VALUE/g_0" in f else np.zeros_like(xa); Sa=f["VALUE/condS_0"][:][ia] if "VALUE/condS_0" in f else np.zeros_like(xa); Tax=f["VALUE/T"][:][ia]
    on=np.where(ga_>1e-4)[0]; onset=(xa[on[0]],Tax[on[0]]) if len(on) else (float("nan"),float("nan"))
    rows.append({k:(float(v) if not isinstance(v,str) else v) for k,v in dict(Md=Md,kind=kind,run=r.name,x=xc/S,M_axis=M[0],T_axis=T[0],P_axis=P[0],S_axis=Sx[0],g_axis=g[0]/Yw,M_avg=Ma,T_avg=Ta,g_avg=ga/Yw,eps=eps,S_max=Sa[xa>0].max(),onset_x=onset[0],onset_T=onset[1],T_min=Tax[xa>0].min()).items()})
print("| M_d | 種別 | run | 評価断面 x/r_t | 軸 M | 軸 T [K] | 軸 P [Pa] | 軸 S | 軸 g/Y_H2O | **質量平均 M** | 質量平均 T | 質量平均 g/Y | コア M 幅 | 軸 S max | onset x (T) |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for w in sorted(rows,key=lambda w:(w["Md"],{"dry":0,"noneq":1,"eq":2}[w["kind"]])):
    on="—" if np.isnan(w["onset_x"]) else f"{w['onset_x']:.1f} ({w['onset_T']:.0f} K)"
    print(f"| {w['Md']} | {w['kind']} | `{w['run']}` | {w['x']:.1f} | {w['M_axis']:.3f} | {w['T_axis']:.1f} | {w['P_axis']:.0f} | {w['S_axis']:.2f} | {w['g_axis']:.2f} | **{w['M_avg']:.3f}** | {w['T_avg']:.1f} | {w['g_avg']:.2f} | {w['eps']:.1f} % | {w['S_max']:.1f} | {on} |")
if a.out: json.dump(rows,open(a.out,"w"),indent=1)
