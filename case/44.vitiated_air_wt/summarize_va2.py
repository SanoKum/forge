"""va2 バッチ (run_0031–0048) の表と図。axis_values.csv (axis_csv_va.py) を全 run に生成してから集計。
cd case/44.vitiated_air_wt && python3 axis_csv_va.py run_00{31..48}_va2_* && design/.venv-opt/bin/python summarize_va2.py"""
import json, glob
from pathlib import Path
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
try: fm.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf"); matplotlib.rcParams["font.family"]=["Noto Sans CJK JP","DejaVu Sans"]
except Exception: pass
CASE=Path(".")
def load(r):
    t=(r/"axis_values.csv").read_text().splitlines(); h=t[1].split(","); d=np.loadtxt(t[2:],delimiter=","); q={k:d[:,i] for i,k in enumerate(h)}
    q["M"]=np.hypot(q["Ux"],q["Uy"])/np.sqrt(np.maximum(q["sonic"]**2 if "sonic" in q else 1,1e-30)) if "sonic" in q else None
    return q
import sys
TAG=sys.argv[1] if len(sys.argv)>1 else "va2"     # va2 | va2fx (va2fx は dry を va2lp から取る)
if TAG=="va2fx":
    runs=sorted(p for p in list(CASE.glob("run_00[6-7][0-9]_va2fx_*"))+list(CASE.glob("run_00[4-6][0-9]_va2lp_*_dry")) if p.is_dir())
else:
    runs=sorted(p for p in CASE.glob("run_00[3-4][0-9]_va2_*") if p.is_dir())
rows=[]; data={}
for r in runs:
    Md=float(r.name.split("_M")[1].split("_")[0]); kind=r.name.split("_")[-1]
    q=load(r); m=json.load(open(r/"metrics.json")); info=json.load(open(r/"prepare_info.json"))
    x=q["x_over_rt"]; xe=info["x_E"]; ie=np.argmin(abs(x-xe)); il=len(x)-1
    g=q.get("g_0",np.zeros_like(x)); S=q.get("condS_0",q["S_post"]); Ts=q.get("condTsat_0",q["Tsat_post"])
    on=np.where(g>1e-4)[0]; onset=(x[on[0]],q["T"][on[0]]) if len(on) else (np.nan,np.nan)
    Yw=0.0286647
    rows.append(dict(Md=Md,kind=kind,run=r.name,T_min=q["T"][x>0].min(),T_exit=q["T"][il],P_exit=q["P"][il],
        S_max=S[x>0].max(),g_exit=g[il],g_frac=g[il]/Yw,g_max=g.max(),onset_x=onset[0],onset_T=onset[1],
        subcool_max=(Ts-q["T"])[x>0].max(),M_exit=m["M_axis_exit"],dM=m["dM_max_rel_Md"],epsM=m["exit_uniformity"].get("eps_M_rms",np.nan)))
    data[(Md,kind)]=q
print("| M_d | 種別 | run | 軸 T min [K] | 出口 T [K] | 出口 P [Pa] | S max | 過冷却 max [K] | g 出口 (H₂O 比) | onset x/r_t (T) | 出口軸 M | ‖ΔM‖∞ | ε_M |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for w in sorted(rows,key=lambda w:(w["Md"],{"dry":0,"noneq":1,"eq":2}[w["kind"]])):
    on = "—" if np.isnan(w["onset_x"]) else f"{w['onset_x']:.1f} ({w['onset_T']:.0f} K)"
    print(f"| {w['Md']} | {w['kind']} | `{w['run']}` | {w['T_min']:.1f} | {w['T_exit']:.1f} | {w['P_exit']:.0f} | {w['S_max']:.2f} | {w['subcool_max']:.1f} | {w['g_exit']:.4f} ({100*w['g_frac']:.0f} %) | {on} | {w['M_exit']:.3f} | {100*w['dM']:.2f} % | {100*w['epsM']:.2f} % |")
json.dump(rows,open(f"study_{TAG}_condensation.json","w"),indent=1)
Mds=sorted({w["Md"] for w in rows}); fig,axs=plt.subplots(3,len(Mds),figsize=(3.3*len(Mds),9),sharex="col")
for j,Md in enumerate(Mds):
    for kind,c in (("dry","k"),("noneq","C0"),("eq","C3")):
        q=data.get((Md,kind)); 
        if q is None: continue
        x=q["x_over_rt"]; m=x>0
        axs[0,j].plot(x[m],q["T"][m],c,lw=1,label=kind)
        if kind!="dry":
            axs[0,j].plot(x[m],q["Tsat_post"][m],c,ls=":",lw=0.8,label=kind+" Tsat" if kind=="noneq" else None)
            axs[1,j].plot(x[m],q["condS_0"][m],c,lw=1,label=kind); axs[2,j].plot(x[m],q["g_0"][m]/0.0286647,c,lw=1,label=kind)
    axs[0,j].set_title(f"M_d {Md}"); axs[1,j].axhline(1,c="gray",lw=.5); axs[1,j].set_yscale("log"); axs[2,j].set_xlabel("x/r_t")
    for a in axs[:,j]: a.grid(alpha=.3)
axs[0,0].set_ylabel("axis T [K] (dotted: T_sat(p_v))"); axs[1,0].set_ylabel("S = p_v/p_sat"); axs[2,0].set_ylabel("g / Y_H2O (凝縮率)")
axs[0,0].legend(fontsize=7); axs[1,0].legend(fontsize=7)
fig.suptitle(f"case/44 {TAG} (Pt 1.137 MPa, Tt 1058 K): dry / 非平衡 (Kw+HK) / 平衡凝縮 — 軸上分布")
fig.tight_layout(); fig.savefig(f"figs/{TAG}_condensation_axis.png",dpi=120); print(f"wrote figs/{TAG}_condensation_axis.png")
