"""壁側に幾何級数で伸縮した直管で、線形場 (ρ=const, u=U0−2ax, v=ar) の 1 step 質量残差を行ごとに測る。
値位置 (centroid / node) × 勾配 (GG / LSQ) × 次数 の組合せで、伸縮メッシュ内部行の O(1) 誤差の有無を見る。"""
import os, subprocess, sys
from pathlib import Path
import numpy as np, h5py
W = Path("/home/sano/work/forge-axisnode"); sys.path.insert(0, str(W/"design"))
from forge_design.meshing.mesh2d import write_msh41_2d, _radial_fracs
BUILD = W/"solver_density_cuda/build"; TOOLS = W/"solver_density_cuda/tools"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial"); HERE = Path(__file__).resolve().parent
L, R, NI, NJ, DT = 1.0, 0.2, 41, 41, 1.0e-3
RHO0, U0, A, P0, GAMMA = 1.2, 1200.0, 500.0, 101325.0, 1.4

def cfg(conv, lsq, van, lim, em=0):
    return f"""mesh: {{meshFormat: "hdf5", discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, nodeAxisDirichlet: 0, nodeValueAtNode: {van}, gradLSQ: {lsq}, nodeReconEdgeMidpoint: {em}, meshFileName: "duct.h5", valueFileName: "duct.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 0, ro: 1.2, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}}
time:
  unsteady: 1
  dualTime: 0
  last: {{control: 0, nStepOuter: 1}}
  deltaT: {{control: 0, dt: {DT}, cfl: 1.0, cfl_pseudo: 1.0, dt_min: 1e-14, dt_max: 1.0, blockDPLUR: 0, lowMachPrecond: 0, detectNaN: 0}}
  outStepStart: 0
  outStepInterval: 1
  timeIntegration: 1
  nStepInner: 1
space: {{convMethod: {conv}, limiter: {lim}}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
"""
BC=("inlet: {physID: 1, kind: slip, outputHDFflg: 0, ints: , floats: }\noutlet: {physID: 2, kind: slip, outputHDFflg: 0, ints: , floats: }\n"
    "wall: {physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }\naxis: {physID: 4, kind: axis, outputHDFflg: 0, ints: , floats: }\n")

def run(first_frac, conv, lsq, van, lim, conv_ic, em=0):
    rd = HERE/f"st_ff{first_frac:g}_c{conv}_l{lsq}_v{van}_lim{lim}_{conv_ic}_em{em}"; rd.mkdir(exist_ok=True)
    xs=np.linspace(0,L,NI); s=_radial_fracs(NJ, first_frac); X,Rr=np.meshgrid(xs,R*s,indexing="ij")
    coords=np.zeros((NI*NJ,3)); coords[:,0]=X.ravel(); coords[:,1]=Rr.ravel()
    nid=lambda i,j:i*NJ+j
    quads=np.array([(nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)) for i in range(NI-1) for j in range(NJ-1)],dtype=np.int64)
    be={"axis":[(nid(i,0),nid(i+1,0)) for i in range(NI-1)],"wall":[(nid(i,NJ-1),nid(i+1,NJ-1)) for i in range(NI-1)],
        "inlet":[(nid(0,j),nid(0,j+1)) for j in range(NJ-1)],"outlet":[(nid(NI-1,j),nid(NI-1,j+1)) for j in range(NJ-1)]}
    write_msh41_2d(rd/"duct.msh",coords,quads,be)
    (rd/"solverConfig.yaml").write_text(cfg(conv,lsq,van,lim,em)); (rd/"bcondConfig.yaml").write_text(BC)
    (rd/"probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    subprocess.run([str(BUILD/"convertGmshToForge"),"duct.msh","duct.h5"],cwd=rd,env=ENV,check=True,capture_output=True,text=True)
    with h5py.File(rd/"duct.h5","r+") as f:
        cc=f["/CELLS/centCoords"][:].reshape(-1,3); nc=f["/MESH/COORD"][:].reshape(-1,3)
        pos = nc if conv_ic=="node" else cc
        x,r=pos[:,0],pos[:,1]; u,v=U0-2*A*x, A*r; n=len(x)
        f["/VALUE/ro"][:]=np.full(n,RHO0,np.float32); f["/VALUE/roUx"][:]=(RHO0*u).astype(np.float32)
        f["/VALUE/roUy"][:]=(RHO0*v).astype(np.float32); f["/VALUE/roUz"][:]=np.zeros(n,np.float32)
        f["/VALUE/roe"][:]=(P0/(GAMMA-1)+0.5*RHO0*(u*u+v*v)).astype(np.float32)
    for q in rd.glob("res_*"): q.unlink()
    subprocess.run([str(TOOLS/"run_case.sh"),str(rd)],env=ENV,capture_output=True,text=True)
    with h5py.File(rd/"res_0.h5") as f0, h5py.File(rd/"res_1.h5") as f1:
        e=(f1["/VALUE/ro"][:].astype(float)-f0["/VALUE/ro"][:].astype(float))/DT/(RHO0*A)
    with h5py.File(rd/"duct.h5") as f: nc=f["/MESH/COORD"][:].reshape(-1,3)
    inner=(nc[:,0]>0.2)&(nc[:,0]<0.8)
    rs=np.unique(np.round(nc[:,1],9))
    # 壁から数えて 1..4 行目 (壁ノード自体は除く) と中央部
    rows={}
    for k in range(1,5):
        m=inner&(np.abs(nc[:,1]-rs[-1-k])<1e-9); rows[f"wall-{k}"]=np.abs(e[m]).max()
    m=inner&(nc[:,1]>0.05)&(nc[:,1]<0.15); rows["mid"]=np.abs(e[m]).max()
    m=inner&(nc[:,1]<1e-12); rows["axis"]=np.abs(e[m]).max()
    return rows

if __name__=="__main__":
    ff=1e-3   # 壁 AR ~ 125, 幾何級数比 q ~ 1.2
    print(f"first_frac={ff}  無次元質量残差 max|e| (線形軸対称場, 厳密値 0)")
    print(f"{'値位置':>8} {'IC規約':>8} {'conv':>4} {'grad':>4} {'lim':>3} | {'wall-1':>9} {'wall-2':>9} {'wall-3':>9} {'wall-4':>9} {'mid':>9} {'axis':>9}")
    import sys
    combos = [(0,"centroid",0),(1,"node",0),(1,"node",1)] if "--em" in sys.argv else [(0,"centroid",0),(1,"node",0)]
    for van,ic,em in combos:
        for conv,lsq,lim in [(0,0,0),(1,0,0),(1,2,0),(1,0,2),(1,2,2)]:
            r=run(ff,conv,lsq,van,lim,ic,em)
            print(f"{('node' if van else 'centroid')+('+em' if em else ''):>10} {ic:>8} {conv:>4} {'LSQ' if lsq else 'GG':>4} {lim:>3} | " + " ".join(f"{r[k]:9.2e}" for k in ["wall-1","wall-2","wall-3","wall-4","mid","axis"]))
