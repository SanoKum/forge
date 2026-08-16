"""自由流保持テスト: 一様圧・静止 (p=const, u=0) で半径運動量残差が出ないか。

r 重み方式の半径運動量は  res_roUy = −Σ_f p·r_f S_f,y + p·A  なので、一様圧では
A = Σ_f r_f S_f,y でないと偽の半径力が残る。A に解析 A_planar を使う (既定) と float32
メトリックの桁落ちぶんだけ残り、離散閉性 (`hoopAreaFromClosure: 1`) では厳密に消える。
測るのは 1 step 陽解法の Δ(ρu_r)/Δt を ρ で割った偽半径加速度 [m/s²]。
"""
import os, subprocess, sys
from pathlib import Path
import numpy as np, h5py

W = Path("/home/sano/work/forge-axisnode"); sys.path.insert(0, str(W/"design"))
from forge_design.meshing.mesh2d import write_msh41_2d, _radial_fracs
BUILD = W/"solver_density_cuda/build"; TOOLS = W/"solver_density_cuda/tools"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
HERE = Path(__file__).resolve().parent
L, R, NI, NJ, DT = 1.0, 0.2, 41, 41, 1.0e-6
P0, RO0, GAMMA, CP = 1.0e6, 4.0, 1.4, 1004.5   # 生産 NS と同じ圧力・密度スケール

CFG = """mesh: {{meshFormat: "hdf5", discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, hoopAreaFromClosure: {clo}, meshFileName: "duct.h5", valueFileName: "duct.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 0, ro: 1.2, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}}
time:
  unsteady: 1
  dualTime: 0
  last: {{control: 0, nStepOuter: 1}}
  deltaT: {{control: 0, dt: {dt}, cfl: 1.0, cfl_pseudo: 1.0, dt_min: 1e-14, dt_max: 1.0, blockDPLUR: 0, lowMachPrecond: 0, detectNaN: 0}}
  outStepStart: 0
  outStepInterval: 1
  timeIntegration: 1
  nStepInner: 1
space: {{convMethod: 0, limiter: 2, pRef: {pref}}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
"""
BC = ("inlet: {physID: 1, kind: slip, outputHDFflg: 0, ints: , floats: }\n"
      "outlet: {physID: 2, kind: slip, outputHDFflg: 0, ints: , floats: }\n"
      "wall: {physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }\n"
      "axis: {physID: 4, kind: axis, outputHDFflg: 0, ints: , floats: }\n")

def run(first_frac, closure, pref=0.0):
    rd = HERE/f"fs_ff{first_frac:g}_clo{closure}_p{pref:g}"; rd.mkdir(exist_ok=True)
    xs = np.linspace(0, L, NI); s = _radial_fracs(NJ, first_frac)
    X, Rr = np.meshgrid(xs, R*s, indexing="ij")
    coords = np.zeros((NI*NJ, 3)); coords[:,0]=X.ravel(); coords[:,1]=Rr.ravel()
    nid = lambda i,j: i*NJ+j
    quads = np.array([(nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)) for i in range(NI-1) for j in range(NJ-1)],dtype=np.int64)
    be = {"axis":[(nid(i,0),nid(i+1,0)) for i in range(NI-1)], "wall":[(nid(i,NJ-1),nid(i+1,NJ-1)) for i in range(NI-1)],
          "inlet":[(nid(0,j),nid(0,j+1)) for j in range(NJ-1)], "outlet":[(nid(NI-1,j),nid(NI-1,j+1)) for j in range(NJ-1)]}
    write_msh41_2d(rd/"duct.msh", coords, quads, be)
    (rd/"solverConfig.yaml").write_text(CFG.format(clo=closure, dt=DT, pref=pref))
    (rd/"bcondConfig.yaml").write_text(BC)
    (rd/"probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    subprocess.run([str(BUILD/"convertGmshToForge"),"duct.msh","duct.h5"],cwd=rd,env=ENV,check=True,capture_output=True,text=True)
    with h5py.File(rd/"duct.h5","r+") as f:                     # 一様静止・一様圧
        n=len(f["/VALUE/ro"]); z=np.zeros(n,dtype=np.float32)
        f["/VALUE/ro"][:]=np.full(n,RO0,dtype=np.float32)
        f["/VALUE/roUx"][:]=z; f["/VALUE/roUy"][:]=z; f["/VALUE/roUz"][:]=z
        f["/VALUE/roe"][:]=np.full(n,P0/(GAMMA-1),dtype=np.float32)
    for x in rd.glob("res_*"): x.unlink()
    subprocess.run([str(TOOLS/"run_case.sh"),str(rd)],env=ENV,capture_output=True,text=True)
    with h5py.File(rd/"res_0.h5") as f0, h5py.File(rd/"res_1.h5") as f1:
        d=(f1["/VALUE/roUy"][:].astype(float)-f0["/VALUE/roUy"][:].astype(float))/DT/RO0
        dx=(f1["/VALUE/roUx"][:].astype(float)-f0["/VALUE/roUx"][:].astype(float))/DT/RO0
    with h5py.File(rd/"duct.h5") as f:
        nc=f["/MESH/COORD"][:].reshape(-1,3)
    wall=nc[:,1].astype(float)>=float(nc[:,1].max())-1e-6
    return np.abs(d).max(), np.abs(d[wall]).max(), np.abs(dx).max()

if __name__=="__main__":
    print(f"{'first_frac':>10} {'AR(壁)':>8} {'A の取り方':>14} {'|a_r| max':>12} {'壁 CV':>12} {'|a_x| max':>12}   [m/s²]")
    import sys
    prefs = [0.0, P0] if "--pref" in sys.argv else [0.0]
    for ff in [0.05, 2e-3, 1e-4]:
        for pref in prefs:
            for clo in (0,1):
                a,aw,ax = run(ff, clo, pref)
                tag = ('解析 A' if clo==0 else '離散閉性 A') + (f' +pRef={pref:g}' if pref else '')
                print(f"{ff:10g} {(L/(NI-1))/(R*ff):8.1f} {tag:>24} {a:12.3e} {aw:12.3e} {ax:12.3e}")
