"""r 重み幾何の離散閉性 Σ_f r_f S_f が壁セル AR でどう崩れるかを測る (幾何のみ)。

厳密には ∮ r n ds = (0, A_planar)。forge は双対面を 1 枚に lump して面重心半径 r_f を掛けるため、
kink した双対面 (median-dual の 2 セグメント) では r_f S_total ≠ Σ_seg r_mid S_seg となり誤差が出る。
薄い壁 CV では目標値 A_planar が個々の項より桁違いに小さいので相対誤差が AR に比例して増える。
"""
import subprocess, sys, os
from pathlib import Path
import numpy as np, h5py

W = Path("/home/sano/work/forge-axisnode"); sys.path.insert(0, str(W/"design"))
from forge_design.meshing.mesh2d import write_msh41_2d, _radial_fracs
BUILD = W/"solver_density_cuda/build"; ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
HERE = Path(__file__).resolve().parent
L, R, NI, NJ = 1.0, 0.2, 41, 41

def build(first_frac, tag, taper=0.0):
    rd = HERE/f"closure_{tag}"; rd.mkdir(exist_ok=True)
    xs = np.linspace(0, L, NI); s = _radial_fracs(NJ, first_frac)
    rw = R*(1.0 - taper*xs/L)                       # taper>0: 壁が傾く (ノズル模擬 → 双対面の kink)
    X, Rr = np.meshgrid(xs, s, indexing="ij"); Rr = Rr*rw[:, None]
    coords = np.zeros((NI*NJ, 3)); coords[:,0]=X.ravel(); coords[:,1]=Rr.ravel()
    nid = lambda i,j: i*NJ+j
    quads = np.array([(nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)) for i in range(NI-1) for j in range(NJ-1)],dtype=np.int64)
    be = {"axis":[(nid(i,0),nid(i+1,0)) for i in range(NI-1)], "wall":[(nid(i,NJ-1),nid(i+1,NJ-1)) for i in range(NI-1)],
          "inlet":[(nid(0,j),nid(0,j+1)) for j in range(NJ-1)], "outlet":[(nid(NI-1,j),nid(NI-1,j+1)) for j in range(NJ-1)]}
    write_msh41_2d(rd/"duct.msh", coords, quads, be)
    (rd/"solverConfig.yaml").write_text(
        'mesh: {meshFormat: "hdf5", discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, meshFileName: "duct.h5", valueFileName: "duct.h5"}\n'
        'gpu: 1\nsolver: "SLAU"\n'
        'physProp: {isCompressible: 1, thermalMethod: 0, viscMethod: 0, ro: 1.2, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}\n'
        'time:\n  unsteady: 1\n  dualTime: 0\n  last: {control: 0, nStepOuter: 1}\n'
        '  deltaT: {control: 0, dt: 1e-6, cfl: 1.0, cfl_pseudo: 1.0, dt_min: 1e-12, dt_max: 1.0, blockDPLUR: 0, lowMachPrecond: 0, detectNaN: 0}\n'
        '  outStepStart: 0\n  outStepInterval: 1\n  timeIntegration: 1\n  nStepInner: 1\n'
        'space: {convMethod: 0, limiter: 2}\nturbulence: {model: "none"}\ninitial: "uniform_p101325_u10"\n')
    (rd/"bcondConfig.yaml").write_text("inlet: {physID: 1, kind: slip, outputHDFflg: 0, ints: , floats: }\noutlet: {physID: 2, kind: slip, outputHDFflg: 0, ints: , floats: }\nwall: {physID: 3, kind: wall, outputHDFflg: 0, ints: , floats: }\naxis: {physID: 4, kind: axis, outputHDFflg: 0, ints: , floats: }\n")
    subprocess.run([str(BUILD/"convertGmshToForge"),"duct.msh","duct.h5"],cwd=rd,env=ENV,check=True,capture_output=True,text=True)
    return rd

def closure(rd, taper=0.0):
    f = h5py.File(rd/"duct.h5")
    st=f["/PLANES/STRUCT"][:]; sv=f["/PLANES/surfVect"][:].reshape(-1,3).astype(np.float64)
    pc=f["/PLANES/centCoords"][:].reshape(-1,3).astype(np.float64); vol=f["/CELLS/volume"][:].astype(np.float64)
    nc=f["/MESH/COORD"][:].reshape(-1,3).astype(np.float64)
    cells=[]; i=0
    while i<len(st):
        nn=st[i]; i+=1+nn; ncl=st[i]; i+=1; cells.append(st[i:i+ncl].tolist()); i+=ncl
    n=len(vol); rSx=np.zeros(n); rSy=np.zeros(n); absS=np.zeros(n)
    for ip,cl in enumerate(cells):
        r=max(pc[ip,1],0.0); tx,ty=r*sv[ip,0],r*sv[ip,1]
        if len(cl)>0 and cl[0]<n: rSx[cl[0]]+=tx; rSy[cl[0]]+=ty; absS[cl[0]]+=abs(ty)
        if len(cl)>1 and cl[1]<n: rSx[cl[1]]-=tx; rSy[cl[1]]-=ty; absS[cl[1]]+=abs(ty)
    err=np.abs(rSy-vol)/vol
    xs=np.unique(np.round(nc[:,0],9)); rw={x: nc[np.abs(nc[:,0]-x)<1e-9,1].max() for x in xs}
    wall=np.array([nc[i,1] >= rw[round(nc[i,0],9)]-1e-12 for i in range(len(nc))])
    return err, absS/vol, wall, nc

if __name__=="__main__":
    print(f"{'first_frac':>10} {'wall Δr [m]':>12} {'AR(壁)':>9} {'max err (全体)':>14} {'wall CV err':>12} {'打消比':>10} {'eps32×打消比':>12}")
    for ff in [0.25, 0.05, 0.01, 2e-3, 4e-4, 1e-4]:
        for taper in (0.0, 0.3):
            rd = build(ff, f"ff{ff:g}_t{taper:g}", taper)
            err, cancel, wall, nc = closure(rd, taper)
            dr = R*ff; ar = (L/(NI-1))/dr
            print(f"{ff:10g} {'taper' if taper else 'plain':>7} {dr:12.3e} {ar:9.1f} {err.max():14.3e} {err[wall].max():12.3e} {cancel[wall].max():10.2e} {cancel[wall].max()*1.2e-7:12.2e}")
