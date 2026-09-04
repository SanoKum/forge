#!/usr/bin/env python3
"""case/48 Cabra H2/N2 浮き上がり火炎: run ディレクトリの config/BC/mesh 変換/IC/入口プロファイルを作る。
  .venv-chem/bin/python setup_cabra_case.py run_dir [--chem 0|1] [--jac 2] [--ji 5] [--tci 0] [--cfl 1] [--conv 1] [--relax 0.7] [--nstep N] [--out N] [--restart bk.h5] [--kfac 1]
ジェット: H2 25 %/N2 75 % (vol), 305 K, バルク 107 m/s (管内 1/7 乗則, 中心 131 m/s); coflow: 1045 K, 3.5 m/s, X O2 .15/H2O .099/N2 .751。1 atm。"""
import argparse, pathlib, shutil, subprocess, os, h5py, numpy as np, yaml, cantera as ct
ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--chem", type=int, default=0); ap.add_argument("--jac", type=int, default=2); ap.add_argument("--ji", type=int, default=5)
ap.add_argument("--tci", type=int, default=0); ap.add_argument("--cfl", type=float, default=1.0); ap.add_argument("--conv", type=int, default=1); ap.add_argument("--relax", type=float, default=0.7)
ap.add_argument("--nstep", type=int, default=20000); ap.add_argument("--out", type=int, default=0); ap.add_argument("--restart", default=None); ap.add_argument("--kfac", type=float, default=1.0)
ap.add_argument("--eps", type=float, default=0.15); ap.add_argument("--precond", type=int, default=2); ap.add_argument("--axisdir", type=int, default=0); ap.add_argument("--nojet", type=int, default=0); ap.add_argument("--single", type=int, default=0); ap.add_argument("--coupling", type=int, default=2); ap.add_argument("--iccol", type=int, default=1); ap.add_argument("--planar", type=int, default=0); ap.add_argument("--jetn2", type=int, default=0); ap.add_argument("--jetcof", type=int, default=0); ap.add_argument("--far", default="slip"); ap.add_argument("--sfr", type=int, default=0); a = ap.parse_args()
HERE = pathlib.Path(__file__).parent; d = pathlib.Path(a.run_dir); d.mkdir(exist_ok=True)
names = ["H2", "O2", "H", "O", "OH", "H2O", "HO2", "H2O2", "N2"]
g = ct.Solution(str(HERE / "mech.yaml")); assert g.species_names == names
g.TPX = 305.0, 101325.0, {"H2": 0.25, "N2": 0.75}; roJ = g.density; YJ = g.Y.copy(); Ub = 107.0; Uc = Ub/0.8167   # 1/7 乗則: mean/center = 49/60
if a.nojet:
    g.TPX = 1045.0, 101325.0, {"O2": 0.15, "H2O": 0.099, "N2": 0.751}; roJ = g.density; YJ = g.Y.copy(); Ub = 3.5; Uc = 3.5
TJ = 1045.0 if a.nojet else 305.0
g.TPX = 1045.0, 101325.0, {"O2": 0.15, "H2O": 0.099, "N2": 0.751}; roC = g.density; YC = g.Y.copy(); UC = 3.5
if a.jetcof:  # ジェット組成 = coflow 組成 (組成ジャンプなし, 9 種輸送は有効)
    g.TPX = TJ, 101325.0, {"O2": 0.15, "H2O": 0.099, "N2": 0.751}; roJ = g.density; YJ = g.Y.copy()
if a.jetn2:   # ジェットを純 N2 に (組成ジャンプはあるが R/cp コントラストは小)
    g.TPX = TJ, 101325.0, {"N2": 1.0}; roJ = g.density; YJ = g.Y.copy()
if a.single:   # 単一種 N2 (組成結合なし): 温度・速度コントラストだけ残す
    names = ["N2"]; g1 = ct.Solution(str(HERE / "mech.yaml")); g1.TPX = TJ, 101325.0, {"N2": 1.0}; roJ = g1.density; YJ = np.array([1.0]); g1.TPX = 1045.0, 101325.0, {"N2": 1.0}; roC = g1.density; YC = np.array([1.0])
print(f"jet: rho={roJ:.4f} Uc={Uc:.1f}; coflow: rho={roC:.4f}; species={names}")
kJ, omJ = a.kfac*1.5*(0.05*Ub)**2, 4.0e4; kC, omC = a.kfac*1.5*(0.05*UC)**2, 100.0
def yf(Y): return ", ".join(f"Y{i}: {Y[i]:.6f}" for i in range(len(names)))
(d / "bcondConfig.yaml").write_text(f"""# Cabra H2/N2 lifted flame (NASA/CR-2004-212887 Table 6.1)
inlet_jet:    {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: {{inletProfile: {0 if a.nojet else 1}}}, floats: {{ro: {roJ:.6f}, Ux: {Ub:.2f}, Uy: 0.0, Uz: 0.0, Ps: 101325.0, k: {kJ:.2f}, omega: {omJ:.1f}, {yf(YJ)}}}}}
inlet_coflow: {{physID: 2, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {{ro: {roC:.6f}, Ux: {UC:.2f}, Uy: 0.0, Uz: 0.0, Ps: 101325.0, k: {kC:.4f}, omega: {omC:.1f}, {yf(YC)}}}}}
outlet:       {{physID: 3, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {{Ps: 101325.0, Pt: 101325.0, Tt: 1045.0}}}}
wall:         {{physID: 4, kind: wall, outputHDFflg: 0, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0}}}}
farfield:     {{physID: 5, kind: {a.far}, outputHDFflg: 0, ints: , floats: {"{Ps: 101325.0, Pt: 101325.0, Tt: 1045.0}" if a.far.startswith("outlet") else ""}}}
axis:         {{physID: 6, kind: {"slip" if a.planar else "axis"}, outputHDFflg: 0, ints: , floats: }}
""")
chem = "" if a.single else f"chemistry: {{enabled: {a.chem}, mechanismFile: \"mech.yaml\", jacobianMode: {a.jac}, jacobianInterval: {a.ji}, tci: {a.tci}, tciTauChem: 1}}"
if a.single: chem = ""
value = a.restart if a.restart else "cabra.h5"
(d / "solverConfig.yaml").write_text(f"""# case/48 Cabra H2/N2: chem={a.chem} cfl_pseudo={a.cfl} conv={a.conv}
mesh: {{meshFormat: "hdf5", discretization: "node", isAxisymmetric: {0 if a.planar else 1}, axisCentroidShift: 1, nodeWallViscGradFlux: 1, nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "cabra.h5", valueFileName: "{value}"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 2, viscMethod: 2, ro: 0.32, visc: 4.0e-5, thermCond: 0.07, cp: 1200.0, gamma: 1.35,
           species: [{", ".join(names)}], speciesDBFile: "species_db.yaml", speciesDiffusionMethod: 1, thermoHrefTemp: 298.15, speciesFaceReconstruction: {a.sfr},
           {chem}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {a.nstep}}}
  deltaT: {{control: 1, dt: 1e-7, cfl: 1.0, cfl_pseudo: {a.cfl}, dt_min: 1e-10, dt_max: 0.01, blockDPLUR: 1, implicitRelax: {a.relax}, lowMachPrecond: {a.precond}, precondEps: {a.eps}, detectNaN: 1, speciesImplicitCoupling: {a.coupling}}}
  outStepStart: 0
  outStepInterval: {a.out if a.out > 0 else min(5000, a.nstep)}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: {a.conv}, limiter: 2}}
turbulence: {{model: "sst", wallTreatmentSST: 1, nodeOmegaWfDirichlet: 1, scalarDiffusion: 1, dilatationCorrection: 0, turbulentSchmidt: 0.7}}
initial: "uniform_p101325_u10"
""")
for f in ("species_db.yaml", "mech.yaml", "probe.yaml"): shutil.copy(HERE / f, d / f)
shutil.copy(HERE / "mesh/cabra.msh", d / "cabra.msh")
# 入口プロファイル (ジェット管入口 x=-Ltube, r∈[0, 2.285 mm]): 1/7 乗則
rj = 0.002285; r = np.concatenate([np.linspace(0, rj*0.9, 30), np.linspace(rj*0.9, rj, 30)[1:]]); u = Uc*np.clip(1 - r/rj, 0, 1)**(1/7)
np.savetxt(d / "inlet_profile_1.csv", np.c_[r, u, 0*u, 0*u], header="y Ux Uy Uz", comments="", fmt="%.7e")
CONV = "/home/sano/work/forge-chem/solver_density_cuda/build/convertGmshToForge"
if not (d / "cabra.h5").exists():
    env = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
    rr = subprocess.run([CONV, "cabra.msh", "cabra.h5"], cwd=d, capture_output=True, text=True, env=env)
    (d / "convert.log").write_text(rr.stdout + rr.stderr); print("convert rc", rr.returncode, [l for l in rr.stdout.splitlines() if "nodeInletCornerWall" in l][:1])
Ru = 8.314462618; Tref = 298.15
db = yaml.safe_load(open(HERE / "species_db.yaml")); W = np.array([db[n]["MW"] for n in names])
def h_molar(sp, T):
    c = sp["nasa9_low"] if T < sp["Tmid"] else sp["nasa9_high"]
    return Ru*T*(-c[0]/T**2 + c[1]*np.log(T)/T + c[2] + c[3]*T/2 + c[4]*T**2/3 + c[5]*T**3/4 + c[6]*T**4/5 + c[7]/T)
W = np.array([db[n]["MW"] for n in names])
def e_mass(Y, T): return sum(Y[i]*(h_molar(db[n], T) - h_molar(db[n], Tref))/W[i] for i, n in enumerate(names)) - (Ru/(1.0/np.sum(Y/W)))*T
if not a.restart:
    with h5py.File(d / "cabra.h5", "r+") as f:
        xyz = f["MESH/COORD"][:].reshape(-1, 3); v = f["VALUE"]; N = v["ro"].shape[0]; x = xyz[:N, 0]; rr_ = xyz[:N, 1]
        jet = (rr_ < rj + 1e-9) & ((x < 0.0) | (a.iccol == 1))   # ジェット柱 (iccol=1: 管内+下流, 0: 管内のみ) を燃料状態に
        wall = v["wall_dist"][:] <= 0.0
        ro = np.where(jet, roJ, roC); U = np.where(jet, Ub, UC); e = np.where(jet, e_mass(YJ, TJ), e_mass(YC, 1045.0))
        U = np.where(wall, 0.0, U)
        v["ro"][:] = ro; v["roUx"][:] = ro*U; v["roUy"][:] = 0.0; v["roUz"][:] = 0.0; v["roe"][:] = ro*(e + 0.5*U**2)
        v["roK"][:] = ro*np.where(jet, kJ, kC); v["roOmega"][:] = ro*np.where(jet, omJ, omC)
        for i in range(len(names)):
            k = f"roY{i}"; val = ro*np.where(jet, YJ[i], YC[i])
            if k in v: v[k][:] = val
            else: v.create_dataset(k, data=val.astype(v["ro"].dtype))
        print(f"IC: N={N} jet nodes={jet.sum()}")
print("done", d)
