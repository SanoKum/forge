#!/usr/bin/env python3
"""case/47 Burrows–Kurkov: run ディレクトリの config/BC/mesh 変換/IC を作る。
  .venv-chem/bin/python setup_bk_case.py run_dir [--chem 0|1] [--tci 0|1] [--cfl 2] [--nstep 20000] [--restart RES.h5]
主流 (vitiated air): M2.44, T 1270 K, p 1 atm, Y(O2 .258, N2 .486, H2O .256); H2: M1, T 254 K, p 1 atm。壁 等温 350 K (VULCAN 2022 に倣う)。"""
import argparse, pathlib, shutil, subprocess, os, h5py, numpy as np, cantera as ct
ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--chem", type=int, default=0); ap.add_argument("--tci", type=int, default=0)
ap.add_argument("--cfl", type=float, default=2.0); ap.add_argument("--nstep", type=int, default=20000); ap.add_argument("--restart", default=None)
ap.add_argument("--jac", type=int, default=2); ap.add_argument("--conv", type=int, default=1); ap.add_argument("--profile", type=int, default=0); ap.add_argument("--out", type=int, default=0); ap.add_argument("--relax", type=float, default=0.7); ap.add_argument("--Twall", type=float, default=350.0); a = ap.parse_args()
HERE = pathlib.Path(__file__).parent; d = pathlib.Path(a.run_dir); d.mkdir(exist_ok=True)
names = ["H2", "O2", "H", "O", "OH", "H2O", "HO2", "H2O2", "N2"]
g = ct.Solution(str(HERE / "mech.yaml")); assert g.species_names == names
g.TPY = 1270.0, 101325.0, {"O2": 0.258, "N2": 0.486, "H2O": 0.256}; aA = g.sound_speed; roA = g.density; UA = 2.44*aA; YA = g.Y.copy()
g.TPY = 254.0, 101325.0, {"H2": 1.0}; aH = g.sound_speed; roH = g.density; UH = 1.0*aH; YH = g.Y.copy()
print(f"air: rho={roA:.5f} U={UA:.1f} (a={aA:.1f}); H2: rho={roH:.5f} U={UH:.1f}")
kA, omA = 1.5*(0.01*UA)**2, 5.0e4; kH, omH = 1.5*(0.02*UH)**2, 5.0e5
def yf(Y): return ", ".join(f"Y{i}: {Y[i]:.6f}" for i in range(len(names)))
(d / "bcondConfig.yaml").write_text(f"""# Burrows-Kurkov (NASA TM X-2828): vitiated air M2.44/1270K/1atm, H2 slot M1/254K/1atm, walls isothermal {a.Twall:.0f} K
inlet_air: {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: {{inletProfile: {a.profile}}}, floats: {{ro: {roA:.6f}, Ux: {UA:.2f}, Uy: 0.0, Uz: 0.0, Ps: 101325.0, k: {kA:.1f}, omega: {omA:.1f}, {yf(YA)}}}}}
inlet_h2:  {{physID: 2, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: {{inletProfile: {a.profile}}}, floats: {{ro: {roH:.6f}, Ux: {UH:.2f}, Uy: 0.0, Uz: 0.0, Ps: 101325.0, k: {kH:.1f}, omega: {omH:.1f}, {yf(YH)}}}}}
outlet:    {{physID: 3, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {{Ps: 101325.0, Pt: 101325.0, Tt: 300.0}}}}
wall:      {{physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: {{wallFunc: 1}}, floats: {{Ts: {a.Twall:.1f}, Ux: 0.0, Uy: 0.0, Uz: 0.0}}}}
""" + ("" if a.profile else "wall_slip: {physID: 5, kind: slip,            outputHDFflg: 0, ints: , floats: }\n"))
chem = f"chemistry: {{enabled: {a.chem}, mechanismFile: \"mech.yaml\", jacobianMode: {a.jac}, tci: {a.tci}}}"
value = a.restart if a.restart else "bk.h5"
(d / "solverConfig.yaml").write_text(f"""# case/47 Burrows-Kurkov: chem={a.chem} tci={a.tci} cfl_pseudo={a.cfl}
mesh: {{meshFormat: "hdf5", discretization: "node", nodeWallViscGradFlux: 1, nodeWallDirichlet: 1, nodeInletCornerWall: {a.profile}, meshFileName: "bk.h5", valueFileName: "{value}"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 2, viscMethod: 2, ro: 0.24, visc: 4.5e-5, thermCond: 0.08, cp: 1300.0, gamma: 1.3,
           species: [{", ".join(names)}], speciesDBFile: "species_db.yaml", speciesDiffusionMethod: 1, thermoHrefTemp: 298.15,
           {chem}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {a.nstep}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: 2.0, cfl_pseudo: {a.cfl}, dt_min: 1e-10, dt_max: 0.001, blockDPLUR: 1, implicitRelax: {a.relax}, lowMachPrecond: 0, detectNaN: 1, speciesImplicitCoupling: 2}}
  outStepStart: 0
  outStepInterval: {a.out if a.out > 0 else min(5000, a.nstep)}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: {a.conv}, limiter: 2}}
turbulence: {{model: "sst", wallTreatmentSST: 1, nodeOmegaWfDirichlet: 1, scalarDiffusion: 1, dilatationCorrection: 0, turbulentSchmidt: 0.5}}
initial: "uniform_p101325_u10"
""")
for f in ("species_db.yaml", "mech.yaml", "probe.yaml"): shutil.copy(HERE / f, d / f)
shutil.copy(HERE / ("mesh/bk_noslip.msh" if a.profile else "mesh/bk.msh"), d / "bk.msh")
CONV = "/home/sano/work/forge-chem/solver_density_cuda/build/convertGmshToForge"
if not (d / "bk.h5").exists():
    env = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
    r = subprocess.run([CONV, "bk.msh", "bk.h5"], cwd=d, capture_output=True, text=True, env=env)
    (d / "convert.log").write_text(r.stdout + r.stderr); print("convert rc", r.returncode, [l for l in r.stdout.splitlines() if "dual" in l or "nodes" in l.lower()][:4])
# IC: 一様主流。スロット (x<0, y<0.004) は H2 状態。sensible datum で roe を組む。
Ru = 8.314462618; Tref = 298.15
import yaml
db = yaml.safe_load(open(HERE / "species_db.yaml")); W = np.array([db[n]["MW"] for n in names])
def h_molar(sp, T):
    c = sp["nasa9_low"] if T < sp["Tmid"] else sp["nasa9_high"]
    return Ru*T*(-c[0]/T**2 + c[1]*np.log(T)/T + c[2] + c[3]*T/2 + c[4]*T**2/3 + c[5]*T**3/4 + c[6]*T**4/5 + c[7]/T)
def e_mass(Y, T): return sum(Y[i]*(h_molar(db[n], T) - h_molar(db[n], Tref))/W[i] for i, n in enumerate(names)) - (Ru/(1.0/np.sum(Y/W)))*T
if not a.restart:
    with h5py.File(d / "bk.h5", "r+") as f:
        xyz = f["MESH/COORD"][:].reshape(-1, 3); v = f["VALUE"]; N = v["ro"].shape[0]
        slot = (xyz[:N, 1] < 0.004 + 1e-9)   # スロット + その下流の帯 (y<4 mm) を H2 ジェット状態に (段階起動: 出口の不連続を避ける)
        ro = np.where(slot, roH, roA); U = np.where(slot, UH, UA)
        e = np.where(slot, e_mass(YH, 254.0), e_mass(YA, 1270.0)) + 0.5*U**2
        # 壁ノード (wall_dist==0): nodeWallDirichlet が u=0 に固定するので、u=0・T=T_wall・p=1 atm で整合させる
        # (流速込み roe のままだと内部エネルギーが下限を割り T=0 → P=0 → NaN: run_0001 初回で発覚)。
        wall = v["wall_dist"][:] <= 0.0
        Rmix = np.where(slot, Ru/(1.0/np.sum(YH/W)), Ru/(1.0/np.sum(YA/W)))
        ro_w = 101325.0/(Rmix*a.Twall); e_w = np.where(slot, e_mass(YH, a.Twall), e_mass(YA, a.Twall))
        ro = np.where(wall, ro_w, ro); U = np.where(wall, 0.0, U); e = np.where(wall, e_w, e)
        v["ro"][:] = ro; v["roUx"][:] = ro*U; v["roUy"][:] = 0.0; v["roUz"][:] = 0.0; v["roe"][:] = ro*e
        v["roK"][:] = ro*np.where(slot, kH, kA); v["roOmega"][:] = ro*np.where(slot, omH, omA)
        for i in range(len(names)):
            k = f"roY{i}"; val = ro*np.where(slot, YH[i], YA[i])
            if k in v: v[k][:] = val
            else: v.create_dataset(k, data=val.astype(v["ro"].dtype))
        print(f"IC written: N={N} slot nodes={slot.sum()}")
if a.profile:
    subprocess.run(["/home/sano/work/forge/.venv-chem/bin/python", str(HERE / "gen_bk_inlet_profiles.py"), str(d), "--UA", f"{UA:.3f}", "--UH", f"{UH:.3f}"], check=True)
print("done", d)
