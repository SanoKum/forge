#!/usr/bin/env python3
"""case/46: 燃焼室 CEA 平衡組成 (A3000_eq, T_c 2788 K, p_c 11.39 bar) の va3 ノズル (case/44 run_0091 メッシュ) に
有限速度化学 (Jachimowski 13 種 + AR, CO2 不活性) を適用する run の IC / BC を作る。
  .venv-chem/bin/python setup_kinetics_case.py run_dir [--chem 0|1] [--jac 2]
IC: run_0091 の収束場 (T_t 1161 K) を T×(T_c/1161), u×sqrt(T_c/1161), p 不変でスケールし、組成は入口組成で一様 (frozen 相当)。"""
import argparse, pathlib, shutil, json, h5py, numpy as np, yaml
ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--chem", type=int, default=1); ap.add_argument("--jac", type=int, default=2)
ap.add_argument("--coupling", type=int, default=2); ap.add_argument("--cfl", type=float, default=2.0); ap.add_argument("--nstep", type=int, default=12000)
a = ap.parse_args()
HERE = pathlib.Path(__file__).parent; SRC = pathlib.Path("/home/sano/work/forge/case/44.vitiated_air_wt/run_0091_va3_M4.19_Lc8_dry")
d = pathlib.Path(a.run_dir); d.mkdir(exist_ok=True)
names = ["H2","O2","H","O","OH","H2O","HO2","H2O2","N2","N","NO","NO2","HNO","AR","CO2"]
cea = json.load(open(HERE/"cea_A3000_eq_chamber.json")); Tc = cea["T_chamber"]; X = cea["X"]
# CO は CO2 へ合算, N2O は落とす (機構外), Ar→AR
Xm = {n: 0.0 for n in names}
for k, v in X.items():
    kk = {"Ar": "AR", "CO": "CO2", "N2O": None}.get(k, k)
    if kk is None or kk not in Xm: continue
    Xm[kk] += v
xs = sum(Xm.values()); Xv = np.array([Xm[n]/xs for n in names])
db = yaml.safe_load(open(HERE/"species_db.yaml")); W = np.array([db[n]["MW"] for n in names]); Wmix = (Xv*W).sum(); Y = Xv*W/Wmix
Ru = 8.314462618; Tref = 298.15
def h_molar(sp, T):
    c = sp["nasa9_low"] if T < sp["Tmid"] else sp["nasa9_high"]
    return Ru*T*(-c[0]/T**2 + c[1]*np.log(T)/T + c[2] + c[3]*T/2 + c[4]*T**2/3 + c[5]*T**3/4 + c[6]*T**4/5 + c[7]/T)
def e_mass(T):  # sensible datum, mixture
    return sum(Y[i]*(h_molar(db[n], T) - h_molar(db[n], Tref))/W[i] for i, n in enumerate(names)) - (Ru/Wmix)*T
# --- IC ---
shutil.copy(SRC/"nozzle.h5", d/"nozzle.h5")
with h5py.File(SRC/"res_12000.h5") as r:
    T0 = r["VALUE/T"][:]; P0 = r["VALUE/P"][:]; Ux0 = r["VALUE/Ux"][:]; Uy0 = r["VALUE/Uy"][:]
sc = Tc/1161.0; T1 = T0*sc; u1x = Ux0*np.sqrt(sc); u1y = Uy0*np.sqrt(sc); P1 = P0
ro1 = P1*Wmix/(Ru*T1)
e1 = np.array([e_mass(t) for t in T1]) + 0.5*(u1x**2 + u1y**2)
with h5py.File(d/"nozzle.h5", "r+") as f:
    v = f["VALUE"]; N = v["ro"].shape[0]
    v["ro"][:] = ro1; v["roUx"][:] = ro1*u1x; v["roUy"][:] = ro1*u1y; v["roUz"][:] = 0.0; v["roe"][:] = ro1*e1
    for i in range(len(names)):
        k = f"roY{i}"
        if k in v: v[k][:] = ro1*Y[i]
        else: v.create_dataset(k, data=(ro1*Y[i]).astype(v["ro"].dtype))
# --- configs ---
yfloats = ", ".join(f"Y{i}: {Y[i]:.8f}" for i in range(len(names)))
(d/"bcondConfig.yaml").write_text(f"""inlet:  {{physID: 1, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {{{yfloats}, Pt: 1139000.0, Tt: {Tc:.2f}, k: 1.0, omega: 18000.0}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 1, ints: , floats: {{Ps: 5000.0, Pt: 5000.0, Tt: 300.0}}}}
wall:   {{physID: 3, kind: slip,             outputHDFflg: 1, ints: , floats: }}
axis:   {{physID: 4, kind: axis,             outputHDFflg: 0, ints: , floats: }}
""")
chem = f"chemistry: {{enabled: {a.chem}, mechanismFile: \"mech13.yaml\", jacobianMode: {a.jac}}}"
(d/"solverConfig.yaml").write_text(f"""# case/46 有限速度化学ノズル: 燃焼室 CEA 平衡 (A3000_eq, T_c {Tc:.0f} K, p_c 11.39 bar) 入口, va3 ノズル (case/44 run_0091 メッシュ)
# chem={a.chem} jacobianMode={a.jac} speciesImplicitCoupling={a.coupling}
mesh: {{meshFormat: "hdf5", discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 2, viscMethod: 0,
           ro: 1.2, visc: 0.0, thermCond: 0.0, cp: 1300.0, gamma: 1.3,
           species: [{", ".join(names)}], speciesDBFile: "species_db.yaml", thermoHrefTemp: 298.15,
           {chem}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {a.nstep}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: 2.0, cfl_pseudo: {a.cfl},
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1, speciesImplicitCoupling: {a.coupling}}}
  outStepStart: 0
  outStepInterval: 4000
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
""")
shutil.copy(HERE/"species_db.yaml", d/"species_db.yaml"); shutil.copy(HERE/"mech13.yaml", d/"mech13.yaml")
shutil.copy(SRC/"probe.yaml", d/"probe.yaml") if (SRC/"probe.yaml").exists() else None
print(f"{d}: Tc={Tc:.1f} X_OH={Xm['OH']/xs:.4f} X_NO={Xm['NO']/xs:.4f} Y=" + " ".join(f"{n}:{y:.4f}" for n, y in zip(names, Y) if y > 1e-4), f"N={N}")
