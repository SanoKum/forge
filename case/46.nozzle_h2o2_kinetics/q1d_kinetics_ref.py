#!/usr/bin/env python3
"""準1次元 (Q1D) 有限速度化学ノズル参照解 (Cantera, mech15.yaml)。
壁半径 r_w(x) はメッシュの壁境界ノードから、スロート状態は CEA 平衡 (cea_A3000_eq_stations.json の throat 列)。
スロート直下流 (A/A*=1.02, 凍結等エントロピーで M≈1.15 に置く) から出口まで
  d(ρuA)=0, ρu du=−dp, h+u²/2=const, dY_s/dx=ω_s/(ρu)
を積分する。--frozen で ω=0。出力 CSV: x, T, p, M, u, Y_*。
  .venv-chem/bin/python q1d_kinetics_ref.py run_dir/nozzle.h5 out.csv [--frozen]"""
import argparse, json, h5py, numpy as np, cantera as ct
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
ap = argparse.ArgumentParser(); ap.add_argument("mesh_h5"); ap.add_argument("out"); ap.add_argument("--frozen", action="store_true")
ap.add_argument("--mech", default="mech15.yaml"); ap.add_argument("--A0", type=float, default=1.02)
ap.add_argument("--start-from-forge", default=None, help="forge res_*.h5: 軸線上で M=--start-M となる点の (T,p,u,Y) を初期状態にする (NO など遅い化学種を CEA 平衡でなく forge の値から始める; 軸の音速線は幾何スロートより下流なので M で選ぶ)")
ap.add_argument("--start-M", type=float, default=1.3)
a = ap.parse_args()
# --- 壁形状 ---
with h5py.File(a.mesh_h5) as f:
    xyz = f["MESH/COORD"][:].reshape(-1, 3); wn = np.unique(f["BCONDS/3/iCells"][:])
xw, rw = xyz[wn, 0], xyz[wn, 1]; o = np.argsort(xw); xw, rw = xw[o], rw[o]
it = np.argmin(rw); x_t, r_t = xw[it], rw[it]; A_t = np.pi*r_t**2
A_of = lambda x: np.pi*np.interp(x, xw, rw)**2
dAdx = lambda x: (A_of(x+1e-4) - A_of(x-1e-4))/2e-4
# --- スロート状態 (CEA 平衡) ---
cea = json.load(open("cea_A3000_eq_stations.json")); Tt_, Pt_ = cea["T"][1], cea["P_bar"][1]*1e5
g = ct.Solution(a.mech); names = g.species_names
X = {}
for k, v in cea["X"].items():
    kk = {"Ar": "AR", "CO": "CO2", "N2O": None}.get(k, k)
    if kk is None or kk not in names: continue
    X[kk] = X.get(kk, 0.0) + v[1]
g.TPX = Tt_, Pt_, X
Y = g.Y.copy(); s_t = g.entropy_mass; h0 = g.enthalpy_mass + 0.5*g.sound_speed**2   # 全エンタルピー (スロート M=1)
rho_t = g.density; a_t = g.sound_speed; mdot = rho_t*a_t*A_t
# 開始点: A/A*=A0 で凍結等エントロピー (組成固定) の超音速解を求める
def state_at_area(Aratio):
    def resid(T):
        g.TPY = T, Pt_, Y; g.SP = s_t, None  # 仮
        return 0
    # 等エントロピー: s(T,p)=s_t, 質量流量 ρuA=mdot, h+u²/2=h0 → p を未知として解く
    def f(p):
        g.SPY = s_t, p, Y; u = np.sqrt(max(2*(h0 - g.enthalpy_mass), 0.0)); return g.density*u*Aratio*A_t - mdot
    p_lo, p_hi = Pt_*0.05, Pt_*0.999   # 超音速枝 (低圧側の根)
    ps = np.linspace(p_lo, p_hi, 400); fs = [f(p) for p in ps]
    roots = [brentq(f, ps[i], ps[i+1]) for i in range(len(ps)-1) if fs[i]*fs[i+1] < 0]
    p = min(roots); g.SPY = s_t, p, Y; u = np.sqrt(2*(h0 - g.enthalpy_mass)); return g.T, p, u
x0 = brentq(lambda x: A_of(x)/A_t - a.A0, x_t, xw[-1]); T0, p0, u0 = state_at_area(a.A0)
if a.start_from_forge:
    with h5py.File(a.start_from_forge) as f:
        v = f["VALUE"]; ax_nodes = np.where(np.abs(xyz[:, 1]) < 1e-9)[0]; ax_nodes = ax_nodes[np.argsort(xyz[ax_nodes, 0])]
        Mf = v["Ux"][:][ax_nodes] / np.maximum(v["sonic"][:][ax_nodes], 1e-30); xa = xyz[ax_nodes, 0]
        j = int(np.argmax((Mf >= a.start_M) & (xa > x_t))); i = ax_nodes[j]; x0 = float(xa[j])
        T0, p0, u0 = float(v["T"][i]), float(v["P"][i]), float(v["Ux"][i])
        yf = np.array([float(v[f"Y{k}"][i]) for k in range(len(names))])   # forge の species 順 = mech15 と同じ順にしておくこと
        Y = yf / yf.sum(); g.TPY = T0, p0, Y; mdot = g.density*u0*A_of(x0)
    print(f"start from forge {a.start_from_forge} at x={xyz[i,0]:.4f}: T={T0:.1f} p={p0:.0f} u={u0:.1f} Y_OH={Y[names.index('OH')]:.3e} Y_NO={Y[names.index('NO')]:.3e}")
print(f"throat x={x_t:.4f} r={r_t:.4f} T*={Tt_:.1f} p*={Pt_:.0f}; start x0={x0:.4f} (A/A*={a.A0}) T0={T0:.1f} p0={p0:.0f} u0={u0:.1f} M0={u0/np.sqrt(ct.gas_constant/g.mean_molecular_weight*g.cp/g.cv*T0):.3f}")
def rhs(x, q):
    u, T = q[0], q[1]; Ys = np.clip(q[2:], 0, None); Ys /= Ys.sum()
    A = A_of(x); rho = mdot/(u*A)
    g.TDY = T, rho, Ys
    p = g.P; R = ct.gas_constant/g.mean_molecular_weight; cp = g.cp_mass
    w = np.zeros_like(Ys) if a.frozen else g.net_production_rates*g.molecular_weights   # kg/m3/s
    dYdx = w/(rho*u)
    hs = g.partial_molar_enthalpies/g.molecular_weights          # J/kg
    Rs = ct.gas_constant/g.molecular_weights
    dRdx = (Rs*dYdx).sum()
    # 運動量: (u − RT/u) u' + R T' = (RT/A) A' − T R'
    # エネルギー: u u' + cp T' = −Σ h_s ω_s/(ρu)
    a11, a12, b1 = u - R*T/u, R, R*T/A*dAdx(x) - T*dRdx
    a21, a22, b2 = u, cp, -(hs*dYdx).sum()
    det = a11*a22 - a12*a21
    du = (b1*a22 - a12*b2)/det; dT = (a11*b2 - a21*b1)/det
    return np.concatenate([[du, dT], dYdx])
q0 = np.concatenate([[u0, T0], Y]); xs = np.linspace(x0, xw[-1], 600)
sol = solve_ivp(rhs, (x0, xw[-1]), q0, method="BDF", t_eval=xs, rtol=1e-8, atol=1e-14, max_step=0.02)
rows = []
for x, q in zip(sol.t, sol.y.T):
    u, T = q[0], q[1]; Ys = np.clip(q[2:], 0, None); Ys /= Ys.sum(); A = A_of(x); rho = mdot/(u*A); g.TDY = T, rho, Ys
    rows.append([x, T, g.P, u/g.sound_speed, u] + list(Ys))
np.savetxt(a.out, np.array(rows), delimiter=",", header="x,T,p,M,u," + ",".join("Y_"+n for n in names), comments="")
r = rows[-1]; print(f"exit x={r[0]:.3f}: T={r[1]:.1f} K p={r[2]:.0f} Pa M={r[3]:.3f} Y_OH={r[5+names.index('OH')]:.2e} Y_NO={r[5+names.index('NO')]:.2e} ({'frozen' if a.frozen else 'finite-rate'})")
