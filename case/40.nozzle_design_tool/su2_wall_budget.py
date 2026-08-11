#!/usr/bin/env python3
"""SU2 run_0042 (flow_fix2) の収束場で、壁ノード半割 CV の粘性エネルギー収支を
SU2 の離散式 (corrected-gradient + AddTauWall + 算術平均) で再構成する。

T_flux = vtu Temperature (壁= Taw 上書き済み primitive)
T_state = 保存量復元 (ρE から) — 勾配はこちらから (SU2 は勾配→上書きの順)
"""
import h5py, numpy as np
from scipy.spatial import cKDTree

sp = '/tmp/claude-1000/-home-sano-work-forge/0281401b-ee96-4551-8093-4eeb884a9b42/scratchpad'
CASE = '/home/sano/work/forge/case/40.nozzle_design_tool'
MESH = f'{CASE}/run_0055_node_yp30_dofonly_tawoff/nozzle.h5'
R = 287.058; gamma = 1.4; cv = R / (gamma - 1); cp_c = gamma * cv
Prt = 0.9; Pr_lam = 0.72
TARGETS = [6564, 6759, 11829]

# ---- forge mesh (トポロジ・面ベクトル共通) ----
fm = h5py.File(MESH, 'r')
nCells = int(fm['MESH'].attrs['nCells'])
nPlanes = int(fm['MESH'].attrs['nPlanes'])
nNormalPlanes = int(fm['MESH'].attrs['nNormalPlanes'])
coordm = fm['MESH/COORD'][:].reshape(-1, 3).astype(np.float64)
strct = fm['PLANES/STRUCT'][:]
surfVect = fm['PLANES/surfVect'][:].reshape(-1, 3).astype(np.float64)
pcent = fm['PLANES/centCoords'][:].reshape(-1, 3).astype(np.float64)
ccent = coordm.copy()   # SU2 vertex-centered: 端点座標=ノード座標

plane_cells = np.full((nPlanes, 2), -1, dtype=np.int64)
ipp = 0
for ip in range(nPlanes):
    nn = strct[ipp]; ipp += 1 + nn
    nc = strct[ipp]; ipp += 1
    for k in range(nc):
        plane_cells[ip, k] = strct[ipp]; ipp += 1
ic0a = plane_cells[:, 0]; ic1a = plane_cells[:, 1]

r_face = np.maximum(np.abs(pcent[:, 1]), 1.0e-20)
S = surfVect * r_face[:, None]
ss = np.linalg.norm(S, axis=1)

# ---- SU2 field → forge ノード順 ----
d = np.load(f'{sp}/su2_flow_fix2.npz')
tree = cKDTree(d['points'][:, :2])
dist, s2f = tree.query(coordm[:, :2])
print(f'[map] max node match dist = {dist.max():.2e}')
rho = d['Density'][s2f].astype(np.float64)
vel = d['Velocity'][s2f, :2].astype(np.float64)
Ux, Uy = vel[:, 0], vel[:, 1]
T_flux = d['Temperature'][s2f].astype(np.float64)   # 壁= Taw overlay 済み
rhoE = d['Energy'][s2f].astype(np.float64)
kT = d['Turb_Kin_Energy'][s2f].astype(np.float64)
T_state = (rhoE / rho - 0.5 * (Ux**2 + Uy**2) - kT) / cv
mu_lam = d['Laminar_Viscosity'][s2f].astype(np.float64)
mu_t = d['Eddy_Viscosity'][s2f].astype(np.float64)
yplus = d['Y_Plus'][s2f].astype(np.float64)
tc_lam = cp_c * mu_lam / Pr_lam

# ---- 壁ノード & Tau_Wall ----
bw = fm['BCONDS/3']
wall_nodes = bw['iCells'][:].astype(np.int64)
is_wall = np.zeros(nCells, bool); is_wall[wall_nodes] = True
# Tau_Wall = rho_w u_tau^2, u_tau from Y_Plus: u_tau = Y+ mu_w/(rho_w y_N)
# y_N = 第一内点距離: 各壁ノードの内部辺の相手ノードまでの壁法線距離で近似
Tau_Wall = np.full(nCells, -1.0)
# 壁法線 = 境界半割面ベクトル (壁ノード所有) 正規化
wall_planes = bw['iPlanes'][:].astype(np.int64)
nrm = surfVect[wall_planes] / np.linalg.norm(surfVect[wall_planes], axis=1)[:, None]
for w, pl, nh in zip(wall_nodes, wall_planes, nrm):
    # 内部辺の相手で non-wall のもの
    cand = []
    for ip in np.where(((ic0a == w) | (ic1a == w)) & (np.arange(nPlanes) < nNormalPlanes))[0]:
        o = int(ic1a[ip] if ic0a[ip] == w else ic0a[ip])
        if not is_wall[o]:
            dn = abs(np.dot(coordm[o, :2] - coordm[w, :2], nh[:2]))
            cand.append(dn)
    yN = min(cand) if cand else np.nan
    utau = yplus[w] * mu_lam[w] / (rho[w] * yN)
    Tau_Wall[w] = rho[w] * utau**2

# ---- GG 勾配 (T_state, 速度) — SU2 は GG/WLS いずれか。GG 近似で評価 ----
vol = np.zeros(nCells)
# forge の volume は VALUE にあるが SU2 場の勾配用に幾何から: forge volume 使用 (同一双対)
ffv = h5py.File(f'{CASE}/run_0055_node_yp30_dofonly_tawoff/res_12000.h5', 'r')['VALUE']
vol = ffv['volume'][:].astype(np.float64)

def green_gauss(phi):
    gx = np.zeros(nCells); gy = np.zeros(nCells)
    i0 = ic0a[:nNormalPlanes]; i1 = ic1a[:nNormalPlanes]
    phif = 0.5 * (phi[i0] + phi[i1])
    np.add.at(gx, i0, phif * S[:nNormalPlanes, 0]); np.add.at(gy, i0, phif * S[:nNormalPlanes, 1])
    np.add.at(gx, i1, -phif * S[:nNormalPlanes, 0]); np.add.at(gy, i1, -phif * S[:nNormalPlanes, 1])
    ib0 = ic0a[nNormalPlanes:]
    np.add.at(gx, ib0, phi[ib0] * S[nNormalPlanes:, 0]); np.add.at(gy, ib0, phi[ib0] * S[nNormalPlanes:, 1])
    return gx / vol, gy / vol

dTdx, dTdy = green_gauss(T_state)
dUxdx, dUxdy = green_gauss(Ux)
dUydx, dUydy = green_gauss(Uy)

# ---- 辺収支 (SU2 式) ----
def su2_edge(ip):
    i0, i1 = plane_cells[ip]
    sxx, syy, sssv = S[ip, 0], S[ip, 1], ss[ip]
    dvec = ccent[i1, :2] - ccent[i0, :2]
    dd = dvec @ dvec
    Um = 0.5 * np.array([Ux[i0] + Ux[i1], Uy[i0] + Uy[i1]])
    mu = 0.5 * (mu_lam[i0] + mu_lam[i1]) + 0.5 * (mu_t[i0] + mu_t[i1])
    # 平均勾配 + corrected (速度も SU2 は correct する)
    def corr(gx0, gy0, phi):
        gx = 0.5 * (gx0[i0] + gx0[i1]); gy = 0.5 * (gy0[i0] + gy0[i1])
        c = ((phi[i1] - phi[i0]) - (gx * dvec[0] + gy * dvec[1])) / dd
        return gx + c * dvec[0], gy + c * dvec[1]
    gux = corr(dUxdx, dUxdy, Ux); guy = corr(dUydx, dUydy, Uy)
    divu = gux[0] + guy[1]   # SU2 axisym は +v/r 項も持つが tau への影響は小さいので近似
    tau_xx = mu * (2 * gux[0]) - 2 / 3 * mu * divu
    tau_yy = mu * (2 * guy[1]) - 2 / 3 * mu * divu
    tau_xy = mu * (gux[1] + guy[0])
    tx = tau_xx * sxx + tau_xy * syy
    ty = tau_xy * sxx + tau_yy * syy
    # AddTauWall (xor)
    tw0, tw1 = Tau_Wall[i0], Tau_Wall[i1]
    if (tw0 > 0) != (tw1 > 0):
        tauw = max(tw0, 0) + max(tw1, 0)
        nx, ny = sxx / sssv, syy / sssv
        Tn = tx * nx + ty * ny
        ttx, tty = tx - Tn * nx, ty - Tn * ny
        sc = tauw * sssv / max(np.hypot(ttx, tty), 1e-30)
        tx *= sc; ty *= sc
    work = tx * Um[0] + ty * Um[1]
    # 熱伝導: corrected-gradient, T_flux (overlay 済み), 算術 k_eff
    keff = 0.5 * (tc_lam[i0] + tc_lam[i1]) + cp_c * 0.5 * (mu_t[i0] + mu_t[i1]) / Prt
    gTx, gTy = corr(dTdx, dTdy, T_flux)
    q = keff * (gTx * sxx + gTy * syy)
    # 比較: T_state 端点 (overlay 無し)
    gTx0, gTy0 = corr(dTdx, dTdy, T_state)
    q0 = keff * (gTx0 * sxx + gTy0 * syy)
    return tx, ty, work, q, q0

rows = []
for w in wall_nodes:
    w = int(w)
    acc = dict(work=0., q=0., q0=0.)
    edges = []
    for ip in np.where(((ic0a == w) | (ic1a == w)) & (np.arange(nPlanes) < nNormalPlanes))[0]:
        sgn = +1.0 if ic0a[ip] == w else -1.0
        tx, ty, work, q, q0 = su2_edge(ip)
        acc['work'] += sgn * work; acc['q'] += sgn * q; acc['q0'] += sgn * q0
        o = int(ic1a[ip] if ic0a[ip] == w else ic0a[ip])
        edges.append((ip, o, sgn, q, work))
    rows.append((w, acc, edges))

qs = np.array([r[1]['q'] for r in rows]); ws = np.array([r[1]['work'] for r in rows])
tot = qs + ws
print(f'\n=== SU2 場での壁ノード粘性エネルギー収支 (W, ノードへ入る向き正, n={len(rows)}) ===')
print(f'q_cond (Taw端点): min {qs.min():.3e} mean {qs.mean():.3e} max {qs.max():.3e}')
print(f'tau·u work      : min {ws.min():.3e} mean {ws.mean():.3e} max {ws.max():.3e}')
print(f'net (q+work)    : min {tot.min():.3e} mean {tot.mean():.3e} max {tot.max():.3e}')

for w in TARGETS:
    r = next(r for r in rows if r[0] == w)
    a = r[1]
    print(f'\n--- node {w}: T_state={T_state[w]:.1f} T_flux(Taw)={T_flux[w]:.1f} mu_t={mu_t[w]:.3e} '
          f'TauW={Tau_Wall[w]:.1f}Pa y+={yplus[w]:.1f}')
    print(f'    q(Taw)={a["q"]:+.4e}  q(state)={a["q0"]:+.4e}  work={a["work"]:+.4e}  '
          f'net={a["q"]+a["work"]:+.4e} W')
    for ip, o, sgn, q, work in r[2]:
        print(f'      ip={ip} -> {o} ({"wall" if is_wall[o] else "int "}): q={sgn*q:+.3e} work={sgn*work:+.3e} '
              f'T_state_o={T_state[o]:.1f} T_flux_o={T_flux[o]:.1f} mu_t_o={mu_t[o]:.3e}')
