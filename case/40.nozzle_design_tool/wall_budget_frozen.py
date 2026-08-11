#!/usr/bin/env python3
"""§11 凍結場収支診断: run_0055 (Taw OFF, 収束) の場を凍結し、壁ノード半割 CV の
エネルギー収支 (対流 SLAU + 粘性仕事 + 伝導) を Python で再構成する。

熱伝導 4 方式:
  A: 現行 forge (T_state, over-relaxed, f 補間 k_eff)  = baseline
  B: Taw overlay + 旧 forge over-relaxed 式 (run_0053 の再現)
  C: Taw overlay + SU2 corrected-gradient 式 (算術平均 k_eff)
  D: C の laminar/turbulent 分離 (成分表示)
補助:
  A2: A と同式で k_eff だけ算術平均 (平均方式差の分離)
  C0: overlay 無しで SU2 corrected-gradient (式差だけの分離)

sign: 「ノード W へ入る向き」を正 (res[ic0]+=F_visc / res[ic1]-=F_visc に整合)。
"""
import h5py, numpy as np

CASE = '/home/sano/work/forge/case/40.nozzle_design_tool'
RUN = f'{CASE}/run_0055_node_yp30_dofonly_tawoff'
RES = f'{RUN}/res_12000.h5'
MESH = f'{RUN}/nozzle.h5'
gamma, cp_c, Prt, Pr_lam = 1.4, 1004.5, 0.9, 0.72
TARGETS = [6564, 6759, 11829]

# ---------------- mesh ----------------
fm = h5py.File(MESH, 'r')
nCells = int(fm['MESH'].attrs['nCells'])
nPlanes = int(fm['MESH'].attrs['nPlanes'])
nNormalPlanes = int(fm['MESH'].attrs['nNormalPlanes'])
coordm = fm['MESH/COORD'][:].reshape(-1, 3).astype(np.float64)
strct = fm['PLANES/STRUCT'][:]
surfVect = fm['PLANES/surfVect'][:].reshape(-1, 3).astype(np.float64)
pcent = fm['PLANES/centCoords'][:].reshape(-1, 3).astype(np.float64)
ccent = fm['CELLS/centCoords'][:].reshape(-1, 3).astype(np.float64)

plane_cells = np.full((nPlanes, 2), -1, dtype=np.int64)
ipp = 0
for ip in range(nPlanes):
    nn = strct[ipp]; ipp += 1 + nn
    nc = strct[ipp]; ipp += 1
    for k in range(nc):
        plane_cells[ip, k] = strct[ipp]; ipp += 1

# 軸対称 r 重み (axisymMethod 0): 面は面重心 r
r_face = np.maximum(np.abs(pcent[:, 1]), 1.0e-20)
S = surfVect * r_face[:, None]
ss = np.linalg.norm(S, axis=1)

# fx (calcStructualVariables_d 再現; r は比で相殺するので unweighted でも同じ)
ic0a = plane_cells[:, 0]; ic1a = plane_cells[:, 1]
nhat = surfVect / np.maximum(np.linalg.norm(surfVect, axis=1), 1e-30)[:, None]
d0p = np.abs(np.einsum('ij,ij->i', nhat, pcent - ccent[ic0a]))
d1p = np.abs(np.einsum('ij,ij->i', nhat, pcent - ccent[np.maximum(ic1a, 0)]))
fx = np.where(ic1a >= 0, d1p / np.maximum(d0p + d1p, 1e-30), 1.0)

# ---------------- field ----------------
ff = h5py.File(RES, 'r')
V = ff['VALUE']
g = lambda k: V[k][:].astype(np.float64)
ro, Ux, Uy, P, T, roe = g('ro'), g('Ux'), g('Uy'), g('P'), g('T'), g('roe')
sonic = g('sonic'); Ht = (roe + P) / ro
vis_lam, vis_turb, thermCond = g('vis_lam'), g('vis_turb'), g('thermCond')
vol = g('volume')
axisym_divU = g('axisym_divU')
Taw_diag = g('Taw_diag')
dUxdx, dUxdy = g('dUxdx'), g('dUxdy')
dUydx, dUydy = g('dUydx'), g('dUydy')
dPdx, dPdy = g('dPdx'), g('dPdy')
drodx, drody = g('drodx'), g('drody')
lim_ro, lim_Ux, lim_Uy, lim_P = g('limiter_ro'), g('limiter_Ux'), g('limiter_Uy'), g('limiter_P')

# ---------------- GG 勾配再現 (owner-state-only) と検証 ----------------
def green_gauss(phi):
    gx = np.zeros(nCells); gy = np.zeros(nCells)
    i0 = ic0a[:nNormalPlanes]; i1 = ic1a[:nNormalPlanes]
    f = fx[:nNormalPlanes]
    phif = f * phi[i0] + (1 - f) * phi[i1]
    np.add.at(gx, i0, phif * S[:nNormalPlanes, 0])
    np.add.at(gy, i0, phif * S[:nNormalPlanes, 1])
    np.add.at(gx, i1, -phif * S[:nNormalPlanes, 0])
    np.add.at(gy, i1, -phif * S[:nNormalPlanes, 1])
    ib0 = ic0a[nNormalPlanes:]
    np.add.at(gx, ib0, phi[ib0] * S[nNormalPlanes:, 0])
    np.add.at(gy, ib0, phi[ib0] * S[nNormalPlanes:, 1])
    return gx / vol, gy / vol

gPx, gPy = green_gauss(P)
scale = np.maximum(np.abs(dPdx), np.abs(dPdx).mean())
errx = np.abs(gPx - dPdx) / scale
print(f'[GG検証 vs 保存 dPdx] max rel err = {errx.max():.3e}  (99%ile {np.percentile(errx,99):.3e})')
dTdx, dTdy = green_gauss(T)

# ---------------- 壁ノード・壁関数量 ----------------
bw = fm['BCONDS/3']
wall_nodes = bw['iCells'][:].astype(np.int64)
wall_planes = bw['iPlanes'][:].astype(np.int64)
fw = h5py.File(RES.replace('res_', 'res_wall_3_'), 'r')['VALUE']
utau_b = fw['utau'][:].astype(np.float64)
ro_b = fw['ro'][:].astype(np.float64)
ypls_b = fw['ypls'][:].astype(np.float64)
Tau_Wall = np.full(nCells, -1.0)
Tau_Wall[wall_nodes] = ro_b * utau_b**2          # ransWallFunction: rho*utau^2 (淀みは 0)
Tau_Wall[wall_nodes[utau_b <= 0]] = 0.0
is_wall = np.zeros(nCells, bool); is_wall[wall_nodes] = True
T_ov = np.where(is_wall & (Taw_diag > 0), Taw_diag, T)   # primitive overlay

# ---------------- 粘性流束 (kernel 再現) ----------------
def visc_edge(ip):
    i0, i1 = plane_cells[ip]
    f = fx[ip]
    sxx, syy, sssv = S[ip, 0], S[ip, 1], ss[ip]
    dcc_v = ccent[i1] - ccent[i0]
    dcc = np.linalg.norm(dcc_v)
    D = max(abs(dcc_v[0] * sxx + dcc_v[1] * syy), 1e-30)
    delta = dcc * sssv**2 / D
    kx = sxx - dcc_v[0] * sssv**2 / D
    ky = syy - dcc_v[1] * sssv**2 / D
    w = lambda a: f * a[i0] + (1 - f) * a[i1]
    Uxf, Uyf = w(Ux), w(Uy)
    dUxdxf, dUxdyf = w(dUxdx), w(dUxdy)
    dUydxf, dUydyf = w(dUydx), w(dUydy)
    dTdxf, dTdyf = w(dTdx), w(dTdy)
    divu = w(axisym_divU)
    mu = w(vis_lam) + w(vis_turb)
    tau_x = mu * ((Ux[i1] - Ux[i0]) / dcc) * delta + mu * (dUxdxf * kx + dUxdyf * ky) \
        + mu * (dUxdxf * sxx + dUydxf * syy) - mu * 2 / 3 * divu * sxx
    tau_y = mu * ((Uy[i1] - Uy[i0]) / dcc) * delta + mu * (dUydxf * kx + dUydyf * ky) \
        + mu * (dUxdyf * sxx + dUydyf * syy) - mu * 2 / 3 * divu * syy
    # AddTauWall (xor)
    tw0, tw1 = Tau_Wall[i0], Tau_Wall[i1]
    rescaled = False
    if (tw0 > 0) != (tw1 > 0):
        tauw = tw0 if tw0 > 0 else tw1
        nx, ny = sxx / sssv, syy / sssv
        Tn = tau_x * nx + tau_y * ny
        ttx, tty = tau_x - Tn * nx, tau_y - Tn * ny
        tmag = np.hypot(ttx, tty)
        sc = tauw * sssv / max(tmag, 1e-30)
        tau_x *= sc; tau_y *= sc; rescaled = True
    work = tau_x * Uxf + tau_y * Uyf
    # 熱伝導 各方式
    tcf = w(thermCond) + cp_c * w(vis_turb) / Prt                 # f 補間 k_eff (forge)
    tca = 0.5 * (thermCond[i0] + thermCond[i1]) \
        + cp_c * 0.5 * (vis_turb[i0] + vis_turb[i1]) / Prt         # 算術 k_eff (SU2)
    tca_lam = 0.5 * (thermCond[i0] + thermCond[i1])
    tca_trb = cp_c * 0.5 * (vis_turb[i0] + vis_turb[i1]) / Prt
    grad_term = dTdxf * kx + dTdyf * ky
    qA = tcf * ((T[i1] - T[i0]) / dcc) * delta + tcf * grad_term
    qA2 = tca * ((T[i1] - T[i0]) / dcc) * delta + tca * grad_term
    qB = tcf * ((T_ov[i1] - T_ov[i0]) / dcc) * delta + tcf * grad_term
    # SU2 corrected-gradient
    gx = 0.5 * (dTdx[i0] + dTdx[i1]); gy = 0.5 * (dTdy[i0] + dTdy[i1])
    dd = dcc_v[0]**2 + dcc_v[1]**2
    def su2q(Tend, tc):
        dT = Tend[i1] - Tend[i0]
        corr = (dT - (gx * dcc_v[0] + gy * dcc_v[1])) / dd
        gcx = gx + corr * dcc_v[0]; gcy = gy + corr * dcc_v[1]
        return tc * (gcx * sxx + gcy * syy)
    qC = su2q(T_ov, tca)
    qC0 = su2q(T, tca)
    qD_lam = su2q(T_ov, tca_lam); qD_trb = su2q(T_ov, tca_trb)
    return dict(tau_x=tau_x, tau_y=tau_y, work=work, rescaled=rescaled,
                qA=qA, qA2=qA2, qB=qB, qC=qC, qC0=qC0, qD_lam=qD_lam, qD_trb=qD_trb)

# ---------------- SLAU 対流 (エネルギー行のみ) ----------------
def recon(phi, dx, dy, lim, ic, dcp):
    return phi[ic] + lim[ic] * (dx[ic] * dcp[0] + dy[ic] * dcp[1])

def slau_Fe(ip):
    icL, icR = plane_cells[ip]
    sxx, syy, sssv = S[ip, 0], S[ip, 1], ss[ip]
    pc = pcent[ip]
    dcpL = pc[:2] - ccent[icL, :2]; dcpR = pc[:2] - ccent[icR, :2]
    ro_L = recon(ro, drodx, drody, lim_ro, icL, dcpL)
    Ux_L = recon(Ux, dUxdx, dUxdy, lim_Ux, icL, dcpL)
    Uy_L = recon(Uy, dUydx, dUydy, lim_Uy, icL, dcpL)
    P_L = recon(P, dPdx, dPdy, lim_P, icL, dcpL)
    ro_R = recon(ro, drodx, drody, lim_ro, icR, dcpR)
    Ux_R = recon(Ux, dUxdx, dUxdy, lim_Ux, icR, dcpR)
    Uy_R = recon(Uy, dUydx, dUydy, lim_Uy, icR, dcpR)
    P_R = recon(P, dPdx, dPdy, lim_P, icR, dcpR)
    v2L = Ux_L**2 + Uy_L**2; v2R = Ux_R**2 + Uy_R**2
    h_p = gamma * P_L / ((gamma - 1) * ro_L) + 0.5 * v2L
    h_m = gamma * P_R / ((gamma - 1) * ro_R) + 0.5 * v2R
    Vn_p = (Ux_L * sxx + Uy_L * syy) / sssv
    Vn_m = (Ux_R * sxx + Uy_R * syy) / sssv
    c_hat = 0.5 * (sonic[icL] + sonic[icR])
    M_p, M_m = Vn_p / c_hat, Vn_m / c_hat
    beta_p = 0.5 * (M_p + abs(M_p)) / M_p if abs(M_p) >= 1 else 0.25 * (M_p + 1)**2 * (2 - M_p)
    beta_m = 0.5 * (M_m - abs(M_m)) / M_m if abs(M_m) >= 1 else 0.25 * (M_m - 1)**2 * (2 + M_m)
    gg = -max(min(M_p, 0.0), -1.0) * min(max(M_m, 0.0), 1.0)
    Vn_hat_abs = (ro_L * abs(Vn_p) + ro_R * abs(Vn_m)) / (ro_L + ro_R)
    Vn_hat_p = (1 - gg) * Vn_hat_abs + gg * abs(Vn_p)
    Vn_hat_m = (1 - gg) * Vn_hat_abs + gg * abs(Vn_m)
    chi = (1 - min(1.0, np.sqrt(0.5 * (v2L + v2R)) / c_hat))**2
    mdot = sssv * 0.5 * ((ro_L * (Vn_p + Vn_hat_p) + ro_R * (Vn_m - Vn_hat_m)) - chi / c_hat * (P_R - P_L))
    return 0.5 * (mdot + abs(mdot)) * h_p + 0.5 * (mdot - abs(mdot)) * h_m

# ---------------- 収支集計 ----------------
inc = {int(w): [] for w in wall_nodes}
for ip in range(nNormalPlanes):
    i0, i1 = plane_cells[ip]
    if is_wall[i0] or is_wall[i1]:
        if is_wall[i0]: inc[int(i0)].append((ip, +1.0))
        if is_wall[i1]: inc[int(i1)].append((ip, -1.0))

edge_cache = {}
rows = []
for w in wall_nodes:
    w = int(w)
    acc = dict(qA=0., qA2=0., qB=0., qC=0., qC0=0., qD_lam=0., qD_trb=0., work=0., conv=0.)
    for ip, sgn in inc[w]:
        if ip not in edge_cache:
            edge_cache[ip] = (visc_edge(ip), slau_Fe(ip))
        e, Fe = edge_cache[ip]
        for k in ('qA', 'qA2', 'qB', 'qC', 'qC0', 'qD_lam', 'qD_trb', 'work'):
            acc[k] += sgn * e[k]
        acc['conv'] += -sgn * Fe    # conv: res[ic0]-=F
    rows.append((w, acc))

arr = {k: np.array([r[1][k] for r in rows]) for k in rows[0][1]}
ids = np.array([r[0] for r in rows])
resA = arr['conv'] + arr['work'] + arr['qA']
volw = vol[ids]; row = ro[ids]
# 温度変化率スケール [K/s]: res_roe/(ro*cp*V)
dTdtA = resA / (row * cp_c * volw)
def dTdt(dq): return (resA + dq - arr['qA']) / (row * cp_c * volw) - dTdtA  # Δ(dT/dt)

print(f'\n=== 壁ノード {len(ids)} 個の凍結場エネルギー収支 (W 単位, ノードへ入る向き正) ===')
print(f'{"":14s} {"min":>12s} {"mean":>12s} {"max":>12s}')
for k in ('conv', 'work', 'qA', 'qA2', 'qC0', 'qB', 'qC', 'qD_lam', 'qD_trb'):
    a = arr[k]
    print(f'{k:14s} {a.min():12.4e} {a.mean():12.4e} {a.max():12.4e}')
print(f'{"res_roe(A)":14s} {resA.min():12.4e} {resA.mean():12.4e} {resA.max():12.4e}')

for tag, key in (('B-A (旧注入)', 'qB'), ('C-A (SU2式)', 'qC'), ('A2-A (平均差)', 'qA2'), ('C0-A (式差のみ)', 'qC0')):
    d = arr[key] - arr['qA']
    ddt = d / (row * cp_c * volw)
    print(f'\nΔ{tag}: flux [W] min {d.min():.4e} mean {d.mean():.4e} max {d.max():.4e}')
    print(f'   → dT/dt [K/s] min {ddt.min():.4e} mean {ddt.mean():.4e} max {ddt.max():.4e}')

print('\n=== 対象ノード詳細 ===')
w2b = {int(w): i for i, w in enumerate(wall_nodes)}
for w in TARGETS:
    i = np.where(ids == w)[0][0]
    b = w2b[w]
    print(f'\n--- node {w}  x={ccent[w,0]:.5f} y={ccent[w,1]:.5f}  V={vol[w]:.4e}')
    print(f'    T[W]={T[w]:.1f}K Taw={Taw_diag[w]:.1f}K  ro={ro[w]:.4f}  mu_t[W]={vis_turb[w]:.3e} '
          f'mu_lam[W]={vis_lam[w]:.3e}  u_tau={utau_b[b]:.2f} y+={ypls_b[b]:.1f} TauW={Tau_Wall[w]:.1f}Pa')
    print(f'    conv={arr["conv"][i]:+.4e} work={arr["work"][i]:+.4e} qA={arr["qA"][i]:+.4e} '
          f'-> res_roe(A)={resA[i]:+.4e} W  (dT/dt={dTdtA[i]:+.3e} K/s)')
    print(f'    qB={arr["qB"][i]:+.4e} qC={arr["qC"][i]:+.4e} qC0={arr["qC0"][i]:+.4e} '
          f'qA2={arr["qA2"][i]:+.4e}  [D: lam={arr["qD_lam"][i]:+.4e} turb={arr["qD_trb"][i]:+.4e}]')
    dB = arr['qB'][i] - arr['qA'][i]; dC = arr['qC'][i] - arr['qA'][i]
    print(f'    ΔB={dB:+.4e} W ({dB/(row[i]*cp_c*volw[i]):+.3e} K/s)   '
          f'ΔC={dC:+.4e} W ({dC/(row[i]*cp_c*volw[i]):+.3e} K/s)')
    for ip, sgn in inc[w]:
        e, Fe = edge_cache[ip]
        o = plane_cells[ip]
        other = int(o[1] if o[0] == w else o[0])
        print(f'      edge ip={ip} -> node {other} ({"wall" if is_wall[other] else "int "}) '
              f'sgn={sgn:+.0f} rescale={e["rescaled"]} | qA={sgn*e["qA"]:+.3e} qB={sgn*e["qB"]:+.3e} '
              f'qC={sgn*e["qC"]:+.3e} work={sgn*e["work"]:+.3e} conv={-sgn*Fe:+.3e} '
              f'| T_other={T[other]:.1f} mu_t_other={vis_turb[other]:.3e}')
