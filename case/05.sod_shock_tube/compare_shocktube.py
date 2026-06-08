#!/usr/bin/env python3
# 衝撃管検証: forge の CPG/TP 計算結果を中心線プロファイルに抽出し、
#   (1) CPG vs 厳密 Sod (Riemann) 解  -> 枠組み検証
#   (2) TP  vs CPG                     -> cp(T) 効果の定量提示
# を 1 枚の図 (rho/P/T/u/Mach) にまとめる。
import sys, numpy as np, h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- 初期条件 (gen_shocktube_ic.py と一致させること) ---
GAM = 1.4; RU = 8.314462618; MW = 0.0280134; R = RU/MW   # N2: R=296.8
X0 = 0.5                       # 隔膜位置
T_FINAL = 0.00025              # 出力 res_1000 の時刻 [s]
PL, TL = 1.0e6, 2000.0
PR, TR = 1.0e5,  400.0
roL = PL/(R*TL); roR = PR/(R*TR)
uL = uR = 0.0

# ----------------------------------------------------------------------
# 厳密 Sod (Riemann) 解  -- Toro, "Riemann Solvers...", Ch.4 (CPG, gamma 一定)
# ----------------------------------------------------------------------
def exact_sod(x, t, g=GAM):
    aL = np.sqrt(g*PL/roL); aR = np.sqrt(g*PR/roR)
    AL = 2.0/((g+1)*roL); BL = (g-1)/(g+1)*PL
    AR = 2.0/((g+1)*roR); BR = (g-1)/(g+1)*PR

    def fK(p, pK, roK, aK, AK, BK):
        if p > pK:   # shock
            return (p-pK)*np.sqrt(AK/(p+BK))
        else:        # rarefaction
            return 2*aK/(g-1)*((p/pK)**((g-1)/(2*g)) - 1.0)

    def dfK(p, pK, roK, aK, AK, BK):
        if p > pK:
            return np.sqrt(AK/(BK+p))*(1.0 - 0.5*(p-pK)/(BK+p))
        else:
            return (1.0/(roK*aK))*(p/pK)**(-(g+1)/(2*g))

    def f(p):
        return (fK(p, PL, roL, aL, AL, BL) + fK(p, PR, roR, aR, AR, BR) + (uR-uL))

    def df(p):
        return dfK(p, PL, roL, aL, AL, BL) + dfK(p, PR, roR, aR, AR, BR)

    # Newton 反復で p_star
    p = 0.5*(PL+PR)
    for _ in range(100):
        dp = -f(p)/df(p)
        p += dp
        if abs(dp) < 1e-8*p:
            break
    p_star = p
    u_star = 0.5*(uL+uR) + 0.5*(fK(p_star, PR, roR, aR, AR, BR)
                               - fK(p_star, PL, roL, aL, AL, BL))

    # 各 x をサンプル (S = (x-X0)/t)
    ro = np.empty_like(x); u = np.empty_like(x); P = np.empty_like(x)
    for i, xi in enumerate(x):
        S = (xi - X0)/t
        if S <= u_star:   # 左側
            if p_star > PL:  # 左 shock
                SL = uL - aL*np.sqrt((g+1)/(2*g)*p_star/PL + (g-1)/(2*g))
                if S < SL:
                    ro[i], u[i], P[i] = roL, uL, PL
                else:
                    ro[i] = roL*((p_star/PL + (g-1)/(g+1))/((g-1)/(g+1)*p_star/PL + 1))
                    u[i], P[i] = u_star, p_star
            else:            # 左 rarefaction
                roLs = roL*(p_star/PL)**(1/g)
                aLs = aL*(p_star/PL)**((g-1)/(2*g))
                SHL = uL - aL; STL = u_star - aLs
                if S < SHL:
                    ro[i], u[i], P[i] = roL, uL, PL
                elif S > STL:
                    ro[i], u[i], P[i] = roLs, u_star, p_star
                else:        # 扇内
                    u[i] = 2/(g+1)*(aL + (g-1)/2*uL + S)
                    c = 2/(g+1)*(aL + (g-1)/2*(uL - S))
                    ro[i] = roL*(c/aL)**(2/(g-1))
                    P[i] = PL*(c/aL)**(2*g/(g-1))
        else:             # 右側
            if p_star > PR:  # 右 shock
                SR = uR + aR*np.sqrt((g+1)/(2*g)*p_star/PR + (g-1)/(2*g))
                if S > SR:
                    ro[i], u[i], P[i] = roR, uR, PR
                else:
                    ro[i] = roR*((p_star/PR + (g-1)/(g+1))/((g-1)/(g+1)*p_star/PR + 1))
                    u[i], P[i] = u_star, p_star
            else:            # 右 rarefaction
                roRs = roR*(p_star/PR)**(1/g)
                aRs = aR*(p_star/PR)**((g-1)/(2*g))
                SHR = uR + aR; STR = u_star + aRs
                if S > SHR:
                    ro[i], u[i], P[i] = roR, uR, PR
                elif S < STR:
                    ro[i], u[i], P[i] = roRs, u_star, p_star
                else:
                    u[i] = 2/(g+1)*(-aR + (g-1)/2*uR + S)
                    c = 2/(g+1)*(aR - (g-1)/2*(uR - S))
                    ro[i] = roR*(c/aR)**(2/(g-1))
                    P[i] = PR*(c/aR)**(2*g/(g-1))
    Temp = P/(ro*R)
    a = np.sqrt(g*P/ro)
    return dict(x=x, ro=ro, u=u, P=P, T=Temp, M=np.abs(u)/a)


def read_profile(h5path):
    with h5py.File(h5path, 'r') as f:
        co = f['MESH/COORD'][:].reshape(-1, 3)
        c = f['MESH/CONNE'][:].reshape(-1, 9)
        x = co[c[:, 1:9], 0].mean(1)
        ro = f['VALUE/ro'][:]; P = f['VALUE/P'][:]; T = f['VALUE/T'][:]
        Ux = f['VALUE/Ux'][:]; sonic = f['VALUE/sonic'][:]
    idx = np.argsort(x)
    M = np.abs(Ux)/np.where(sonic > 1e-9, sonic, 1e-9)
    return dict(x=x[idx], ro=ro[idx], P=P[idx], T=T[idx], u=Ux[idx], M=M[idx])


tp = read_profile(sys.argv[1] if len(sys.argv) > 1 else 'run_0001_slau_tp_n2/res_1000.h5')
cpg = read_profile(sys.argv[2] if len(sys.argv) > 2 else 'run_0002_slau_cpg/res_1000.h5')
ex = exact_sod(cpg['x'], T_FINAL)

fields = [('ro', r'$\rho$ [kg/m$^3$]'), ('P', 'P [Pa]'), ('T', 'T [K]'),
          ('u', 'u [m/s]'), ('M', 'Mach')]
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
axes = axes.ravel()
for ax, (k, lab) in zip(axes, fields):
    ax.plot(ex['x'], ex[k], 'k-', lw=1.5, label='exact Sod (CPG)')
    ax.plot(cpg['x'], cpg[k], 'C0o', ms=2.5, label='forge CPG (SLAU)')
    ax.plot(tp['x'], tp[k], 'C3.', ms=3, label='forge TP N2')
    ax.set_xlabel('x [m]'); ax.set_ylabel(lab); ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
axes[-1].axis('off')
fig.suptitle(f'N2 shock tube  (PL/PR=10, TL/TR=2000/400K)  t={T_FINAL*1e6:.0f} us\n'
             'CPG vs exact = framework check ; TP vs CPG = cp(T) effect', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('shocktube_comparison.png', dpi=130)
print('wrote shocktube_comparison.png')

# --- 定量誤差: CPG vs 厳密 (L2 相対) ---
def l2rel(a, b):
    return np.sqrt(np.mean((a-b)**2))/np.sqrt(np.mean(b**2))
print('--- CPG vs exact Sod (L2 relative error) ---')
for k in ('ro', 'P', 'u'):
    print(f'  {k:3s}: {l2rel(cpg[k], ex[k])*100:6.2f} %')
# --- TP vs CPG 差分 (波速・状態の cp(T) シフト) ---
print('--- TP vs CPG (max abs diff over domain) ---')
for k, u in (('T', 'K'), ('u', 'm/s'), ('P', 'Pa'), ('ro', 'kg/m3')):
    print(f'  {k:3s}: {np.max(np.abs(tp[k]-cpg[k])):10.3g} {u}')
print(f'TP  centerline: T[{tp["T"].min():.1f},{tp["T"].max():.1f}]K  '
      f'P[{tp["P"].min():.3e},{tp["P"].max():.3e}]  u_max={np.max(np.abs(tp["u"])):.1f}')
print(f'CPG centerline: T[{cpg["T"].min():.1f},{cpg["T"].max():.1f}]K  '
      f'P[{cpg["P"].min():.3e},{cpg["P"].max():.3e}]  u_max={np.max(np.abs(cpg["u"])):.1f}')
print(f'exact p_star/u_star region: P_contact={ex["P"][np.argmin(np.abs(ex["x"]-X0))]:.4e}')
