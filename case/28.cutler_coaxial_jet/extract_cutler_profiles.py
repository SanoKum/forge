#!/usr/bin/env python3
# Cutler 超音速同軸 He/空気ジェット (AIAA-99-3588) の定量照合用プロファイル抽出。
#   res_*.h5 (forge 出力) + 内蔵メッシュ (MESH/COORD,CONNE) から:
#     - 中心軸上の He モル分率 X_He(x) 減衰と He ポテンシャルコア長
#     - x/D ステーションでの径方向プロファイル (X_He, 軸速度 Ux, 全温 Tt)
#   を抽出し CSV + PNG に出力する。実験データ (cutler_exp_centerline.csv 等) が
#   あれば重ねる (列: x_over_D,X_He / r_over_D,X_He)。
#
# 使い方: python3 extract_cutler_profiles.py <res_*.h5> [out_prefix]
#   軸対称 2D: x=軸方向, y=半径方向 r。ジェット出口径 D=10.24mm (r_jet=5.12mm)。
import sys, numpy as np, h5py
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

RES = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else RES.replace('.h5','')
D_JET = 10.24e-3                  # 中心ジェット出口径 [m]
MW = {'He':0.0040026,'N2':0.0280134,'O2':0.0319988}
# 化学種順は solverConfig の species: と同じ。binary=[He,N2], ternary=[He,O2,N2]。
def species_list(V):
    return ['He','N2'] if 'Y2' not in V else ['He','O2','N2']

h = h5py.File(RES,'r'); V = h['VALUE']
sp = species_list(V)
Y = [V['Y%d'%i][()] for i in range(len(sp))]
# モル分率 X_He = (Y_He/MW_He)/Σ(Y_s/MW_s)
invM = sum(Y[i]/MW[sp[i]] for i in range(len(sp)))
X_He = (Y[0]/MW['He'])/invM
ro = V['ro'][()]; T = V['T'][()]
Ux = V['Ux'][()]; Uy = V['Uy'][()] if 'Uy' in V else np.zeros_like(Ux)
Umag = np.sqrt(Ux*Ux+Uy*Uy)
# 全温 Tt ≈ T + |u|²/(2 cp)。cp は混合 (簡易に T 依存無視の代表値) — ここでは観測量として T も出す。
cp_rep = {'He':5193.0,'N2':1040.0,'O2':918.0}
cp_mix = sum(Y[i]*cp_rep[sp[i]] for i in range(len(sp)))
Tt = T + Umag*Umag/(2.0*cp_mix)

# セル重心 (XDMF Mixed quad: CONNE = [5,n0,n1,n2,n3]×ncell)
co = h['MESH']['COORD'][()].reshape(-1,3)
CN = h['MESH']['CONNE'][()].astype(np.int64).reshape(-1,5)[:,1:5]
xc = co[CN,0].mean(1); rc = co[CN,1].mean(1)

# --- 中心軸プロファイル (各 x 位置で最も軸に近いセルを採用しバンド散らばりを除く) ---
axis = rc < 0.0030
xa_all = xc[axis]; rc_a = rc[axis]
XHe_all=X_He[axis]; U_all=Umag[axis]; Tt_all=Tt[axis]
# x をビン化し各ビンで最小 r のセルを選ぶ
nb = 400
edges = np.linspace(xa_all.min(), xa_all.max(), nb+1)
bi = np.clip(np.digitize(xa_all, edges)-1, 0, nb-1)
xa=[]; XHe_a=[]; U_a=[]; Tt_a=[]
for b in range(nb):
    m = bi==b
    if not m.any(): continue
    j = np.argmin(rc_a[m])
    xa.append(xa_all[m][j]); XHe_a.append(XHe_all[m][j]); U_a.append(U_all[m][j]); Tt_a.append(Tt_all[m][j])
xa=np.array(xa); XHe_a=np.array(XHe_a); U_a=np.array(U_a); Tt_a=np.array(Tt_a)
o=np.argsort(xa); xa,XHe_a,U_a,Tt_a = xa[o],XHe_a[o],U_a[o],Tt_a[o]
xD = xa/D_JET
# ポテンシャルコア長: 中心軸 X_He が出口値 (=ジェットコアのモル分率, binary 1.0 / ternary 0.95)
# の 99%/95% を下回る最初の x。出口値で正規化し binary/ternary 双方に対応。
X_jet = float(np.nanmax(XHe_a[xD < 2.0])) if (xD < 2.0).any() else float(XHe_a[0])
def core_len(frac):
    below = np.where(XHe_a < frac*X_jet)[0]
    return xa[below[0]]/D_JET if len(below) else np.nan
core99, core95 = core_len(0.99), core_len(0.95)

np.savetxt(OUT+'_centerline.csv',
           np.column_stack([xD, XHe_a, U_a, Tt_a]),
           header='x_over_D,X_He,Umag[m/s],Tt[K]', delimiter=',', comments='')

# --- 径方向プロファイル (x/D ステーション) ---
stations = [2,5,10,15,20]
rad = {}
for s in stations:
    xt = s*D_JET; m = np.abs(xc-xt) < 0.0025
    if m.sum() < 4: continue
    oo = np.argsort(rc[m]); rr = rc[m][oo]
    rad[s] = (rr/D_JET, X_He[m][oo], Umag[m][oo], Tt[m][oo])

# --- プロット ---
fig, ax = plt.subplots(2,2, figsize=(12,9))
ax[0,0].plot(xD, XHe_a, '-', lw=1.5)
ax[0,0].axhline(0.99*X_jet,ls=':',c='gray'); ax[0,0].set(xlabel='x/D', ylabel='centerline $X_{He}$', title=f'He core (X_jet={X_jet:.2f}): x/D(99%)={core99:.1f}, x/D(95%)={core95:.1f}')
ax[0,0].grid(alpha=0.3)
ax[0,1].plot(xD, U_a, '-', lw=1.5, color='C1'); ax[0,1].set(xlabel='x/D', ylabel='centerline |U| [m/s]'); ax[0,1].grid(alpha=0.3)
for s,(r,xh,u,tt) in rad.items():
    ax[1,0].plot(r, xh, label=f'x/D={s}')
    ax[1,1].plot(r, u, label=f'x/D={s}')
ax[1,0].set(xlabel='r/D', ylabel='$X_{He}$', title='radial He mole fraction'); ax[1,0].legend(fontsize=8); ax[1,0].grid(alpha=0.3); ax[1,0].set_xlim(0,3)
ax[1,1].set(xlabel='r/D', ylabel='|U| [m/s]', title='radial velocity'); ax[1,1].legend(fontsize=8); ax[1,1].grid(alpha=0.3); ax[1,1].set_xlim(0,3)

# 実験データ重ね (あれば)
import os
exp = 'cutler_exp_centerline.csv'
if os.path.exists(exp):
    E = np.genfromtxt(exp, delimiter=',', names=True)
    ax[0,0].plot(E['x_over_D'], E['X_He'], 'ks', ms=4, label='Cutler exp'); ax[0,0].legend(fontsize=8)

fig.suptitle(f'Cutler coaxial jet — {os.path.basename(RES)}  ({"/".join(sp)})')
fig.tight_layout(); fig.savefig(OUT+'_profiles.png', dpi=110)
print(f'species={sp}')
print(f'He potential core: x/D(X=0.99)={core99:.2f}, x/D(X=0.95)={core95:.2f}  (D={D_JET*1e3:.2f}mm)')
print(f'centerline X_He: x/D=0->{XHe_a[0]:.3f}, =5->{np.interp(5,xD,XHe_a):.3f}, =10->{np.interp(10,xD,XHe_a):.3f}, =20->{np.interp(20,xD,XHe_a):.3f}')
print(f'centerline |U|:  x/D=0->{U_a[0]:.0f}, =10->{np.interp(10,xD,U_a):.0f}, =20->{np.interp(20,xD,U_a):.0f} m/s')
print(f'wrote {OUT}_centerline.csv, {OUT}_profiles.png')
