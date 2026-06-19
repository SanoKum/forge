#!/usr/bin/env python3
"""
各化学種の cp(T) と エンタルピー h(T) を温度軸でプロットする。
- cp は NASA-9 (CEA) を Tlo=200K でクランプ外挿、h は線形外挿 (solver thermo_d.cuh と同一規約)。
- エンタルピー基準オフセット (thermoHrefTemp) の効果を可視化:
    左下: 絶対 (生成込み) エンタルピー  h_abs(T)         … 種ごとに桁違い (H2O≈-13.4MJ/kg)
    右下: オフセット後 sensible エンタルピー h(T)-h(Tref)  … 全種 Tref で 0、桁が揃う
出力: thermo_href_compare.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RU = 8.314462618
TMID = 1000.0
TLO, THI = 200.0, 6000.0
TREF = 298.15

MW = {"N2": 0.0280134, "H2O": 0.0180153, "O2": 0.0319988}
LO = {
 "N2":  np.array([2.210371497e4,-3.818461820e2,6.082738360,-8.530914410e-3,1.384646189e-5,-9.625793620e-9,2.519705809e-12,7.108460860e2,-1.076003744e1]),
 "H2O": np.array([-3.947960830e4,5.755731020e2,9.317826530e-1,7.222712860e-3,-7.342557370e-6,4.955043490e-9,-1.336933246e-12,-3.303974310e4,1.724205775e1]),
 "O2":  np.array([-3.425563420e4,4.847000970e2,1.119010961e0,4.293889240e-3,-6.836300520e-7,-2.023372700e-9,1.039040018e-12,-3.391454870e3,1.849699470e1]),
}
HI = {
 "N2":  np.array([5.877124060e5,-2.239249073e3,6.066949220,-6.139685500e-4,1.491806679e-7,-1.923105485e-11,1.061954386e-15,1.283210415e4,-1.586640027e1]),
 "H2O": np.array([1.034972096e6,-2.412698562e3,4.646110780,2.291998307e-3,-6.836830480e-7,9.426468930e-11,-4.822380530e-15,-1.384286509e4,-7.978148510e0]),
 "O2":  np.array([-1.037939022e6,2.344830282e3,1.819732036e0,1.267847582e-3,-2.188067988e-7,2.053719572e-11,-8.193467050e-16,-1.689010929e4,1.738716506e1]),
}

def coeffs(sp, T):
    return np.where(T[...,None] < TMID, LO[sp], HI[sp])

def cpR(sp, T):
    a = coeffs(sp, T); Ti = 1.0/T
    return a[...,0]*Ti*Ti + a[...,1]*Ti + a[...,2] + a[...,3]*T + a[...,4]*T**2 + a[...,5]*T**3 + a[...,6]*T**4

def hRT(sp, T):
    a = coeffs(sp, T); Ti = 1.0/T; lnT = np.log(T)
    return (-a[...,0]*Ti*Ti + a[...,1]*lnT*Ti + a[...,2] + a[...,3]*T/2 + a[...,4]*T**2/3
            + a[...,5]*T**3/4 + a[...,6]*T**4/5 + a[...,7]*Ti)

def cp_mass(sp, T):                       # クランプ外挿 (solver と同一)
    R = RU/MW[sp]; Tc = np.clip(T, TLO, THI)
    return R*cpR(sp, Tc)

def h_mass(sp, T):                        # 線形外挿 (solver と同一)
    R = RU/MW[sp]; Tc = np.clip(T, TLO, THI)
    h = R*Tc*hRT(sp, Tc)
    lo = T < TLO
    if np.any(lo):
        cp_lo = R*cpR(sp, np.full_like(T, TLO))
        h = np.where(lo, R*TLO*hRT(sp, np.full_like(T, TLO)) + cp_lo*(T-TLO), h)
    return h

T = np.linspace(120.0, 1200.0, 600)
species = ["N2", "O2", "H2O"]
colors = {"N2":"tab:blue", "O2":"tab:green", "H2O":"tab:red"}

fig, ax = plt.subplots(1, 3, figsize=(16, 5))

# (1) cp(T)
for sp in species:
    ax[0].plot(T, cp_mass(sp, T), color=colors[sp], label=sp)
ax[0].axvline(TLO, ls=":", color="gray", lw=1)
ax[0].text(TLO+5, ax[0].get_ylim()[0], "Tlo=200K\n(cp clamp)", fontsize=8, color="gray", va="bottom")
ax[0].set_xlabel("T [K]"); ax[0].set_ylabel("cp [J/(kg·K)]")
ax[0].set_title("(1) Specific heat cp(T)  (clamped below 200K)")
ax[0].legend(); ax[0].grid(alpha=0.3)

# (2) 絶対エンタルピー h_abs(T) — 桁違い
for sp in species:
    ax[1].plot(T, h_mass(sp, T)/1e6, color=colors[sp], label=f"{sp}")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].axvline(TREF, ls="--", color="gray", lw=1)
ax[1].text(TREF+5, ax[1].get_ylim()[1]*0.9, f"Tref={TREF}K", fontsize=8, color="gray")
ax[1].set_xlabel("T [K]"); ax[1].set_ylabel("h_abs [MJ/kg]")
ax[1].set_title("(2) Absolute enthalpy (with formation)\nH2O sits at -13.4 MJ/kg (orders apart)")
ax[1].legend(); ax[1].grid(alpha=0.3)

# (3) オフセット後 sensible エンタルピー h(T)-h(Tref)
for sp in species:
    href = h_mass(sp, np.array([TREF]))[0]
    ax[2].plot(T, (h_mass(sp, T)-href)/1e6, color=colors[sp], label=f"{sp}  (h_ref={href/1e6:+.3f})")
ax[2].axhline(0, color="k", lw=0.6)
ax[2].axvline(TREF, ls="--", color="gray", lw=1)
ax[2].plot([TREF],[0],"ko",ms=4)
ax[2].set_xlabel("T [K]"); ax[2].set_ylabel("h - h(Tref) [MJ/kg]")
ax[2].set_title(f"(3) Offset (sensible) enthalpy h-h(Tref)\nthermoHrefTemp={TREF}: all species 0 at Tref")
ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)

plt.tight_layout()
plt.savefig("thermo_href_compare.png", dpi=130)
print("wrote thermo_href_compare.png")
# 参考値
for sp in species:
    print(f"  {sp}: h_abs(298.15)={h_mass(sp,np.array([TREF]))[0]/1e6:+.4f} MJ/kg, "
          f"cp(300)={cp_mass(sp,np.array([300.0]))[0]:.1f}, cp(150,clamp)={cp_mass(sp,np.array([150.0]))[0]:.1f} J/kgK")
