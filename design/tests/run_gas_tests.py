#!/usr/bin/env python3
"""semi-perfect (NASA-9, frozen) ガスモデルの単体検証。

1. 一定 cp 擬似種で CPG に機械精度退化 (ν, A/A*, T(M))
2. NASA-9 の cp が文献値 (N2/O2/CO2/H2O) に一致
3. 混合擬似種 (mixture_pseudo_species) が多成分混合と一致
4. MOC カーネルの pm_nu/pm_mach がガスモデルで往復
5. 有効 γ 一定の CPG が燃焼ガスの A/A* を系統的に誤る (semi-perfect の存在意義)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forge_design.gas.semiperfect as spm  # noqa: E402
from forge_design.gas import GasCPG, GasSemiPerfect  # noqa: E402
from forge_design.geometry.moc_kernel import pm_mach, pm_mach_vec, pm_nu  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


# 1. CPG 退化
spm.SPECIES_NASA9["AIRX"] = dict(MW=0.0289647, low=[0, 0, 3.5, 0, 0, 0, 0, -1043.1, 3.0],
                                 high=[0, 0, 3.5, 0, 0, 0, 0, -1043.1, 3.0])
g = GasSemiPerfect({"AIRX": 1.0}, Tt=1000.0, n_tab=20000)
cpg = GasCPG(gamma=1.4, cp=float(g.cp_mass(500.0)[0]))
Ms = np.array([1.5, 2.0, 3.0, 4.0, 5.0])
check(f"CPG 退化: ν 一致 (max|Δ| {np.max(np.abs(g.nu(Ms)-cpg.nu(Ms))):.1e})",
      float(np.max(np.abs(g.nu(Ms) - cpg.nu(Ms)))) < 1e-6)
check(f"CPG 退化: A/A* 一致 (max rel {np.max(np.abs(g.area_ratio(Ms)/cpg.area_ratio(Ms)-1)):.1e})",
      float(np.max(np.abs(g.area_ratio(Ms) / cpg.area_ratio(Ms) - 1))) < 1e-6)
check(f"CPG 退化: γ* = {g.gamma_star:.6f}", abs(g.gamma_star - 1.4) < 1e-6)
check(f"CPG 退化: T* = {g.T_star:.3f} (理論 833.333)", abs(g.T_star - 1000 / 1.2) < 0.01)

# 2. NASA-9 文献値
for sp, T, ref, tol in (("N2", 300, 1040, 2), ("CO2", 1000, 1234, 2), ("H2O", 1000, 2290, 5), ("O2", 300, 918, 2)):
    v = float(GasSemiPerfect({sp: 1.0}, Tt=1500.0).cp_mass(T)[0])
    check(f"NASA-9 cp {sp}@{T}K = {v:.1f} (文献 ≈{ref})", abs(v - ref) < tol)

# 3. 混合擬似種
Y = {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.022, "N2": 0.7245}
mix = spm.mixture_pseudo_species(Y, "MIXT")
spm.SPECIES_NASA9["MIXT"] = dict(MW=mix["MIXT"]["MW"], low=mix["MIXT"]["nasa9_low"],
                                 high=mix["MIXT"]["nasa9_high"])
gm = GasSemiPerfect(Y, Tt=1000.0); g1 = GasSemiPerfect({"MIXT": 1.0}, Tt=1000.0)
Ts = np.array([300., 900., 1000., 1500.])
check(f"混合擬似種: cp 一致 (max|Δ| {np.max(np.abs(gm.cp_mass(Ts)-g1.cp_mass(Ts))):.1e})",
      float(np.max(np.abs(gm.cp_mass(Ts) - g1.cp_mass(Ts)))) < 1e-8)
check(f"混合擬似種: h 一致 (max|Δ| {np.max(np.abs(gm.h_mass(Ts)-g1.h_mass(Ts))):.1e})",
      float(np.max(np.abs(gm.h_mass(Ts) - g1.h_mass(Ts)))) < 1e-4)
check(f"混合擬似種: R 一致 ({gm.R:.4f} vs {g1.R:.4f})", abs(gm.R - g1.R) < 1e-9)

# 4. MOC カーネル往復
for M in (1.2, 2.5, 4.0):
    nu = float(pm_nu(M, gm)); Mb = pm_mach(nu, gm)
    check(f"pm_nu/pm_mach ガス往復 M={M}: {Mb:.6f}", abs(Mb - M) < 1e-4)
Mv = pm_mach_vec(pm_nu(Ms, gm), gm)
check(f"pm_mach_vec ガス往復 (max|Δ| {np.max(np.abs(Mv-Ms)):.1e})", float(np.max(np.abs(Mv - Ms))) < 1e-4)

# 5. 有効 γ 一定 CPG の系統誤差
ar_sp = float(gm.area_ratio(4.0))
ar_star = float(GasCPG(gamma=gm.gamma_star).area_ratio(4.0))
ar_m4 = float(GasCPG(gamma=float(gm.gamma(gm.T_of_M(4.0))[0])).area_ratio(4.0))
print(f"info 燃焼ガス A/A*(M4): semi-perfect {ar_sp:.3f} / CPG(γ*={gm.gamma_star:.3f}) {ar_star:.3f} / CPG(γ_M4) {ar_m4:.3f}")
check("semi-perfect の A/A* は γ* と γ_M4 の CPG の間にある (単一 γ では表現不能)",
      min(ar_star, ar_m4) < ar_sp < max(ar_star, ar_m4))
check(f"CPG(γ*) は出口半径を >5% 誤る ({100*(np.sqrt(ar_star/ar_sp)-1):+.1f}%)",
      abs(np.sqrt(ar_star / ar_sp) - 1) > 0.05)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
