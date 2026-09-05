#!/usr/bin/env python3
"""Phase A (混合分率インフラ) の検証: res_*.h5 の xi / xiVar / chi を点検する。
  python3 check_mixfrac.py run_dir step
- 反応 OFF なら Bilger ξ は Y_H2 / Y_H2,fuel と厳密一致する (H は H2 と H2O にしか無く、H2O は coflow 由来で線形)。
- 0 ≤ ξ''² ≤ ξ̃(1-ξ̃) (実現可能性)、入口で ξ''²≈0、せん断層で最大。
- χ̃ = cChi β* ω ξ''² の桁 (Cabra せん断層で O(10–1000) 1/s が文献値)。
"""
import sys, h5py, numpy as np
run, step = sys.argv[1], sys.argv[2]
with h5py.File(f"{run}/cabra.h5") as m: xyz = m["MESH/COORD"][:].reshape(-1, 3)
with h5py.File(f"{run}/res_{step}.h5") as f:
    V = f["VALUE"]; N = V["T"].shape[0]
    xi, xv, chi = V["xi"][:], V["xiVar"][:], V["chi"][:]
    YH2 = V["Y0"][:]; T = V["T"][:]
x, y = xyz[:N, 0], xyz[:N, 1]; D = 0.00457
Y_H2_fuel = 0.25 * 2.01588 / (0.25 * 2.01588 + 0.75 * 28.0134)
xi_ref = np.clip(YH2 / Y_H2_fuel, 0, 1)
err = np.abs(xi - xi_ref)
print(f"xi: range {xi.min():.4f}..{xi.max():.4f}; |xi - Y_H2/Y_H2,fuel| max {err.max():.2e} mean {err.mean():.2e} (反応 OFF なら ~1e-6)")
viol = xv - xi * (1 - xi)
print(f"xiVar: range {xv.min():.3e}..{xv.max():.3e}; realizability max(xiVar - xi(1-xi)) = {viol.max():.2e}; negative count {np.sum(xv < 0)}")
print(f"chi: range {chi.min():.3e}..{chi.max():.3e} 1/s")
for zd in (5, 10, 20, 30):
    m = np.where((np.abs(x - zd * D) < 0.5 * D) & (y < 0.03))[0]
    i = m[np.argmax(xv[m])]
    print(f"  z/d={zd:2d}: max xiVar={xv[m].max():.3e} at r={y[i]*1e3:.2f} mm (xi={xi[i]:.3f}, chi={chi[i]:.3g} 1/s, sqrt(var)/xi={np.sqrt(xv[i])/max(xi[i],1e-9):.2f})")
ax = np.where(np.abs(y) < 1e-9)[0]; o = np.argsort(x[ax])
print("axis xi @ z/d 10/20/30/40:", [f"{np.interp(zd*D, x[ax][o], xi[ax][o]):.3f}" for zd in (10, 20, 30, 40)])
