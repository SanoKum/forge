#!/usr/bin/env python3
"""ユーザ指定の上壁プロファイル (2026-08-19, 単位 mm) と現行 Fig.3 メッシュ形状 (make_nozzle_fig3.py) の比較。

ユーザ式 (半高 y [mm], x [mm], throat x=0):
  12.7                                        (-60 <= x <= -38)
  0.16 - 0.33 x                               (-38 < x < -21.2727)
  2.5 + 0.015512821 x^2 + 0.000243078 x^3     (-21.2727 < x < 0)
  2.5 + 0.00194105 x^2 - 7.99459e-5 x^3       (0 < x < 8.093175942)
  2.457620744 + 0.015709255 x                 (8.093175942 < x < 95)
出力: geometry_user_vs_current.png, 標準出力に差分表。
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "mesh"))
from make_nozzle_fig3 import y_upper as y_current_cm   # 現行 (cm)
from nozzle_user_profile import y_user_mm              # ユーザ式 (mm)

GAM = 1.4

def area_mach(ar, sup):
    """A/A* -> M (1D 等エントロピー)"""
    from scipy.optimize import brentq
    f = lambda M: (1/M)*((2/(GAM+1))*(1+(GAM-1)/2*M*M))**((GAM+1)/(2*(GAM-1))) - ar
    if ar <= 1.0 + 1e-12: return 1.0
    return brentq(f, 1.0, 10.0) if sup else brentq(f, 1e-6, 1.0)

x = np.linspace(-38.0, 95.0, 2661)                    # 現行メッシュの範囲 (mm)
yu = np.array([y_user_mm(v) for v in x])
yc = np.array([y_current_cm(v/10.0)*10.0 for v in x])
d = yu - yc
print("x range compared: [-38, 95] mm (current mesh has no straight inlet; user has -60..-38 at 12.7)")
print(f"throat half [mm]  user {y_user_mm(0.0):.4f}  current {y_current_cm(0.0)*10:.4f}")
print(f"exit   half [mm]  user {y_user_mm(95.0):.4f}  current {y_current_cm(9.5)*10:.4f}")
print(f"A/A*_exit         user {y_user_mm(95.0)/2.5:.4f}  current {y_current_cm(9.5)*10/2.5:.4f}")
print(f"max |dy| [mm] converging: {np.abs(d[x<0]).max():.3f} at x={x[x<0][np.abs(d[x<0]).argmax()]:.1f}")
print(f"max |dy| [mm] diverging : {np.abs(d[x>0]).max():.3f} at x={x[x>0][np.abs(d[x>0]).argmax()]:.1f}")
# throat curvature (2nd deriv): user 2*0.015512821 (conv) / 2*0.00194105 (div)
h=1e-3
def d2(f, x0): return (f(x0+h)-2*f(x0)+f(x0-h))/h/h
print(f"y'' at throat [1/mm]: user conv {2*0.015512821:.5f} div {2*0.00194105:.5f} ; current conv {d2(lambda v:y_current_cm(v/10)*10,-1e-3):.5f} div {d2(lambda v:y_current_cm(v/10)*10,1e-3):.5f}")
print(f"radius of curvature at throat R=1/y'' [mm]: user conv {1/(2*0.015512821):.1f}, div {1/(2*0.00194105):.1f} ; current conv {1/max(d2(lambda v:y_current_cm(v/10)*10,-1e-3),1e-9):.1f} div inf (linear)")
# 1D isentropic
xs = x[x>0]
Mu = np.array([area_mach(y_user_mm(v)/2.5, True) for v in xs])
Mc = np.array([area_mach(y_current_cm(v/10)*10/2.5, True) for v in xs])
pu = (1+(GAM-1)/2*Mu**2)**(-GAM/(GAM-1)); pc = (1+(GAM-1)/2*Mc**2)**(-GAM/(GAM-1))
for xq in [5, 10, 20, 30, 50, 70, 95]:
    i = np.argmin(np.abs(xs-xq))
    print(f"x={xq:3d} mm  M user {Mu[i]:.3f} cur {Mc[i]:.3f} (d {Mu[i]-Mc[i]:+.3f})  p/p0 user {pu[i]:.4f} cur {pc[i]:.4f} (d {100*(pu[i]/pc[i]-1):+.1f}%)")

fig, ax = plt.subplots(3, 1, figsize=(9, 10), sharex=False)
xf = np.linspace(-60, 95, 3101)
ax[0].plot(xf, [y_user_mm(v) for v in xf], "r-", label="user (2026-08-19)")
ax[0].plot(x, yc, "b--", label="current mesh (make_nozzle_fig3.py)")
ax[0].plot(xf, [-y_user_mm(v) for v in xf], "r-"); ax[0].plot(x, -yc, "b--")
ax[0].set_ylabel("y [mm]"); ax[0].legend(); ax[0].grid(alpha=.3); ax[0].set_title("Wyslouzil Fig.3 nozzle: user profile vs current mesh")
ax[1].plot(x, d, "k-"); ax[1].axhline(0, color="gray", lw=.5)
ax[1].set_ylabel("y_user - y_current [mm]"); ax[1].grid(alpha=.3)
ax[2].plot(xs, Mu, "r-", label="M user"); ax[2].plot(xs, Mc, "b--", label="M current")
ax2 = ax[2].twinx(); ax2.plot(xs, pu, "r:", label="p/p0 user"); ax2.plot(xs, pc, "b:", label="p/p0 current"); ax2.set_ylabel("p/p0")
ax[2].set_ylabel("M (1D isentropic, γ=1.4)"); ax[2].set_xlabel("x [mm]"); ax[2].legend(loc="center right"); ax[2].grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "geometry_user_vs_current.png"), dpi=130)
print("wrote geometry_user_vs_current.png")
