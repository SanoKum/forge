"""摂動オーダー仮説の検証: R を上げると Sauer/円弧の壁角不整合は縮むか。"""
import sys
sys.path.insert(0, "design")
import numpy as np
from forge_design.geometry.transonic import SauerThroat
from forge_design.geometry.bezier import MachBezier
from forge_design.geometry.moc_inverse import inverse_design

G, Md = 1.4, 4.0
print(f"{'R':>5} {'ubar(壁)':>9} {'Sauer θ_w':>10} {'円弧 θ':>8} {'不整合':>8} "
      f"{'設計壁 tanθ[1]':>14} {'円弧 tanθ':>10} {'比':>6}")
for R in (2.0, 5.0, 10.0, 20.0):
    st = SauerThroat(R=R, gamma=G)
    x0, rr, MM, tt = st.starting_line(M_start=1.05, n=41)
    ub = float(st.ubar(x0, rr[-1]))
    th_arc = float(np.arcsin(min(x0 / R, 1.0)))
    mism = np.degrees(th_arc - tt[-1])
    h = 1e-4
    M0 = float(st.mach(x0, 0.0)); dM0 = (st.mach(x0+h,0.)-st.mach(x0-h,0.))/(2*h)
    xd = 14.0
    eps = (1/Md)*((2+(G-1)*Md*Md)/(G+1))**((G+1)/(2*(G-1)))
    x_end = xd + 2.3*np.sqrt(eps)/np.tan(np.arcsin(1/Md))
    bz = MachBezier.from_constraints(x0, xd, start=(M0, dM0), free_cp=[2.0,2.9,3.6],
                                     end=(Md,0.,0.))
    tgt = lambda x: Md if x >= xd else float(np.atleast_1d(bz(float(x)))[0])  # noqa
    res = inverse_design(st, tgt, x_axis_end=float(x_end), n_axis=400, n_start=41,
                         gamma=G, th_wall0=th_arc)
    w = res["wall"]
    t1, ta = float(np.tan(w[1,2])), float(np.tan(th_arc))
    print(f"{R:5.1f} {ub:9.4f} {np.degrees(tt[-1]):10.3f} {np.degrees(th_arc):8.3f} "
          f"{mism:8.3f} {t1:14.4f} {ta:10.4f} {t1/ta:6.3f}")
