#!/usr/bin/env python3
"""BK 入口の壁法則 BL プロファイル (1/7 乗則) を inlet_profile_<physID>.csv に書く。
  空気入口 (physID 1, y∈[0.00476, 0.09376]): 上下壁に δ の BL。H2 スロット (physID 2, y∈[0,0.004]): チャネル形。
  .venv-chem/bin/python gen_bk_inlet_profiles.py run_dir --UA 1784.1 --UH 1217.0 [--delta 0.006]"""
import argparse, numpy as np, pathlib
ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--UA", type=float, required=True); ap.add_argument("--UH", type=float, required=True)
ap.add_argument("--delta", type=float, default=0.006); ap.add_argument("--n", type=float, default=7.0)
ap.add_argument("--Tw", type=float, default=350.0); ap.add_argument("--TA", type=float, default=1270.0); ap.add_argument("--TH", type=float, default=254.0)
ap.add_argument("--RA", type=float, default=329.7); ap.add_argument("--RH", type=float, default=4124.0); ap.add_argument("--cpA", type=float, default=1400.0); ap.add_argument("--cpH", type=float, default=14500.0)
ap.add_argument("--p", type=float, default=101325.0); a = ap.parse_args()
d = pathlib.Path(a.run_dir)
def bl(dist, delta): return np.clip(dist/delta, 0, 1)**(1.0/a.n)
y0, y1 = 0.00476, 0.09376; y = np.concatenate([np.linspace(y0, y0+a.delta, 60), np.linspace(y0+a.delta, y1-a.delta, 40)[1:-1], np.linspace(y1-a.delta, y1, 60)])
u = a.UA*bl(y-y0, a.delta)*bl(y1-y, a.delta)
# 温度は Crocco–Busemann (r=0.89): T = Tw + (Te−Tw) u/Ue + r Ue²/(2cp) (u/Ue)(1−u/Ue), ρ = p/(R T)。壁節点で T=Tw に整合。
def crocco(u, Ue, Te, Tw, cp): f = u/Ue; return Tw + (Te-Tw)*f + 0.89*Ue**2/(2*cp)*f*(1-f)
T = crocco(u, a.UA, a.TA, a.Tw, a.cpA); ro = a.p/(a.RA*T)
np.savetxt(d/"inlet_profile_1.csv", np.c_[y, u, 0*u, 0*u, ro], header="y Ux Uy Uz ro", comments="", fmt="%.7e")
h = 0.004; ys = np.linspace(0, h, 81); us = a.UH*(1.0-np.abs(2*ys/h-1.0))**(1.0/a.n)
Ts = crocco(us, a.UH, a.TH, a.Tw, a.cpH); ros = a.p/(a.RH*Ts)
np.savetxt(d/"inlet_profile_2.csv", np.c_[ys, us, 0*us, 0*us, ros], header="y Ux Uy Uz ro", comments="", fmt="%.7e")
print("wrote", d/"inlet_profile_1.csv", d/"inlet_profile_2.csv", "u_air max", u.max(), "mean/U", u.mean()/a.UA)
