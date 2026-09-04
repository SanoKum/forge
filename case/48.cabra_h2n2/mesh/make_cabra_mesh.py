#!/usr/bin/env python3
"""Cabra H2/N2 浮き上がり火炎 (Berkeley vitiated coflow burner) の軸対称平面構造メッシュ (node 用)。
x: 軸方向 (ノズル出口 x=0 から x=Lx), y: 半径。ジェット管 ID 4.57 mm (r_j=2.285), OD 6.35 mm (リップ 0.89 mm),
coflow 半径 105 mm (計算域は R=100 mm まで、外周は slip)。ジェット管は上流に L_tube だけ延長 (管内 no-slip で発達)。
  .venv-mesh/bin/python make_cabra_mesh.py [--Lx 0.25] [--R 0.10] [--Ltube 0.02] [--out cabra.msh]"""
import argparse, gmsh
ap = argparse.ArgumentParser(); ap.add_argument("--Lx", type=float, default=0.25); ap.add_argument("--R", type=float, default=0.10)
ap.add_argument("--Ltube", type=float, default=0.02); ap.add_argument("--nx", type=int, default=300); ap.add_argument("--nxt", type=int, default=40)
ap.add_argument("--ny_jet", type=int, default=30); ap.add_argument("--ny_lip", type=int, default=6); ap.add_argument("--ny_co1", type=int, default=90); ap.add_argument("--ny_co2", type=int, default=50)
ap.add_argument("--out", default="cabra.msh"); a = ap.parse_args()
rj, ro_ = 0.002285, 0.003175; R1 = 0.030   # 内側 coflow ブロック境界 (r=30 mm まで細かく)
gmsh.initialize(); gmsh.model.add("cabra")
P = lambda x, y: gmsh.model.geo.addPoint(x, y, 0.0); L = gmsh.model.geo.addLine
# 点: 管上流 (x=-Ltube), 出口面 (x=0), 出口 (x=Lx)
pT0 = P(-a.Ltube, 0); pT1 = P(-a.Ltube, rj)
pE0 = P(0, 0); pE1 = P(0, rj); pE2 = P(0, ro_); pE3 = P(0, R1); pE4 = P(0, a.R)
pC0 = P(a.Lx, 0); pC1 = P(a.Lx, rj); pC2 = P(a.Lx, ro_); pC3 = P(a.Lx, R1); pC4 = P(a.Lx, a.R)
pA2 = P(-a.Ltube, ro_); pA3 = P(-a.Ltube, R1); pA4 = P(-a.Ltube, a.R)   # coflow も上流 Ltube から (管外壁沿い)
# 管内ブロック T: (pT0,pE0,pE1,pT1); 管壁ブロックは固体 (無し): coflow 上流ブロック A1: (pA2,pE2,pE3,pA3), A2: (pA3,pE3,pE4,pA4)
t_axis = L(pT0, pE0); e_jet = L(pE0, pE1); t_top = L(pE1, pT1); t_in = L(pT1, pT0)
a_bot = L(pA2, pE2); e_lip = L(pE1, pE2); e_co1 = L(pE2, pE3); a_mid = L(pE3, pA3); a_in1 = L(pA3, pA2)
e_co2 = L(pE3, pE4); a_top = L(pE4, pA4); a_in2 = L(pA4, pA3)
c_axis = L(pE0, pC0); c_o0 = L(pC0, pC1); c_m1 = L(pC1, pE1); c_o1 = L(pC1, pC2); c_m2 = L(pC2, pE2)
c_o2 = L(pC2, pC3); c_m3 = L(pC3, pE3); c_o3 = L(pC3, pC4); c_top = L(pC4, pE4)
def surf(ls):
    cl = gmsh.model.geo.addCurveLoop(ls); return gmsh.model.geo.addPlaneSurface([cl])
S_T = surf([t_axis, e_jet, t_top, t_in]); S_A1 = surf([a_bot, e_co1, a_mid, a_in1]); S_A2 = surf([-a_mid, e_co2, a_top, a_in2])
S_C0 = surf([c_axis, c_o0, c_m1, -e_jet]); S_C1 = surf([-c_m1, c_o1, c_m2, -e_lip]); S_C2 = surf([-c_m2, c_o2, c_m3, -e_co1]); S_C3 = surf([-c_m3, c_o3, c_top, -e_co2])
gmsh.model.geo.synchronize(); T = gmsh.model.geo.mesh.setTransfiniteCurve
for l in (t_axis, t_top, a_bot, a_mid, a_top): T(l, a.nxt + 1, "Progression", 1.0)
T(c_axis, a.nx + 1, "Progression", 1.008)
for l in (c_m1, c_m2, c_m3, c_top): T(l, a.nx + 1, "Progression", 1.0/1.008)
# 半径: ジェット (壁側に寄せ), リップ (等間隔), coflow 内側 (リップ側に寄せ、外へ幾何級数), 外側 (粗)
for l in (e_jet, t_in, c_o0): T(l, a.ny_jet + 1, "Bump", 0.25)
for l in (e_lip, c_o1): T(l, a.ny_lip + 1, "Progression", 1.0)
# co1: pE2→pE3 (内→外) で内側細かく: 幾何級数 r
M = a.ny_co1; dy1 = 1.5e-4; lo, hi = 1.0001, 1.5
for _ in range(80):
    r = 0.5*(lo+hi); s = dy1*(r**M - 1)/(r - 1)
    if s > (R1 - ro_): hi = r
    else: lo = r
T(e_co1, M + 1, "Progression", r); T(c_o2, M + 1, "Progression", r); T(a_in1, M + 1, "Progression", 1.0/r)
for l in (e_co2, c_o3): T(l, a.ny_co2 + 1, "Progression", 1.02)
T(a_in2, a.ny_co2 + 1, "Progression", 1.0/1.02)
for sfc in (S_T, S_A1, S_A2, S_C0, S_C1, S_C2, S_C3):
    gmsh.model.geo.mesh.setTransfiniteSurface(sfc); gmsh.model.geo.mesh.setRecombine(2, sfc)
gmsh.model.geo.synchronize()
gmsh.model.addPhysicalGroup(1, [t_in], 1, "inlet_jet"); gmsh.model.addPhysicalGroup(1, [a_in1, a_in2], 2, "inlet_coflow")
gmsh.model.addPhysicalGroup(1, [c_o0, c_o1, c_o2, c_o3], 3, "outlet")
gmsh.model.addPhysicalGroup(1, [t_top, a_bot, e_lip], 4, "wall")          # 管内壁・管外壁・リップ端面
gmsh.model.addPhysicalGroup(1, [a_top, c_top], 5, "farfield")            # 外周 (slip)
gmsh.model.addPhysicalGroup(1, [t_axis, c_axis], 6, "axis")
gmsh.model.addPhysicalGroup(2, [S_T, S_A1, S_A2, S_C0, S_C1, S_C2, S_C3], 7, "fluid")
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1); gmsh.model.mesh.generate(2)
print(f"nodes={len(gmsh.model.mesh.getNodes()[0])} co1 ratio={r:.4f} outer dy≈{dy1*r**(M-1):.2e}"); gmsh.write(a.out); gmsh.finalize()
