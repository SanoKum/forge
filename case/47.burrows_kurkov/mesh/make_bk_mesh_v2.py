#!/usr/bin/env python3
"""Burrows–Kurkov 平面構造メッシュ v2: 主流ブロックを上下 2 分割し、上壁側にも第 1 間隔 dy1 のクラスタ (壁関数 y+~30-50)。
形状は make_bk_mesh.py と同一。全壁 no-slip (wall=4)、入口 BL プロファイル + nodeInletCornerWall 併用前提。
  .venv-mesh/bin/python make_bk_mesh_v2.py [--nx 360] [--dy1 5e-5] [--ny_half 100] [--out bk_v2_noslip.msh]"""
import argparse, gmsh
ap = argparse.ArgumentParser(); ap.add_argument("--nx", type=int, default=360); ap.add_argument("--nxu", type=int, default=50)
ap.add_argument("--nxs", type=int, default=30); ap.add_argument("--ny_slot", type=int, default=24); ap.add_argument("--ny_lip", type=int, default=4)
ap.add_argument("--ny_half", type=int, default=100); ap.add_argument("--dy1", type=float, default=5e-5); ap.add_argument("--out", default="bk_v2_noslip.msh"); a = ap.parse_args()
Lu, Ls, Lc = 0.05, 0.02, 0.356; hs, hl, Ha = 0.004, 0.00076, 0.089
y1 = hs; y2 = hs + hl; ytop0 = y2 + Ha; ytop1 = 0.105; ym0 = y2 + Ha/2; ym1 = (y2 + ytop1)/2
gmsh.initialize(); gmsh.model.add("bk2")
P = lambda x, y: gmsh.model.geo.addPoint(x, y, 0.0); L = gmsh.model.geo.addLine
pS0 = P(-Ls, 0); pS1 = P(0, 0); pS2 = P(0, y1); pS3 = P(-Ls, y1); pL2 = P(0, y2)
pA0 = P(-Lu, y2); pAm = P(-Lu, ym0); pA1 = P(-Lu, ytop0); pLm = P(0, ym0); pA2 = P(0, ytop0)
pC0 = P(Lc, 0); pC1 = P(Lc, y1); pC2 = P(Lc, y2); pCm = P(Lc, ym1); pC3 = P(Lc, ytop1)
s_bot = L(pS0, pS1); s_face = L(pS1, pS2); s_top = L(pS2, pS3); s_in = L(pS3, pS0); lip_face = L(pS2, pL2)
a_bot = L(pA0, pL2); a_face_lo = L(pL2, pLm); a_mid = L(pLm, pAm); a_in_lo = L(pAm, pA0)
a_face_hi = L(pLm, pA2); a_top = L(pA2, pA1); a_in_hi = L(pA1, pAm)
c_bot = L(pS1, pC0); c_out1 = L(pC0, pC1); c_mid1 = L(pC1, pS2); c_out2 = L(pC1, pC2); c_mid2 = L(pC2, pL2)
c_out3a = L(pC2, pCm); c_mid3 = L(pCm, pLm); c_out3b = L(pCm, pC3); c_top = L(pC3, pA2)
def surf(lines):
    cl = gmsh.model.geo.addCurveLoop(lines); return gmsh.model.geo.addPlaneSurface([cl])
S_slot = surf([s_bot, s_face, s_top, s_in]); S_alo = surf([a_bot, a_face_lo, a_mid, a_in_lo]); S_ahi = surf([-a_mid, a_face_hi, a_top, a_in_hi])
S_C1 = surf([c_bot, c_out1, c_mid1, -s_face]); S_C2 = surf([-c_mid1, c_out2, c_mid2, -lip_face])
S_C3a = surf([-c_mid2, c_out3a, c_mid3, -a_face_lo]); S_C3b = surf([-c_mid3, c_out3b, c_top, -a_face_hi])
gmsh.model.geo.synchronize(); T = gmsh.model.geo.mesh.setTransfiniteCurve
for l in (s_bot, s_top): T(l, a.nxs + 1, "Progression", 1.0)
for l in (a_bot, a_mid, a_top): T(l, a.nxu + 1, "Bump", 0.3)
T(c_bot, a.nx + 1, "Progression", 1.004)
for l in (c_mid1, c_mid2, c_mid3, c_top): T(l, a.nx + 1, "Progression", 1.0/1.004)
for l in (s_face, s_in, c_out1): T(l, a.ny_slot + 1, "Bump", 0.15)
for l in (lip_face, c_out2): T(l, a.ny_lip + 1, "Progression", 1.0)
# 半区間 Ha/2 を dy1 から幾何級数: 比 r を二分法で
M = a.ny_half; lo, hi = 1.0001, 1.5
for _ in range(80):
    r = 0.5*(lo+hi); s = a.dy1*(r**M - 1)/(r - 1)
    if s > Ha/2: hi = r
    else: lo = r
# 下半分: 壁 (下) → 中線: 下→上向きの線は r, 上→下向きは 1/r。上半分: 中線 → 上壁: 上壁側を細かく = 上→下向きが r
T(a_face_lo, M + 1, "Progression", r)        # pL2→pLm (下→上)
T(a_in_lo,   M + 1, "Progression", 1.0/r)    # pAm→pA0 (上→下)
T(c_out3a,   M + 1, "Progression", r)        # pC2→pCm (下→上)
T(a_face_hi, M + 1, "Progression", 1.0/r)    # pLm→pA2 (下→上): 上端で細かく
T(a_in_hi,   M + 1, "Progression", r)        # pA1→pAm (上→下): 上端で細かく
T(c_out3b,   M + 1, "Progression", 1.0/r)    # pCm→pC3 (下→上)
for sfc in (S_slot, S_alo, S_ahi, S_C1, S_C2, S_C3a, S_C3b):
    gmsh.model.geo.mesh.setTransfiniteSurface(sfc); gmsh.model.geo.mesh.setRecombine(2, sfc)
gmsh.model.geo.synchronize()
gmsh.model.addPhysicalGroup(1, [a_in_lo, a_in_hi], 1, "inlet_air"); gmsh.model.addPhysicalGroup(1, [s_in], 2, "inlet_h2")
gmsh.model.addPhysicalGroup(1, [c_out1, c_out2, c_out3a, c_out3b], 3, "outlet")
gmsh.model.addPhysicalGroup(1, [c_bot, lip_face, c_top, s_bot, s_top, a_bot, a_top], 4, "wall")
gmsh.model.addPhysicalGroup(2, [S_slot, S_alo, S_ahi, S_C1, S_C2, S_C3a, S_C3b], 7, "fluid")
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1); gmsh.model.mesh.generate(2)
print(f"nodes={len(gmsh.model.mesh.getNodes()[0])} r={r:.4f} mid dy≈{a.dy1*r**(M-1):.2e}"); gmsh.write(a.out); gmsh.finalize()
