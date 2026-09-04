#!/usr/bin/env python3
"""Burrows–Kurkov (NASA TM X-2828) 2D 平面構造メッシュ (gmsh Python API, node 用)。
形状 [m]: 主流ダクト高 0.089 (ステップ上), スロット高 0.004, リップ 0.00076 → ステップ高 0.00476。
ステップ後 (x=0) の燃焼器は上壁が y=0.09376 → 0.105 (x=0.356) へ直線拡大 (BL 補償)。上流ダクト長 0.05, スロット長 0.02。
  .venv-mesh/bin/python make_bk_mesh.py [--nx 360] [--dy1 5e-5] [--out bk.msh]"""
import argparse, gmsh, math
ap = argparse.ArgumentParser(); ap.add_argument("--nx", type=int, default=360); ap.add_argument("--nxu", type=int, default=50)
ap.add_argument("--nxs", type=int, default=30); ap.add_argument("--ny_slot", type=int, default=24); ap.add_argument("--ny_lip", type=int, default=4)
ap.add_argument("--ny_air", type=int, default=140); ap.add_argument("--dy1", type=float, default=5e-5); ap.add_argument("--out", default="bk.msh"); a = ap.parse_args()
Lu, Ls, Lc = 0.05, 0.02, 0.356; hs, hl, Ha = 0.004, 0.00076, 0.089
y1 = hs; y2 = hs + hl; ytop0 = y2 + Ha; ytop1 = 0.105
gmsh.initialize(); gmsh.model.add("bk")
P = lambda x, y: gmsh.model.geo.addPoint(x, y, 0.0)
# points
pS0 = P(-Ls, 0); pS1 = P(0, 0); pS2 = P(0, y1); pS3 = P(-Ls, y1)                        # slot
pL2 = P(0, y2)                                                                            # step top (lip top corner)
pA0 = P(-Lu, y2); pA1 = P(-Lu, ytop0); pA2 = P(0, ytop0)                                  # air duct
pC0 = P(Lc, 0); pC1 = P(Lc, y1); pC2 = P(Lc, y2); pC3 = P(Lc, ytop1)                      # exit
L = gmsh.model.geo.addLine
# slot block: pS0-pS1-pS2-pS3
s_bot = L(pS0, pS1); s_face = L(pS1, pS2); s_top = L(pS2, pS3); s_in = L(pS3, pS0)
# lip step face: pS2-pL2 (wall), air duct: pL2-pA0? no: air duct bottom is pA0-pL2
a_bot = L(pA0, pL2); a_face = L(pL2, pA2); a_top = L(pA2, pA1); a_in = L(pA1, pA0)
lip_face = L(pS2, pL2)
# combustor: 3 horizontal blocks sharing x-lines
c_bot = L(pS1, pC0); c_out1 = L(pC0, pC1); c_mid1 = L(pC1, pS2)                           # C1: slot wake
c_out2 = L(pC1, pC2); c_mid2 = L(pC2, pL2)                                               # C2: lip wake
c_out3 = L(pC2, pC3); c_top = L(pC3, pA2)                                                # C3: air
def surf(lines):
    cl = gmsh.model.geo.addCurveLoop(lines); s = gmsh.model.geo.addPlaneSurface([cl]); return s
S_slot = surf([s_bot, s_face, s_top, s_in])
S_air  = surf([a_bot, a_face, a_top, a_in])
S_C1 = surf([c_bot, c_out1, c_mid1, -s_face])
S_C2 = surf([-c_mid1, c_out2, c_mid2, -lip_face])
S_C3 = surf([-c_mid2, c_out3, c_top, -a_face])
gmsh.model.geo.synchronize()
T = gmsh.model.geo.mesh.setTransfiniteCurve
# x 方向
for l in (s_bot, s_top): T(l, a.nxs + 1, "Progression", 1.0)
for l in (a_bot, a_top): T(l, a.nxu + 1, "Bump", 0.3)
# 進行比は曲線の向きに沿う: c_bot は入口→出口 (1.004), c_mid1/c_mid2/c_top は出口→入口向きなので逆数 (列を揃えないと平行四辺形セルで skew 0.98)
T(c_bot, a.nx + 1, "Progression", 1.004)
for l in (c_mid1, c_mid2, c_top): T(l, a.nx + 1, "Progression", 1.0/1.004)
# y 方向: スロット (両壁に寄せる), リップ, 主流 (下壁+リップ近傍に寄せ, 上壁にも寄せる)
for l in (s_face, s_in, c_out1): T(l, a.ny_slot + 1, "Bump", 0.15)
for l in (lip_face, c_out2): T(l, a.ny_lip + 1, "Progression", 1.0)
# 主流: 幾何級数で下側を細かく (最初の間隔 dy1)。Bump だと両端が細かい: 下壁側 dy1, 上壁側は緩め
r = 1.0
# 幾何級数比 r を dy1 と長さ・分割から解く: Ha = dy1 (r^N - 1)/(r - 1)
N = a.ny_air
lo, hi = 1.0001, 1.5
for _ in range(80):
    r = 0.5*(lo+hi); s = a.dy1*(r**N - 1)/(r - 1)
    if s > Ha: hi = r
    else: lo = r
for l in (a_face, a_in, c_out3): T(l, N + 1, "Progression", r)
# 主流の向き: a_face は pL2→pA2 (下→上) で進行比 r>1 は下が細かい。a_in は pA1→pA0 (上→下) なので 1/r。c_out3 は pC2→pC3 (下→上)。
T(a_in, N + 1, "Progression", 1.0/r)
for s in (S_slot, S_air, S_C1, S_C2, S_C3):
    gmsh.model.geo.mesh.setTransfiniteSurface(s); gmsh.model.geo.mesh.setRecombine(2, s)
gmsh.model.geo.synchronize()
gmsh.model.addPhysicalGroup(1, [a_in], 1, "inlet_air")
gmsh.model.addPhysicalGroup(1, [s_in], 2, "inlet_h2")
gmsh.model.addPhysicalGroup(1, [c_out1, c_out2, c_out3], 3, "outlet")
gmsh.model.addPhysicalGroup(1, [s_bot, c_bot, s_top, lip_face, a_bot, a_top, c_top], 4, "wall")
gmsh.model.addPhysicalGroup(2, [S_slot, S_air, S_C1, S_C2, S_C3], 7, "fluid")
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.model.mesh.generate(2)
nn = len(gmsh.model.mesh.getNodes()[0]); print(f"nodes={nn} dy1={a.dy1} r={r:.4f} top dy≈{a.dy1*r**(N-1):.2e}")
gmsh.write(a.out); gmsh.finalize()
