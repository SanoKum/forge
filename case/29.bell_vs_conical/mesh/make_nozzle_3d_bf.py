#!/usr/bin/env python3
r"""ノズル 90° セクタ 3D を **バタフライ(O-grid)** で生成 (Gmsh 4.1 .msh 直接出力).

make_nozzle_3d.py (単一 squircle) では壁第一層の厚さが周方向で不均一だった
(対称面で粗く、対角で細かい)。ここでは:

  - **core**: 半径 r_inner = fcore*R(x) の小さい squircle ブロック (軸近傍, 非特異, BL 無し)。
    squircle の外縁 (u=1 ∪ v=1) はちょうど半径 r_inner の円弧上 (y^2+z^2=r_inner^2)。
  - **ring**: その外縁ノードを **半径方向に** r_inner→R(x) へ延ばす環状ブロック。
    radial 分布は壁(n=Nrad)で固定第一層厚 dwall, 内側へ幾何級数成長。
    → **全周で壁第一層厚が dwall に均一** (ご要望)。radial 線は真に半径方向=壁直交。

境界 physID: inlet=1, outlet=2, wall=3, sym_y0=4(y=0面), sym_z0=5(z=0面)。fluid=8。
座標 m。
"""
from __future__ import annotations
import argparse
import math
import numpy as np


def axial_stations(x0, xth, x1, nconv, ndiv, b):
    tc = np.linspace(0, 1, nconv + 1); sc = np.sinh(b * tc) / math.sinh(b)
    xc = x0 + (xth - x0) * sc
    td = np.linspace(0, 1, ndiv + 1); sd = 1.0 - np.sinh(b * (1 - td)) / math.sinh(b)
    xd = xth + (x1 - xth) * sd
    return np.concatenate([xc[:-1], xd])


def radial_nodes(r_inner, R, Nrad, dwall):
    """r[0]=r_inner .. r[Nrad]=R, 壁側(末尾)セル厚=dwall, 内側へ幾何成長。"""
    L = R - r_inner
    if dwall * Nrad >= L:               # dwall が大きすぎる→一様
        return np.linspace(r_inner, R, Nrad + 1)
    # solve growth g: dwall*(g^Nrad - 1)/(g-1) = L
    lo, hi = 1.0 + 1e-6, 3.0
    for _ in range(100):
        g = 0.5 * (lo + hi)
        s = dwall * (g**Nrad - 1.0) / (g - 1.0)
        if s > L: hi = g
        else: lo = g
    g = 0.5 * (lo + hi)
    # cell sizes from wall inward: dwall, dwall*g, ...
    cells = dwall * g**np.arange(Nrad)          # [Nrad], cells[0] adjacent to wall
    r = np.empty(Nrad + 1); r[Nrad] = R
    for n in range(Nrad - 1, -1, -1):
        r[n] = r[n + 1] - cells[Nrad - 1 - n]
    r[0] = r_inner
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contour", default="contour_conical.csv")
    ap.add_argument("--rt", type=float, default=10.0)
    ap.add_argument("--nc", type=int, default=12, help="core squircle 1 辺分割 (周方向 ring = 2*nc)")
    ap.add_argument("--nrad", type=int, default=34, help="ring 半径方向分割")
    ap.add_argument("--fcore", type=float, default=0.5, help="r_inner/R")
    ap.add_argument("--dwall", type=float, default=0.008, help="壁第一層厚 [mm] (全周一定)")
    ap.add_argument("--nconv", type=int, default=46)
    ap.add_argument("--ndiv", type=int, default=84)
    ap.add_argument("--bax", type=float, default=2.0)
    ap.add_argument("--out", default="nozzle3d_bf.msh")
    args = ap.parse_args()

    cont = np.loadtxt(args.contour, delimiter=",", skiprows=1, usecols=(0, 1))
    cx, cy = cont[:, 0], cont[:, 1]; o = np.argsort(cx); cx, cy = cx[o], cy[o]
    xa = axial_stations(cx.min(), 0.0, cx.max(), args.nconv, args.ndiv, args.bax)
    Rw = np.interp(xa, cx, cy); Rw[np.argmin(np.abs(xa))] = args.rt
    Nx = len(xa) - 1
    Nc, Nrad = args.nc, args.nrad
    Mcirc = 2 * Nc                      # 周方向セグメント数
    sc = np.linspace(0.0, 1.0, Nc + 1)  # core squircle param (uniform)

    n_core = (Nc + 1) * (Nc + 1)
    n_ring = (Mcirc + 1) * Nrad         # n=1..Nrad
    n_percs = n_core + n_ring

    def cidx(i, j): return i * (Nc + 1) + j            # core local (0-based)

    # core outer-edge local id for ring circumferential index m (n=0 蓋)
    def outer_core_local(m):
        if m <= Nc: return cidx(Nc, m)                 # u=1 edge, v=sc[m]
        return cidx(Mcirc - m, Nc)                     # v=1 edge, u=sc[2Nc-m]

    def ring_local(m, n):
        if n == 0: return outer_core_local(m)
        return n_core + (n - 1) * (Mcirc + 1) + m

    # circumferential unit directions (cosθ,sinθ) from core outer edge (radius-independent)
    dirs = np.empty((Mcirc + 1, 2))
    for m in range(Mcirc + 1):
        if m <= Nc:
            v = sc[m]; dirs[m] = (math.sqrt(1 - 0.5 * v * v), v / math.sqrt(2))
        else:
            u = sc[Mcirc - m]; dirs[m] = (u / math.sqrt(2), math.sqrt(1 - 0.5 * u * u))

    # ---- nodes ----
    coords = np.empty(((Nx + 1) * n_percs, 3))
    for a in range(Nx + 1):
        R = Rw[a] * 1e-3; x = xa[a] * 1e-3; ri = args.fcore * R
        base = a * n_percs
        # core
        for i in range(Nc + 1):
            u = sc[i]
            for j in range(Nc + 1):
                v = sc[j]
                y = ri * u * math.sqrt(1 - 0.5 * v * v)
                z = ri * v * math.sqrt(1 - 0.5 * u * u)
                coords[base + cidx(i, j)] = (x, y, z)
        # ring (n=1..Nrad)
        rr = radial_nodes(ri, R, Nrad, args.dwall * 1e-3)
        for n in range(1, Nrad + 1):
            for m in range(Mcirc + 1):
                yz = rr[n] * dirs[m]
                coords[base + n_core + (n - 1) * (Mcirc + 1) + m] = (x, yz[0], yz[1])
    Nnode = coords.shape[0]

    def gid(a, loc): return a * n_percs + loc + 1     # 1-based

    # ---- hexes ----
    hexes = []
    for a in range(Nx):
        for i in range(Nc):
            for j in range(Nc):
                hexes.append([gid(a, cidx(i, j)), gid(a + 1, cidx(i, j)), gid(a + 1, cidx(i + 1, j)), gid(a, cidx(i + 1, j)),
                              gid(a, cidx(i, j + 1)), gid(a + 1, cidx(i, j + 1)), gid(a + 1, cidx(i + 1, j + 1)), gid(a, cidx(i + 1, j + 1))])
        for m in range(Mcirc):
            for n in range(Nrad):
                hexes.append([gid(a, ring_local(m, n)), gid(a + 1, ring_local(m, n)), gid(a + 1, ring_local(m + 1, n)), gid(a, ring_local(m + 1, n)),
                              gid(a, ring_local(m, n + 1)), gid(a + 1, ring_local(m, n + 1)), gid(a + 1, ring_local(m + 1, n + 1)), gid(a, ring_local(m + 1, n + 1))])

    # ---- boundary quads ----
    inlet, outlet, wall, sym_y0, sym_z0 = [], [], [], [], []
    for i in range(Nc):
        for j in range(Nc):
            inlet.append([gid(0, cidx(i, j)), gid(0, cidx(i + 1, j)), gid(0, cidx(i + 1, j + 1)), gid(0, cidx(i, j + 1))])
            outlet.append([gid(Nx, cidx(i, j)), gid(Nx, cidx(i, j + 1)), gid(Nx, cidx(i + 1, j + 1)), gid(Nx, cidx(i + 1, j))])
    for m in range(Mcirc):
        for n in range(Nrad):
            inlet.append([gid(0, ring_local(m, n)), gid(0, ring_local(m + 1, n)), gid(0, ring_local(m + 1, n + 1)), gid(0, ring_local(m, n + 1))])
            outlet.append([gid(Nx, ring_local(m, n)), gid(Nx, ring_local(m, n + 1)), gid(Nx, ring_local(m + 1, n + 1)), gid(Nx, ring_local(m + 1, n))])
    for a in range(Nx):
        for m in range(Mcirc):     # wall = ring outer (n=Nrad)
            wall.append([gid(a, ring_local(m, Nrad)), gid(a + 1, ring_local(m, Nrad)), gid(a + 1, ring_local(m + 1, Nrad)), gid(a, ring_local(m + 1, Nrad))])
        for i in range(Nc):        # sym y0 (u=0, core i=0 edge)
            sym_y0.append([gid(a, cidx(0, i)), gid(a, cidx(0, i + 1)), gid(a + 1, cidx(0, i + 1)), gid(a + 1, cidx(0, i))])
        for n in range(Nrad):      # sym y0 (ring m=Mcirc)
            sym_y0.append([gid(a, ring_local(Mcirc, n)), gid(a, ring_local(Mcirc, n + 1)), gid(a + 1, ring_local(Mcirc, n + 1)), gid(a + 1, ring_local(Mcirc, n))])
        for i in range(Nc):        # sym z0 (v=0, core j=0 edge)
            sym_z0.append([gid(a, cidx(i, 0)), gid(a + 1, cidx(i, 0)), gid(a + 1, cidx(i + 1, 0)), gid(a, cidx(i + 1, 0))])
        for n in range(Nrad):      # sym z0 (ring m=0)
            sym_z0.append([gid(a, ring_local(0, n)), gid(a + 1, ring_local(0, n)), gid(a + 1, ring_local(0, n + 1)), gid(a, ring_local(0, n + 1))])

    groups = [(11, 1, "inlet", inlet), (12, 2, "outlet", outlet), (13, 3, "wall", wall),
              (14, 4, "sym_y0", sym_y0), (15, 5, "sym_z0", sym_z0)]
    nquad = sum(len(g[3]) for g in groups)

    def bbox(idl):
        p = coords[np.array(idl).ravel() - 1]; return p.min(0), p.max(0)

    with open(args.out, "w") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n6\n2 1 \"inlet\"\n2 2 \"outlet\"\n2 3 \"wall\"\n2 4 \"sym_y0\"\n2 5 \"sym_z0\"\n3 8 \"fluid\"\n$EndPhysicalNames\n")
        f.write("$Entities\n0 0 5 1\n")
        for ent, phys, name, quads in groups:
            lo, hi = bbox(quads); f.write(f"{ent} {lo[0]:.9g} {lo[1]:.9g} {lo[2]:.9g} {hi[0]:.9g} {hi[1]:.9g} {hi[2]:.9g} 1 {phys} 0\n")
        lo, hi = coords.min(0), coords.max(0)
        f.write(f"1 {lo[0]:.9g} {lo[1]:.9g} {lo[2]:.9g} {hi[0]:.9g} {hi[1]:.9g} {hi[2]:.9g} 1 8 0\n$EndEntities\n")
        f.write(f"$Nodes\n1 {Nnode} 1 {Nnode}\n3 1 0 {Nnode}\n")
        for t in range(1, Nnode + 1): f.write(f"{t}\n")
        for c in coords: f.write(f"{c[0]:.9g} {c[1]:.9g} {c[2]:.9g}\n")
        f.write("$EndNodes\n")
        nele = nquad + len(hexes)
        f.write(f"$Elements\n{len(groups)+1} {nele} 1 {nele}\n")
        eid = 1
        for ent, phys, name, quads in groups:
            f.write(f"2 {ent} 3 {len(quads)}\n")
            for q in quads: f.write(f"{eid} {q[0]} {q[1]} {q[2]} {q[3]}\n"); eid += 1
        f.write(f"3 1 5 {len(hexes)}\n")
        for h in hexes: f.write(f"{eid} " + " ".join(map(str, h)) + "\n"); eid += 1
        f.write("$EndElements\n")
    print(f"wrote {args.out}: {Nnode} nodes, {len(hexes)} hexes, {nquad} bquads | Nx={Nx} Nc={Nc} Mcirc={Mcirc} Nrad={Nrad}")
    print(f"  wall first-layer dwall={args.dwall}mm (uniform circumferentially), r_inner=fcore*R, fcore={args.fcore}")


if __name__ == "__main__":
    main()
