#!/usr/bin/env python3
r"""ノズルの 90° セクタ 3D 構造ヘックスメッシュを生成 (Gmsh 4.1 .msh 直接出力).

目的: 軸中心 k の過大が「軸対称定式化のアーティファクトか / モデル(ひずみ生産)由来か」を
切り分けるため、同じノズルを **真の 3D** (isAxisymmetric=0) で解く。軸特異(ポーラ崩れ)を
避けるため断面は **非特異 squircle 写像** を使う:

  (u,v) in [0,1]^2  ->  y = R(x)*u*sqrt(1 - v^2/2),  z = R(x)*v*sqrt(1 - u^2/2)

これは原点(=軸)で非特異(ヤコビアン>0)、u=1 と v=1 が円弧(壁)に乗る。軸はメッシュの
正規な「辺」になり、スロートでもセルが潰れない。

境界 physID: inlet=1, outlet=2, wall=3(円弧), sym_y0=4(y=0面), sym_z0=5(z=0面)。fluid=8。
壁付近(u,v->1)へクラスタリング。座標は m (contour mm を 0.001 倍)。
"""
from __future__ import annotations
import argparse
import math
import numpy as np


def wall_stretch(n, b):
    """[0,1] を n+1 点へ、s=1(壁)側に密にクラスタ (sinh 引き伸ばし)。"""
    t = np.linspace(0.0, 1.0, n + 1)
    if b <= 1e-6:
        return t
    return 1.0 - np.sinh(b * (1.0 - t)) / math.sinh(b)


def axial_stations(x0, xth, x1, nconv, ndiv, b):
    """x0..xth(throat)..x1 を、throat 近傍に密にクラスタした軸ステーション列。"""
    # convergent: x0 -> 0, throat 側(end)を密に
    tc = np.linspace(0, 1, nconv + 1)
    sc = np.sinh(b * tc) / math.sinh(b)          # 0..1, end(=throat) 側密
    xc = x0 + (xth - x0) * sc
    # divergent: 0 -> x1, throat 側(start)を密に
    td = np.linspace(0, 1, ndiv + 1)
    sd = 1.0 - np.sinh(b * (1.0 - td)) / math.sinh(b)  # start(=throat) 側密
    xd = xth + (x1 - xth) * sd
    return np.concatenate([xc[:-1], xd])         # throat を 1 個だけ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contour", default="contour_conical.csv")
    ap.add_argument("--rt", type=float, default=10.0)
    ap.add_argument("--nr", type=int, default=32, help="断面 1 辺の分割")
    ap.add_argument("--nconv", type=int, default=50)
    ap.add_argument("--ndiv", type=int, default=90)
    ap.add_argument("--bwall", type=float, default=2.2, help="壁クラスタ強さ")
    ap.add_argument("--bax", type=float, default=2.0, help="軸方向 throat クラスタ強さ")
    ap.add_argument("--out", default="nozzle3d.msh")
    args = ap.parse_args()

    cont = np.loadtxt(args.contour, delimiter=",", skiprows=1, usecols=(0, 1))
    cx, cy = cont[:, 0], cont[:, 1]
    o = np.argsort(cx); cx, cy = cx[o], cy[o]
    x0, x1 = cx.min(), cx.max()

    xa = axial_stations(x0, 0.0, x1, args.nconv, args.ndiv, args.bax)  # mm
    Rw = np.interp(xa, cx, cy)                    # 壁半径(mm)
    Rw[np.argmin(np.abs(xa))] = args.rt          # throat を正確に Rt

    Nx = len(xa) - 1
    Nr = args.nr
    s = wall_stretch(Nr, args.bwall)             # [0,1], 壁側密
    npl = (Nr + 1) * (Nr + 1)                    # 断面ノード数

    def nid(a, i, j):
        return a * npl + i * (Nr + 1) + j + 1    # 1-based, 順序=位置

    # ---- nodes ----
    coords = np.empty(((Nx + 1) * npl, 3))
    for a in range(Nx + 1):
        R = Rw[a] * 1e-3
        x = xa[a] * 1e-3
        for i in range(Nr + 1):
            u = s[i]
            for j in range(Nr + 1):
                v = s[j]
                y = R * u * math.sqrt(1.0 - 0.5 * v * v)
                z = R * v * math.sqrt(1.0 - 0.5 * u * u)
                coords[nid(a, i, j) - 1] = (x, y, z)
    Nnode = coords.shape[0]

    # ---- hexes (Gmsh order) ----
    hexes = []
    for a in range(Nx):
        for i in range(Nr):
            for j in range(Nr):
                n = [nid(a, i, j), nid(a + 1, i, j), nid(a + 1, i + 1, j), nid(a, i + 1, j),
                     nid(a, i, j + 1), nid(a + 1, i, j + 1), nid(a + 1, i + 1, j + 1), nid(a, i + 1, j + 1)]
                hexes.append(n)

    # ---- boundary quads per BC ----
    inlet, outlet, wall, sym_y0, sym_z0 = [], [], [], [], []
    for i in range(Nr):
        for j in range(Nr):
            inlet.append([nid(0, i, j), nid(0, i + 1, j), nid(0, i + 1, j + 1), nid(0, i, j + 1)])
            outlet.append([nid(Nx, i, j), nid(Nx, i, j + 1), nid(Nx, i + 1, j + 1), nid(Nx, i + 1, j)])
    for a in range(Nx):
        for j in range(Nr):  # i=Nr face (arc, wall)
            wall.append([nid(a, Nr, j), nid(a + 1, Nr, j), nid(a + 1, Nr, j + 1), nid(a, Nr, j + 1)])
        for i in range(Nr):  # j=Nr face (arc, wall)
            wall.append([nid(a, i, Nr), nid(a, i + 1, Nr), nid(a + 1, i + 1, Nr), nid(a + 1, i, Nr)])
        for j in range(Nr):  # i=0 face (y=0 symmetry)
            sym_y0.append([nid(a, 0, j), nid(a, 0, j + 1), nid(a + 1, 0, j + 1), nid(a + 1, 0, j)])
        for i in range(Nr):  # j=0 face (z=0 symmetry)
            sym_z0.append([nid(a, i, 0), nid(a + 1, i, 0), nid(a + 1, i + 1, 0), nid(a, i + 1, 0)])

    surf_groups = [(11, 1, "inlet", inlet), (12, 2, "outlet", outlet), (13, 3, "wall", wall),
                   (14, 4, "sym_y0", sym_y0), (15, 5, "sym_z0", sym_z0)]
    nquad = sum(len(g[3]) for g in surf_groups)

    # bbox helper
    def bbox(idlist):
        pts = coords[np.array(idlist).ravel() - 1]
        return pts.min(0), pts.max(0)

    # ---- write Gmsh 4.1 ----
    with open(args.out, "w") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n6\n")
        f.write('2 1 "inlet"\n2 2 "outlet"\n2 3 "wall"\n2 4 "sym_y0"\n2 5 "sym_z0"\n3 8 "fluid"\n')
        f.write("$EndPhysicalNames\n")
        # entities: 0 points, 0 curves, 5 surfaces, 1 volume
        f.write("$Entities\n0 0 5 1\n")
        for ent, phys, name, quads in surf_groups:
            lo, hi = bbox(quads)
            f.write(f"{ent} {lo[0]:.9g} {lo[1]:.9g} {lo[2]:.9g} {hi[0]:.9g} {hi[1]:.9g} {hi[2]:.9g} 1 {phys} 0\n")
        lo, hi = coords.min(0), coords.max(0)
        f.write(f"1 {lo[0]:.9g} {lo[1]:.9g} {lo[2]:.9g} {hi[0]:.9g} {hi[1]:.9g} {hi[2]:.9g} 1 8 0\n")
        f.write("$EndEntities\n")
        # nodes (single block, dim=3 entTag=1, tags 1..N in order)
        f.write("$Nodes\n")
        f.write(f"1 {Nnode} 1 {Nnode}\n")
        f.write(f"3 1 0 {Nnode}\n")
        for t in range(1, Nnode + 1):
            f.write(f"{t}\n")
        for c in coords:
            f.write(f"{c[0]:.9g} {c[1]:.9g} {c[2]:.9g}\n")
        f.write("$EndNodes\n")
        # elements
        nele = nquad + len(hexes)
        f.write("$Elements\n")
        f.write(f"{len(surf_groups) + 1} {nele} 1 {nele}\n")
        eid = 1
        for ent, phys, name, quads in surf_groups:
            f.write(f"2 {ent} 3 {len(quads)}\n")
            for q in quads:
                f.write(f"{eid} {q[0]} {q[1]} {q[2]} {q[3]}\n"); eid += 1
        f.write(f"3 1 5 {len(hexes)}\n")
        for h in hexes:
            f.write(f"{eid} " + " ".join(map(str, h)) + "\n"); eid += 1
        f.write("$EndElements\n")

    print(f"wrote {args.out}: {Nnode} nodes, {len(hexes)} hexes, {nquad} boundary quads")
    print(f"  Nx={Nx} Nr={Nr}  x=[{x0:.2f},{x1:.2f}]mm  throat R={args.rt}mm")
    # quick quality: min/max cell edge near throat
    ath = int(np.argmin(np.abs(xa)))
    print(f"  throat station a={ath} x={xa[ath]:.3f}mm R={Rw[ath]:.3f}mm")


if __name__ == "__main__":
    main()
