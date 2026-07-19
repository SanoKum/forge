#!/usr/bin/env python3
"""wavy.h5 (forge cell メッシュ) から同一節点・同一 hex 結線の gmsh .msh (4.1) を再生成する。

目的: node (median-dual) 変換は convertGmshToForge + solverConfig `discretization: node` 経由が
正規ルートだが、wavy は合成 Fluent 由来で .msh が存在しない (生成スクリプトも散逸済み)。
そこで forge h5 の MESH/COORD (11^3 節点) と MESH/CONNE (XDMF hex: [9, n0..n7] x 1000) から
.msh を復元し、cell 検証 (run_0005/0011) と**厳密に同一の節点座標**で node 検証を可能にする。

境界 quad は「4 節点が同一境界平面上にある hex 面」から抽出し、physID は mapping.yaml と同じ
1:xmin 2:xmax 3:ymin 4:ymax 5:zmin 6:zmax を割り当てる。

使い方: python3 make_wavy_node_msh.py [wavy.h5] [wavy_node.msh]
"""
import sys

import h5py
import numpy as np


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "wavy.h5"
    dst = sys.argv[2] if len(sys.argv) > 2 else "wavy_node.msh"

    with h5py.File(src, "r") as f:
        co = f["MESH/COORD"][:].reshape(-1, 3).astype(np.float64)
        cn = f["MESH/CONNE"][:].reshape(-1, 9)
    assert (cn[:, 0] == 9).all(), "XDMF hex (type 9) 前提"
    hexes = cn[:, 1:]
    nN = co.shape[0]
    nC = hexes.shape[0]

    # 境界 quad 抽出: hex の 6 面 (VTK/gmsh hex の局所面定義) で全 4 節点が境界平面上のもの
    faces_local = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                   (3, 2, 6, 7), (0, 3, 7, 4), (1, 2, 6, 5)]
    lo, hi = co.min(0), co.max(0)
    tol = 1e-9
    onb = []  # (axis, side) ごとの節点マスク
    for ax in range(3):
        onb.append(np.abs(co[:, ax] - lo[ax]) < tol)
        onb.append(np.abs(co[:, ax] - hi[ax]) < tol)
    quads = {p: [] for p in range(1, 7)}  # physID -> [quad nodes]
    for h in hexes:
        for fl in faces_local:
            q = h[list(fl)]
            for p in range(6):
                if onb[p][q].all():
                    # physID: 1:xmin 2:xmax 3:ymin 4:ymax 5:zmin 6:zmax = (ax*2+side)+1
                    quads[p + 1].append(q)
                    break
    counts = {p: len(v) for p, v in quads.items()}
    print("boundary quad counts:", counts)
    assert all(c == counts[1] for c in counts.values()), "各面の quad 数が不一致"

    nQ = sum(counts.values())
    names = {1: "xmin", 2: "xmax", 3: "ymin", 4: "ymax", 5: "zmin", 6: "zmax"}
    with open(dst, "w") as w:
        w.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        w.write("$PhysicalNames\n7\n")
        for p in range(1, 7):
            w.write(f'2 {p} "{names[p]}"\n')
        w.write('3 7 "fluid"\n$EndPhysicalNames\n')
        # Entities: 点/線 0、面 6 (phys 1-6)、体 1 (phys 7)。bbox は参考値 (reader は先頭 9 欄のみ読む)
        w.write("$Entities\n0 0 6 1\n")
        for p in range(1, 7):
            w.write(f"{p} {lo[0]} {lo[1]} {lo[2]} {hi[0]} {hi[1]} {hi[2]} 1 {p} 0\n")
        w.write(f"7 {lo[0]} {lo[1]} {lo[2]} {hi[0]} {hi[1]} {hi[2]} 1 7 0\n$EndEntities\n")
        # Nodes: 1 ブロック (dim3, entTag7)。reader は tag を読み飛ばし出現順に格納する
        w.write(f"$Nodes\n1 {nN} 1 {nN}\n")
        w.write(f"3 7 0 {nN}\n")
        for i in range(nN):
            w.write(f"{i+1}\n")
        for i in range(nN):
            w.write(f"{co[i,0]:.17g} {co[i,1]:.17g} {co[i,2]:.17g}\n")
        w.write("$EndNodes\n")
        # Elements: 面 6 ブロック (quad=type3) + 体 1 ブロック (hex=type5)
        w.write(f"$Elements\n7 {nQ + nC} 1 {nQ + nC}\n")
        eid = 1
        for p in range(1, 7):
            w.write(f"2 {p} 3 {counts[p]}\n")
            for q in quads[p]:
                w.write(f"{eid} {q[0]+1} {q[1]+1} {q[2]+1} {q[3]+1}\n")
                eid += 1
        w.write(f"3 7 5 {nC}\n")
        for h in hexes:
            w.write(f"{eid} " + " ".join(str(n + 1) for n in h) + "\n")
            eid += 1
        w.write("$EndElements\n")
    print(f"wrote {dst}: {nN} nodes, {nC} hex, {nQ} boundary quads")

    # periodic 妥当性: xmin/xmax の (y,z) 節点集合が厳密一致するか
    yz0 = np.sort(co[onb[0]][:, 1:].round(12).view("f8,f8"), order=["f0", "f1"], axis=0)
    yz1 = np.sort(co[onb[1]][:, 1:].round(12).view("f8,f8"), order=["f0", "f1"], axis=0)
    print("xmin/xmax (y,z) node match:", bool((yz0 == yz1).all()))


if __name__ == "__main__":
    main()
