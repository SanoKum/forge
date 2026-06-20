#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Salome 出力 (.med / .unv) を forge 取り込み用の gmsh msh4 (.msh, version 4.1) に変換する。

注意: forge は msh4 (4.1) を読む。msh2 は不可 (case/36 README 参照)。よって 4.1 で書き出す。

使い方 (Salome 不要・通常の python + meshio):
    python3 med_to_gmsh.py passive_porous.med
    python3 med_to_gmsh.py passive_solid.unv   passive_solid.msh

forge 規約に合わせて以下を保証する:
  * Physical Curve(dim=1): inlet=1, outlet=2, wall_slip=3, wall=4
  * Physical Surface(dim=2): fluid=10
    -> 同ディレクトリ gmsh 版 make_mesh.py と physID 一致 (同じ bcondConfig.yaml が使える)。
  * どの境界グループにも属さない内部線要素 (スロット口・ブロック界面など) は
    **書き出さない** (forge は dim=1 要素を全て境界として読むため、physID=0 の
    内部線が残ると 'couldn't find physID 0' で落ちる)。
  * 2D 面要素は quad(=3)/tri(=2)、境界は line(=1)、z=0。

依存: pip install meshio h5py   (MED 読み込みに h5py が要る)
"""
import sys
import os
import numpy as np
import meshio

# forge と突き合わせる physID (gmsh 版 make_mesh.py と一致させること)
BND_ID = {"inlet": 1, "outlet": 2, "wall_slip": 3, "wall": 4}   # dim=1
SURF_ID = {"fluid": 10}                                          # dim=2
DIM = dict({k: 1 for k in BND_ID}, **{k: 2 for k in SURF_ID})
ALL_ID = dict(BND_ID, **SURF_ID)

LINE_TYPES = {"line", "line3"}
SURF_TYPES = {"triangle", "triangle6", "quad", "quad8", "quad9"}


def norm(name):
    return name.strip().lower().lstrip("_").replace(" ", "")


def build_gmsh_mesh(m):
    """meshio で読んだメッシュ (cell_sets にグループを持つ) を、
    forge 用 physID を付与し内部線を除いた gmsh 出力用 meshio.Mesh に変換する。"""
    # --- 各 cell ブロックの各要素に physID を割り当てる ---
    # cell_sets[name] = [idx_block0, idx_block1, ...] (各ブロック内の要素 index)
    nblk = len(m.cells)
    phys = [np.zeros(len(cb.data), dtype=np.int64) for cb in m.cells]

    # グループ名 -> physID をマッチ (大文字小文字/前置 _ を吸収)
    matched = {}
    for gname, sets in m.cell_sets.items():
        key = norm(gname)
        if key not in ALL_ID:
            continue
        pid = ALL_ID[key]
        matched[key] = pid
        for ib in range(nblk):
            idx = sets[ib] if ib < len(sets) else None
            if idx is None:
                continue
            idx = np.asarray(idx, dtype=np.int64)
            if idx.size:
                phys[ib][idx] = pid

    # fluid グループが cell_sets に無い場合の保険: 全 surface 要素を fluid=10 に
    if "fluid" not in matched:
        print("  [warn] 'fluid' group not found in cell_sets; tagging all surface cells as fluid=10")
        for ib, cb in enumerate(m.cells):
            if cb.type in SURF_TYPES:
                phys[ib][:] = SURF_ID["fluid"]
        matched["fluid"] = SURF_ID["fluid"]

    # --- 内部線要素 (physID==0 の line) を除外して再構築 ---
    new_cells, new_phys = [], []
    dropped_lines = 0
    for ib, cb in enumerate(m.cells):
        p = phys[ib]
        if cb.type in LINE_TYPES:
            keep = p > 0
            dropped_lines += int(np.count_nonzero(~keep))
            if np.count_nonzero(keep) == 0:
                continue
            new_cells.append(meshio.CellBlock(cb.type, cb.data[keep]))
            new_phys.append(p[keep])
        elif cb.type in SURF_TYPES:
            # 念のため未タグの surface も fluid に
            p = np.where(p > 0, p, SURF_ID["fluid"])
            new_cells.append(meshio.CellBlock(cb.type, cb.data))
            new_phys.append(p)
        else:
            # 他次元 (vertex 等) は捨てる
            continue

    field_data = {name: np.array([pid, DIM[name]]) for name, pid in matched.items()}

    out = meshio.Mesh(
        points=m.points,
        cells=new_cells,
        cell_data={"gmsh:physical": new_phys,
                   "gmsh:geometrical": [pp.copy() for pp in new_phys]},
        field_data=field_data,
    )
    return out, dropped_lines


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(src)[0] + ".msh"

    m = meshio.read(src)
    print(f"[read] {src}")
    print(f"  cell blocks : {[ (cb.type, len(cb.data)) for cb in m.cells ]}")
    print(f"  cell_sets   : {list(m.cell_sets.keys())}")

    out, dropped_lines = build_gmsh_mesh(m)

    # meshio の gmsh4.1 ライタは多セルタイプ時に node 単位の gmsh:dim_tags を要求するため、
    # まず msh2 で書き出し、インストール済み gmsh で msh4.1 に変換する (物理グループ保持)。
    import subprocess, shutil
    tmp2 = os.path.splitext(dst)[0] + ".v2.msh"
    meshio.write(tmp2, out, file_format="gmsh22", binary=False)
    print(f"[write] {tmp2} (msh2, dropped {dropped_lines} internal line elems)")

    gmsh_bin = shutil.which("gmsh")
    if gmsh_bin:
        subprocess.run([gmsh_bin, tmp2, "-save", "-format", "msh4", "-o", dst],
                       check=True, stdout=subprocess.DEVNULL)
        os.remove(tmp2)
        print(f"[write] {dst} (msh4.1; forge はこちらを読む)")
    else:
        os.replace(tmp2, dst)
        print("  [WARN] gmsh が見つからないため msh2 のまま出力。forge は msh4 必須:")
        print("         gmsh %s -save -format msh4 -o %s" % (dst, dst))

    verify(dst)


def verify(path):
    """書き出した .msh を読み直して physical group を確認。"""
    m = meshio.read(path)
    print("-" * 60)
    print(f"[verify] {path}")
    print(f"  physical names: {m.field_data}")
    # physID ごとの要素数
    counts = {}
    for cb, pdat in zip(m.cells, m.cell_data.get("gmsh:physical", [])):
        for pid in np.unique(pdat):
            counts[(cb.type, int(pid))] = int(np.count_nonzero(pdat == pid))
    for (ctype, pid), n in sorted(counts.items()):
        names = [k for k, v in m.field_data.items() if int(v[0]) == pid]
        print(f"    {ctype:9s} physID={pid:3d} ({','.join(names) or '?'}): {n}")
    # 健全性チェック
    line_pids = set()
    surf_pids = set()
    for cb, pdat in zip(m.cells, m.cell_data.get("gmsh:physical", [])):
        if cb.type in LINE_TYPES:
            line_pids |= set(int(x) for x in np.unique(pdat))
        elif cb.type in SURF_TYPES:
            surf_pids |= set(int(x) for x in np.unique(pdat))
    ok = True
    if 0 in line_pids:
        print("  [FAIL] physID=0 の境界線が残っている"); ok = False
    for name, pid in BND_ID.items():
        if pid not in line_pids:
            print(f"  [WARN] 境界 '{name}' (physID={pid}) が見つからない")
    if SURF_ID["fluid"] not in surf_pids:
        print("  [FAIL] fluid surface が無い"); ok = False
    if len(surf_pids - {0}) != 1:
        print(f"  [WARN] fluid surface physID が複数: {surf_pids}")
    print(f"  => {'OK' if ok else 'CHECK FAILED'}")


if __name__ == "__main__":
    main()
