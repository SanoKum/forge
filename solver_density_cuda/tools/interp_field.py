#!/usr/bin/env python3
"""過去 h5 (任意メッシュ) の場を最近傍interpolateで新メッシュ入力 h5 に貼る (cross-mesh restart)。

メッシュを変えた (quad↔tri↔構造, 解像度変更) ときに、収束済みの場を初期値として移植する
標準ツール。**メッシュ変更後は uniform 初期値から始めず本ツールで restart する**
(超音速/衝撃波/SST は uniform IC から発散しやすい)。同一メッシュの restart は restart_field.py。

- SRC: res_*.h5 (primitives P,T,Ux,Uy,Uz,ro[,k,omega]) でも input h5 (conserved) でも可。
- DST: convertGmshToForge 直後の新メッシュ input h5。
- 各 DST セル重心に最も近い SRC セル重心の値を採用 (scipy cKDTree, 最近傍)。
- 移植する /VALUE/: 保存量 ro,roUx,roUy,roUz,roe・乱流 roK,roOmega・スカラー輸送 roY* (あれば)。
- **wall_dist は移植しない** (新メッシュで convert 時に計算済みの値を使う)。

usage: interp_field.py SRC.h5 DST_input.h5 [--gamma 1.4]
"""
import argparse, sys, os
import numpy as np, h5py
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from res_h5_to_vtu import parse_conne


def centroids(f):
    # 入力 h5 は /CELLS/centCoords を持つのでそれを優先 (node/median-dual の
    # 双対 CONNE は parse_conne で解析できないため必須。res_*.h5 は従来経路)
    if "CELLS/centCoords" in f:
        return np.array(f["CELLS/centCoords"]).reshape(-1, 3)[:, :2]
    coord = np.array(f["MESH/COORD"]).reshape(-1, 3)[:, :2]
    nc = f["VALUE/ro"].shape[0]
    if nc == coord.shape[0]:
        # node-centered res (median-dual): 値の位置はノード座標そのもの
        # (MESH/CONNE は可視化用 primal トポロジで parse できない)
        return coord
    conn, offs, _ = parse_conne(np.array(f["MESH/CONNE"]), nc)
    c = np.zeros((nc, 2)); s = 0
    for i, o in enumerate(offs):
        c[i] = coord[conn[s:o]].mean(axis=0); s = o
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("dst")
    ap.add_argument("--gamma", type=float, default=1.4)
    a = ap.parse_args(); g = a.gamma

    with h5py.File(a.src, "r") as s:
        cs = centroids(s); V = s["VALUE"]
        if "P" in V and "Ux" in V:            # res (primitives)
            ro = np.array(V["ro"]); P = np.array(V["P"])
            Ux = np.array(V["Ux"]); Uy = np.array(V["Uy"]); Uz = np.array(V["Uz"])
            fields = {"ro": ro, "roUx": ro*Ux, "roUy": ro*Uy, "roUz": ro*Uz,
                      "roe": P/(g-1.0) + 0.5*ro*(Ux**2+Uy**2+Uz**2)}
            if "k" in V and "omega" in V:
                fields["roK"] = ro*np.array(V["k"]); fields["roOmega"] = ro*np.array(V["omega"])
            for key in V:                      # scalar transport Y* -> roY*
                if key.startswith("Y") and key[1:].isdigit():
                    fields["ro"+key] = ro*np.array(V[key])
        else:                                  # input (conserved)
            fields = {n: np.array(V[n]) for n in
                      ["ro","roUx","roUy","roUz","roe","roK","roOmega"] if n in V}
            for key in V:
                if key.startswith("roY") and key[3:].isdigit():
                    fields[key] = np.array(V[key])

    tree = cKDTree(cs)
    with h5py.File(a.dst, "r+") as d:
        cd = centroids(d)
        _, idx = tree.query(cd)
        moved = []
        for name, arr in fields.items():
            ds = "VALUE/"+name
            if ds in d and name != "wall_dist":
                d[ds][...] = arr[idx].astype(d[ds].dtype); moved.append(name)
        print(f"interp {a.src} -> {a.dst}: {len(cd)} dst cells, moved {moved} (wall_dist kept)")


if __name__ == "__main__":
    main()
