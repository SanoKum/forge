#!/usr/bin/env python3
"""Fluent (CFF HDF5) メッシュ -> forge 入力 HDF5 への変換ツール。

世の中で配布されている ANSYS Fluent の HDF5 メッシュ (`*.msh.h5` / `*.cas.h5`,
いわゆる CFF = Common Fluids Format) を読み込み、forge の density solver が
`mesh.readMesh()` / `variables.readValueHDF5()` で読める HDF5 (+ XDMF) に変換する。

forge の内部メッシュは「面 (plane) ベース」で、Fluent も「面 (face) ベース」なので
両者はほぼ 1:1 で対応する。本ツールは Gmsh を経由せず、Fluent の面情報から
直接 forge HDF5 を書き出す。

----------------------------------------------------------------------------
forge 側の前提 (mesh/gmshReader.hpp::writeInputH5, mesh/mesh.cpp::readMesh より)
----------------------------------------------------------------------------
- /MESH 属性 : nNodes, nCells, nPlanes, nNormalPlanes, nBPlanes, nBconds (int32)
- /MESH/COORD  : [x,y,z]*nNodes (float)
- /MESH/CONNE  : XDMF mixed topology 用 [typeCode, node...]*nCells (int32)
- /PLANES/STRUCT    : [nNodes, node..., nCells, cell...]*nPlanes  (内部面が先, 境界面が後)
- /PLANES/surfVect  : [sx,sy,sz]*nPlanes (大きさ=面積, 向きは iCells[0]->iCells[1])
- /PLANES/surfArea  : [area]*nPlanes
- /PLANES/centCoords: [x,y,z]*nPlanes
- /CELLS/STRUCT     : [nNodes, node..., nPlanes, plane..., nPlanesDir, dir..., ieleType]*nCells
- /CELLS/volume     : [vol]*nCells
- /CELLS/centCoords : [x,y,z]*nCells
- /CELLS/regionId   : [region]*nCells (int32)
- /VALUE/{ro,roUx,roUy,roUz,roe,wall_dist} (必須) , roK,roOmega (任意)
- /BCONDS/<physID>/{iPlanes,iBPlanes,iCells} (int32) + 属性 bcondKind(文字列)

境界条件は **physID (整数)** で実行時 bcondConfig.yaml と突き合わされる
(boundaryCond.cpp::readBcondConfig)。H5 の bcondKind 属性は yaml により上書き
されるため、本ツールでは「Fluent ゾーン名 -> physID」の対応表さえあればよい。

----------------------------------------------------------------------------
Fluent CFF レイアウトについて (重要)
----------------------------------------------------------------------------
CFF の HDF5 内部レイアウトは Fluent のバージョンで差がある。本ツールは代表的な
レイアウトを自動探索するが、合わない場合は `dump` サブコマンドで構造を確認し、
`--nodes-path` 等で明示できる。ゾーン名が取れない場合でも、`dump` が出すゾーン
ID で対応表 (mapping) を書けば変換できる。

使い方:
    # 1) まず構造とゾーン一覧を確認
    python3 fluent_h5_to_forge.py dump  mesh.cas.h5

    # 2) ゾーン -> physID 対応表と bcondConfig.yaml の雛型を出力
    python3 fluent_h5_to_forge.py template mesh.cas.h5 -o mapping.yaml --bcond bcondConfig.yaml

    # 3) 変換 (mapping.yaml を編集後)
    python3 fluent_h5_to_forge.py convert mesh.cas.h5 mesh.h5 --map mapping.yaml \
            --ro 1.225 --p 101325 --ux 0 --bcond bcondConfig.yaml
"""

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("h5py が必要です: pip install h5py")


# --------------------------------------------------------------------------
# 要素タイプコード
# --------------------------------------------------------------------------
# Gmsh element type id (forge: cell.ieleType, elementType.hpp の mapElementFromGmshID)
GMSH_TRI, GMSH_QUAD, GMSH_TET, GMSH_HEX, GMSH_PRISM, GMSH_PYRAMID = 2, 3, 4, 5, 6, 7
# XDMF mixed topology code (forge: /MESH/CONNE 先頭, gmshReader.hpp の CONNE0)
XDMF_CODE = {GMSH_TRI: 4, GMSH_QUAD: 5, GMSH_TET: 6, GMSH_PYRAMID: 7,
             GMSH_PRISM: 8, GMSH_HEX: 9}
GMSH_NAME = {GMSH_TRI: "triangle", GMSH_QUAD: "quad", GMSH_TET: "tetra",
             GMSH_HEX: "hex", GMSH_PRISM: "prism", GMSH_PYRAMID: "pyramid"}


# ==========================================================================
# Fluent CFF HDF5 リーダ
# ==========================================================================
class FluentMesh:
    """Fluent CFF HDF5 から nodes / faces / cells / zones を抽出する。

    属性:
        dim          : 2 or 3
        coords       : (nNodes, 3) float64 (2D は z=0)
        nCells       : セル数
        cell_region  : (nCells,) int  各セルの region 番号 (cell zone 由来)
        faces        : list of dict {nodes:[..], c0:int, c1:int(or -1), zone:int}
                       (node/cell id はすべて 0-based)
        zones        : dict zone_id -> {name, ztype, is_boundary, nfaces}
    """

    def __init__(self, path, nodes_path=None, scale=1.0, verbose=True):
        self.path = path
        self.verbose = verbose
        self.scale = scale
        self.f = h5py.File(path, "r")
        self.mesh_grp = self._find_mesh_group()
        self._check_units()
        self._read_nodes(nodes_path)
        self._read_cells()
        self._read_faces()

    def _check_units(self):
        u = self.mesh_grp.attrs.get("units")
        if u is not None:
            u = u.decode() if isinstance(u, bytes) else str(u)
            if u not in ("m", "meter", "meters") and self.scale == 1.0:
                self._log(f"  注意: メッシュ単位は '{u}'。--scale で SI(m) へ換算してください "
                          f"(例 in->m: --scale 0.0254, mm->m: --scale 0.001)")

    def close(self):
        self.f.close()

    def _log(self, *a):
        if self.verbose:
            print(*a)

    # --- mesh group 探索 ---
    def _find_mesh_group(self):
        if "meshes" in self.f:
            g = self.f["meshes"]
            keys = list(g.keys())
            if not keys:
                sys.exit("/meshes 配下が空です。dump で構造を確認してください。")
            self._log(f"mesh group: /meshes/{keys[0]} (ほか {len(keys)-1} 個)")
            return g[keys[0]]
        # フォールバック: nodes / faces を直接含むグループを探す
        for name in ("Mesh", "mesh", "MESH"):
            if name in self.f and "nodes" in self.f[name]:
                return self.f[name]
        sys.exit("mesh group (/meshes/<id>) が見つかりません。dump で構造を確認してください。")

    # --- nodes ---
    def _read_nodes(self, nodes_path):
        if nodes_path:
            ds = self.f[nodes_path]
            coords = np.asarray(ds)
        else:
            grp = self._require(self.mesh_grp, "nodes")
            cgrp = self._require(grp, "coords")

            # coords 配下のデータセットを minId (= 大域 node id 開始) 順に連結
            def _start(k):
                a = cgrp[k]
                if "minId" in a.attrs:
                    return int(np.array(a.attrs["minId"]).ravel()[0])
                return _intkey(k)

            keys = [k for k in cgrp.keys() if isinstance(cgrp[k], h5py.Dataset)]
            parts = []
            for key in sorted(keys, key=_start):
                arr = np.asarray(cgrp[key])
                parts.append(arr.reshape(-1, 1) if arr.ndim == 1 else arr)
            coords = np.vstack(parts) if len(parts) > 1 else parts[0]
        coords = np.atleast_2d(coords)
        if coords.shape[1] == 1 and coords.shape[0] % 3 == 0:
            coords = coords.reshape(-1, 3)
        self.dim = coords.shape[1]
        if self.dim == 2:
            coords = np.column_stack([coords, np.zeros(len(coords))])
        elif self.dim != 3:
            sys.exit(f"node 座標の次元が不正です: shape={coords.shape}")
        self.coords = coords.astype(np.float64) * float(self.scale)
        self._log(f"nodes: {len(self.coords)}  (dim={self.dim}, scale={self.scale})")

    # --- cells ---
    def _read_cells(self):
        ncells = 0
        if "cellCount" in self.mesh_grp.attrs:
            ncells = int(np.array(self.mesh_grp.attrs["cellCount"]).ravel()[0])

        grp = self.mesh_grp.get("cells")
        zt = []
        if grp is not None and "zoneTopology" in grp:
            zt = _read_zone_topology(grp["zoneTopology"])
            for z in zt:
                if "maxId" in z:
                    ncells = max(ncells, int(z["maxId"]))
        # ctype からセル総数を補完 (cells/ctype/<sec>/cell-types)
        if ncells == 0 and grp is not None and "ctype" in grp:
            tot = 0
            for k in grp["ctype"].keys():
                sec = grp["ctype"][k]
                ds = sec["cell-types"] if isinstance(sec, h5py.Group) else sec
                tot += np.asarray(ds).size
            ncells = tot
        if ncells == 0:
            sys.exit("cell 数を決定できません。dump で /cells を確認してください。")

        self.nCells = int(ncells)
        self.cell_region = np.zeros(self.nCells, dtype=np.int32)
        self.region_names = {}
        for region_idx, z in enumerate(zt):
            if "minId" in z and "maxId" in z:
                lo, hi = int(z["minId"]), int(z["maxId"])
                self.cell_region[lo - 1:hi] = region_idx
                self.region_names[region_idx] = z.get("name", f"region_{region_idx}")
        self._log(f"cells: {self.nCells}  (regions={len(set(self.cell_region.tolist()))})")

    # --- faces ---
    def _read_faces(self):
        grp = self._require(self.mesh_grp, "faces")
        zt = _read_zone_topology(grp["zoneTopology"]) if "zoneTopology" in grp else []
        # zone を 大域 face id 範囲 [minId, maxId] で引けるようにする
        zranges = [(int(z["minId"]), int(z["maxId"]), z)
                   for z in zt if "minId" in z and "maxId" in z]

        self.zones = {}
        self.faces = []

        nodes_grp = self._require(grp, "nodes")
        c0_grp = self._require(grp, "c0")
        c1_grp = grp.get("c1")

        sec_keys = sorted(nodes_grp.keys(), key=_intkey)
        for sk in sec_keys:
            face_nodes, counts = _read_face_nodes(nodes_grp[sk])
            c0ds = c0_grp[sk]
            c0 = np.asarray(c0ds).ravel().astype(np.int64)
            if c1_grp is not None and sk in c1_grp:
                c1 = np.asarray(c1_grp[sk]).ravel().astype(np.int64)
            else:
                c1 = np.zeros_like(c0)

            # この section の 大域 face id 開始 (zone 対応付け用)
            gstart = None
            if "minId" in c0ds.attrs:
                gstart = int(np.array(c0ds.attrs["minId"]).ravel()[0])

            zinfo = None
            if gstart is not None:
                for lo, hi, z in zranges:
                    if lo <= gstart <= hi:
                        zinfo = z
                        break
            if zinfo is None and len(zranges) == len(sec_keys):
                # フォールバック: minId 昇順で section と zone を順対応
                zinfo = sorted(zranges, key=lambda r: r[0])[sec_keys.index(sk)][2]

            zid = int(zinfo["id"]) if zinfo and "id" in zinfo else _intkey(sk)
            zname = zinfo.get("name") if zinfo else None
            ztype = zinfo.get("zoneType") if zinfo else None

            is_boundary = bool(np.all(c1 <= 0))
            z = self.zones.setdefault(zid, {
                "name": zname or f"zone_{zid}",
                "zoneType": ztype,
                "is_boundary": is_boundary,
                "nfaces": 0,
            })
            z["nfaces"] += len(c0)
            if not is_boundary:
                z["is_boundary"] = False

            off = 0
            for i, n in enumerate(counts):
                nd = face_nodes[off:off + n] - 1   # 1-based -> 0-based
                off += n
                a = int(c0[i]) - 1                 # 1-based -> 0-based
                b = int(c1[i]) - 1
                self.faces.append({
                    "nodes": nd.astype(np.int64),
                    "c0": a,
                    "c1": b if b >= 0 else -1,
                    "zone": zid,
                })
        self._log(f"faces: {len(self.faces)}  (zones={len(self.zones)})")

    # --- helper ---
    @staticmethod
    def _require(grp, name):
        if name not in grp:
            sys.exit(f"必須グループ '{name}' が {grp.name} 配下にありません。"
                     f" dump で構造を確認してください。")
        return grp[name]


# --------------------------------------------------------------------------
# CFF レイアウト補助
# --------------------------------------------------------------------------
def _intkey(k):
    try:
        return int(k)
    except (ValueError, TypeError):
        return k


def _read_face_nodes(ds_or_grp):
    """face-node 接続を (flat_nodes, counts) で返す。

    レイアウトの揺れに対応:
      - グループで {nodes, nnodes} を持つ
      - 2D データセット (nfaces, k)  (均一サイズ)
      - 1D データセット (ragged, 先頭に各 face の節点数を含む場合もある)
    """
    if isinstance(ds_or_grp, h5py.Group):
        nodes = np.asarray(ds_or_grp["nodes"]).ravel().astype(np.int64)
        nn = np.asarray(ds_or_grp["nnodes"]).ravel().astype(np.int64)
        return nodes, nn
    arr = np.asarray(ds_or_grp)
    if arr.ndim == 2:
        k = arr.shape[1]
        return arr.ravel().astype(np.int64), np.full(arr.shape[0], k, dtype=np.int64)
    # 1D: 均一サイズ想定が難しいので、そのまま 1 列とは扱えない -> エラー
    sys.exit("face nodes が 1D で nnodes が併設されていません。"
             " --nodes-path 等で構造を確認してください (dump 参照)。")


def _split_names(arr):
    """name データセットを per-zone のリストにする。

    CFF では複数ゾーンの名前が 1 つの文字列に ';' 連結されている
    (例: shape (1,) の b'wall;inlet;outlet')。1 ゾーン 1 要素の形式にも対応。
    """
    vals = []
    for x in np.atleast_1d(np.asarray(arr)).ravel().tolist():
        if isinstance(x, (bytes, np.bytes_)):
            x = x.decode("utf-8", "replace")
        vals.append(str(x))
    if len(vals) == 1 and ";" in vals[0]:
        return [s for s in vals[0].split(";")]
    return vals


def _read_zone_topology(node):
    """zoneTopology を [{id,minId,maxId,zoneType,name,...}, ...] にして返す。

    対応形式:
      - グループ形式 (各フィールドが個別データセット; CFF 標準)。name は ';' 連結文字列。
      - compound dtype 1 本
      - 列順が決まった 2D int 配列
    """
    # グループ形式 (各フィールドが個別データセット)
    if isinstance(node, h5py.Group):
        nz_attr = node.attrs.get("nZones")
        per = {}        # per-zone の数値配列
        names = None
        for k in node.keys():
            kl = k.lower()
            if kl == "name":
                names = _split_names(node[k])
            elif kl == "fields":
                continue  # フィールド名一覧であり zone データではない
            else:
                per[k] = np.asarray(node[k]).ravel()
        # ゾーン数を決定 (per-zone 配列長を優先)
        nz = 0
        for k in ("id", "minId", "minid", "zoneType", "zonetype", "maxId", "maxid"):
            if k in per:
                nz = len(per[k])
                break
        if nz == 0:
            nz = int(np.array(nz_attr).ravel()[0]) if nz_attr is not None else \
                 (len(names) if names else 0)
        out = []
        for i in range(nz):
            rec = {}
            for k, v in per.items():
                if i < len(v):
                    rec[k] = v[i]
            if names and i < len(names):
                rec["name"] = names[i]
            out.append(_normalize_zone_rec(rec))
        return out

    arr = np.asarray(node)
    # compound
    if arr.dtype.names:
        out = []
        for rec in arr:
            d = {name: rec[name] for name in arr.dtype.names}
            out.append(_normalize_zone_rec(d))
        return out
    # plain 2D int 配列: 列順は Fluent CFF の慣例 (id, childId, ... minId, maxId)
    # 確実でないため、最小限 id/minId/maxId を推定する。
    arr = np.atleast_2d(arr)
    out = []
    for row in arr:
        row = row.tolist()
        # 慣例: [id, childZoneId, ..., minId, maxId] のことが多い。
        # 末尾 2 つを (minId, maxId) とみなし、先頭を id とする。
        rec = {"id": int(row[0])}
        if len(row) >= 2:
            rec["minId"] = int(row[-2])
            rec["maxId"] = int(row[-1])
        out.append(_normalize_zone_rec(rec))
    return out


def _normalize_zone_rec(rec):
    out = {}
    for k, v in rec.items():
        kl = k.lower()
        if isinstance(v, bytes):
            v = v.decode("utf-8", "replace")
        if isinstance(v, np.generic):
            v = v.item()
        if kl in ("id", "zoneid"):
            out["id"] = int(v)
        elif kl in ("minid", "min_id", "min"):
            out["minId"] = int(v)
        elif kl in ("maxid", "max_id", "max"):
            out["maxId"] = int(v)
        elif kl in ("zonetype", "type", "bctype"):
            out["zoneType"] = int(v) if not isinstance(v, str) else v
        elif kl in ("name", "zonename"):
            out["name"] = str(v)
        else:
            out[k] = v
    out.setdefault("id", out.get("id", -1))
    return out


# Fluent zoneType -> forge bcond kind の推定マップ
FLUENT_ZTYPE_KIND = {
    3: "wall",
    4: "inlet_uniformVelocity",   # pressure-inlet
    5: "outflow",                 # pressure-outlet
    7: "slip",                    # symmetry
    8: "slip",                    # axis (軸対称軸; 必要なら別途調整)
    9: "inlet_uniformVelocity",   # pressure-far-field
    10: "inlet_uniformVelocity",  # velocity-inlet
    12: "periodic",
    14: "slip",                   # porous-jump など (暫定)
    20: "inlet_uniformVelocity",  # mass-flow-inlet
    36: "outflow",                # outflow
    37: "outflow",                # target-mass-flow
}


# ==========================================================================
# forge メッシュ構築
# ==========================================================================
class ForgeMesh:
    def __init__(self, fl: FluentMesh, triangulate=False, verbose=True):
        self.fl = fl
        self.verbose = verbose
        self.dim = fl.dim
        self.coords = fl.coords
        self.nNodes = len(fl.coords)
        self.nCells = fl.nCells
        self.triangulate = triangulate
        self._build()

    def _log(self, *a):
        if self.verbose:
            print(*a)

    def _build(self):
        fl = self.fl
        # cell -> faces 隣接
        cell_faces = [[] for _ in range(self.nCells)]
        for fi, fc in enumerate(fl.faces):
            cell_faces[fc["c0"]].append(fi)
            if fc["c1"] >= 0:
                cell_faces[fc["c1"]].append(fi)

        # --- cell の節点集合・要素タイプ・順序付き接続 (元の面で判定; viz/iNodes 用) ---
        self.cell_nodes = []   # 順序付き (可能なら) 節点
        self.cell_etype = []   # gmsh id
        for ic in range(self.nCells):
            faces = [fl.faces[fi]["nodes"] for fi in cell_faces[ic]]
            etype, ordered = _reconstruct_cell(faces, self.coords, self.dim)
            self.cell_etype.append(etype)
            self.cell_nodes.append(ordered)

        # --- 面リスト (任意で四角/多角形を三角形に分割) ---
        # 3D の非平面四角面 (Fluent の境界層 prism 側面など) は単一法線で表せず
        # forge のスキームが不安定化する。--triangulate で全面を三角形化し平面化する。
        if self.triangulate and self.dim == 3:
            wfaces = _triangulate_faces(fl.faces)
            self._log(f"  triangulate: faces {len(fl.faces)} -> {len(wfaces)}")
        else:
            wfaces = fl.faces
        self._wfaces = wfaces

        # --- 面 (plane) を内部->境界の順に並べ替え ---
        internal = [fi for fi, fc in enumerate(wfaces) if fc["c1"] >= 0]
        boundary = [fi for fi, fc in enumerate(wfaces) if fc["c1"] < 0]
        order = internal + boundary
        self.nNormalPlanes = len(internal)
        self.nBPlanes = len(boundary)
        self.nPlanes = len(order)

        # 旧 index -> 新 index
        new_of = np.empty(len(wfaces), dtype=np.int64)
        for new_i, old_i in enumerate(order):
            new_of[old_i] = new_i

        # plane 配列
        self.pl_nodes = []
        self.pl_cells = []
        self.surfVect = np.zeros((self.nPlanes, 3))
        self.surfArea = np.zeros(self.nPlanes)
        self.pl_cent = np.zeros((self.nPlanes, 3))

        # cell の centroid (node 平均, forge 準拠)
        self.cell_cent = np.zeros((self.nCells, 3))
        for ic in range(self.nCells):
            nd = self.cell_nodes[ic]
            self.cell_cent[ic] = self.coords[nd].mean(axis=0)

        # cell -> planes
        cell_planes = [[] for _ in range(self.nCells)]
        cell_pdir = [[] for _ in range(self.nCells)]
        self.cell_vol = np.zeros(self.nCells)

        for new_i, old_i in enumerate(order):
            fc = wfaces[old_i]
            nd = fc["nodes"]
            c0 = fc["c0"]
            c1 = fc["c1"]
            pts = self.coords[nd]
            sv = _area_vector(pts, self.dim)
            cent = pts.mean(axis=0)
            # 体積寄与 (sv と同じ未反転向きに対応する厳密モーメント)
            moment = _face_volume_moment(pts, self.dim)

            # 向きをそろえる: 内部は c0->c1, 境界は c0 から外向き
            if c1 >= 0:
                d = self.cell_cent[c1] - self.cell_cent[c0]
            else:
                d = cent - self.cell_cent[c0]
            if np.dot(d, sv) < 0.0:
                sv = -sv
                moment = -moment

            area = float(np.linalg.norm(sv))
            self.surfVect[new_i] = sv
            self.surfArea[new_i] = area
            self.pl_cent[new_i] = cent
            self.pl_nodes.append(nd)
            if c1 >= 0:
                self.pl_cells.append([c0, c1])
            else:
                self.pl_cells.append([c0])

            # iPlanesDir と体積 (発散定理, 面三角分割で厳密) を蓄積
            cell_planes[c0].append(new_i)
            cell_pdir[c0].append(1)
            self.cell_vol[c0] += moment            # owner: 外向き +sv
            if c1 >= 0:
                cell_planes[c1].append(new_i)
                cell_pdir[c1].append(-1)
                self.cell_vol[c1] -= moment        # neighbor: 外向き -sv

        self.cell_vol /= float(self.dim)
        self.cell_vol = np.abs(self.cell_vol)
        self.cell_planes = cell_planes
        self.cell_pdir = cell_pdir

        nbad = int(np.sum(self.cell_vol < 1e-20))
        if nbad:
            self._log(f"  WARNING: 体積がほぼ 0 のセルが {nbad} 個あります")

        # --- 境界条件 (zone -> plane) ---
        # zone id -> [新 plane index]
        self.bzone_planes = {}
        for old_i in boundary:
            zid = wfaces[old_i]["zone"]
            self.bzone_planes.setdefault(zid, []).append(int(new_of[old_i]))

        self._log(f"planes: total={self.nPlanes} internal={self.nNormalPlanes} "
                  f"boundary={self.nBPlanes}")
        etypes = {}
        for e in self.cell_etype:
            etypes[GMSH_NAME.get(e, f"id{e}")] = etypes.get(GMSH_NAME.get(e, f"id{e}"), 0) + 1
        self._log(f"  cell types: {etypes}")
        n_poly = sum(1 for e in self.cell_etype if e == 0)
        if n_poly:
            self._log(f"  注意: 標準要素でないセル(多面体等)が {n_poly} 個。"
                      f"体積/接続は発散定理で正しく計算され解は回るが、forge の"
                      f"結果可視化(output CONNE)はこれらのセルで正しく描画されない。")

    def validate(self, strict=True):
        """生成した plane/cell の幾何・トポロジーを厳密に検査する。

        検査項目 (forge が前提とする規約):
          1. 各面 surfArea == |surfVect|。
          2. 退化面 (面積≈0) が無いこと。
          3. 内部面の向き: surfVect が iCells[0]->iCells[1] (D=(cc1-cc0)·sv>0)。
          4. 境界面の向き: surfVect が owner から外向き ((pc-cc0)·sv>0)。
          5. iPlanesDir 整合: owner=+1, neighbor=-1。
          6. 隣接整合: 面の節点が両隣セルの節点集合に含まれること (誤接続検出)。
          7. 各セルの閉性: Σ dir·sv ≈ 0 (free-stream/GCL の必要条件)。
          8. 体積 > 0、面のセル数 (内部=2, 境界=1)。

        重大な不整合 (向き/隣接/体積/面積) があれば strict 時に例外で停止する。
        """
        sv = self.surfVect
        area = self.surfArea
        nP = self.nPlanes
        eps = 1e-12
        problems = []
        warns = []

        # 1. area == |surfVect|
        svmag = np.linalg.norm(sv, axis=1)
        amax = float(area.max()) if nP else 1.0
        d_area = float(np.max(np.abs(area - svmag))) if nP else 0.0
        if d_area > 1e-5 * amax:
            problems.append(f"surfArea != |surfVect| (max diff {d_area:.2e})")

        # 2. 退化面
        n_degen = int(np.sum(area < 1e-12 * amax))
        if n_degen:
            problems.append(f"退化面 (面積≈0) が {n_degen} 枚")

        # 3/4. 向き
        nNP = self.nNormalPlanes
        bad_int = bad_bnd = 0
        for ip in range(nP):
            cls = self.pl_cells[ip]
            if len(cls) == 2:
                d = self.cell_cent[cls[1]] - self.cell_cent[cls[0]]
                if np.dot(d, sv[ip]) <= 0.0:
                    bad_int += 1
            else:
                d = self.pl_cent[ip] - self.cell_cent[cls[0]]
                if np.dot(d, sv[ip]) <= 0.0:
                    bad_bnd += 1
        if bad_int:
            problems.append(f"内部面の向き不整合 (c0->c1 でない) {bad_int} 枚")
        if bad_bnd:
            problems.append(f"境界面の向き不整合 (外向きでない) {bad_bnd} 枚")

        # 5. iPlanesDir 整合 + 6. 隣接整合 + 7. 閉性
        cell_node_sets = [set(int(n) for n in nd) for nd in self.cell_nodes]
        bad_dir = 0
        bad_adj = 0
        closure = np.zeros((self.nCells, 3))
        for ic in range(self.nCells):
            for pl, dr in zip(self.cell_planes[ic], self.cell_pdir[ic]):
                closure[ic] += dr * sv[pl]
                # owner(iCells[0]) は +1, neighbor は -1 のはず
                want = 1 if self.pl_cells[pl][0] == ic else -1
                if dr != want:
                    bad_dir += 1
                # 隣接整合: 面節点 ⊂ セル節点
                if not set(int(n) for n in self.pl_nodes[pl]).issubset(cell_node_sets[ic]):
                    bad_adj += 1
        if bad_dir:
            problems.append(f"iPlanesDir 不整合 {bad_dir} 件")
        if bad_adj:
            problems.append(f"隣接整合違反 (面節点が隣接セルに無い) {bad_adj} 件")

        # 7. 閉性
        mean_area_cell = np.array([np.mean([area[pl] for pl in self.cell_planes[ic]])
                                   for ic in range(self.nCells)])
        clo = np.linalg.norm(closure, axis=1) / np.maximum(mean_area_cell, eps)
        clo_max = float(clo.max()) if self.nCells else 0.0
        if clo_max > 1e-3:
            problems.append(f"セル閉性 |Σdir·sv|/面積 が大 (max {clo_max:.2e})")
        elif clo_max > 1e-5:
            warns.append(f"セル閉性 max {clo_max:.2e} (float32 丸めとしてはやや大)")

        # 8b. 薄さゼロ壁 (baffle) 検出: 同位置 (重心一致) の境界面が複数あるか。
        #     forge の境界面は「1 セル + ゲスト」前提で、内部薄板壁 (両側に実セル) は
        #     想定外。fan のハブ等がこれに当たり、後流が準アーティファクトになりうる。
        n_baffle = 0
        bpl = list(range(nNP, nP))
        if bpl:
            bc = np.array([self.pl_cent[i] for i in bpl])
            try:
                from scipy.spatial import cKDTree
                tree = cKDTree(bc)
                pairs = tree.query_pairs(r=1e-9 * max(amax ** 0.5, 1.0))
                n_baffle = len(pairs)
            except Exception:
                n_baffle = -1  # scipy 無し: 未判定

        # 8. 体積・面のセル数
        n_badvol = int(np.sum(self.cell_vol <= 0.0))
        if n_badvol:
            problems.append(f"体積<=0 のセル {n_badvol} 個")
        n_badpc = sum(1 for ip in range(nP)
                      if len(self.pl_cells[ip]) not in (1, 2))
        if n_badpc:
            problems.append(f"面のセル数が 1/2 でない面 {n_badpc} 枚")

        # レポート
        self._log("  --- 幾何/トポロジー検証 ---")
        self._log(f"    surfArea==|surfVect| max diff {d_area:.2e}; 退化面 {n_degen}")
        self._log(f"    向き不整合 internal={bad_int} boundary={bad_bnd}; iPlanesDir 不整合 {bad_dir}")
        self._log(f"    隣接整合違反 {bad_adj}; セル閉性 |Σdir·sv|/面積 max {clo_max:.2e}")
        self._log(f"    体積<=0 {n_badvol}; 面のセル数異常 {n_badpc}")
        if n_baffle > 0:
            warns.append(f"薄さゼロ壁(baffle)候補 {n_baffle} 対 (同位置の境界面)。"
                         f"forge は内部薄板壁を想定設計しておらず、後流が準アーティファクトに"
                         f"なりうる。透過させたい場合は内部面化が必要。")
        for w in warns:
            self._log(f"    [WARN] {w}")
        if problems:
            for p in problems:
                self._log(f"    [FAIL] {p}")
            if strict:
                raise RuntimeError("メッシュ検証に失敗: " + " / ".join(problems))
        else:
            self._log("    検証 OK")
        return not problems


# --------------------------------------------------------------------------
# 幾何ユーティリティ
# --------------------------------------------------------------------------
def _triangulate_faces(faces):
    """四角/多角形の面を扇状に三角形へ分割する (3D)。

    共有面は 1 度だけ分割され、両側セルが同じ三角形を参照する (保存性を保つ)。
    c0/c1/zone は引き継ぐ。三角形は常に平面なので、非平面四角面による
    単一法線の不整合を解消できる。
    """
    out = []
    for fc in faces:
        nd = fc["nodes"]
        if len(nd) <= 3:
            out.append(fc)
            continue
        for k in range(1, len(nd) - 1):
            out.append({
                "nodes": np.array([nd[0], nd[k], nd[k + 1]], dtype=np.int64),
                "c0": fc["c0"], "c1": fc["c1"], "zone": fc["zone"],
            })
    return out


def _area_vector(pts, dim):
    """面の面積ベクトル (大きさ=面積) を返す。向きは未確定 (後で補正)。"""
    if dim == 2:
        # 辺要素 (2 点): 法線 = (dy,-dx,0), 大きさ=長さ
        p0, p1 = pts[0], pts[1]
        t = p1 - p0
        return np.array([t[1], -t[0], 0.0])
    # 多角形: 0.5 * sum(v_i x v_{i+1})
    n = len(pts)
    sv = np.zeros(3)
    for i in range(n):
        sv += np.cross(pts[i], pts[(i + 1) % n])
    return 0.5 * sv


def _face_volume_moment(pts, dim):
    """発散定理による体積の面寄与 ∫_face (x·n) dA を厳密に返す。

    返り値 M は (未反転の面ベクトル sv0 = _area_vector(pts) に対応する向きでの)
    ∫_face x·n dA。セル体積は V = (1/dim) Σ_face (±M) で、± は面ベクトルを
    そのセルの外向きへそろえた符号。

    平面三角形では ∫ x·n dA = (三角形重心)·(三角形面ベクトル) が厳密。
    非平面多角形は扇状に三角形分割して三角形ごとに厳密積分し総和する
    (節点平均重心を使う単一面近似より厳密。tet/平面では従来と一致)。
    """
    if dim == 2:  # 辺要素: 中点·辺法線 が厳密
        c = 0.5 * (pts[0] + pts[1])
        t = pts[1] - pts[0]
        sv = np.array([t[1], -t[0], 0.0])
        return float(np.dot(c, sv))
    n = len(pts)
    if n == 3:
        c = pts.mean(axis=0)
        av = 0.5 * np.cross(pts[1] - pts[0], pts[2] - pts[0])
        return float(np.dot(c, av))
    # 多角形: pts[0] を扇の要に三角形分割
    M = 0.0
    for k in range(1, n - 1):
        av = 0.5 * np.cross(pts[k] - pts[0], pts[k + 1] - pts[0])
        c = (pts[0] + pts[k] + pts[k + 1]) / 3.0
        M += float(np.dot(c, av))
    return M


def _signed_volume_tet(a, b, c, d):
    return np.dot(np.cross(b - a, c - a), d - a) / 6.0


def _reconstruct_cell(faces, coords, dim):
    """セルの面リストから (gmsh要素タイプ, 順序付き節点) を返す。

    順序の再構成は標準要素 (tri/quad/tet/hex/prism/pyramid) のみ対応。
    体積は本ツールでは別途発散定理で計算するため、順序が完全でなくても
    解は正しい (順序は主に XDMF 可視化と ieleType 判定に効く)。
    未対応形状は (0, ユニーク節点) を返す。
    """
    if dim == 2:
        return _reconstruct_polygon(faces, coords)
    nf = len(faces)
    sizes = sorted(len(f) for f in faces)
    uniq = _unique_nodes(faces)
    nv = len(uniq)
    try:
        if nf == 4 and sizes == [3, 3, 3, 3] and nv == 4:
            return GMSH_TET, _order_tet(faces, coords)
        if nf == 6 and all(s == 4 for s in sizes) and nv == 8:
            return GMSH_HEX, _order_hex(faces, coords)
        if nf == 5 and sizes == [3, 3, 4, 4, 4] and nv == 6:
            return GMSH_PRISM, _order_prism(faces, coords)
        if nf == 5 and sizes == [3, 3, 3, 3, 4] and nv == 5:
            return GMSH_PYRAMID, _order_pyramid(faces, coords)
    except Exception as e:  # 再構成失敗時はユニーク節点で続行
        print(f"  WARNING: cell 再構成失敗 ({e}); ユニーク節点で続行")
    return 0, np.array(uniq, dtype=np.int64)


def _unique_nodes(faces):
    seen = []
    s = set()
    for f in faces:
        for n in f:
            n = int(n)
            if n not in s:
                s.add(n)
                seen.append(n)
    return seen


def _cell_edges(faces):
    """セルの全辺 (節点ペア集合) を返す。"""
    edges = set()
    for f in faces:
        m = len(f)
        for i in range(m):
            a, b = int(f[i]), int(f[(i + 1) % m])
            edges.add((min(a, b), max(a, b)))
    return edges


def _reconstruct_polygon(edges_faces, coords):
    """2D: 辺リストから多角形ループを作り CCW にそろえる。"""
    # edges_faces: 各要素は 2 節点
    adj = {}
    for e in edges_faces:
        a, b = int(e[0]), int(e[1])
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    start = next(iter(adj))
    loop = [start]
    prev = None
    cur = start
    while True:
        nbrs = adj[cur]
        nxt = nbrs[0] if nbrs[0] != prev else nbrs[1]
        if nxt == start:
            break
        loop.append(nxt)
        prev, cur = cur, nxt
        if len(loop) > len(adj):
            break
    loop = np.array(loop, dtype=np.int64)
    # CCW 化 (shoelace > 0)
    p = coords[loop]
    area2 = np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1])
    if area2 < 0:
        loop = loop[::-1]
    et = GMSH_TRI if len(loop) == 3 else (GMSH_QUAD if len(loop) == 4 else 0)
    return et, loop


def _order_tet(faces, coords):
    base = [int(x) for x in faces[0]]
    alln = set(_unique_nodes(faces))
    apex = (alln - set(base)).pop()
    a, b, c = base
    if _signed_volume_tet(coords[a], coords[b], coords[c], coords[apex]) < 0:
        b, c = c, b
    return np.array([a, b, c, apex], dtype=np.int64)


def _order_pyramid(faces, coords):
    quad = next(f for f in faces if len(f) == 4)
    base = [int(x) for x in quad]
    alln = set(_unique_nodes(faces))
    apex = (alln - set(base)).pop()
    a, b, c, d = base
    # 体積符号で base 向きを補正
    vol = (_signed_volume_tet(coords[a], coords[b], coords[c], coords[apex])
           + _signed_volume_tet(coords[a], coords[c], coords[d], coords[apex]))
    if vol < 0:
        a, b, c, d = a, d, c, b
    return np.array([a, b, c, d, apex], dtype=np.int64)


def _order_prism(faces, coords):
    tris = [[int(x) for x in f] for f in faces if len(f) == 3]
    edges = _cell_edges(faces)
    bot = tris[0]
    top_set = set(_unique_nodes(faces)) - set(bot)

    def partner(x):
        for y in top_set:
            if (min(x, y), max(x, y)) in edges:
                return y
        raise RuntimeError("prism partner not found")

    top = [partner(x) for x in bot]
    a, b, c = bot
    if _signed_volume_tet(coords[a], coords[b], coords[c], coords[top[0]]) < 0:
        bot = [a, c, b]
        top = [partner(x) for x in bot]
    return np.array(bot + top, dtype=np.int64)


def _order_hex(faces, coords):
    quads = [[int(x) for x in f] for f in faces]
    edges = _cell_edges(faces)
    bot = quads[0]
    bot_set = set(bot)
    # 底面と節点を共有しない面 = 対面 (天面)
    top_face = next(q for q in quads[1:] if not (set(q) & bot_set))
    top_set = set(top_face)

    def partner(x):
        for y in top_set:
            if (min(x, y), max(x, y)) in edges:
                return y
        raise RuntimeError("hex partner not found")

    top = [partner(x) for x in bot]
    a, b, c, d = bot
    if _signed_volume_tet(coords[a], coords[b], coords[c], coords[top[0]]) < 0:
        bot = [a, d, c, b]
        top = [partner(x) for x in bot]
    return np.array(bot + top, dtype=np.int64)


# ==========================================================================
# 初期値・壁距離
# ==========================================================================
def make_initial(fm: ForgeMesh, args):
    n = fm.nCells
    ro = float(args.ro)
    ux, uy, uz = float(args.ux), float(args.uy), float(args.uz)
    gam = float(args.gamma)
    if args.p is not None:
        P = float(args.p)
    elif args.t is not None:
        # P = ro * R * T, R = cp*(gam-1)/gam, cp 不明なら空気 R=287
        R = float(args.R)
        P = ro * R * float(args.t)
    else:
        P = 101325.0
    roe = P / (gam - 1.0) + 0.5 * ro * (ux * ux + uy * uy + uz * uz)
    val = {
        "ro": np.full(n, ro),
        "roUx": np.full(n, ro * ux),
        "roUy": np.full(n, ro * uy),
        "roUz": np.full(n, ro * uz),
        "roe": np.full(n, roe),
    }
    if args.rok is not None:
        val["roK"] = np.full(n, ro * float(args.rok))
    if args.roomega is not None:
        val["roOmega"] = np.full(n, ro * float(args.roomega))
    return val


def compute_wall_distance(fm: ForgeMesh, wall_zone_ids):
    """壁ゾーンの境界面重心までの最短距離をセル重心ごとに計算。"""
    n = fm.nCells
    if not wall_zone_ids:
        print("  wall_dist: 壁ゾーン未指定のため 0 で埋めます (粘性/RANS では再計算が必要)")
        return np.zeros(n)
    pts = []
    for zid in wall_zone_ids:
        for p in fm.bzone_planes.get(zid, []):
            pts.append(fm.pl_cent[p])
    if not pts:
        print("  wall_dist: 該当する壁面が見つからず 0 で埋めます")
        return np.zeros(n)
    pts = np.array(pts)
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(pts)
        d, _ = tree.query(fm.cell_cent)
        return d
    except ImportError:
        print("  wall_dist: scipy 無し -> brute force")
        d = np.empty(n)
        for i in range(n):
            d[i] = np.min(np.linalg.norm(pts - fm.cell_cent[i], axis=1))
        return d


# ==========================================================================
# forge HDF5 書き出し
# ==========================================================================
def write_forge_h5(out_path, fm: ForgeMesh, values, wall_dist, zone_phys,
                   zone_kind, fdtype=np.float32):
    idt = np.int32
    n_phys = len(zone_phys)

    with h5py.File(out_path, "w") as f:
        # /MESH
        g = f.create_group("MESH")
        for k, v in (("nCells", fm.nCells), ("nNodes", fm.nNodes),
                     ("nPlanes", fm.nPlanes), ("nNormalPlanes", fm.nNormalPlanes),
                     ("nBPlanes", fm.nBPlanes), ("nBconds", n_phys)):
            g.attrs.create(k, np.array(v, dtype=idt))

        f.create_dataset("MESH/COORD",
                         data=fm.coords.reshape(-1).astype(fdtype))

        # CONNE (XDMF mixed)
        conne = []
        for ic in range(fm.nCells):
            code = XDMF_CODE.get(fm.cell_etype[ic], 2)  # 不明は polyline(2)
            conne.append(code)
            conne.extend(int(x) for x in fm.cell_nodes[ic])
        f.create_dataset("MESH/CONNE", data=np.array(conne, dtype=idt))

        # /PLANES
        pstruct = []
        for ip in range(fm.nPlanes):
            nd = fm.pl_nodes[ip]
            pstruct.append(len(nd))
            pstruct.extend(int(x) for x in nd)
            cl = fm.pl_cells[ip]
            pstruct.append(len(cl))
            pstruct.extend(int(x) for x in cl)
        f.create_dataset("PLANES/STRUCT", data=np.array(pstruct, dtype=idt))
        f.create_dataset("PLANES/surfVect",
                         data=fm.surfVect.reshape(-1).astype(fdtype))
        f.create_dataset("PLANES/surfArea", data=fm.surfArea.astype(fdtype))
        f.create_dataset("PLANES/centCoords",
                         data=fm.pl_cent.reshape(-1).astype(fdtype))

        # /CELLS
        cstruct = []
        for ic in range(fm.nCells):
            nd = fm.cell_nodes[ic]
            cstruct.append(len(nd))
            cstruct.extend(int(x) for x in nd)
            pls = fm.cell_planes[ic]
            cstruct.append(len(pls))
            cstruct.extend(int(x) for x in pls)
            dirs = fm.cell_pdir[ic]
            cstruct.append(len(dirs))
            cstruct.extend(int(x) for x in dirs)
            cstruct.append(int(fm.cell_etype[ic]))
        f.create_dataset("CELLS/STRUCT", data=np.array(cstruct, dtype=idt))
        f.create_dataset("CELLS/volume", data=fm.cell_vol.astype(fdtype))
        f.create_dataset("CELLS/centCoords",
                         data=fm.cell_cent.reshape(-1).astype(fdtype))
        f.create_dataset("CELLS/regionId",
                         data=fm.fl.cell_region.astype(idt))

        # /VALUE
        for name, arr in values.items():
            f.create_dataset("VALUE/" + name, data=arr.astype(fdtype))
        f.create_dataset("VALUE/wall_dist", data=wall_dist.astype(fdtype))

        # /BCONDS  (physID ごと)
        for zid, planes in fm.bzone_planes.items():
            phys = zone_phys[zid]
            grp = f.create_group(f"BCONDS/{phys}")
            iplanes = np.array(sorted(planes), dtype=idt)
            ibplanes = iplanes - fm.nNormalPlanes
            icells = np.array([fm.pl_cells[p][0] for p in sorted(planes)], dtype=idt)
            grp.create_dataset("iPlanes", data=iplanes)
            grp.create_dataset("iBPlanes", data=ibplanes)
            grp.create_dataset("iCells", data=icells)
            kind = zone_kind.get(zid, "wall")
            # 元コンバータ (HighFive createAttribute<std::string>) と同じ
            # 可変長 UTF-8 文字列属性で書く
            grp.attrs.create("bcondKind", kind,
                             dtype=h5py.string_dtype("utf-8"))

    # XDMF
    _write_xdmf(out_path, fm, list(values.keys()) + ["wall_dist"])
    print(f"-> {out_path} を書き出しました")


def _write_xdmf(out_path, fm, value_names):
    out_path = Path(out_path)
    base = out_path.name
    conne_dim = sum(1 + len(fm.cell_nodes[ic]) for ic in range(fm.nCells))
    lines = [
        "<?xml version='1.0' ?>",
        "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>",
        "<Xdmf>", "  <Domain>",
        "    <Grid GridType='Collection' CollectionType='Spatial' Name='Mixed'>",
        "      <Grid Name='mesh'>",
        f"        <Topology Type='Mixed' NumberOfElements='{fm.nCells}'>",
        f"          <DataItem Format='HDF' DataType='Int' Dimensions='{conne_dim}'>",
        f"           {base}:MESH/CONNE",
        "          </DataItem>", "        </Topology>",
        "        <Geometry Type='XYZ'>",
        f"          <DataItem Format='HDF' DataType='Float' Dimensions='{fm.nNodes*3}'>",
        f"           {base}:MESH/COORD",
        "          </DataItem>", "        </Geometry>",
    ]
    for name in value_names:
        lines += [
            f"        <Attribute Name='{name}' Center='Cell'>",
            f"          <DataItem Format='HDF' DataType='Float' Dimensions='{fm.nCells}'>",
            f"            {base}:VALUE/{name}",
            "          </DataItem>", "        </Attribute>",
        ]
    lines += ["      </Grid>", "    </Grid>", "  </Domain>", "</Xdmf>"]
    out_path.with_suffix(".xmf").write_text("\n".join(lines) + "\n")


# ==========================================================================
# mapping / bcondConfig 入出力
# ==========================================================================
def load_mapping(path):
    """zone -> physID 対応表を読む (yaml/json)。

    形式 (どちらか):
      by_name: { wall: 3, inlet: 1, ... }
      by_id:   { 5: 3, 7: 1, ... }
    また kind を併記してもよい:
      wall:  { physID: 3, kind: wall }
    """
    import yaml
    data = yaml.safe_load(Path(path).read_text())
    by_name, by_id, kinds = {}, {}, {}
    if "by_name" in data or "by_id" in data:
        for k, v in (data.get("by_name") or {}).items():
            by_name[str(k)] = int(v)
        for k, v in (data.get("by_id") or {}).items():
            by_id[int(k)] = int(v)
    else:
        for k, v in data.items():
            if isinstance(v, dict):
                pid = int(v["physID"])
                if "kind" in v:
                    kinds[str(k)] = str(v["kind"])
            else:
                pid = int(v)
            by_name[str(k)] = pid
    return by_name, by_id, kinds


def resolve_zone_phys(fl: FluentMesh, by_name, by_id):
    """各境界 zone id -> physID を解決する。"""
    zone_phys = {}
    unresolved = []
    for zid, z in fl.zones.items():
        if not z["is_boundary"]:
            continue
        if zid in by_id:
            zone_phys[zid] = by_id[zid]
        elif z["name"] in by_name:
            zone_phys[zid] = by_name[z["name"]]
        else:
            unresolved.append((zid, z["name"]))
    return zone_phys, unresolved


def load_bcond_kinds(path):
    """bcondConfig.yaml から physID -> kind を読む (wall_dist 判定/属性用)。"""
    import yaml
    out = {}
    try:
        cfg = yaml.safe_load(Path(path).read_text())
    except FileNotFoundError:
        return out
    for name, d in (cfg or {}).items():
        if isinstance(d, dict) and "physID" in d:
            out[int(d["physID"])] = str(d.get("kind", "wall"))
    return out


# ==========================================================================
# サブコマンド
# ==========================================================================
def cmd_dump(args):
    with h5py.File(args.input, "r") as f:
        print(f"=== HDF5 tree: {args.input} ===")

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"  [D] {name}  shape={obj.shape} dtype={obj.dtype}")
            else:
                print(f"  [G] {name}")
            for ak, av in obj.attrs.items():
                print(f"        @{ak} = {av!r}")
        f.visititems(visit)
    print("\n=== zones ===")
    fl = FluentMesh(args.input, nodes_path=args.nodes_path, verbose=True)
    for zid, z in sorted(fl.zones.items()):
        kind = "boundary" if z["is_boundary"] else "interior"
        print(f"  id={zid:<6} name={z['name']:<20} type={z['zoneType']} "
              f"{kind:<9} nfaces={z['nfaces']}")
    fl.close()


def cmd_template(args):
    fl = FluentMesh(args.input, nodes_path=args.nodes_path, verbose=True)
    bzones = [(zid, z) for zid, z in sorted(fl.zones.items()) if z["is_boundary"]]
    lines = ["# Fluent zone -> forge physID 対応表 (本ツール template 出力)",
             "# 編集して convert --map に渡す。kind は任意 (wall_dist 判定に使用)。"]
    pid = 1
    bc_lines = ["# bcondConfig.yaml 雛型 (本ツール template 出力)",
                "# kind/floats を物理条件に合わせて編集すること。"]
    for zid, z in bzones:
        name = z["name"]
        guess = _guess_kind(z)
        lines.append(f"{name}: {{ physID: {pid}, kind: {guess} }}  # zone id={zid}")
        bc_lines += [
            f"{name}:",
            f"  physID: {pid}",
            f"  kind: {guess}",
            "  outputHDFflg: 0",
            "  ints:",
            "  floats:",
        ]
        if guess.startswith("inlet"):
            bc_lines += ["    Ux: 0.0", "    Uy: 0.0", "    Uz: 0.0",
                         "    Ps: 101325.0", "    Ts: 300.0"]
        elif guess == "outflow":
            bc_lines += ["    Ps: 101325.0"]
        pid += 1
    out = args.output or "mapping.yaml"
    Path(out).write_text("\n".join(lines) + "\n")
    print(f"-> {out} を書き出しました")
    if args.bcond:
        Path(args.bcond).write_text("\n".join(bc_lines) + "\n")
        print(f"-> {args.bcond} (雛型) を書き出しました")
    fl.close()


def _guess_kind(z):
    # 名前のキーワードを優先する。mesh のみの *.msh.h5 では BC 未設定で
    # zoneType が wall(3) 既定になりがちで当てにならないため (例: 入口 in1 が type3)。
    name = (z.get("name") or "").lower()
    if "inlet" in name or "inflow" in name or name.startswith("in"):
        return "inlet_uniformVelocity"
    if "outlet" in name or "outflow" in name or "exit" in name or name.startswith("out"):
        return "outflow"
    if "sym" in name or "axis" in name:
        return "slip"
    if "far" in name:
        return "inlet_uniformVelocity"
    if "wall" in name:
        return "wall"
    # 名前で判らなければ Fluent zoneType から推定 (*.cas.h5 では信頼できる)
    zt = z.get("zoneType")
    if isinstance(zt, (int, np.integer)) and int(zt) in FLUENT_ZTYPE_KIND:
        return FLUENT_ZTYPE_KIND[int(zt)]
    return "wall"


def cmd_convert(args):
    fl = FluentMesh(args.input, nodes_path=args.nodes_path,
                    scale=float(args.scale), verbose=True)

    if not args.map:
        sys.exit("--map (zone->physID 対応表) が必要です。"
                 " 先に `template` で雛型を作ってください。")
    by_name, by_id, map_kinds = load_mapping(args.map)
    zone_phys, unresolved = resolve_zone_phys(fl, by_name, by_id)
    if unresolved:
        print("ERROR: 以下の境界 zone が mapping で解決できません:")
        for zid, name in unresolved:
            print(f"   id={zid} name={name}")
        sys.exit(1)

    # physID -> kind (bcondConfig 優先, 無ければ mapping の kind, 既定 wall)
    phys_kind = {}
    if args.bcond:
        phys_kind = load_bcond_kinds(args.bcond)
    zone_kind = {}
    for zid in zone_phys:
        phys = zone_phys[zid]
        if phys in phys_kind:
            zone_kind[zid] = phys_kind[phys]
        elif fl.zones[zid]["name"] in map_kinds:
            zone_kind[zid] = map_kinds[fl.zones[zid]["name"]]
        else:
            zone_kind[zid] = _guess_kind(fl.zones[zid])

    fm = ForgeMesh(fl, triangulate=args.triangulate, verbose=True)

    # 生成した plane/cell の幾何・トポロジー (面ベクトル向き・隣接・閉性・体積) を検証。
    # 重大な不整合があれば例外で停止 (--no-validate で抑止)。
    fm.validate(strict=not args.no_validate)

    # 壁ゾーン = kind が wall で始まるもの
    wall_zone_ids = [zid for zid, k in zone_kind.items() if k.startswith("wall")]
    values = make_initial(fm, args)
    wall_dist = compute_wall_distance(fm, wall_zone_ids)

    fdtype = np.float64 if args.double else np.float32
    write_forge_h5(args.output, fm, values, wall_dist, zone_phys, zone_kind,
                   fdtype=fdtype)

    # 解決状況のサマリ
    print("\n=== boundary summary ===")
    for zid in sorted(zone_phys):
        z = fl.zones[zid]
        print(f"  zone id={zid:<5} name={z['name']:<18} -> physID={zone_phys[zid]:<3} "
              f"kind={zone_kind[zid]:<22} nplanes={len(fm.bzone_planes.get(zid, []))}")
    print("\n注意: bcondConfig.yaml の physID が上の physID と一致している必要があります。")
    fl.close()


def build_argparser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pd = sub.add_parser("dump", help="HDF5 構造とゾーン一覧を表示")
    pd.add_argument("input")
    pd.add_argument("--nodes-path", default=None)
    pd.set_defaults(func=cmd_dump)

    pt = sub.add_parser("template", help="zone->physID 対応表と bcondConfig 雛型を出力")
    pt.add_argument("input")
    pt.add_argument("-o", "--output", default="mapping.yaml")
    pt.add_argument("--bcond", default=None, help="bcondConfig.yaml 雛型の出力先")
    pt.add_argument("--nodes-path", default=None)
    pt.set_defaults(func=cmd_template)

    pc = sub.add_parser("convert", help="Fluent .h5 -> forge .h5 変換")
    pc.add_argument("input")
    pc.add_argument("output")
    pc.add_argument("--map", required=True, help="zone->physID 対応表 (yaml/json)")
    pc.add_argument("--bcond", default=None,
                    help="bcondConfig.yaml (kind 判定/壁距離に使用)")
    pc.add_argument("--nodes-path", default=None)
    pc.add_argument("--scale", default=1.0,
                    help="座標スケール係数 (in->m: 0.0254, mm->m: 0.001)")
    pc.add_argument("--double", action="store_true",
                    help="geom_float=double ビルド向けに float64 で書く")
    pc.add_argument("--triangulate", action="store_true",
                    help="3D の四角/多角形面を三角形に分割 (非平面面の不整合回避)")
    pc.add_argument("--no-validate", action="store_true",
                    help="生成メッシュの幾何/トポロジー検証で停止しない (既定は検証して不整合なら停止)")
    # 初期値
    pc.add_argument("--ro", default=1.225, help="初期密度")
    pc.add_argument("--ux", default=0.0)
    pc.add_argument("--uy", default=0.0)
    pc.add_argument("--uz", default=0.0)
    pc.add_argument("--p", default=None, help="初期静圧 (指定時 roe を計算)")
    pc.add_argument("--t", default=None, help="初期温度 (P 未指定時 P=ro*R*T)")
    pc.add_argument("--R", default=287.0, help="気体定数 (T 指定時)")
    pc.add_argument("--gamma", default=1.4)
    pc.add_argument("--rok", default=None, help="初期 k (RANS, roK=ro*k)")
    pc.add_argument("--roomega", default=None, help="初期 omega (RANS)")
    pc.set_defaults(func=cmd_convert)

    return p


def main():
    args = build_argparser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
