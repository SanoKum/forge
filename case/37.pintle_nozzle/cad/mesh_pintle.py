#!/usr/bin/env python3
"""pintle_fluid_half.step を Netgen で非構造メッシュ化し、forge 互換 gmsh 4.1 (forge.msh) を書く。

実行 (forge の mesh venv で):
    source /home/sano/work/forge/.venv-mesh/bin/activate
    cd case/37.pintle_nozzle/cad
    python3 mesh_pintle.py            # -> forge.msh
    # 続けて:
    #   convertGmshToForge forge.msh forge.h5   (要 solverConfig.yaml + bcondConfig.yaml)
    #   check_mesh_quality.py forge.h5          (AR<=1000 / skew<=0.9)

背景 (なぜ Netgen を直接使うか):
- Salome 公式バイナリは登録必須でCIに乗らない。Salome が内部で使う NETGEN エンジン
  (STEP取込 + viscous layers) は pip の netgen-mesher で同等に得られる。
- Netgen の "Gmsh2 Format" 出力は $PhysicalNames を書かないため convertGmshToForge が
  読めない。そこで本スクリプトが Netgen mesh から forge 必須要件を満たす msh4.1 を直接書く:
  (1) gmsh 4.1 形式  (2) 1エンティティ=物理タグ1個  (3) 流体ボリュームに 3D 物理グループ
  (4) 外表面を全て面グループで覆う  (5) 線形要素のみ。

面命名は STEP の面重心の幾何条件で自動分類 (対話ピック不要・再現可能)。
"""
import argparse
import math
from collections import defaultdict
from netgen.occ import OCCGeometry
from netgen.meshing import MeshingParameters, BoundaryLayerParameters

# ===== 設定 (build_geom.py の PARAMS と整合させる。品質 FAIL 時はここを詰める) =====
STEP = "pintle_fluid_half.step"
OUT  = "forge.msh"

# 形状定数 (面分類のしきい値に使用。build_geom.py の PARAMS と一致させる)
G = dict(x_in=0.0, x_exit=35.0, x_tip=15.0, z_bot=-20.0, Rp=1.5, Rt=2.5)

# メッシュサイズ (throat 環状すきま Rt-Rp=1mm を割るため maxh<gap が必須。
#                coarse だと表面自己交差・品質 FAIL になる)
MAXH            = 0.4
GRADING         = 0.3
CURVATURESAFETY = 2.0

# 出力座標スケール。CAD/メッシュは mm で作るが forge は SI (m) 前提のため既定 0.001。
# (これを 1.0 で書くと 35mm ノズルが 35m になる。旧 run_0001/0002 はこのバグ持ち。)
SCALE = 0.001

# --- 近壁解像度 ---
# (A) 壁面近傍の tet 細分 (推奨・堅牢): 壁 face に局所 maxh を与え近壁を細かくする。
#     prism 無し・クラッシュ無し。壁関数 (y+~30-80) と相性が良い。None で無効。
#     例: WALL_MAXH=0.12 + MAXH=0.5 で ~50万 tet (小型形状のため細分は cell 数増に注意)。
WALL_MAXH = None          # 壁面の局所メッシュ長 [mm] (None=無効)

# (B) prism 境界層 (netgen boundary_layers)。
#     ** 現状この形状 (T字給気・スロート・ピントル基部の凹角) では BL 生成が失敗する **
#     project_boundaries 無し=segfault、有りでも "too many attempts" で体積充填が失敗。
#     使うには凹角エッジを fillet 化する必要がある。USE_BL=True は fillet 済み形状でのみ。
#     半割では project_boundaries で層を対称面/入口/出口に沿わせるのが必須。
USE_BL     = False
BL_THICK   = [0.02, 0.035, 0.06, 0.10]  # 各層厚 [mm] (壁から外向き)
BL_BOUND   = "wall.*"                   # 正規表現で wall, wall_pintle にマッチ
BL_PROJECT = "symmetry|inlet|outlet"    # 層をこれらの面に沿って終端させる (半割で必須)

# 物理 ID 割り当て (bcondConfig.yaml と一致させること)
PHYS = {"fluid": (3, 1), "wall": (2, 10), "inlet": (2, 11),
        "outlet": (2, 12), "symmetry": (2, 13), "wall_pintle": (2, 14)}
SURF_ENT = {"wall": 1, "inlet": 2, "outlet": 3, "symmetry": 4, "wall_pintle": 5}


def face_center(f):
    c = f.center
    try:
        return c.x, c.y, c.z
    except AttributeError:
        return c[0], c[1], c[2]


def classify(cx, cy, cz):
    r = math.hypot(cy, cz)
    if abs(cz - G['z_bot']) < 1.0:
        return "inlet"
    if abs(cx - G['x_exit']) < 1.0:
        return "outlet"
    # wall_pintle を symmetry より先に判定 (順序が重要): ピントル先端キャップ (底面) は
    # 重心 y~0 のため |y| 判定を先にすると symmetry に誤分類され BC が slip になる。
    if (G['x_in'] - 2 <= cx <= G['x_tip'] + 2) and r < 0.5 * (G['Rp'] + G['Rt']):
        return "wall_pintle"                            # 軸近傍 x<=tip -> ピントル
    if abs(cy) < 0.5:                                    # y=0 対称面 (複数枚)
        return "symmetry"
    return "wall"


def write_msh41(mesh, path):
    """Netgen mesh -> forge 互換 gmsh 4.1 を書く (節点 1ブロック, 要素は面グループ別ブロック)。"""
    N = len(mesh.Points())
    pts = [mesh.Points()[i].p for i in range(1, N + 1)]
    tets = [[v.nr for v in e.vertices] for e in mesh.Elements3D()]
    # 現状 tet のみ対応。prism/pyramid (境界層 USE_BL 有効時) が出たら node 順マップが必要
    # なので silent 不正を避けて明示的に失敗させる。
    if any(len(t) != 4 for t in tets):
        raise SystemExit("ERROR: 非 tet の 3D 要素あり (prism/pyramid)。write_msh41 は tet のみ対応。"
                         " 境界層を使うなら gmsh type 6/7 の書き出しと node 順マップを追加すること。")
    tris = defaultdict(list)
    for e in mesh.Elements2D():
        tris[mesh.GetBCName(e.index - 1)].append([v.nr for v in e.vertices])
    nTri = sum(len(v) for v in tris.values())
    nTet = len(tets)
    nElem = nTri + nTet

    with open(path, "w") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n%d\n" % len(PHYS))
        for nm, (d, pid) in PHYS.items():
            f.write('%d %d "%s"\n' % (d, pid, nm))
        f.write("$EndPhysicalNames\n")
        # Entities: 0 点, 0 曲線, len(SURF_ENT) 面, 1 体 (各 1 物理タグ)
        f.write("$Entities\n0 0 %d 1\n" % len(SURF_ENT))
        for nm, et in sorted(SURF_ENT.items(), key=lambda kv: kv[1]):
            f.write("%d 0 0 0 0 0 0 1 %d 0\n" % (et, PHYS[nm][1]))
        f.write("1 0 0 0 0 0 0 1 %d 0\n" % PHYS["fluid"][1])
        f.write("$EndEntities\n")
        # Nodes: 1 ブロック (entityDim=3, entityTag=1)
        f.write("$Nodes\n1 %d 1 %d\n" % (N, N))
        f.write("3 1 0 %d\n" % N)
        for i in range(1, N + 1):
            f.write("%d\n" % i)
        for x, y, z in pts:
            f.write("%.10g %.10g %.10g\n" % (x * SCALE, y * SCALE, z * SCALE))
        f.write("$EndNodes\n")
        # Elements: 面グループ別の三角形ブロック + 四面体ブロック (タグは 1..nElem 連番)
        f.write("$Elements\n%d %d 1 %d\n" % (len(SURF_ENT) + 1, nElem, nElem))
        tag = 1
        for nm, et in sorted(SURF_ENT.items(), key=lambda kv: kv[1]):
            tl = tris[nm]
            f.write("2 %d 2 %d\n" % (et, len(tl)))
            for v in tl:
                f.write("%d %d %d %d\n" % (tag, v[0], v[1], v[2]))
                tag += 1
        f.write("3 1 4 %d\n" % nTet)
        for v in tets:
            f.write("%d %d %d %d %d\n" % (tag, v[0], v[1], v[2], v[3]))
            tag += 1
        f.write("$EndElements\n")
    return N, nTet, {k: len(v) for k, v in tris.items()}


def main():
    global STEP, OUT, MAXH, WALL_MAXH, SCALE
    ap = argparse.ArgumentParser(description="STEP -> Netgen tet mesh -> forge msh4.1")
    ap.add_argument("--step", default=STEP)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--maxh", type=float, default=MAXH)
    ap.add_argument("--wall-maxh", type=float, default=WALL_MAXH,
                    help="壁 face の局所メッシュ長 [mm] (省略=無効)")
    ap.add_argument("--scale", type=float, default=SCALE,
                    help="出力座標スケール (既定 0.001 = mm->m, forge は SI)")
    a = ap.parse_args()
    STEP, OUT, MAXH, WALL_MAXH, SCALE = a.step, a.out, a.maxh, a.wall_maxh, a.scale

    shape = OCCGeometry(STEP).shape
    cnt = defaultdict(int)
    for f in shape.faces:
        nm = classify(*face_center(f))
        f.name = nm
        if WALL_MAXH is not None and nm in ("wall", "wall_pintle"):
            f.maxh = WALL_MAXH                       # (A) 近壁 tet 細分
        cnt[nm] += 1
    for s in shape.solids:
        s.name = "fluid"
    print("face groups:", dict(cnt))
    missing = set(SURF_ENT) - set(cnt)
    if missing:
        raise SystemExit("ERROR: 未検出の面グループ %s (classify しきい値を見直す)" % missing)

    mp = MeshingParameters(maxh=MAXH, grading=GRADING, curvaturesafety=CURVATURESAFETY)
    if USE_BL:  # (B) prism 境界層 (現状この形状では失敗する。凹角 fillet 化が前提)
        blp = BoundaryLayerParameters(boundary=BL_BOUND, thickness=BL_THICK,
                                      new_material=None, domain="fluid", outside=False,
                                      project_boundaries=BL_PROJECT)
        mesh = OCCGeometry(shape).GenerateMesh(mp=mp, boundary_layers=[blp])
    else:
        mesh = OCCGeometry(shape).GenerateMesh(mp=mp)
    print("mesh: 3D", len(mesh.Elements3D()), "2D", len(mesh.Elements2D()),
          "pts", len(mesh.Points()))

    N, nTet, tri = write_msh41(mesh, OUT)
    print("wrote %s  nodes=%d tets=%d tris=%s" % (OUT, N, nTet, tri))


if __name__ == "__main__":
    main()
