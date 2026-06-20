#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""case/36.passive_pseudoshock_control 用 Salome (GEOM/SMESH) 2D 平面メッシュ生成。

Matsuo et al. (1988) のパッシブコントロール試験部 (多孔壁+キャビティ) を 2D 平面で再現する。
同ディレクトリの gmsh 版 `make_mesh.py` と **同一形状・同一グループ名・同一 physID** を作る
ことを目的とし、forge へは meshio で gmsh msh2 に変換して取り込む (`med_to_gmsh.py` 参照)。

実行 (Salome TUI):
    salome -t make_mesh_salome.py                # POROUS=True
    salome -t make_mesh_salome.py args:solid     # POROUS=False (固体壁・比較基準)
    salome -t make_mesh_salome.py args:porous,netgen   # 構造ダクトを使わず全 NETGEN

出力 (このスクリプトのあるディレクトリに生成):
    passive_porous.med / .unv   (POROUS=True)
    passive_solid.med  / .unv   (POROUS=False)

座標系: x=流れ方向(右が下流), y=上下(中心線 y=0, 上下対称), z 無し(平面)。原点 x=0 は入口。
全寸法は mm で書き、出力時に SCALE=0.001 を掛けて [m] にする (forge は SI)。

--- forge 取り込み規約 (gmshReader.hpp より) ---
  * 2D 面要素(quad/tri) + 1D 境界線要素(line) を持つ平面メッシュ。押し出し不要。
  * 境界 = Physical Curve(dim=1), 流体 = Physical Surface(dim=2)。physID は整数で
    bcondConfig.yaml と突き合わせる。スロット口/ブロック界面などの内部エッジは
    どのグループにも入れない (→ 線要素として書き出されず conformal に連続する)。
  * z=0、巻き方向は forge 側で CW/CCW 自動判定するため不問。
"""
import sys
import math

# ============================================================================
# 1. パラメータ (全てここに集約)
# ============================================================================
# --- 実行モード (args で上書き可) ---
POROUS         = True    # True: 多孔壁+キャビティ, False: 連続固体壁
STRUCTURED_DUCT = True   # True: 直管部を Quadrangle 写像(構造)。False: 全 NETGEN quad-dominant

# args:solid / args:porous / args:netgen / args:struct を解釈
for _a in sys.argv[1:]:
    for _tok in _a.replace("args:", "").split(","):
        t = _tok.strip().lower()
        if t == "solid":   POROUS = False
        elif t == "porous": POROUS = True
        elif t == "netgen": STRUCTURED_DUCT = False
        elif t == "struct": STRUCTURED_DUCT = True

NAME    = "passive_porous" if POROUS else "passive_solid"

# --- 形状 [mm] ---
HIN_HALF = 10.7545       # 入口半高さ (入口高さ 21.509mm, M_in=1.689)
TAN1     = math.tan(math.radians(1.0))   # 片側 1deg こう配
X_IN     = 0.0
X_RUN    = 30.0          # 助走領域端 (slip -> no-slip 切替)
X_C0     = 130.0         # キャビティ部 上流端 (多孔壁上流端)
X_C1     = 162.0         # キャビティ部 下流端 (=130+32)
X_OUT    = 690.0         # 出口 (=162+528)

PLATE_T  = 10.0          # 多孔板厚
CAV_D    = 20.0          # キャビティ深さ
N_SLOT   = 8             # スロット本数 (流れ方向)
SLOT_PITCH = 4.0         # スロットピッチ
SLOT_W   = 0.8           # スロット幅 (開口率 8*0.8/32 = 0.2)
SLOT_C0  = X_C0 + 2.0    # 最初のスロット中心 (=132) -> 132,136,...,160

# --- メッシュ密度 (構造ダクト) ---
NY       = 64            # ダクト高さ方向セル数
NX_RUN   = 50            # 助走領域 (x 0->30)
NX_DIFF  = 120           # 拡大ダクト (x 30->130)
NX_DOWN  = 300           # 下流直管 (x 162->690)
CLUSTER_C = 8.0          # 壁近傍クラスタリング強度 (density f(t)=1+C*(2t-1)^2; 大きいほど壁に密)

# --- 特性長 [m] (非構造領域 / NETGEN・local size) ---
CL_CORE  = 0.00035       # コア(キャビティ部ダクト本体)
CL_SLOT  = 0.00015       # スロット幅方向 (0.8mm を ~5 セル -> 法線方向 >=4 セル)
CL_CAV   = 0.0006        # キャビティ空洞

SCALE    = 0.001         # mm -> m
TOL_MM   = 1.0e-4        # 座標一致判定許容 [mm]

# 出力先 = このスクリプトと同じディレクトリ
import os
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
MED_PATH = os.path.join(HERE, NAME + ".med")
UNV_PATH = os.path.join(HERE, NAME + ".unv")


def yw(x):
    """中心線からの半高さ [mm]。"""
    return HIN_HALF if x <= X_RUN else HIN_HALF + (x - X_RUN) * TAN1


slot_centers = [SLOT_C0 + i * SLOT_PITCH for i in range(N_SLOT)]
slots = [(c - SLOT_W / 2.0, c + SLOT_W / 2.0) for c in slot_centers]

print("=" * 70)
print(f"[make_mesh_salome] NAME={NAME}  POROUS={POROUS}  STRUCTURED_DUCT={STRUCTURED_DUCT}")
print(f"  yw: in={yw(0):.4f} c0={yw(X_C0):.4f} c1={yw(X_C1):.4f} out={yw(X_OUT):.4f} mm")
print(f"  slots(center mm) = {slot_centers}")
print("=" * 70)


# ============================================================================
# 2. Salome 初期化
# ============================================================================
import salome
salome.salome_init()

from salome.geom import geomBuilder
geompy = geomBuilder.New()

from salome.smesh import smeshBuilder
import SMESH
smesh = smeshBuilder.New()


# ============================================================================
# 3. GEOM: 流体面を構築 -> Partition で conformal な単一コンパウンドに
# ============================================================================
def V(x_mm, y_mm):
    return geompy.MakeVertex(x_mm * SCALE, y_mm * SCALE, 0.0)


def face_from_pts(pts_mm):
    """(x_mm,y_mm) の列から閉ポリラインで平面 face を作る。"""
    verts = [V(x, y) for (x, y) in pts_mm]
    wire = geompy.MakePolyline(verts, True)   # isClosed=True
    return geompy.MakeFaceWires([wire], True)


# --- ダクト本体ブロック (フル高さ; 上壁 +yw, 下壁 -yw) ---
def duct_rect(xa, xb):
    """xa->xb のフル高さ矩形(台形) face。反時計回り: 左下->右下->右上->左上。"""
    return face_from_pts([(xa, -yw(xa)), (xb, -yw(xb)), (xb, yw(xb)), (xa, yw(xa))])


faces = []
B1 = duct_rect(X_IN,  X_RUN);  faces.append(B1)   # 助走 (slip)
B2 = duct_rect(X_RUN, X_C0);   faces.append(B2)   # 拡大ダクト (no-slip)
B4 = duct_rect(X_C1,  X_OUT);  faces.append(B4)   # 下流直管 (no-slip)


# --- キャビティ部コア (x[130,162], フル高さ) ---
def core_wall_pts(sign, left_to_right=True):
    """コアの上(sign=+1)/下(sign=-1)壁線の点列(左->右)。
    POROUS のときスロット口で分割した点を含む。"""
    y = lambda xx: sign * yw(xx)
    pts = [(X_C0, y(X_C0))]
    if POROUS:
        for (sl, sr) in slots:
            pts.append((sl, y(sl)))   # 壁セグメント終端 = スロット左
            pts.append((sr, y(sr)))   # スロット口 (内部境界)
        pts.append((X_C1, y(X_C1)))
    else:
        pts.append((X_C1, y(X_C1)))
    return pts if left_to_right else list(reversed(pts))


# コア多角形 (反時計回り): 下壁(左->右) -> 右辺(下->上) -> 上壁(右->左) -> 左辺(上->下)
core_pts = []
core_pts += core_wall_pts(-1, left_to_right=True)            # 下壁 左->右
core_pts += core_wall_pts(+1, left_to_right=False)           # 上壁 右->左
B3 = face_from_pts(core_pts);  faces.append(B3)              # コア


# --- キャビティ + スロット (comb 形状, 上下) ---
def comb_pts(sign):
    """sign=+1 上側 / -1 下側。床(Yf)->スロット指->天井(Yc) の閉ループ点列。"""
    Yf = lambda xx: sign * (yw(xx) + PLATE_T)            # キャビティ床 (=多孔板上面)
    Yc = lambda xx: sign * (yw(xx) + PLATE_T + CAV_D)    # キャビティ天井
    yw_ = lambda xx: sign * yw(xx)                        # ダクト壁線
    pts = [(X_C0, Yf(X_C0))]
    for (sl, sr) in slots:
        pts.append((sl, Yf(sl)))     # 床: 前のセグメント -> スロット左床
        pts.append((sl, yw_(sl)))    # スロット左側壁 (床->壁線)
        pts.append((sr, yw_(sr)))    # スロット口 (内部境界, core と共有)
        pts.append((sr, Yf(sr)))     # スロット右側壁 (壁線->床)
    pts.append((X_C1, Yf(X_C1)))     # 床 -> 右端
    pts.append((X_C1, Yc(X_C1)))     # 右側壁 (床->天井)
    pts.append((X_C0, Yc(X_C0)))     # 天井 (右->左)
    return pts


if POROUS:
    comb_top = face_from_pts(comb_pts(+1));  faces.append(comb_top)
    comb_bot = face_from_pts(comb_pts(-1));  faces.append(comb_bot)


# --- Partition: 全 face を imprint/glue して conformal な単一形状に ---
part = geompy.MakePartition(faces, [], [], [], geompy.ShapeType["FACE"], 0, [], 0)
geompy.addToStudy(part, NAME + "_geom")


# ============================================================================
# 4. サブ形状の分類 (座標ベース; 名前付きグループを作る)
# ============================================================================
def cdg_mm(shape):
    x, y, z = geompy.PointCoordinates(geompy.MakeCDG(shape))
    return x / SCALE, y / SCALE


def edge_mid_mm(e):
    p = geompy.MakeVertexOnCurve(e, 0.5)
    x, y, z = geompy.PointCoordinates(p)
    return x / SCALE, y / SCALE


def edge_len_mm(e):
    return geompy.BasicProperties(e)[0] / SCALE


def edge_is_vertical(e):
    vs = geompy.SubShapeAll(e, geompy.ShapeType["VERTEX"])
    if len(vs) != 2:
        return False
    x0, y0, _ = geompy.PointCoordinates(vs[0])
    x1, y1, _ = geompy.PointCoordinates(vs[1])
    return abs(x1 - x0) / SCALE < 1.0e-3


all_faces = geompy.SubShapeAll(part, geompy.ShapeType["FACE"])
all_edges = geompy.SubShapeAll(part, geompy.ShapeType["EDGE"])

# --- 各エッジが何枚の face に属するか (内部=2, 境界=1) を中点キーで数える ---
def key_of(xy):
    return (round(xy[0], 3), round(xy[1], 3))

face_count = {}
for f in all_faces:
    for e in geompy.SubShapeAll(f, geompy.ShapeType["EDGE"]):
        k = key_of(edge_mid_mm(e))
        face_count[k] = face_count.get(k, 0) + 1


def is_internal(e):
    return face_count.get(key_of(edge_mid_mm(e)), 1) >= 2


def near(a, b, tol=0.05):
    return abs(a - b) <= tol


# --- face 分類 (centroid で) ---
def classify_face(f):
    cx, cy = cdg_mm(f)
    if X_IN + 0.1 < cx < X_RUN - 0.1:
        return "B1"
    if X_RUN + 0.1 < cx < X_C0 - 0.1:
        return "B2"
    if X_C1 + 0.1 < cx < X_OUT - 0.1:
        return "B4"
    if X_C0 - 0.1 < cx < X_C1 + 0.1:
        if abs(cy) < yw(cx) + 0.1:
            return "CORE"
        return "COMB_TOP" if cy > 0 else "COMB_BOT"
    return "OTHER"


fmap = {"B1": [], "B2": [], "B4": [], "CORE": [], "COMB_TOP": [], "COMB_BOT": [], "OTHER": []}
for f in all_faces:
    fmap[classify_face(f)].append(f)

# --- edge 分類 ---
g_inlet, g_outlet, g_slip, g_wall = [], [], [], []
# 1D 制御用グループ
e_height = []        # 高さ方向 縦エッジ (inlet/outlet/x30/x130/x162, |y|<yw) -> NY
e_xrun, e_xdiff, e_xdown = [], [], []  # 直管壁 (上下)
e_coreline = []      # コア壁線 (x[130,162], |y|~yw): セグメント+スロット口
e_combfine, e_combcoarse = [], []      # キャビティ/スロット


for e in all_edges:
    mx, my = edge_mid_mm(e)
    L = edge_len_mm(e)
    vert = edge_is_vertical(e)
    internal = is_internal(e)
    aty = abs(my)

    # --- 境界グループ (count==1 のみ) ---
    if not internal:
        if near(mx, X_IN):
            g_inlet.append(e)
        elif near(mx, X_OUT):
            g_outlet.append(e)
        elif (X_IN < mx < X_RUN) and near(aty, yw(mx), 0.2):
            g_slip.append(e)          # 助走領域 上下壁 = slip
        else:
            g_wall.append(e)          # それ以外の固体壁 = no-slip

    # --- 1D 制御グループ (内部/境界問わず) ---
    on_wallline = near(aty, yw(mx), 0.25)   # ダクト壁線 (上下) 付近か
    if vert and aty < 0.5:
        # 中心線をまたぐフル高さ縦エッジ
        e_height.append(e)
    elif on_wallline and (X_IN < mx < X_RUN):
        e_xrun.append(e)
    elif on_wallline and (X_RUN < mx < X_C0):
        e_xdiff.append(e)
    elif on_wallline and (X_C1 < mx < X_OUT):
        e_xdown.append(e)
    elif on_wallline and (X_C0 - 0.05 <= mx <= X_C1 + 0.05):
        e_coreline.append(e)          # コア壁線 (セグメント + スロット口)
    elif aty > yw(mx) + 0.3:
        # comb (キャビティ + スロット指)
        if L < 1.2:
            e_combfine.append(e)      # スロット幅 0.8mm 系
        else:
            e_combcoarse.append(e)    # 床/側壁/天井
    else:
        e_combcoarse.append(e)        # 取りこぼし救済


def make_group(name, shape_type, sublist):
    grp = geompy.CreateGroup(part, shape_type)
    if sublist:
        geompy.UnionList(grp, sublist)
    geompy.addToStudyInFather(part, grp, name)
    return grp


ST_EDGE = geompy.ShapeType["EDGE"]
ST_FACE = geompy.ShapeType["FACE"]

grp_inlet  = make_group("inlet",     ST_EDGE, g_inlet)
grp_outlet = make_group("outlet",    ST_EDGE, g_outlet)
grp_slip   = make_group("wall_slip", ST_EDGE, g_slip)
grp_wall   = make_group("wall",      ST_EDGE, g_wall)
grp_fluid  = make_group("fluid",     ST_FACE, all_faces)

print(f"[geom] faces={len(all_faces)} edges={len(all_edges)}")
print(f"[geom] B1={len(fmap['B1'])} B2={len(fmap['B2'])} B4={len(fmap['B4'])} "
      f"CORE={len(fmap['CORE'])} COMB_TOP={len(fmap['COMB_TOP'])} COMB_BOT={len(fmap['COMB_BOT'])}")
print(f"[geom] boundary edges: inlet={len(g_inlet)} outlet={len(g_outlet)} "
      f"wall_slip={len(g_slip)} wall={len(g_wall)}")


# ============================================================================
# 5. SMESH
# ============================================================================
Mesh = smesh.Mesh(part, NAME)


def add_1d_submesh(geom_list, name, n_seg=None, local_len=None, density_expr=None):
    """geom_list の各エッジを束ねたグループに Regular_1D + 1D 仮説を割り当てる。"""
    if not geom_list:
        return None
    grp = make_group("_1d_" + name, ST_EDGE, geom_list)
    algo = Mesh.Segment(geom=grp)
    if n_seg is not None:
        hyp = algo.NumberOfSegments(n_seg)
        if density_expr is not None:
            try:
                hyp.SetDistrType(2)               # 2 = distribution with analytic density
                hyp.SetConversionMode(1)
                hyp.SetExpressionFunction(density_expr)
            except Exception as ex:
                print(f"  [warn] density on {name} failed ({ex}); uniform にフォールバック")
    elif local_len is not None:
        algo.LocalLength(local_len)
    return algo


if STRUCTURED_DUCT:
    # ---- 構造ダクト: 直管ブロックは Quadrangle 写像、コア/comb は NETGEN quad-dominant ----
    # 1D: 高さ方向は両壁クラスタリング (対称 density), 流れ方向はブロックごと等分割
    dens = "1+%g*(2*t-1)^2" % CLUSTER_C
    add_1d_submesh(e_height, "height", n_seg=NY, density_expr=dens)
    add_1d_submesh(e_xrun,  "xrun",  n_seg=NX_RUN)
    add_1d_submesh(e_xdiff, "xdiff", n_seg=NX_DIFF)
    add_1d_submesh(e_xdown, "xdown", n_seg=NX_DOWN)
    add_1d_submesh(e_coreline,    "coreline",   local_len=CL_SLOT)
    add_1d_submesh(e_combfine,    "combfine",   local_len=CL_SLOT)
    add_1d_submesh(e_combcoarse,  "combcoarse", local_len=CL_CAV)

    # 取りこぼしエッジ用の global 1D フォールバック
    Mesh.Segment().LocalLength(CL_CORE)

    # global 2D = NETGEN 2D (quad-dominant): コア/comb をカバー
    n2d = Mesh.Triangle(algo=smeshBuilder.NETGEN_2D)
    try:
        p2d = n2d.Parameters()
        p2d.SetMaxSize(CL_CORE)
        p2d.SetMinSize(CL_SLOT * 0.5)
        p2d.SetOptimize(1)
        p2d.SetFineness(2)
        p2d.SetQuadAllowed(1)
    except Exception as ex:
        print(f"  [warn] NETGEN_2D params: {ex}")

    # 直管ブロックを Quadrangle (Mapping) で上書き (構造)
    struct_faces = list(fmap["B1"]) + list(fmap["B2"]) + list(fmap["B4"])
    if not POROUS:
        struct_faces += list(fmap["CORE"])
    for i, f in enumerate(struct_faces):
        grp = make_group("_q_%d" % i, ST_FACE, [f])
        Mesh.Quadrangle(geom=grp)

else:
    # ---- 全 NETGEN 1D-2D quad-dominant (堅牢なフォールバック) ----
    algo = Mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
    params = algo.Parameters()
    params.SetMaxSize(CL_CORE)
    params.SetMinSize(CL_SLOT * 0.5)
    params.SetOptimize(1)
    params.SetFineness(2)
    try:
        params.SetQuadAllowed(1)
    except Exception as ex:
        print(f"  [warn] SetQuadAllowed: {ex}")
    # local size でスロット/コア/キャビティを解像
    for e in e_coreline + e_combfine:
        try: params.SetLocalSizeOnShape(e, CL_SLOT)
        except Exception: pass
    for e in e_combcoarse:
        try: params.SetLocalSizeOnShape(e, CL_CAV)
        except Exception: pass


# ---- 計算 ----
ok = Mesh.Compute()
print(f"[smesh] Compute -> {ok}")


# ---- 名前付きグループを mesh 側に作成 (MED/UNV に保存される) ----
def mesh_group_on_geom(geom_grp, name, elem_type):
    g = Mesh.GroupOnGeom(geom_grp, name, elem_type)
    return g


mg_inlet  = mesh_group_on_geom(grp_inlet,  "inlet",     SMESH.EDGE)
mg_outlet = mesh_group_on_geom(grp_outlet, "outlet",    SMESH.EDGE)
mg_slip   = mesh_group_on_geom(grp_slip,   "wall_slip", SMESH.EDGE)
mg_wall   = mesh_group_on_geom(grp_wall,   "wall",      SMESH.EDGE)
mg_fluid  = mesh_group_on_geom(grp_fluid,  "fluid",     SMESH.FACE)


# ============================================================================
# 6. エクスポート (MED + UNV)
# ============================================================================
try:
    Mesh.ExportMED(MED_PATH, auto_groups=0)
except TypeError:
    Mesh.ExportMED(MED_PATH)
Mesh.ExportUNV(UNV_PATH)
print(f"[export] MED -> {MED_PATH}")
print(f"[export] UNV -> {UNV_PATH}")


# ============================================================================
# 7. 検証用統計
# ============================================================================
nq = Mesh.NbQuadrangles()
nt = Mesh.NbTriangles()
nfaces = Mesh.NbFaces()
nedges = Mesh.NbEdges()
print("-" * 70)
print("[stats]")
print(f"  2D cells: quad={nq}  tri={nt}  total={nfaces}  "
      f"(quad ratio={100.0*nq/max(nfaces,1):.1f}%)")
print(f"  1D boundary segments (all): {nedges}")
for name, g in [("inlet", mg_inlet), ("outlet", mg_outlet),
                ("wall_slip", mg_slip), ("wall", mg_wall), ("fluid", mg_fluid)]:
    print(f"  group '{name}': {g.Size()} elems")
print(f"  slot width target cells: {SLOT_W*SCALE/CL_SLOT:.1f} "
      f"(SLOT_W={SLOT_W}mm / CL_SLOT={CL_SLOT*1e3:.3f}mm)")

# 最小セルサイズ (面積 control の最小辺長換算の目安)
try:
    area_min, area_max = Mesh.GetMinMax(SMESH.FT_Area)
    print(f"  cell area: min={area_min:.3e} max={area_max:.3e} m^2 "
          f"(min sqrt={math.sqrt(max(area_min,0)):.3e} m)")
except Exception:
    pass
print("-" * 70)
print(f"[done] {NAME}: 次に  python3 med_to_gmsh.py {NAME}.med  で gmsh msh2 に変換")

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()
