# mesh_salome.py -- Salome (SMESH) で STEP -> NETGEN tet + Viscous Layers (prism BL) -> MED
#
# 実行 (ヘッドレス):
#   /home/sano/opt/salome/SALOME-9.14.0-native-UB24.04-SRC/salome -t <this file> args:<step>,<out.med>
#   (無引数なら pintle_fluid_half.step -> pintle_salome.med)
#
# netgen 生 API の boundary_layers はこの形状の凹角 (T字給気・スロート・ピントル基部) で
# 破綻するため、SMESH の StdMeshers_ViscousLayers (凹角の層縮退処理を持つ) を使う。
# 面グループは mesh_pintle.py と同じ幾何条件で自動分類 (対話ピック不要)。
import math
import sys

import salome

salome.salome_init()
import GEOM  # noqa: E402
import SMESH  # noqa: E402
from salome.geom import geomBuilder  # noqa: E402
from salome.smesh import smeshBuilder  # noqa: E402

# ===== 設定 =====
# 既定値。cwd に mesh_salome_config.json があれば上書き (テンプレは cad/ に同梱)。
# 実際に使った値は OUT と同じ場所に mesh_settings_used.json として保存する。
import json
import os

CFG = dict(
    step="pintle_fluid_half.step",
    out="pintle_salome.med",
    # 設計寸法 [mm] (面分類しきい値; build_geom.py の PARAMS と整合させる)
    # plume_L>0 でプルーム領域あり (plume_far/plume_out/base グループが増える)
    geometry=dict(x_in=0.0, x_exit=35.0, x_tip=15.0, z_bot=-20.0, Rp=1.5, Rt=2.5,
                  tip_len=3.0, plume_L=60.0, plume_R=25.0, Re=4.33),
    # 大域メッシュ [mm] (プルームありのときは粗く: プルーム内は maxh_mm、
    #                    ノズル系の面は nozzle_maxh_mm の局所サイズで細かく)
    maxh_mm=3.0,
    minh_mm=0.05,
    fineness=3,                # NETGEN Fineness (0-5, 3=Moderate)
    nozzle_maxh_mm=0.5,        # ノズル系の面 (wall/wall_pintle/inlet) の局所サイズ
    # ピントル先端 (ノーズ+キャップ) の局所細分 [mm]。強曲率 (r_tip~0.4mm) で
    # Viscous Layers が層を落とすのを防ぐ。None で無効。
    tip_maxh_mm=0.12,
    lip_maxh_mm=0.3,           # 出口リップ円エッジの局所細分 (プルームあり時の skew 対策)
    bl_on_base=True,           # base (リップ接続壁板) にも viscous layer を張る
    # Viscous Layers: 総厚 / 層数 / stretch (第一層厚 ~ TOTAL*(s-1)/(s^N-1))
    vl_total_mm=0.25,
    vl_nlayers=4,
    vl_stretch=1.3,
)
if os.path.exists("mesh_salome_config.json"):
    with open("mesh_salome_config.json") as _f:
        _user = json.load(_f)
    for _k, _v in _user.items():
        if _k == "geometry" and isinstance(_v, dict):
            CFG["geometry"].update(_v)
        else:
            CFG[_k] = _v
    print("loaded mesh_salome_config.json:", _user)

argv = sys.argv[1:]
if argv:
    if len(argv) >= 1 and argv[0]:
        CFG["step"] = argv[0]
    if len(argv) >= 2 and argv[1]:
        CFG["out"] = argv[1]

STEP = CFG["step"]
OUT = CFG["out"]
G_MM = CFG["geometry"]
MAXH_MM = CFG["maxh_mm"]
MINH_MM = CFG["minh_mm"]
VL_TOTAL_MM = CFG["vl_total_mm"]
VL_NLAYER = CFG["vl_nlayers"]
VL_STRETCH = CFG["vl_stretch"]

geompy = geomBuilder.New()
smesh = smeshBuilder.New()

shape = geompy.ImportSTEP(STEP, False, True)
faces = geompy.ExtractShapes(shape, geompy.ShapeType["FACE"], True)

# 単位スケール検出: 設計 x 全長 (プルーム込み) が bbox で幾つに見えるか (1.0=mm, 0.001=m)
bb = geompy.BoundingBox(shape)   # (xmin,xmax,ymin,ymax,zmin,zmax)
x_len_mm = (G_MM["x_exit"] + G_MM.get("plume_L", 0)) - G_MM["x_in"]
unit_scale = (bb[1] - bb[0]) / x_len_mm
print("unit_scale (GEOM units per mm):", unit_scale)
G = {k: v * unit_scale for k, v in G_MM.items()}
MAXH = MAXH_MM * unit_scale
MINH = MINH_MM * unit_scale
VL_TOTAL = VL_TOTAL_MM * unit_scale
TOL_BIG = 1.0 * unit_scale       # 面分類の粗しきい値 (1mm 相当)
TOL_SYM = 0.5 * unit_scale       # 対称面判定 (0.5mm 相当)


def classify(x, y, z):
    r = math.hypot(y, z)
    has_plume = G.get("plume_L", 0) > 0
    x_end = G["x_exit"] + (G.get("plume_L", 0) if has_plume else 0)
    if abs(z - G["z_bot"]) < TOL_BIG:
        return "inlet"
    if has_plume:
        if abs(x - x_end) < TOL_BIG:
            return "plume_out"                    # プルーム下流円盤
        if abs(x - G["x_exit"]) < TOL_BIG and r > G["Re"]:
            return "base"                         # 出口面のノズル外環状板
        if x > G["x_exit"] + TOL_BIG and r > 0.5 * G["plume_R"]:
            return "plume_far"                    # プルーム外周円筒
    else:
        if abs(x - G["x_exit"]) < TOL_BIG:
            return "outlet"
    # wall_pintle を symmetry より先に判定する (順序が重要)。
    # ピントル先端キャップ (底面) は重心 y~0 のため、|y| 判定を先にすると
    # symmetry に誤分類され viscous layer の除外リストに入る + BC も slip になる。
    if (G["x_in"] - 2 * TOL_BIG <= x <= G["x_tip"] + 2 * TOL_BIG) \
            and r < 0.5 * (G["Rp"] + G["Rt"]):
        return "wall_pintle"
    if abs(y) < TOL_SYM:
        return "symmetry"
    return "wall"


bins = {}
face_cdg = {}
for f in faces:
    cm = geompy.MakeCDG(f)
    x, y, z = geompy.PointCoordinates(cm)
    bins.setdefault(classify(x, y, z), []).append(f)
    face_cdg[f] = (x, y, z)
print("face groups:", {k: len(v) for k, v in bins.items()})
HAS_PLUME = G_MM.get("plume_L", 0) > 0
REQ = {"inlet", "symmetry", "wall", "wall_pintle"} | \
      ({"plume_out", "plume_far", "base"} if HAS_PLUME else {"outlet"})
missing = REQ - set(bins)
if missing:
    raise SystemExit("ERROR: 未検出の面グループ %s" % missing)

ggrp = {}
for name, fl in bins.items():
    g = geompy.CreateGroup(shape, geompy.ShapeType["FACE"])
    geompy.UnionList(g, fl)
    g.SetName(name)
    geompy.addToStudyInFather(shape, g, name)
    ggrp[name] = g

mesh = smesh.Mesh(shape, "pintle")
algo = mesh.Tetrahedron(smeshBuilder.NETGEN_1D2D3D)
par = algo.Parameters()
par.SetMaxSize(MAXH)
par.SetMinSize(MINH)
par.SetFineness(CFG["fineness"])
par.SetSecondOrder(0)       # forge は線形要素のみ
par.SetOptimize(1)

# ピントル先端 (ノーズ+キャップ) の局所細分。強曲率 (r_tip) 面で Viscous Layers が
# 層を落とすのを防ぐ (実測: 局所細分なしだと先端の層あり率 78%)。
if CFG.get("tip_maxh_mm"):
    tip_maxh = CFG["tip_maxh_mm"] * unit_scale
    x_nose = (G_MM["x_tip"] - G_MM["tip_len"]) * unit_scale
    n_tip = 0
    for f in bins["wall_pintle"]:
        if face_cdg[f][0] >= x_nose:
            # SetLocalSizeOnShape は study 登録済みオブジェクトが必要 -> 先に publish
            geompy.addToStudyInFather(shape, f, "tip_face_%d" % n_tip)
            par.SetLocalSizeOnShape(f, tip_maxh)
            n_tip += 1
    print("tip local size %.4g on %d faces" % (tip_maxh, n_tip))

# Viscous Layers: 層を張らない面 = inlet/outlet/symmetry (isFacesToIgnore=True)。
# bl_on_base=True (既定) なら base (出口リップに接続する外側の壁板) にも層を張る:
# リップ角を prism が両側から巻き、リップ近傍のセル品質と ω の挙動が改善する。
no_layer = ["inlet", "symmetry"] + \
           (["plume_out", "plume_far"] if HAS_PLUME else ["outlet"])
if HAS_PLUME and not CFG.get("bl_on_base", True):
    no_layer.append("base")
ignore_ids = []
for n in no_layer:
    for f in bins.get(n, []):
        ignore_ids.append(geompy.GetSubShapeID(shape, f))
algo.ViscousLayers(VL_TOTAL, VL_NLAYER, VL_STRETCH, ignore_ids, True)

# プルームありのとき: 出口リップ (base 面の内円エッジ, r=Re) に局所細分。
# 拡大部壁の終端×base×プルームの 90°角で skew>0.9 のセルが出るのを防ぐ。
if HAS_PLUME and CFG.get("lip_maxh_mm"):
    lip = CFG["lip_maxh_mm"] * unit_scale
    Re_s = G_MM["Re"] * unit_scale
    n_lip = 0
    for f in bins.get("base", []):
        for e in geompy.ExtractShapes(f, geompy.ShapeType["EDGE"], True):
            cm = geompy.MakeCDG(e)
            ex, ey, ez = geompy.PointCoordinates(cm)
            if abs(math.hypot(ey, ez) - Re_s * 2 / math.pi) < Re_s:  # 内円 (半円 CDG r=2R/π)
                r_e = math.hypot(ey, ez)
                if r_e < 1.5 * Re_s:      # 外円 (2*plume_R/π) は除外
                    geompy.addToStudyInFather(shape, e, "lip_edge_%d" % n_lip)
                    par.SetLocalSizeOnShape(e, lip)
                    n_lip += 1
    print("lip local size %.4g on %d edges" % (lip, n_lip))

# プルームありのとき: ノズル系の面 (wall/wall_pintle/inlet) に局所細サイズ
# (大域 maxh はプルーム用に粗いので、ノズル側を nozzle_maxh_mm で維持する)
if HAS_PLUME and CFG.get("nozzle_maxh_mm"):
    nsz = CFG["nozzle_maxh_mm"] * unit_scale
    n_loc = 0
    for ngrp in ("wall", "wall_pintle", "inlet"):
        for f in bins.get(ngrp, []):
            geompy.addToStudyInFather(shape, f, "noz_face_%d" % n_loc)
            par.SetLocalSizeOnShape(f, nsz)
            n_loc += 1
    print("nozzle local size %.4g on %d faces" % (nsz, n_loc))

ok = mesh.Compute()
print("Compute:", ok)
if not ok:
    raise SystemExit("ERROR: mesh.Compute() failed")

for name, g in ggrp.items():
    mesh.GroupOnGeom(g, name, SMESH.FACE)

# 要素型の内訳 (prism が入ったか確認)
print("nodes:", mesh.NbNodes(), "tet:", mesh.NbTetras(), "prism:", mesh.NbPrisms(),
      "pyram:", mesh.NbPyramids(), "hex:", mesh.NbHexas())
mesh.ExportMED(OUT)
print("wrote", OUT)

# 使った設定を成果物と並べて保存 (トレーサビリティ; run ディレクトリに残る)
used = dict(CFG)
used["unit_scale"] = unit_scale
used_path = os.path.join(os.path.dirname(os.path.abspath(OUT)), "mesh_settings_used.json")
with open(used_path, "w") as f:
    json.dump(used, f, indent=2, ensure_ascii=False)
print("wrote", used_path)
