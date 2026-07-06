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

# ===== 設定 (mesh_pintle.py と整合) =====
STEP = "pintle_fluid_half.step"
OUT = "pintle_salome.med"
if len(sys.argv) > 1 and "," in ",".join(sys.argv):  # salome -t script args:a,b
    pass
argv = sys.argv[1:]
if argv:
    if len(argv) >= 1 and argv[0]:
        STEP = argv[0]
    if len(argv) >= 2 and argv[1]:
        OUT = argv[1]

# 設計寸法 [mm]。Salome は STEP の mm 単位を m へ変換して取り込むことがあるため、
# 取り込み後の bbox から unit_scale を自動検出し、しきい値・メッシュ長をスケールする。
G_MM = dict(x_in=0.0, x_exit=35.0, x_tip=15.0, z_bot=-20.0, Rp=1.5, Rt=2.5)
MAXH_MM = 0.5
MINH_MM = 0.05
# Viscous Layers: 総厚 / 層数 / stretch (第一層厚 ~ TOTAL*(s-1)/(s^N-1))
VL_TOTAL_MM = 0.25
VL_NLAYER = 4
VL_STRETCH = 1.3

geompy = geomBuilder.New()
smesh = smeshBuilder.New()

shape = geompy.ImportSTEP(STEP, False, True)
faces = geompy.ExtractShapes(shape, geompy.ShapeType["FACE"], True)

# 単位スケール検出: 設計 x 全長 35mm が bbox で幾つに見えるか (1.0=mm, 0.001=m)
bb = geompy.BoundingBox(shape)   # (xmin,xmax,ymin,ymax,zmin,zmax)
unit_scale = (bb[1] - bb[0]) / (G_MM["x_exit"] - G_MM["x_in"])
print("unit_scale (GEOM units per mm):", unit_scale)
G = {k: v * unit_scale for k, v in G_MM.items()}
MAXH = MAXH_MM * unit_scale
MINH = MINH_MM * unit_scale
VL_TOTAL = VL_TOTAL_MM * unit_scale
TOL_BIG = 1.0 * unit_scale       # 面分類の粗しきい値 (1mm 相当)
TOL_SYM = 0.5 * unit_scale       # 対称面判定 (0.5mm 相当)


def classify(x, y, z):
    r = math.hypot(y, z)
    if abs(z - G["z_bot"]) < TOL_BIG:
        return "inlet"
    if abs(x - G["x_exit"]) < TOL_BIG:
        return "outlet"
    if abs(y) < TOL_SYM:
        return "symmetry"
    if (G["x_in"] - 2 * TOL_BIG <= x <= G["x_tip"] + 2 * TOL_BIG) \
            and r < 0.5 * (G["Rp"] + G["Rt"]):
        return "wall_pintle"
    return "wall"


bins = {}
for f in faces:
    cm = geompy.MakeCDG(f)
    x, y, z = geompy.PointCoordinates(cm)
    bins.setdefault(classify(x, y, z), []).append(f)
print("face groups:", {k: len(v) for k, v in bins.items()})
missing = {"inlet", "outlet", "symmetry", "wall", "wall_pintle"} - set(bins)
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
par.SetFineness(3)          # Moderate
par.SetSecondOrder(0)       # forge は線形要素のみ
par.SetOptimize(1)

# Viscous Layers: 層を張らない面 = inlet/outlet/symmetry (isFacesToIgnore=True)
ignore_ids = []
for n in ("inlet", "outlet", "symmetry"):
    for f in bins[n]:
        ignore_ids.append(geompy.GetSubShapeID(shape, f))
algo.ViscousLayers(VL_TOTAL, VL_NLAYER, VL_STRETCH, ignore_ids, True)

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
