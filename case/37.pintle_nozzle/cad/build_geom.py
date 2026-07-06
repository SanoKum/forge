# build_geom.py  --  pintle nozzle 流体ドメイン (半割モデル) を生成する [単位: mm]
#
# 形状: 下方から立ち上がる円形給気管 -> 90°ベンド -> 収縮拡大ノズル。
#       ノズル軸上にピントルが左端壁から +x に刺さる。x-z 平面 (y=0) 対称の半割。
#
# 実行:
#   freecadcmd build_geom.py            # GUI 不要 (出力: pintle_fluid_half.step)
#   (GUI の場合は Python コンソールに貼り付け)
#
# 構築順序 (ブーリアンを壊さないコツ):
#   全体を作る -> ピントルを cut -> 最後に y>=0 で half-cut
#
# 座標系:
#   x = ノズル軸 (流れは +x に排出) / z = 鉛直 (給気管は下から +z で立ち上がる)
#   対称面 = y=0 (x-z 平面)、半割は y>=0 を残す
#
# 既定はノズル輪郭をコーン (直線) 近似。曲線輪郭に拡張するときは
# make_chamber_nozzle の収縮・スロート部の LineSegment を Part.Arc に差し替える。

import FreeCAD as App
import Part, math
from FreeCAD import Vector

# ===== パラメータ (実設計値: 論文 Fig.2 のノズル寸法 + ベンド/パイプの実寸で埋める) =====
# 論文スケール仮組み (Pan/Chen/Ye 2017 Fig.2)。確実な値は L=35mm / ER=3.0 / 2.5=Rt のみ。
# 図は径方向が非スケールのため Rc・収縮部は推定。給気管/ベンドは論文に無く暫定値。
P = dict(
    # --- ノズル+チャンバー (x 軸まわりの回転体, r は z 方向) ---
    x_in     = 0.0,                 # 左端壁 (ピントル取付面) の x
    Lc       = 5.0,                 # チャンバー直管長 (x_in -> x_cyl_end)
    Rc       = 6.0,                 # チャンバー半径 (推定; 収縮比 Rc/Rt=2.4)
    x_throat = 15.0,                # スロート x (全長 35 の ~42%)
    Rt       = 2.5,                 # スロート半径 (= 論文 Fig.2 の 2.5mm)
    x_exit   = 35.0,                # 出口 x (= 論文 L=35mm)
    Re       = 2.5 * math.sqrt(3.0),# 出口半径 (膨張比 3.0 -> Re = Rt*sqrt(3) = 4.33)
    # --- 給気管 (垂直円柱; 論文に無いため暫定。Rbend/Zjoin/stub は現行 vcyl-only では未使用) ---
    Rpipe    = 3.0,                 # 給気管半径 (暫定; ~Rc/2)
    Xfeed    = 4.0,                 # 給気管中心の x (チャンバー直管部の下)
    Rbend    = 6.0,                 # (予約: 将来の滑らかエルボ用)
    Zjoin    = -3.0,                # (予約: 同上)
    z_bot    = -20.0,               # 入口 (パイプ下端) の z (下方 20mm)
    stub     = 3.0,                 # (予約: 同上)
    # --- ピントル (x 軸まわりの回転体) ---
    Rp       = 1.5,                 # ピントル胴半径 (2.92 vs 2.5 のスロート閉塞から ~1.5)
    x_tip    = 15.0,                # 先端 x (=スロート。位置を振るときはここを変える)
    tip_len  = 3.0,                 # 先端ノーズ長
    r_tip    = 0.4,                 # 先端の小半径 (>0)。尖らせると軸上に退化頂点ができ
                                    # Netgen 表面メッシュが破綻するため blunt 化する
    eps      = 0.5,                 # 取付壁を貫通させる余白 (壁裏に薄膜を残さない)
)


def make_chamber_nozzle(P):
    """x-z 平面の閉プロファイル (x, 0, r) を x 軸まわりに 360°回転 -> チャンバー+ノズル流体"""
    x_cyl_end = P['x_in'] + P['Lc']
    pts = [
        Vector(P['x_in'],     0, 0),
        Vector(P['x_in'],     0, P['Rc']),    # 左端壁を上へ
        Vector(x_cyl_end,     0, P['Rc']),    # 直管
        Vector(P['x_throat'], 0, P['Rt']),    # 収縮 (コーン)
        Vector(P['x_exit'],   0, P['Re']),    # 拡大 (コーン)
        Vector(P['x_exit'],   0, 0),          # 出口を軸へ
    ]
    edges = [Part.LineSegment(pts[i], pts[i + 1]).toShape() for i in range(len(pts) - 1)]
    edges.append(Part.LineSegment(pts[-1], pts[0]).toShape())  # 軸に沿って閉じる
    face = Part.Face(Part.Wire(edges))
    return face.revolve(Vector(0, 0, 0), Vector(1, 0, 0), 360.0)


def make_feed(P):
    """垂直円柱の給気管部品を [list] で返す -> 下から立ち上がりチャンバー下面へ入る。
    流れはチャンバー内で 90°向きを変えてノズルへ向かう。inlet は管下端 (z=z_bot) の円盤。
    list で返すのは chamber と一括 multi-fuse する方が OCC の operand 取りこぼしを避けられるため。

    NOTE: 当初は垂直+水平円柱の sharp エルボにしたが、膝部の円柱×円柱×chamber 三重交差が
    sliver 面 (面積~25) を生み Netgen 表面メッシュを破綻させたため垂直円柱のみに簡素化した。
    滑らかな 90°ベンド (fillet / torus) は将来の改良 (boolean 安定化が要)。"""
    Xf, Rp, zb = P['Xfeed'], P['Rpipe'], P['z_bot']
    # 垂直円柱: 入口下端 z_bot -> 軸 (z=0)。chamber (z=[-Rc,Rc]) と十分重ねて融合・sliver 回避。
    vcyl = Part.makeCylinder(Rp, 0.0 - zb, Vector(Xf, 0, zb), Vector(0, 0, 1))
    return [vcyl]


def make_pintle(P):
    """取付壁裏 (x=x_in-eps) から先端 (x_tip) まで。胴=円筒, 先端=円弧 (small arc)。
    コーン先端にするなら最後の Part.Arc を LineSegment に差し替える。"""
    x0  = P['x_in'] - P['eps']
    xbe = P['x_tip'] - P['tip_len']         # 胴とノーズの境
    rt  = P['r_tip']
    pa  = Vector(xbe, 0, P['Rp'])
    pt  = Vector(P['x_tip'], 0, rt)         # ノーズ先端 (小半径 rt; 尖らせない)
    pc  = Vector(P['x_tip'], 0, 0)          # 先端フラットキャップの軸側
    pb  = Vector(P['x_tip'] - P['tip_len'] * 0.3, 0, P['Rp'] * 0.7 + rt * 0.3)  # ノーズの膨らみ
    head = [Vector(x0, 0, 0), Vector(x0, 0, P['Rp']), pa]            # 壁裏 -> 上 -> 胴端
    edges = [Part.LineSegment(head[i], head[i + 1]).toShape() for i in range(len(head) - 1)]
    edges.append(Part.Arc(pa, pb, pt).toShape())                    # ノーズ円弧 (small arc)
    edges.append(Part.LineSegment(pt, pc).toShape())                # 先端フラットキャップ (半径 rt)
    edges.append(Part.LineSegment(pc, head[0]).toShape())           # 軸で閉じる
    face = Part.Face(Part.Wire(edges))
    return face.revolve(Vector(0, 0, 0), Vector(1, 0, 0), 360.0)


def main():
    # OCC boolean の operand 取りこぼし対策: chamber + 給気管部品を fuzzy で一括 multi-fuse。
    FUZ = 1.0e-2
    big = 1000.0
    yneg = Part.makeBox(big, big, big, Vector(-big / 2, -big, -big / 2))  # y<0 領域

    def assemble():
        f = make_chamber_nozzle(P).fuse(make_feed(P), FUZ)  # 流体 = チャンバー + パイプ (一括)
        f = f.cut(make_pintle(P))                           # ピントルを引く
        f = f.cut(yneg)                                     # 半割: y<0 を削る (y>=0 を残す)
        return f.removeSplitter()                           # 同一面の継ぎ目を統合 (メッシュに優しい)

    # OCC boolean は稀に run 間で取りこぼす (operand 消失 / 部分clip)。
    # 期待形状 (x[0,70], y[0,Rc], z[z_bot,Rc], 単一閉ソリッド) を満たすまで作り直す。
    fluid = None
    for attempt in range(5):
        fluid = assemble()
        b = fluid.BoundBox
        ok = (len(fluid.Solids) == 1 and bool(fluid.Shells) and fluid.Shells[0].isClosed()
              and abs(b.YMin) < 1e-2 and abs(b.XMax - P['x_exit']) < 1e-2
              and abs(b.ZMin - P['z_bot']) < 1e-2)
        if ok:
            break
        print(f"  retry {attempt}: bad bbox x[{b.XMin:.1f},{b.XMax:.1f}] "
              f"y[{b.YMin:.2f},{b.YMax:.1f}] z[{b.ZMin:.1f},{b.ZMax:.1f}] solids {len(fluid.Solids)}")

    b = fluid.BoundBox
    closed = bool(fluid.Shells) and fluid.Shells[0].isClosed()
    print(f"valid: {fluid.isValid()}  solids: {len(fluid.Solids)}  closed: {closed}")
    print(f"bbox x[{b.XMin:.2f},{b.XMax:.2f}] y[{b.YMin:.3f},{b.YMax:.2f}] z[{b.ZMin:.2f},{b.ZMax:.2f}]")
    fluid.exportStep("pintle_fluid_half.step")
    print("wrote pintle_fluid_half.step")


main()
