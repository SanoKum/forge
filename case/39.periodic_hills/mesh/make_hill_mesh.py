#!/usr/bin/env python3
"""周期丘 (periodic hills, Re=10595) 用の構造 hex メッシュを gmsh 4.1 形式で直接生成する。

幾何: ERCOFTAC/Almeida 標準丘形状。多項式定義は NASA TMR 配布の hill-geometry.dat
(https://tmbwg.github.io/turbmodels/Other_LES_Data/2Dhill_periodic/hill-geometry.dat)
の正本をそのまま使用 (単位 mm, 丘高さ h=28mm, 半丘 x∈[0,54mm])。
領域は 9h × 3.035h × 4.5h、山頂が x=0 と x=9h (x/z 周期)。
上壁は平坦 y=3.035h。風上斜面は x=9h−54mm〜9h に鏡映で配置。

格子: DES 標準構成 (ESAIM Proc. 16 (2007) 133-145 の 160×100×60, Δx+,Δz+<35,
Δy1+≲1 / ATAAC D3.2-36 mandatory 160×160×60 相当) に合わせた 160×100×60。
- x: 一様 (Δx=0.056h)
- y: 各カラムを壁 y_w(x)〜3.035h で両側幾何ストレッチ (50+50)。正規化分布は
  全カラム共通で、平坦部で y1=1.5e-3 h (→y1+≲1) となるよう公比を決める。
  山頂カラムは局所高さが小さいぶん y1 も比例して薄くなる (u_tau 最大部で有利)。
- z: 一様 (Δz=0.075h, Δz+~35)

msh4.1 ライタは case/37.pintle_nozzle/cad/med_to_msh41.py の形式に準拠
(1エンティティ=物理タグ1 / 体に3D物理グループ / 線形要素)。
physID は case/38 チャネルと同一規約: xlo=1 xhi=2 ylo=3 yhi=4 zlo=5 zhi=6 fluid=7。
周期対 (xlo↔xhi, zlo↔zhi) の節点座標は同一コードパスで生成するため厳密一致する。

使い方:
    python3 make_hill_mesh.py [出力.msh] [--nx 160] [--ny 100] [--nz 60] [--y1 1.5e-3]
    (--y1 は平坦部カラムの第一層厚 [h 単位]。粗メッシュ試行用に解像度を引数化)
"""
import argparse
import sys

import numpy as np

# ---------------- 幾何パラメータ ----------------
H_MM = 28.0                # 丘高さ [mm] (正本の単位)
H = H_MM * 1.0e-3          # 丘高さ [m]
LX = 9.0 * H               # 周期長
LY_TOP = 3.035 * H         # 上壁高さ
LZ = 4.5 * H               # スパン
X_FOOT_MM = 54.0           # 丘すその x [mm]

# ---------------- 格子パラメータ (既定; CLI で上書き可) ----------------
NX, NY, NZ = 160, 100, 60
NYH = NY // 2              # 片側の要素数
Y1_FLAT = 1.5e-3 * H       # 平坦部カラム (高さ 3.035h) での第一層厚

# NASA TMR hill-geometry.dat の多項式 (x[mm] -> y[mm])。係数は正本の転記。
_SEG = [
    (0.0, 9.0, (2.800000000000e+01, 0.0, 6.775070969851e-03, -2.124527775800e-03)),
    (9.0, 14.0, (2.507355893131e+01, 9.754803562315e-01, -1.016116352781e-01, 1.889794677828e-03)),
    (14.0, 20.0, (2.579601052357e+01, 8.206693007457e-01, -9.055370274339e-02, 1.626510569859e-03)),
    (20.0, 30.0, (4.046435022819e+01, -1.379581654948e+00, 1.945884504128e-02, -2.070318932190e-04)),
    (30.0, 40.0, (1.792461334664e+01, 8.743920332081e-01, -5.567361123058e-02, 6.277731764683e-04)),
    (40.0, 54.0, (5.639011190988e+01, -2.010520359035e+00, 1.644919857549e-02, 2.674976141766e-05)),
]


NOCLIP = False  # main() が CLI から設定

def hill_half_mm(x_mm):
    """半丘 (山頂 x=0 → すそ x=54mm) の高さ [mm]。区間外は 0。

    NOCLIP=True では区間 1 の min(28,·) 頭打ちを外す (多項式そのまま)。山頂 ±3.19mm の
    C1 折れ (壁勾配 1.2° 不連続) が除去され、逸脱は最大 +0.0047mm (=1.7e-4 h)。"""
    if x_mm >= X_FOOT_MM:
        return 0.0
    for seg_i, (lo, hi, c) in enumerate(_SEG):
        if lo <= x_mm <= hi:
            y = c[0] + c[1] * x_mm + c[2] * x_mm**2 + c[3] * x_mm**3
            if NOCLIP and seg_i == 0:
                return max(0.0, y)
            return min(H_MM, max(0.0, y))
    raise ValueError(x_mm)


def wall_y(x):
    """下壁 y_w(x) [m]。山頂 x=0、鏡映山頂 x=LX。"""
    x_mm = x * 1.0e3
    lx_mm = LX * 1.0e3
    if x_mm <= lx_mm * 0.5:
        return hill_half_mm(x_mm) * 1.0e-3
    return hill_half_mm(lx_mm - x_mm) * 1.0e-3


def solve_progression(q1, n):
    """(r-1)/(r^n-1) = q1 となる公比 r を二分法で解く (q1 = y1/半高さ)。"""
    lo, hi = 1.0 + 1e-12, 2.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        val = (r - 1.0) / (r**n - 1.0)
        if val > q1:
            lo = r
        else:
            hi = r
    return 0.5 * (lo + hi)


def main():
    global NX, NY, NZ, NYH, Y1_FLAT
    ap = argparse.ArgumentParser()
    ap.add_argument("output", nargs="?", default=None)
    ap.add_argument("--nx", type=int, default=160)
    ap.add_argument("--ny", type=int, default=100)
    ap.add_argument("--nz", type=int, default=60)
    ap.add_argument("--y1", type=float, default=1.5e-3, help="平坦部第一層厚 [h 単位]")
    ap.add_argument("--shift", action="store_true",
                    help="丘を半周期ずらす (山頂を x=4.5h に、周期継ぎ目を平坦床に置く切り分け用)")
    ap.add_argument("--noclip", action="store_true",
                    help="山頂の min(28,poly) クリップを外す (C1 折れ除去; 公式形状から最大 +0.0047mm)")
    ap.add_argument("--xcluster", type=float, default=0.0,
                    help="山頂 (x=0/Lx) への x クラスタリング強度 a (0=一様)。x(ξ)=ξ−(a/2π)sin(2πξ) で"
                         " 山頂 Δx が一様比 (1−a) 倍に細かく、中央が (1+a) 倍に粗くなる (例 0.6 → 2.5:1)")
    args = ap.parse_args()
    global NOCLIP
    NOCLIP = args.noclip
    NX, NY, NZ = args.nx, args.ny, args.nz
    assert NY % 2 == 0, "ny は偶数 (両側ストレッチ)"
    NYH = NY // 2
    Y1_FLAT = args.y1 * H
    dst = args.output or f"hill_des_{NX}x{NY}x{NZ}.msh"

    # --- 多項式の健全性チェック (区間接続の連続性・端点値) ---
    for (lo, hi, _), (lo2, _, _) in zip(_SEG[:-1], _SEG[1:]):
        assert abs(hi - lo2) < 1e-12
        a = hill_half_mm(hi - 1e-9)
        b = hill_half_mm(hi + 1e-9)
        assert abs(a - b) < 0.05, (hi, a, b)  # 接続段差 < 0.05mm (=1.8e-3 h)
    assert abs(hill_half_mm(0.0) - H_MM) < (0.01 if NOCLIP else 1e-9)
    assert hill_half_mm(X_FOOT_MM - 1e-9) < 0.05
    # 周期対の厳密一致 (両端は同一コードパス)
    assert wall_y(0.0) == wall_y(LX)

    # --- y 正規化分布 (全カラム共通): 両側幾何ストレッチ ---
    q1 = Y1_FLAT / (0.5 * (LY_TOP - 0.0))      # 平坦部カラムの片側正規化第一層
    r = solve_progression(q1, NYH)
    g = (r ** np.arange(NYH + 1) - 1.0) / (r**NYH - 1.0)   # 0..1 片側
    frac = np.empty(NY + 1)
    frac[: NYH + 1] = 0.5 * g
    frac[NYH:] = 1.0 - 0.5 * g[::-1]
    assert abs(frac[0]) < 1e-15 and abs(frac[-1] - 1.0) < 1e-14
    assert np.all(np.diff(frac) > 0.0)

    xi = np.linspace(0.0, 1.0, NX + 1)
    if args.xcluster > 0.0:
        a = args.xcluster
        assert a < 1.0
        x = LX * (xi - (a / (2.0 * np.pi)) * np.sin(2.0 * np.pi * xi))
    else:
        x = LX * xi
    z = np.linspace(0.0, LZ, NZ + 1)
    if args.shift:
        yw = np.array([wall_y((xi + 0.5 * LX) % LX) for xi in x])
    else:
        yw = np.array([wall_y(xi) for xi in x])

    # --- 統計出力 (投入判断用) ---
    dx = LX / NX
    dz = LZ / NZ
    y1_flat = frac[1] * (LY_TOP - 0.0)
    y1_crest = frac[1] * (LY_TOP - H)
    dy_mid = (frac[NYH] - frac[NYH - 1]) * LY_TOP
    print(f"progression r = {r:.6f}")
    print(f"dx = {dx/H:.4f} h, dz = {dz/H:.4f} h")
    print(f"y1(flat) = {y1_flat/H:.3e} h, y1(crest) = {y1_crest/H:.3e} h")
    print(f"dy(mid,flat column) = {dy_mid/H:.4f} h")
    print(f"AR(first layer) ~ {dx/y1_crest:.0f}")
    print(f"cells = {NX*NY*NZ}, nodes = {(NX+1)*(NY+1)*(NZ+1)}")

    # --- 節点座標: idx(i,j,k) = ((i*(NY+1))+j)*(NZ+1)+k ---
    ny1, nz1 = NY + 1, NZ + 1
    npts = (NX + 1) * ny1 * nz1
    pts = np.empty((npts, 3))
    for i in range(NX + 1):
        ycol = yw[i] + frac * (LY_TOP - yw[i])
        base = i * ny1 * nz1
        for j in range(ny1):
            b2 = base + j * nz1
            pts[b2:b2 + nz1, 0] = x[i]
            pts[b2:b2 + nz1, 1] = ycol[j]
            pts[b2:b2 + nz1, 2] = z

    def idx(i, j, k):
        return (i * ny1 + j) * nz1 + k + 1   # gmsh は 1 始まり

    # --- hex 要素 (VTK/gmsh 線形 hexahedron 順) ---
    hexes = []
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                hexes.append((idx(i, j, k), idx(i + 1, j, k),
                              idx(i + 1, j + 1, k), idx(i, j + 1, k),
                              idx(i, j, k + 1), idx(i + 1, j, k + 1),
                              idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1)))

    # --- 境界 quad ---
    quads = {}
    quads["xlo"] = [(idx(0, j, k), idx(0, j + 1, k), idx(0, j + 1, k + 1), idx(0, j, k + 1))
                    for j in range(NY) for k in range(NZ)]
    quads["xhi"] = [(idx(NX, j, k), idx(NX, j + 1, k), idx(NX, j + 1, k + 1), idx(NX, j, k + 1))
                    for j in range(NY) for k in range(NZ)]
    quads["ylo"] = [(idx(i, 0, k), idx(i + 1, 0, k), idx(i + 1, 0, k + 1), idx(i, 0, k + 1))
                    for i in range(NX) for k in range(NZ)]
    quads["yhi"] = [(idx(i, NY, k), idx(i + 1, NY, k), idx(i + 1, NY, k + 1), idx(i, NY, k + 1))
                    for i in range(NX) for k in range(NZ)]
    quads["zlo"] = [(idx(i, j, 0), idx(i + 1, j, 0), idx(i + 1, j + 1, 0), idx(i, j + 1, 0))
                    for i in range(NX) for j in range(NY)]
    quads["zhi"] = [(idx(i, j, NZ), idx(i + 1, j, NZ), idx(i + 1, j + 1, NZ), idx(i, j + 1, NZ))
                    for i in range(NX) for j in range(NY)]

    PHYS = {"xlo": 1, "xhi": 2, "ylo": 3, "yhi": 4, "zlo": 5, "zhi": 6, "fluid": 7}
    names = ["xlo", "xhi", "ylo", "yhi", "zlo", "zhi"]

    nsurf = sum(len(v) for v in quads.values())
    with open(dst, "w") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n%d\n" % (len(names) + 1))
        for nm in names:
            f.write('2 %d "%s"\n' % (PHYS[nm], nm))
        f.write('3 %d "fluid"\n' % PHYS["fluid"])
        f.write("$EndPhysicalNames\n")
        f.write("$Entities\n0 0 %d 1\n" % len(names))
        for nm in names:
            f.write("%d 0 0 0 0 0 0 1 %d 0\n" % (PHYS[nm], PHYS[nm]))
        f.write("1 0 0 0 0 0 0 1 %d 0\n" % PHYS["fluid"])
        f.write("$EndEntities\n")
        f.write("$Nodes\n1 %d 1 %d\n" % (npts, npts))
        f.write("3 1 0 %d\n" % npts)
        f.write("".join("%d\n" % i for i in range(1, npts + 1)))
        f.write("".join("%.16g %.16g %.16g\n" % (p[0], p[1], p[2]) for p in pts))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d %d 1 %d\n" % (len(names) + 1, nsurf + len(hexes), nsurf + len(hexes)))
        tag = 1
        for nm in names:
            f.write("2 %d 3 %d\n" % (PHYS[nm], len(quads[nm])))
            for q in quads[nm]:
                f.write("%d %d %d %d %d\n" % (tag, q[0], q[1], q[2], q[3]))
                tag += 1
        f.write("3 1 5 %d\n" % len(hexes))
        for hx in hexes:
            f.write("%d %d %d %d %d %d %d %d %d\n" % ((tag,) + hx))
            tag += 1
        f.write("$EndElements\n")
    print("wrote", dst)


if __name__ == "__main__":
    main()
