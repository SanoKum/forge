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

直交品質改善 (--wall-block):
    単一ブロックの y 分布 ycol=yw(x)+frac*(LY_TOP-yw(x)) は、どの高さでも局所勾配が
    (1-frac)*yw'(x) になる。壁近傍 (y+~1 用に frac≈0 が続く領域) では (1-frac)≈1 なので
    丘の勾配がほぼそのまま何層も残り非直交になる。--wall-block は壁直上に
    「中間断面 y_mid(x) = wallBlockBase*h + wallBlockRelax*yw(x)」を挟む 2 ブロックにし、
    ブロック A (壁→中間断面) の上端勾配を wallBlockRelax 倍に圧縮する。
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
    """(r-1)/(r^n-1) = q1 となる公比 r を二分法で解く (q1 = y1/半高さ)。

    val(r)=(r-1)/(r^n-1) は r>0 で単調減少、r=1 (一様分割) で val=1/n。
    q1<=1/n (第一層が平均より小=クラスタ) は r>=1 側、q1>1/n (第一層が平均より大=
    継ぎ目でのセルサイズ継承など) は r<1 側で解が存在する — 片側 [1,2] 固定だと
    q1>1/n のとき解なしで二分法が静かに r≈1 (縮退・要求した q1 を無視) へ収束するバグが
    あった (wall-block 継ぎ目のセルサイズ不整合で発覚)。"""
    if q1 <= 1.0 / n:
        lo, hi = 1.0 + 1e-12, 2.0
        while (hi - 1.0) / (hi**n - 1.0) > q1 and hi < 1.0e6:
            hi *= 2.0
    else:
        lo, hi = 1.0e-9, 1.0 - 1e-12
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
    ap.add_argument("--wall-block", action="store_true",
                    help="壁直上に緩勾配の中間断面を挟む 2 ブロック y 分布 (直交品質改善)")
    ap.add_argument("--wall-block-relax", type=float, default=0.3,
                    help="中間断面の勾配緩和係数 (y_mid'=relax*yw', 既定 0.3)")
    ap.add_argument("--wall-block-base", type=float, default=1.0,
                    help="中間断面の基準高さ [h 単位] (y_mid=base*h+relax*yw(x)、既定 1.0)")
    ap.add_argument("--ny-low", type=int, default=0,
                    help="ブロック A (壁→中間断面) の分割数 (--wall-block 時、既定 NY の 30%% 目安)")
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

    # --- 2 ブロック y 分布 (--wall-block): 壁直上に緩勾配の中間断面を挟む ---
    if args.wall_block:
        RELAX = args.wall_block_relax
        Y_BASE = args.wall_block_base * H
        assert 0.0 < RELAX < 1.0
        y_mid = Y_BASE + RELAX * yw                      # 中間断面 (勾配 = RELAX*yw')
        thickA = y_mid - yw
        thickB = LY_TOP - y_mid
        margin = 0.05 * H
        assert thickA.min() > margin, (
            f"wall-block-base too small: block A thickness min={thickA.min()/H:.3f}h "
            f"(need > {margin/H:.2f}h; increase --wall-block-base or decrease --wall-block-relax)")
        assert thickB.min() > margin, "wall-block-base too large: block B collapses near top wall"

        NY_LOW = args.ny_low or max(10, round(NY * 0.3))
        NY_UP = NY - NY_LOW
        if NY_UP % 2 != 0:          # NY_UP は両側ストレッチ用に偶数が必要 → ブロック A 側で 1 吸収
            NY_LOW += 1
            NY_UP -= 1
        assert NY_UP >= 2 and NY_LOW >= 2, \
            f"NY_LOW={NY_LOW} NY_UP={NY_UP} too small (adjust --ny-low or --ny)"
        NY_UPH = NY_UP // 2

        # ブロック A: 片側幾何ストレッチ (壁で細かく、中間断面で粗く)。壁第一層の絶対厚みは
        # 元設計と同じ目標 y1_target(x)=Y1_FLAT*(LY_TOP-yw(x))/LY_TOP (元の y1_flat/y1_crest を再現)
        # を狙う。ブロック A の厚み thickA(x) は元の局所列高さ (LY_TOP-yw) よりずっと薄いため、
        # q1 (正規化第一層比) は列ごとに解き直す必要がある (定数 q1 を流用すると壁第一層が
        # 桁違いに薄くなり AR が爆発する — 実測で確認済み)。
        y1_target = Y1_FLAT * (LY_TOP - yw) / LY_TOP
        q1A_arr = y1_target / thickA
        fracA_cols = np.empty((NX + 1, NY_LOW + 1))
        for i in range(NX + 1):
            rAi = solve_progression(q1A_arr[i], NY_LOW)
            fracA_cols[i] = (rAi ** np.arange(NY_LOW + 1) - 1.0) / (rAi**NY_LOW - 1.0)

        # ブロック B: 非対称両側ストレッチ。中間断面側の第一層厚をブロック A の最終セル厚に
        # 一致させ (継ぎ目のセルサイズ連続性)、上壁側の第一層厚は元設計と同じ y1_target を狙う。
        # 単純に元の対称ストレッチ (q1 一定) を流用すると、中間断面直上がいきなり Y1_FLAT 相当の
        # 薄い層に戻ってしまい (ブロック A 最終セルの ~300 倍のギャップ)、境界に歪なセルができる
        # (実測で確認・要修正点)。
        dyA_last = thickA * (1.0 - fracA_cols[:, -2])         # ブロック A 最終セル厚 (列ごと)
        halfB = 0.5 * thickB
        q1_lo = dyA_last / halfB
        q1_hi = y1_target / halfB
        assert q1_lo.max() < 1.0 and q1_hi.max() < 1.0, "wall-block: NY_UP または RELAX を見直す"
        fracB = np.empty((NX + 1, NY_UP + 1))
        for i in range(NX + 1):
            r_lo = solve_progression(q1_lo[i], NY_UPH)
            g_lo = (r_lo ** np.arange(NY_UPH + 1) - 1.0) / (r_lo**NY_UPH - 1.0)
            r_hi = solve_progression(q1_hi[i], NY_UPH)
            g_hi = (r_hi ** np.arange(NY_UPH + 1) - 1.0) / (r_hi**NY_UPH - 1.0)
            fracB[i, : NY_UPH + 1] = 0.5 * g_lo
            fracB[i, NY_UPH:] = 1.0 - 0.5 * g_hi[::-1]

        slope_ratio = np.abs(np.gradient(y_mid, x)).max() / max(np.abs(np.gradient(yw, x)).max(), 1e-30)
        print(f"[wall-block] NY_LOW={NY_LOW} NY_UP={NY_UP}  RELAX={RELAX} Y_BASE={Y_BASE/H:.2f}h")
        print(f"[wall-block] thickA: {thickA.min()/H:.3f}-{thickA.max()/H:.3f}h  "
              f"thickB: {thickB.min()/H:.3f}-{thickB.max()/H:.3f}h")
        print(f"[wall-block] max|dy_mid/dx| / max|dyw/dx| = {slope_ratio:.3f} (目標 ~{RELAX})")

    # --- 節点座標: idx(i,j,k) = ((i*(NY+1))+j)*(NZ+1)+k ---
    ny1, nz1 = NY + 1, NZ + 1
    npts = (NX + 1) * ny1 * nz1
    pts = np.empty((npts, 3))
    for i in range(NX + 1):
        if args.wall_block:
            ycolA = yw[i] + thickA[i] * fracA_cols[i, :-1]  # 壁→中間断面 (中間断面点は B 側で重複生成)
            ycolB = y_mid[i] + thickB[i] * fracB[i]          # 中間断面→上壁
            ycol = np.concatenate([ycolA, ycolB])
        else:
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
