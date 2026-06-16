#!/usr/bin/env python3
"""nozzle_3d_wallres.geo を生成する。

in-plane (x,y) メッシュは nozzle_2d_visc.geo と完全同一 (同じ spline・nx・
縦線の Bump 0.02 → 輪郭壁 y1≈2µm を保持)。変えるのは z 押し出しのみ:

  - 押し出し幅 = 1.27 cm (= 12.7 mm, Wyslouzil 論文の断面幅)
  - z 方向を「対称両側ジオメトリック集中」で ~36 層に分割し、front/back 側壁の
    第一層が 2D の壁第一層 (~2µm) を模すようにする (実 no-slip 側壁を解像)。

Gmsh の Extrude は `Layers{ {n_1,..,n_k}, {h_1,..,h_k=1.0} }` で
非一様層 (n_i 要素 / 累積正規化高さ h_i) を指定できる。ここでは 1 層 = 1 要素
(n_i すべて 1) とし、両端に集中させた累積高さ列を算出して書き出す。

Physical Surface の割り当て (inlet/outlet/wall/frontback) は nozzle_2d_visc.geo の
e[] 捕捉マッピングをそのまま流用する (in-plane トポロジが同一なので妥当)。
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "nozzle_2d_visc.geo")
OUT = os.path.join(HERE, "nozzle_3d_wallres.geo")

# --- 押し出しパラメータ -------------------------------------------------
SPAN_CM = 1.27           # 押し出し幅 [cm] (= 12.7 mm, Wyslouzil)
SPAN_M = SPAN_CM * 0.01  # ScalingFactor 0.01 適用後の実幅 [m]
FIRST_LAYER_M = 2.0e-6   # 側壁第一層厚さ目標 [m] (2D の y1≈1.6〜9µm を模す)
NHALF = 18               # 片側の層数 (全 2*NHALF = 36 層)


def solve_ratio(a, nhalf, target=0.5):
    """初項 a・項数 nhalf の等比級数和が target になる公比 r を二分法で解く。"""
    def series_sum(r):
        if abs(r - 1.0) < 1e-12:
            return a * nhalf
        return a * (r ** nhalf - 1.0) / (r - 1.0)
    lo, hi = 1.0, 10.0
    # a*nhalf が target 未満であることを前提 (集中が必要)
    assert series_sum(lo) < target, "first layer too thick for uniform fit"
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if series_sum(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def build_heights():
    a = FIRST_LAYER_M / SPAN_M           # 正規化第一層厚さ
    r = solve_ratio(a, NHALF, target=0.5)
    # 前半 (z=0 側) の層厚と累積位置
    thick = [a * r ** k for k in range(NHALF)]
    pos = [0.0]
    for t in thick:
        pos.append(pos[-1] + t)
    pos[-1] = 0.5  # 数値誤差を吸収して中央へスナップ
    # 後半は前半を中央対称にミラー
    heights = pos[1:]  # 前半の層上面 (0 を除く) : NHALF 個, 末尾 0.5
    for k in range(1, NHALF + 1):
        heights.append(1.0 - pos[NHALF - k])
    heights[-1] = 1.0  # 末尾を厳密に 1.0
    return r, a, heights


def main():
    with open(TEMPLATE) as f:
        lines = f.readlines()
    # テンプレートの in-plane 部分 (Recombine Surface(4) まで) を切り出す。
    # 行 222 (Recombine Surface(4);) の直後で切る。
    cut = None
    for i, ln in enumerate(lines):
        if "Recombine Surface(4)" in ln:
            cut = i + 1
            break
    assert cut is not None, "Recombine Surface(4) not found in template"
    head = "".join(lines[:cut])

    r, a, heights = build_heights()
    nlayers = len(heights)
    n_list = ", ".join(["1"] * nlayers)
    h_list = ", ".join(f"{h:.10f}" for h in heights)

    # 検証用: 第一層・中央層の実厚 [µm]
    first_um = (heights[0]) * SPAN_M * 1e6
    mid_idx = nlayers // 2
    mid_um = (heights[mid_idx] - heights[mid_idx - 1]) * SPAN_M * 1e6

    extrude = f"""

// ====================================================================
// z 押し出し (壁解像クラスタ) — make_nozzle_3d_wallres.py が生成
//   幅 {SPAN_CM} cm (= {SPAN_CM*10:.1f} mm, Wyslouzil 断面幅)
//   {nlayers} 層 (対称両側集中, 公比 r={r:.4f})
//   第一層 ≈ {first_um:.2f} µm (2D 壁第一層を模す), 中央層 ≈ {mid_um:.1f} µm
// 返り値 e[] は入力面ごとに 6 要素: [top(back), volume, lat0, lat1, lat2, lat3]
// ====================================================================
e[] = Extrude {{0, 0, {SPAN_CM}}} {{
  Surface{{1}}; Surface{{2}}; Surface{{3}}; Surface{{4}};
  Layers{{ {{{n_list}}}, {{{h_list}}} }}; Recombine;
}};
// S1 loop {{l_in, l1, l_mid1, l_m1}}  -> [inlet, wall_up, interior, wall_lo]
// S2 loop {{l_mid1,-l_m2,-l_mid2,-l2}}-> [interior, wall_lo, interior, wall_up]
// S3 loop {{l_mid2,-l_m3,-l_mid3,-l3}}-> [interior, wall_lo, interior, wall_up]
// S4 loop {{l_mid3,-l_m4,-l5,-l4}}    -> [interior, wall_lo, outlet, wall_up]

Physical Surface("inlet", 1)  = {{e[2]}};
Physical Surface("outlet", 2) = {{e[22]}};
// 輪郭壁 (上下) = no-slip wall
Physical Surface("wall", 3)   = {{e[3], e[11], e[17], e[23],  e[5], e[9], e[15], e[21]}};
// front(z=0)=元面1..4, back(z=span)=e[0],e[6],e[12],e[18]。3D では bcond で no-slip 壁にする
Physical Surface("frontback", 4) = {{1, 2, 3, 4,  e[0], e[6], e[12], e[18]}};
Physical Volume("fluid", 5) = {{e[1], e[7], e[13], e[19]}};
"""

    with open(OUT, "w") as f:
        f.write(head)
        f.write(extrude)

    print(f"wrote {OUT}")
    print(f"  span      = {SPAN_CM} cm ({SPAN_M*1e3:.2f} mm)")
    print(f"  nlayers   = {nlayers} (NHALF={NHALF} per side)")
    print(f"  ratio r   = {r:.4f}")
    print(f"  1st layer = {first_um:.3f} µm (target {FIRST_LAYER_M*1e6:.1f} µm)")
    print(f"  mid layer = {mid_um:.2f} µm")


if __name__ == "__main__":
    main()
