#!/usr/bin/env python3
"""Wyslouzil JCP 113,7317 (2000) Fig.3 ノズルの .geo を生成する (修正ジオメトリ)。

旧 `nozzle_H.geo` は別 Wyslouzil 論文の "H" ノズル (スロート半高 2.25mm, A/A*=1.69) で
Fig.3 実験ノズルと不一致だった。本ファイルは Fig.3 論文記載の諸元で作り直す:

  - スロート: 高さ 5 mm (半高 2.5 mm) × 幅 12.7 mm
  - 収束部: 38 mm
  - 発散部: 95 mm, 直線壁, A/A*_exit = 1.58  (壁間 ~1.75deg ≈ 1.8deg)
  - 押し出し幅 z = 12.7 mm

usage:
  make_nozzle_fig3.py 2d      -> nozzle_fig3_2d.geo   (1 層 pseudo-2D, frontback=別physID)
  make_nozzle_fig3.py 3d      -> nozzle_fig3_3d.geo   (z 36 層 壁解像クラスタ)
単位 cm (Mesh.ScalingFactor 0.01)。throat を x=0 に置く。
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# --- 幾何パラメータ [cm] ---  解釈(A): throat 全高 5mm, 壁間角 1.8° (片側 0.9°)
Y_TH   = 0.25          # throat 半高 (2.5 mm = 全高 5mm)
LC     = 3.8           # 収束部長 (38 mm)
LD     = 9.5           # 発散部長 (95 mm)
HALF_ANGLE_DEG = 0.9   # 発散直線壁の片側角 (壁間 1.8° = 論文値)
Y_IN   = 1.27          # 入口半高 (full 25.4 mm = flow straightener)
SLOPE  = np.tan(np.radians(HALF_ANGLE_DEG))  # 発散直線壁の傾き
Y_EXIT = Y_TH + SLOPE * LD       # 出口半高 (~0.399 cm)
AR     = Y_EXIT / Y_TH           # 出口/スロート面積比 (~1.597 ≈ 論文 1.58)

# --- メッシュ分解能 ---
NX_CONV = 70           # 収束部 streamwise 分割
NX_DIV  = 230          # 発散部 streamwise 分割 (凝縮域を細かく)
NY      = 120          # 壁法線方向 (BL クラスタ)
BUMP    = 0.004        # 壁集中 (両端集中)

# --- z 押し出し (3D 用) ---
SPAN   = 1.27          # 12.7 mm
SPAN_M = SPAN * 0.01
FIRST_LAYER_M = 2.0e-6
NHALF  = 18


def y_upper(x):
    """上壁半高 y(x)。throat x=0。"""
    if x >= 0.0:
        return Y_TH + SLOPE * x                      # 発散: 直線
    t = -x / LC                                      # 収束: smoothstep (throat で slope0)
    return Y_TH + (Y_IN - Y_TH) * (3*t*t - 2*t*t*t)


def contour_points():
    xs_c = np.linspace(-LC, 0.0, NX_CONV+1)[:-1]     # 収束 (throat 含めない)
    xs_d = np.linspace(0.0, LD, NX_DIV+1)            # 発散 (throat 含む)
    xs = np.concatenate([xs_c, xs_d])
    ys = np.array([y_upper(x) for x in xs])
    return xs, ys


def z_layer_heights():
    a = FIRST_LAYER_M / SPAN_M
    def ssum(r):
        return a*NHALF if abs(r-1) < 1e-12 else a*(r**NHALF - 1)/(r-1)
    lo, hi = 1.0, 10.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if ssum(mid) < 0.5: lo = mid
        else: hi = mid
    r = 0.5*(lo+hi)
    thick = [a*r**k for k in range(NHALF)]
    pos = [0.0]
    for t in thick: pos.append(pos[-1]+t)
    pos[-1] = 0.5
    h = pos[1:]
    for k in range(1, NHALF+1): h.append(1.0 - pos[NHALF-k])
    h[-1] = 1.0
    return r, h


HEADER = """Geometry.LineNumbers = 1;
Mesh.ScalingFactor = 0.01; // cm -> m

// Wyslouzil JCP 113,7317 (2000) Fig.3 nozzle (corrected geometry)
//   throat half-height {yth} cm, converging {lc} cm, diverging {ld} cm (linear),
//   A/A*_exit = {ar} (exit half {yex:.3f} cm), wall half-angle {ang:.3f} deg
"""


def build_geo(mode):
    xs, ys = contour_points()
    n = len(xs)
    lines = []
    lines.append(HEADER.format(yth=Y_TH, lc=LC, ld=LD, ar=AR, yex=Y_EXIT,
                               ang=np.degrees(np.arctan(SLOPE))))
    # points: upper wall (inlet->outlet), then lower wall (inlet->outlet)
    pid = 1
    up_ids = []
    for x, y in zip(xs, ys):
        lines.append(f"Point({pid}) = {{{x:.6f}, {y:.6f}, 0.0, 0.1}};")
        up_ids.append(pid); pid += 1
    lo_ids = []
    for x, y in zip(xs, ys):
        lines.append(f"Point({pid}) = {{{x:.6f}, {-y:.6f}, 0.0, 0.1}};")
        lo_ids.append(pid); pid += 1

    p_in_top, p_out_top = up_ids[0], up_ids[-1]
    p_in_bot, p_out_bot = lo_ids[0], lo_ids[-1]

    # curves
    l_top = 1; lines.append(f"Spline({l_top}) = {{{','.join(map(str, up_ids))}}};")
    l_bot = 2; lines.append(f"Spline({l_bot}) = {{{','.join(map(str, lo_ids))}}};")
    l_out = 3; lines.append(f"Line({l_out}) = {{{p_out_top}, {p_out_bot}}};")
    l_in  = 4; lines.append(f"Line({l_in}) = {{{p_in_top}, {p_in_bot}}};")

    # transfinite
    nx = n
    lines.append(f"Transfinite Line {{{l_top}, {l_bot}}} = {nx};")
    lines.append(f"Transfinite Line {{{l_in}, {l_out}}} = {NY} Using Bump {BUMP};")

    # surface  (loop: top, out, -bot, -in)
    lines.append(f"Curve Loop(1) = {{{l_top}, {l_out}, -{l_bot}, -{l_in}}};")
    lines.append("Plane Surface(1) = {1};")
    lines.append(f"Transfinite Surface {{1}} = {{{p_in_top}, {p_out_top}, {p_out_bot}, {p_in_bot}}};")
    lines.append("Recombine Surface(1);")

    # extrude
    if mode == "2d":
        zspec = "Layers{1}"
        cmt = "// pseudo-2D (1 層), frontback=slip 用に別 physID"
    else:
        r, h = z_layer_heights()
        nlay = len(h)
        nlist = ", ".join(["1"]*nlay)
        hlist = ", ".join(f"{v:.10f}" for v in h)
        zspec = f"Layers{{ {{{nlist}}}, {{{hlist}}} }}"
        first_um = h[0]*SPAN_M*1e6
        cmt = f"// 3D 壁解像: {nlay} 層 対称集中 (r={r:.4f}), 第一層 {first_um:.2f} µm"
    lines.append("")
    lines.append(cmt)
    lines.append(f"e[] = Extrude {{0, 0, {SPAN}}} {{ Surface{{1}}; {zspec}; Recombine; }};")
    # 1 面 Extrude の返り値: e[0]=back(top), e[1]=volume, e[2..5]=lateral(loop 順)
    #   loop = {l_top, l_out, -l_bot, -l_in} -> e[2]=upper wall, e[3]=outlet, e[4]=lower wall, e[5]=inlet
    lines.append('Physical Surface("inlet", 1)  = {e[5]};')
    lines.append('Physical Surface("outlet", 2) = {e[3]};')
    lines.append('Physical Surface("wall", 3)   = {e[2], e[4]};')
    lines.append('Physical Surface("frontback", 4) = {1, e[0]};')
    lines.append('Physical Volume("fluid", 5) = {e[1]};')
    return "\n".join(lines) + "\n"


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "2d"
    assert mode in ("2d", "3d")
    geo = build_geo(mode)
    out = os.path.join(HERE, f"nozzle_fig3_{mode}.geo")
    with open(out, "w") as f:
        f.write(geo)
    print(f"wrote {out}")
    print(f"  throat half={Y_TH}cm  exit half={Y_EXIT:.3f}cm  A/A*={AR}  "
          f"half-angle={np.degrees(np.arctan(SLOPE)):.3f}deg (total {2*np.degrees(np.arctan(SLOPE)):.2f})")
    xs, ys = contour_points()
    print(f"  streamwise pts={len(xs)}  NY={NY}  x range [{xs.min():.2f},{xs.max():.2f}]cm")
    if mode == "3d":
        r, h = z_layer_heights()
        print(f"  z layers={len(h)} r={r:.4f} first={h[0]*SPAN_M*1e6:.2f}µm")


if __name__ == "__main__":
    main()
