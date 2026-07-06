#!/usr/bin/env python3
"""forge 入力 H5 メッシュの品質指標 (アスペクト比・スキューネス) を計算・判定する。

使い方:
    python3 check_mesh_quality.py mesh.h5 [--ar-max 1000] [--skew-max 0.9] [--mode auto|2d|3d]

- アスペクト比 (AR): 各セルの 最長辺/最短辺。薄い境界層セルの監視用。
- スキューネス (equiangle skew): 内角 θ から max[(θmax-ideal)/(180-ideal),(ideal-θmin)/ideal]
  (0=直交/正三角、1=退化)。四角形 ideal=90、三角形 ideal=60。
- 既定しきい値: AR ≤ 1000, skew ≤ 0.9 (本リポジトリの計算前メッシュ品質ルール、
  AGENTS.md / forge-calculation-workflow.md)。VERDICT を PASS/FAIL で返す。

2D / 3D の扱い (重要):
- **2D モード** (準2D・z 押し出し・軸対称): セルを x-y 平面に射影し多角形として測る (従来挙動)。
- **3D モード** (完全非構造 3D: tet/prism/hex/pyramid が z 方向にも構造を持つ): x-y 射影は
  無意味 (z 平行エッジが一点に潰れ AR=inf 等の偽値) なため、**各セルの面ごとに equiangle skew
  を測りその最大**を skew、**全辺の最長/最短**を AR とする (真の 3D 形状で判定)。
- `--mode auto` (既定): 3D セル型を含み z 平面が 3 枚以上なら 3D、それ以外は 2D と自動判定。

`solver_density_cuda/tools/res_h5_to_vtu.py` の parse_conne を流用して CONNE を読む。
"""
import sys, os, argparse
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from res_h5_to_vtu import parse_conne

# VTK cell type -> 面 (セル局所節点 index)。parse_conne の節点順 (gmsh/XDMF/VTK) 準拠。
# TETRA=10, HEXAHEDRON=12, WEDGE(prism)=13, PYRAMID=14
FACES_3D = {
    10: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
    12: [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
    13: [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    14: [(0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
}
VTK_3D = set(FACES_3D)


def poly_skew(pts):
    """多角形 (n,2 or n,3) の equiangle skew を返す (cell_metrics の skew 部と同一定義)。"""
    n = len(pts)
    angs = []
    for i in range(n):
        a = pts[(i - 1) % n] - pts[i]
        b = pts[(i + 1) % n] - pts[i]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if na == 0 or nb == 0:
            angs.append(0.0); continue
        c = np.clip(np.dot(a, b) / (na * nb), -1.0, 1.0)
        angs.append(np.degrees(np.arccos(c)))
    angs = np.array(angs)
    ideal = 90.0 if n == 4 else (60.0 if n == 3 else 180.0 * (n - 2) / n)
    return max((angs.max() - ideal) / (180.0 - ideal), (ideal - angs.min()) / ideal)


def cell_metrics(pts):
    """pts: (n,2) 多角形頂点 (CCW or CW)。(aspect_ratio, skew) を返す (2D モード用)。"""
    n = len(pts)
    edges = [np.linalg.norm(pts[(i + 1) % n] - pts[i]) for i in range(n)]
    emin, emax = min(edges), max(edges)
    ar = emax / emin if emin > 0 else float("inf")
    return ar, poly_skew(pts)


def cell_metrics_3d(pts, vtk):
    """pts: (m,3) セル節点。面ごとの equiangle skew の最大と、全辺の 最長/最短 を返す。"""
    faces = FACES_3D[vtk]
    skew = 0.0
    edges = []
    for fc in faces:
        fp = pts[list(fc)]
        skew = max(skew, poly_skew(fp))
        for i in range(len(fc)):
            edges.append(np.linalg.norm(fp[(i + 1) % len(fc)] - fp[i]))
    emin, emax = min(edges), max(edges)
    ar = emax / emin if emin > 0 else float("inf")
    return ar, skew


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--ar-max", type=float, default=1000.0)
    ap.add_argument("--skew-max", type=float, default=0.9)
    ap.add_argument("--mode", choices=["auto", "2d", "3d"], default="auto")
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        coord3 = np.array(f["MESH/COORD"]).reshape(-1, 3)
        conne = np.array(f["MESH/CONNE"])
        ncells = f["VALUE/ro"].shape[0] if "VALUE/ro" in f else None
    if ncells is None:
        # count from CONNE
        from res_h5_to_vtu import _count_cells
        ncells = _count_cells(conne)
    conn, offs, types = parse_conne(conne, ncells)

    # モード判定: 3D セル型を含み z 平面が 3 枚以上なら 3D (x-y 射影が破綻するため)
    if a.mode == "auto":
        nz = len(np.unique(np.round(coord3[:, 2], 6)))
        is3d = bool(np.isin(types, list(VTK_3D)).any()) and nz > 2
    else:
        is3d = (a.mode == "3d")
    print(f"mode: {'3D (per-face skew / edge AR)' if is3d else '2D (x-y projection)'}")

    coord = coord3 if is3d else coord3[:, :2]
    ars = np.empty(ncells); skews = np.empty(ncells)
    s = 0
    for i, o in enumerate(offs):
        pts = coord[conn[s:o]]; s = o
        ars[i], skews[i] = (cell_metrics_3d(pts, types[i]) if is3d else cell_metrics(pts))

    ar_bad = int((ars > a.ar_max).sum())
    sk_bad = int((skews > a.skew_max).sum())
    print(f"cells: {ncells}")
    print(f"aspect ratio : max={ars.max():.1f}  p99={np.percentile(ars,99):.1f}  "
          f"mean={ars.mean():.1f}  >{a.ar_max:.0f}: {ar_bad} ({100*ar_bad/ncells:.2f}%)")
    print(f"skewness     : max={skews.max():.3f}  p99={np.percentile(skews,99):.3f}  "
          f"mean={skews.mean():.3f}  >{a.skew_max:.2f}: {sk_bad} ({100*sk_bad/ncells:.2f}%)")
    ok = (ars.max() <= a.ar_max) and (skews.max() <= a.skew_max)
    # 少数の外れ値は許容するソフト判定も併記
    soft = (ar_bad <= 0.001 * ncells) and (sk_bad <= 0.001 * ncells)
    print(f"VERDICT: {'PASS' if ok else ('SOFT-PASS (<0.1% outliers)' if soft else 'FAIL')} "
          f"(AR<= {a.ar_max:.0f}, skew<= {a.skew_max:.2f})")
    sys.exit(0 if (ok or soft) else 1)


if __name__ == "__main__":
    main()
