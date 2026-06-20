#!/usr/bin/env python3
"""forge 入力 H5 メッシュの品質指標 (アスペクト比・スキューネス) を計算・判定する。

使い方:
    python3 check_mesh_quality.py mesh.h5 [--ar-max 1000] [--skew-max 0.9]

- アスペクト比 (AR): 各セルの 最長辺/最短辺。薄い境界層セルの監視用。
- スキューネス (equiangle skew): 各セルの内角 θ から max[(θmax-90)/90,(90-θmin)/90]
  (四角形基準。0=直交、1=退化)。三角形は 60 基準。
- 既定しきい値: AR ≤ 1000, skew ≤ 0.9 (本リポジトリの計算前メッシュ品質ルール、
  AGENTS.md / forge-calculation-workflow.md)。VERDICT を PASS/FAIL で返す。

`solver_density_cuda/tools/res_h5_to_vtu.py` の parse_conne を流用して CONNE を読む。
"""
import sys, os, argparse
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from res_h5_to_vtu import parse_conne


def cell_metrics(pts):
    """pts: (n,2) 多角形頂点 (CCW or CW)。(aspect_ratio, skew) を返す。"""
    n = len(pts)
    edges = [np.linalg.norm(pts[(i + 1) % n] - pts[i]) for i in range(n)]
    emin, emax = min(edges), max(edges)
    ar = emax / emin if emin > 0 else float("inf")
    # 内角
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
    skew = max((angs.max() - ideal) / (180.0 - ideal), (ideal - angs.min()) / ideal)
    return ar, skew


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--ar-max", type=float, default=1000.0)
    ap.add_argument("--skew-max", type=float, default=0.9)
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        coord = np.array(f["MESH/COORD"]).reshape(-1, 3)[:, :2]
        conne = np.array(f["MESH/CONNE"])
        ncells = f["VALUE/ro"].shape[0] if "VALUE/ro" in f else None
    if ncells is None:
        # count from CONNE
        from res_h5_to_vtu import _count_cells
        ncells = _count_cells(conne)
    conn, offs, types = parse_conne(conne, ncells)

    ars = np.empty(ncells); skews = np.empty(ncells)
    s = 0
    for i, o in enumerate(offs):
        pts = coord[conn[s:o]]; s = o
        ars[i], skews[i] = cell_metrics(pts)

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
