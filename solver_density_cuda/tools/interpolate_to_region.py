"""
forge H5 region 補間スクリプト

既存の解析結果 H5 から、プリューム延長メッシュなど新しいメッシュの
特定の流体領域 (physId) へ物理量を最近傍補間 (nearest-neighbor) する。

補間するフィールド:
    保存量: ro, roUx, roUy, roUz, roe
    乱流量: roK, roOmega (ソースに存在する場合のみ)
    壁距離: wall_dist (ソースに存在する場合のみ)

使い方:
    python3 interpolate_to_region.py \\
        --src-mesh  <source_mesh.h5>   \\
        --src-result <res_10000.h5>    \\
        --dst       <target.h5>        \\
        --region    8                  \\
        [--k-neighbors 1]

引数:
    --src-mesh    ソースのメッシュ H5 (/CELLS/centCoords を含む)
    --src-result  ソースの結果 H5 (/VALUE/ に保存量・乱流量)
    --dst         パッチ先の H5 (/CELLS/regionId + /VALUE/ を含む)
    --region      対象の physId (デフォルト: 8)
    --k-neighbors 補間に使う最近傍数 (1=最近傍, >1=逆距離加重, デフォルト: 1)
"""

import argparse
import sys
import numpy as np
import h5py
from scipy.spatial import cKDTree

INTERP_FIELDS = ['ro', 'roUx', 'roUy', 'roUz', 'roe']
OPTIONAL_FIELDS = ['roK', 'roOmega', 'wall_dist']


def load_src_centroids(mesh_h5_path: str) -> np.ndarray:
    with h5py.File(mesh_h5_path, 'r') as f:
        cc = f['/CELLS/centCoords'][:]
    return cc.reshape(-1, 3)


def load_src_fields(result_h5_path: str) -> dict:
    fields = {}
    with h5py.File(result_h5_path, 'r') as f:
        val = f['/VALUE']
        for key in INTERP_FIELDS + OPTIONAL_FIELDS:
            if key in val:
                fields[key] = val[key][:]
    return fields


def idw_weights(distances: np.ndarray) -> np.ndarray:
    """逆距離加重 (距離ゼロは完全一致扱い)。"""
    eps = 1e-30
    w = 1.0 / np.maximum(distances, eps)
    return w / w.sum(axis=1, keepdims=True)


def interpolate(
    src_centroids: np.ndarray,
    src_fields: dict,
    dst_centroids: np.ndarray,
    dst_mask: np.ndarray,
    k: int,
) -> dict:
    tree = cKDTree(src_centroids)
    dst_pts = dst_centroids[dst_mask]

    if k == 1:
        _, idx = tree.query(dst_pts, k=1)
        result = {}
        for key, arr in src_fields.items():
            result[key] = arr[idx]
    else:
        distances, idx = tree.query(dst_pts, k=k)
        weights = idw_weights(distances)            # (n_dst, k)
        result = {}
        for key, arr in src_fields.items():
            vals = arr[idx]                         # (n_dst, k)
            result[key] = (weights * vals).sum(axis=1)

    return result


def patch_dst(dst_h5_path: str, dst_mask: np.ndarray, interp_vals: dict) -> None:
    with h5py.File(dst_h5_path, 'r+') as f:
        val = f['/VALUE']
        for key, vals in interp_vals.items():
            if key not in val:
                print(f"  スキップ: /VALUE/{key} が dst に存在しません")
                continue
            arr = val[key][:]
            arr[dst_mask] = vals.astype(arr.dtype)
            val[key][:] = arr


def main() -> None:
    parser = argparse.ArgumentParser(
        description='forge H5: 既存結果から特定 region へ nearest-neighbor 補間'
    )
    parser.add_argument('--src-mesh',    required=True, help='ソース メッシュ H5')
    parser.add_argument('--src-result',  required=True, help='ソース 結果 H5')
    parser.add_argument('--dst',         required=True, help='パッチ先 H5')
    parser.add_argument('--region',      type=int, default=8, help='対象 physId (デフォルト: 8)')
    parser.add_argument('--k-neighbors', type=int, default=1,
                        help='補間に使う最近傍数 (1=NN, >1=逆距離加重, デフォルト: 1)')
    args = parser.parse_args()

    print(f"ソース メッシュ : {args.src_mesh}")
    print(f"ソース 結果    : {args.src_result}")
    print(f"パッチ先       : {args.dst}")
    print(f"対象 region    : physId={args.region}")
    print(f"補間方式       : k={args.k_neighbors} ({'最近傍' if args.k_neighbors == 1 else '逆距離加重'})")
    print()

    print("1. ソース セントロイド読み込み...")
    src_centroids = load_src_centroids(args.src_mesh)
    print(f"   {len(src_centroids)} セル")

    print("2. ソース フィールド読み込み...")
    src_fields = load_src_fields(args.src_result)
    print(f"   フィールド: {list(src_fields.keys())}")

    print("3. パッチ先メッシュ読み込み...")
    with h5py.File(args.dst, 'r') as f:
        region_ids  = f['/CELLS/regionId'][:]
        dst_cc      = f['/CELLS/centCoords'][:]
    dst_centroids = dst_cc.reshape(-1, 3)
    dst_mask = region_ids == args.region
    n_target = int(dst_mask.sum())
    print(f"   対象セル数: {n_target} / {len(region_ids)}")

    if n_target == 0:
        print(f"ERROR: physId={args.region} のセルが見つかりません")
        sys.exit(1)

    print("4. 補間中...")
    interp_vals = interpolate(
        src_centroids, src_fields, dst_centroids, dst_mask, args.k_neighbors
    )

    print("5. パッチ先 H5 へ書き込み...")
    patch_dst(args.dst, dst_mask, interp_vals)

    print()
    print("=== 補間後の値サマリ (physId={}) ===".format(args.region))
    with h5py.File(args.dst, 'r') as f:
        for key in INTERP_FIELDS + OPTIONAL_FIELDS:
            if f'/VALUE/{key}' in f:
                arr = f[f'/VALUE/{key}'][:][dst_mask]
                print(f"  {key:12s}: min={arr.min():+.4e}  max={arr.max():+.4e}")

    print("\n補間完了")


if __name__ == '__main__':
    main()
