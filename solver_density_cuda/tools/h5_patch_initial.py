"""
per-region 初期値を forge の H5 ファイルに書き込む共通ユーティリティ。

使い方:
    from h5_patch_initial import patch_initial_by_region

    patch_initial_by_region(
        "case.h5",
        regions=[
            {"physId": 9, "ro": 1.225, "p": 101325.0, "Ux": 0.0, "Uy": 0.0, "Uz": 0.0},
        ],
        gamma=1.4,
        cp=1004.5,
    )

regions の各要素:
    physId  : Gmsh の Physical Surface/Volume タグ番号
    ro      : 密度 [kg/m^3]
    p       : 静圧 [Pa]
    Ux/Uy/Uz: 速度成分 [m/s]
"""

import h5py
import numpy as np


def patch_initial_by_region(
    h5_path: str,
    regions: list,
    gamma: float,
    cp: float,
) -> None:
    with h5py.File(h5_path, "r+") as f:
        if "/CELLS/regionId" not in f:
            raise KeyError(
                f"{h5_path} に /CELLS/regionId が存在しません。"
                " convertGmshToForge を最新ビルドで再実行してください。"
            )

        region_ids = f["/CELLS/regionId"][:]
        ro_arr   = f["/VALUE/ro"][:]
        roUx_arr = f["/VALUE/roUx"][:]
        roUy_arr = f["/VALUE/roUy"][:]
        roUz_arr = f["/VALUE/roUz"][:]
        roe_arr  = f["/VALUE/roe"][:]

        for ric in regions:
            phys_id = ric["physId"]
            mask = (region_ids == phys_id)
            n_cells = int(mask.sum())
            if n_cells == 0:
                print(f"  警告: physId={phys_id} に対応するセルが見つかりません")
                continue

            ro = float(ric["ro"])
            p  = float(ric["p"])
            Ux = float(ric.get("Ux", 0.0))
            Uy = float(ric.get("Uy", 0.0))
            Uz = float(ric.get("Uz", 0.0))
            ke = 0.5 * (Ux**2 + Uy**2 + Uz**2)
            e  = p / (gamma - 1.0) + ro * ke

            ro_arr[mask]   = ro
            roUx_arr[mask] = ro * Ux
            roUy_arr[mask] = ro * Uy
            roUz_arr[mask] = ro * Uz
            roe_arr[mask]  = e

            print(f"  physId={phys_id}: {n_cells} セル → ro={ro:.4g}, p={p:.4g}, U=({Ux},{Uy},{Uz})")

        f["/VALUE/ro"][:]   = ro_arr
        f["/VALUE/roUx"][:] = roUx_arr
        f["/VALUE/roUy"][:] = roUy_arr
        f["/VALUE/roUz"][:] = roUz_arr
        f["/VALUE/roe"][:]  = roe_arr

    print(f"パッチ完了: {h5_path}")
