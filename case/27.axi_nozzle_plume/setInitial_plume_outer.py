"""
axi_nozzle プリューム延長メッシュ用 初期値パッチスクリプト。

outer 領域 (physId=9: s_plume_in + s_plume_out) に大気静止条件を設定する。
ノズル内部 (physId=8: s_conv + s_div) はそのまま (convertGmshToForge のベース IC を維持)。

使い方 (run ディレクトリから):
    python ../../../../case/23.axi_nozzle/setInitial_plume_outer.py axi_nozzle.h5

または直接 H5 パスを指定:
    python /path/to/setInitial_plume_outer.py /path/to/axi_nozzle.h5
"""

import sys
import os

# 共通ユーティリティのパスを追加
_repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_repo_root, "solver_density_cuda", "tools"))

from h5_patch_initial import patch_initial_by_region

# ---- ガス定数 ----
GAMMA = 1.4
CP    = 1004.5  # [J/(kg·K)]

# ---- outer 領域初期条件 ----
P_OUT  = 1000.0   # [Pa]
T_OUT  = 300.0    # [K]
R_gas  = CP * (GAMMA - 1.0) / GAMMA  # 287.06 J/(kg·K)
RO_OUT = P_OUT / (R_gas * T_OUT)

if __name__ == "__main__":
    h5_path = sys.argv[1] if len(sys.argv) > 1 else "axi_nozzle.h5"

    print(f"初期値パッチ: {h5_path}")
    print(f"  outer 領域 (physId=9): p={P_OUT} Pa, T={T_OUT} K, ro={RO_OUT:.5f} kg/m^3, Ux=5 m/s")

    patch_initial_by_region(
        h5_path,
        regions=[
            {
                "physId": 9,   # "outer": s_plume_in + s_plume_out
                "ro": RO_OUT,
                "p":  P_OUT,
                "Ux": 5.0,
                "Uy": 0.0,
                "Uz": 0.0,
            },
        ],
        gamma=GAMMA,
        cp=CP,
    )
