#!/usr/bin/env python3
"""W4 実行 (暫定近似): run_0043 の CFD 場から D 発データ線を抽出し、M_d=4 を
目標とした march を実行する。

**注意 (2026-08-15, 外部レビューで確認)**: `march_wall_from_dataline` の
$J^-_W=\\nu(M_{\\rm target})$ 一定方式は軸対称からの真の逆 MOC ではなく、
平面 simple-wave を模した暫定近似 (`moc_cancel.py` docstring 参照)。ここでの
実行結果は「近似が実データでどう振る舞うか」の記録であり、最終設計として
使える壁座標ではない。$D$ の局所 $M$ (~2.26) と $M_d$=4 の $J^-$ ギャップが
大きく (~20°)、march が発散するケースを確認している。
"""
import json
import sys
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE.parents[1] / "design"))

from forge_design.geometry.moc_cancel import (  # noqa: E402
    extract_dataline_cfd, march_wall_from_dataline)

RUN = CASE / "run_0043_conecheck_staged"
info = json.loads((RUN / "prepare_info.json").read_text())
GAMMA, SCALE, M_D = 1.4, 0.01, 4.0
x_D, r_D, theta_D_deg = info["x_D"], info["r_D"], info["theta_D_deg"]

print(f"D=(x={x_D:.4f}, r={r_D:.4f}), theta_D={theta_D_deg:.4f}deg, M_target={M_D}")

dl = extract_dataline_cfd(RUN, x_D, r_D, gamma=GAMMA, scale=SCALE, n_line=31,
                          wall_theta_geom=np.deg2rad(theta_D_deg))
print(f"dataline: x_axis={dl.x_axis:.4f}, P0 diag: {dl.p0_diag}")
print(f"dataline theta range: {np.rad2deg(dl.theta[0]):.4f} .. "
      f"{np.rad2deg(dl.theta[-1]):.4f} deg")
print(f"dataline M range: {dl.M[0]:.4f} .. {dl.M[-1]:.4f}")

summary = {"x_D": x_D, "r_D": r_D, "theta_D_deg": theta_D_deg, "M_target": M_D,
          "dataline_p0_diag": dl.p0_diag}
try:
    res = march_wall_from_dataline(dl, M_target=M_D, theta_tol_deg=0.5,
                                   eps_M_rel=0.02, max_iter=2000)
    wall = res["wall"]
    print(f"terminated_by={res['terminated_by']} n_iter={res['n_iter']}")
    print(f"wall span: x=[{wall[0,0]:.4f}, {wall[-1,0]:.4f}] "
          f"r=[{wall[0,1]:.4f}, {wall[-1,1]:.4f}]")
    print(f"final: theta={np.rad2deg(wall[-1,2]):.4f}deg M={wall[-1,3]:.4f}")
    dx = np.diff(wall[:, 0])
    dr = np.diff(wall[:, 1])
    th_avg = 0.5 * (wall[:-1, 2] + wall[1:, 2])
    mask = np.abs(dx) > 1e-9
    stream_err = np.max(np.abs(dr[mask] / dx[mask] - np.tan(th_avg[mask])))
    print(f"streamline self-consistency max|Delta|: {stream_err:.3e}")
    print(f"theta monotonic non-increasing: {np.all(np.diff(wall[:,2]) <= 1e-6)}")
    np.savetxt(RUN / "w4_cancellation_wall.csv", wall, delimiter=",",
               header="x_rt,r_rt,theta_rad,M", comments="")
    summary.update({
        "status": "completed", "terminated_by": res["terminated_by"],
        "n_iter": res["n_iter"],
        "wall_x_span": [float(wall[0, 0]), float(wall[-1, 0])],
        "wall_r_span": [float(wall[0, 1]), float(wall[-1, 1])],
        "final_theta_deg": float(np.rad2deg(wall[-1, 2])),
        "final_M": float(wall[-1, 3]), "streamline_max_delta": float(stream_err)})
except RuntimeError as exc:
    # 既知の限界 (モジュール docstring 参照): D の局所 M と M_target の J⁻
    # ギャップが大きいと march が発散する。ここでは記録として残し、クラッシュ
    # させない (この近似の失敗モードそのものが調査結果の一部)。
    print(f"march 失敗 (既知の限界 — moc_cancel.py docstring 参照): {exc}")
    summary.update({"status": "failed", "failure_reason": str(exc)})

(RUN / "w4_cancellation_summary.json").write_text(json.dumps(summary, indent=1))
print("saved:", RUN / "w4_cancellation_summary.json")
