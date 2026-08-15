r"""RANS 場からの排除厚 $\delta^*(x)$ 抽出 (A12 = axis-Mach 粘性補正)。

計画: plans/active/tooling-nozzle-axismach-viscous-deltastar.md §4。

定義 (圧縮性・軸対称の平面近似 — $\delta^* \ll r_w$ なので曲率補正は省く):

$$\delta^*(x) = \int_0^{y_e} \left(1 - \frac{\rho u}{(\rho u)_e}\right) dy,
  \qquad y = \text{壁からの法線距離 (鉛直近似)}$$

**教訓 (run_0030 で実測、繰り返し禁止)**:
- 探索窓はフリーストリームに届くまで取る (窓 0.35 r* は下流で BL 縁の手前で切れ、
  「δ* が下流ほど減る」という非物理な結論を出した — 撤回済み)。既定 1.5 r*。
- スロート近傍 (x ≲ 8 r_t) はコア流が半径方向に未一様で「ρu 極大 = BL 縁」の
  判定が破綻する。呼び出し側で相関へフォールバックすること。
"""
from __future__ import annotations

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator


def deltastar_from_run(mesh_h5, res_h5, wall_xy, scale: float,
                       x_stations, window: float = 1.5, n_prof: int = 400,
                       edge_frac: float = 0.995) -> dict:
    """node RANS run から各 x ステーションの δ* を測る。

    wall_xy: (n,2) 物理壁 [m 単位 or r_t 単位 — scale と整合していれば可]。
    x_stations: r_t 単位。戻り: dict(x, dstar, ue_rel, y_edge, ok) 全て r_t 単位。
    """
    with h5py.File(mesh_h5) as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(res_h5) as f:
        ro = f["/VALUE/ro"][:]
        Ux = f["/VALUE/Ux"][:]
    if len(ro) != len(nc):
        raise ValueError("node run でない (VALUE 長 != 節点数)")
    x_n, r_n = nc[:, 0] / scale, nc[:, 1] / scale
    itp = LinearNDInterpolator(np.c_[x_n, r_n], ro * Ux)
    wx, wr = np.asarray(wall_xy)[:, 0], np.asarray(wall_xy)[:, 1]
    out = {"x": [], "dstar": [], "y_edge": [], "ok": []}
    for xs in np.atleast_1d(x_stations):
        rw = float(np.interp(xs, wx, wr))
        y = np.linspace(0.0, min(window, 0.95 * rw), n_prof)   # 壁 → 内側
        q = itp(np.full(n_prof, xs), rw - y)
        good = np.isfinite(q)
        if good.sum() < n_prof // 2:
            out["x"].append(float(xs)); out["dstar"].append(np.nan)
            out["y_edge"].append(np.nan); out["ok"].append(False)
            continue
        q = np.where(good, q, 0.0)
        qe = float(np.max(q))
        # 縁 = 壁から見て最初に edge_frac·(ρu)_max へ達する y
        idx = np.argmax(q >= edge_frac * qe)
        y_e = float(y[idx])
        m = y <= y_e
        integ = float(np.trapezoid(1.0 - q[m] / max(qe, 1e-30), y[m]))
        out["x"].append(float(xs)); out["dstar"].append(integ)
        out["y_edge"].append(y_e); out["ok"].append(bool(idx > 3))
    return {k: np.asarray(v) for k, v in out.items()}
