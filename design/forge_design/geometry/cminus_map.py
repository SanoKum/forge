"""v1 MOC 設計場の凍結 C⁻ マップ (壁点 → 軸着地点)。

`feedback/euler_loop.py` (v2 帰還ループ) と `evaluate/runner_wt.py`
(`design_chain` の B7 ブートストラップ, $x_{\\mathrm{reach}}$ 抽出) の両方が使うため
独立モジュールに置く (`euler_loop` は `runner_wt` を import するので、後者に置くと
循環 import になる)。
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import LinearNDInterpolator

from .moc_kernel import pm_mach


def build_cminus_map(pts, wall_inv, gamma: float) -> np.ndarray:
    """**v2**: v1 の MOC 設計場で各壁点から C⁻ を軸まで前進積分し (x_w → x_axis)
    マップを作る (以後のパスで凍結して使う)。"""
    xy = np.array([[p.x, p.r] for p in pts])
    th = np.array([p.th for p in pts])
    nu = np.array([p.nu for p in pts])
    itp_th = LinearNDInterpolator(xy, th)
    itp_nu = LinearNDInterpolator(xy, nu)
    rows = []
    for xw, rw in wall_inv[::4, :2]:
        x, r = float(xw), float(rw)
        ok = True
        for _ in range(4000):
            t = itp_th(x, r)
            v = itp_nu(x, r)
            if not (np.isfinite(t) and np.isfinite(v)):
                ok = False
                break
            M = pm_mach(max(float(v), 1e-9))
            mu = np.arcsin(1.0 / max(M, 1.0 + 1e-9))
            dr = -min(0.02, r)
            x += dr / np.tan(float(t) - mu)
            r += dr
            if r <= 1e-9:
                break
        if ok and r <= 1e-6:
            rows.append((float(xw), x))
    return np.asarray(rows)
