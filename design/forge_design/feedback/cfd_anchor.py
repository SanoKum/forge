r"""CFD アンカー遷音速データ (①モード F, plan §9.2 B6)。

pass CFD 場からスロート近傍の Cauchy データを抽出し、`SauerThroat` 互換の
インターフェース (`starting_line` / `mach`) で `inverse_design` と `design_chain`
に供給する。動機 (B5 実測): Sauer 一次解は実 CFD に対し軸上 $M'$ が +25%、
$M''$ が 2.3 倍ずれ、starting line 半径方向も壁に向かって ΔM −0.09 / Δθ −0.5°
まで開く — 接合波の主因は「設計 Cauchy データと実流の不整合」だった。

方針:
- **軸**: 軸帯セルの $M(x)$ を $x_0$ 近傍の 4 次多項式で最小二乗フィットし、
  $x_0$ は「フィット上で $M=M_{start}$」の根として決める (Sauer の
  `x_axis_of_mach` と同じ意味論)。$M,M',M''(x_0)$ もフィット微分で返す
  (窓 ±0.15/0.25/0.40 で ±0.02 以内の安定を B5 で確認済み)。
- **starting line**: $x=x_0$ 縦線上の $M(r),\theta(r)$ をセル場の線形補間で
  サンプルし、軸対称性を持つ基底 (M: $r$ の偶多項式 / θ: 奇多項式) で最小二乗
  フィット。$r=0$ (セルが無い) と壁足 (セル外) はフィットの評価で自然に埋まる。
  壁足の θ は `design_chain` 側が従来どおり円弧接線で上書きする (壁 BC の厳密値)。
- **凍結**: 遷音速場は下流壁にほぼ依存しない (実測 ΔM=0.0024) ため、抽出元 run
  1 つに固定して使い回す。パスごとに再抽出しない (目標が動くと ΔM 評価の
  ゴールポストが移動するため)。

使い方: problem YAML の `geometry.throat_anchor_run: <run_dir>` を指定すると
`runner_wt.design_chain` が本モジュールで throat を構築して Sauer を置き換える。
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np


class CFDThroat:
    """CFD 場から抽出した遷音速アンカー (SauerThroat 互換の必要最小)。"""

    def __init__(self, mesh_h5, res_h5, scale: float, R: float,
                 gamma: float = 1.4, fit_halfwidth: float = 0.4,
                 axis_band_frac: float = 0.09, x0_guess: float = 0.26) -> None:
        self.R, self.g, self.scale = float(R), float(gamma), float(scale)
        with h5py.File(mesh_h5) as nz:
            cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
        with h5py.File(res_h5) as f:
            Ux, Uy, son = f["/VALUE/Ux"][:], f["/VALUE/Uy"][:], f["/VALUE/sonic"][:]
        x, r = cc[:, 0] / scale, cc[:, 1] / scale     # r* 単位
        M = np.hypot(Ux, Uy) / np.maximum(son, 1e-9)
        th = np.arctan2(Uy, Ux)
        # --- 軸 M(x) の局所 4 次フィット ---
        ax = (r < axis_band_frac) & (np.abs(x - x0_guess) < fit_halfwidth)
        if ax.sum() < 12:
            raise ValueError(f"軸帯セル不足 ({int(ax.sum())})")
        self._axis_c = np.polyfit(x[ax] - x0_guess, M[ax], 4)
        self._axis_x0g = float(x0_guess)
        self._axis_win = (float(x[ax].min()), float(x[ax].max()))
        # --- 縦線サンプル用の近傍データ (starting_line で利用) ---
        near = (np.abs(x - x0_guess) < 1.0) & (r < 1.3)
        from scipy.interpolate import LinearNDInterpolator
        self._iM = LinearNDInterpolator(np.c_[x[near], r[near]], M[near])
        self._iTh = LinearNDInterpolator(np.c_[x[near], r[near]], th[near])
        self._r_cell_min = float(r[near].min())
        self._line_fit_diag: dict = {}

    # -- SauerThroat 互換 ------------------------------------------------------
    def mach(self, x, r):
        """軸上 (r=0) 専用: design_chain のアンカー差分が呼ぶ。"""
        if np.any(np.asarray(r) != 0.0):
            raise NotImplementedError("CFDThroat.mach は軸上 (r=0) のみ")
        xq = np.asarray(x, dtype=float) - self._axis_x0g
        return np.polyval(self._axis_c, xq)

    def x_axis_of_mach(self, M_start: float) -> float:
        xs = np.linspace(self._axis_win[0], self._axis_win[1], 2000)
        Ms = np.polyval(self._axis_c, xs - self._axis_x0g)
        if not (Ms.min() < M_start < Ms.max()):
            raise ValueError(f"M_start={M_start} が軸フィット範囲 "
                             f"[{Ms.min():.3f},{Ms.max():.3f}] 外")
        return float(np.interp(M_start, Ms, xs))   # 近スロートで単調増加

    def starting_line(self, M_start: float = 1.05, n: int = 25):
        """x=x0 縦線 (軸→円弧壁足) の (x0, r, M, θ)。x0 は CFD 軸 M=M_start の根。"""
        x0 = self.x_axis_of_mach(M_start)
        r_w = 1.0 + self.R - float(np.sqrt(max(self.R ** 2 - x0 ** 2, 0.0)))  # 円弧厳密
        # セルのある帯でサンプル → 偶/奇多項式で最小二乗 (r=0・壁足へは基底で外挿)
        rs = np.linspace(self._r_cell_min * 1.05, r_w * 0.985, 160)
        Mv = np.array([float(self._iM(x0, ri)) for ri in rs])
        Tv = np.array([float(self._iTh(x0, ri)) for ri in rs])
        ok = np.isfinite(Mv) & np.isfinite(Tv)
        rs, Mv, Tv = rs[ok], Mv[ok], Tv[ok]
        if len(rs) < 20:
            raise ValueError(f"縦線サンプル不足 ({len(rs)})")
        BM = np.stack([rs ** (2 * k) for k in range(5)], axis=1)   # 偶: 1..r^8
        BT = np.stack([rs ** (2 * k + 1) for k in range(4)], axis=1)  # 奇: r..r^7
        cM, *_ = np.linalg.lstsq(BM, Mv, rcond=None)
        cT, *_ = np.linalg.lstsq(BT, Tv, rcond=None)
        self._line_fit_diag = {
            "x0": float(x0), "r_wall": float(r_w),
            "rms_M": float(np.sqrt(np.mean((BM @ cM - Mv) ** 2))),
            "rms_theta_deg": float(np.rad2deg(np.sqrt(np.mean((BT @ cT - Tv) ** 2)))),
            "n_samples": int(len(rs)),
        }
        r = np.linspace(0.0, r_w, n)
        Mo = np.stack([r ** (2 * k) for k in range(5)], axis=1) @ cM
        To = np.stack([r ** (2 * k + 1) for k in range(4)], axis=1) @ cT
        return float(x0), r, Mo, To


def from_run(run_dir, scale: float, R: float, gamma: float = 1.4) -> CFDThroat:
    """run ディレクトリ (nozzle.h5 + 最新 res_*.h5) から CFDThroat を作る。"""
    rd = Path(run_dir)
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if not res:
        raise FileNotFoundError(f"{rd} に res_*.h5 が無い")
    return CFDThroat(rd / "nozzle.h5", res[-1], scale, R, gamma)
