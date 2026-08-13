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
- **到達不能域の目標 (B7)**: `local_anchor(x_c)` は $x_0$ の点アンカーと同じ狭窓
  局所フィットを任意の中心点で行う (`design_chain` が $x_{\mathrm{reach}}$ で
  呼ぶ)。$[x_0,x_{\mathrm{reach}})$ 全体は**両端の $(M,M',M'')$ だけで決まる
  5 次 Hermite** (`MachBezier.from_constraints` に free_cp なしで渡す) で繋ぐ —
  この区間全体を 1 本の曲線として直接フィットする方式は、近スロート域のセル列が
  疎 (cell モードに軸上 DOF が無い) なため高次多項式が Runge 振動・bin 平均は
  非単調・等調回帰は勾配を潰す、と軒並み破綻したため不採用 (経緯は
  `local_anchor` の docstring 参照)。

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
        """res_h5: 単一ファイルまたはリスト。リストなら M,θ をスナップショット平均
        してから抽出する (中域リミットサイクルの位相ノイズを落とす — 軸 M 測定を
        末尾平均で行うのと同じ理由)。"""
        self.R, self.g, self.scale = float(R), float(gamma), float(scale)
        with h5py.File(mesh_h5) as nz:
            cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
        res_list = [res_h5] if isinstance(res_h5, (str, Path)) else list(res_h5)
        Ms, ths = [], []
        for rf in res_list:
            with h5py.File(rf) as f:
                Ux, Uy, son = f["/VALUE/Ux"][:], f["/VALUE/Uy"][:], f["/VALUE/sonic"][:]
            Ms.append(np.hypot(Ux, Uy) / np.maximum(son, 1e-9))
            ths.append(np.arctan2(Uy, Ux))
        x, r = cc[:, 0] / scale, cc[:, 1] / scale     # r* 単位
        M = np.mean(Ms, axis=0)
        th = np.mean(ths, axis=0)
        # --- 軸 M(x) の局所 4 次フィット (点アンカー M,M',M'' 用の狭窓) ---
        ax = (r < axis_band_frac) & (np.abs(x - x0_guess) < fit_halfwidth)
        if ax.sum() < 12:
            raise ValueError(f"軸帯セル不足 ({int(ax.sum())})")
        self._axis_c = np.polyfit(x[ax] - x0_guess, M[ax], 4)
        self._axis_x0g = float(x0_guess)
        self._axis_win = (float(x[ax].min()), float(x[ax].max()))
        # 広窓 (x 制限なし) の軸帯生データも保持 — local_anchor() が使う
        ax_wide = r < axis_band_frac
        self._axis_x_wide = x[ax_wide].copy()
        self._axis_M_wide = M[ax_wide].copy()
        self._axis_band_frac = float(axis_band_frac)
        # 近スロートの全セル (r 制限緩め) — axis_segment_curve() の列ごと外挿が使う
        seg = (x < 6.0) & (r < 0.6)
        self._seg_x, self._seg_r, self._seg_M = x[seg].copy(), r[seg].copy(), M[seg].copy()
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

    def axis_segment_curve(self, x_lo: float, x_hi: float, band: float | None = None,
                           pad: float = 0.25, lam: float | None = 1e-5):
        r"""$[x_{lo},x_{hi}]$ の軸 $M(x)$ を **C2 平滑化スプライン**で返す (B7 改)。

        **単調性は仮定しない** — 到達不能域には接合波の軸着地による実在の局所こぶ
        (実測 run_0017: x/rt≈1.4 に山、3 スナップ平均でも再現) があり、以前の
        等調回帰・両端 Hermite はこれを潰して**人工残差 ±0.05〜0.07 を作った**
        (run_0020 で実証)。到達不能域の目標は「実流がやること」そのもの — 山ごと写す。

        **測定規約との整合が最重要**: ΔM 評価 (`metrics.extract.axis_mach`) は
        「r < axis_band のセルの値」をそのまま使う**帯規約**で、近スロートの
        径方向勾配では真の軸値より +0.02〜0.04 大きい。ループは測定値を目標へ
        寄せるので、**目標も同じ帯規約で抽出する** (r→0 の偶基底外挿 [真値] を
        使うと規約差がそのまま系統残差になる — 実測 0.04)。抽出は列ごと
        (構造化メッシュの x 列を丸めで同定) に r<band セルを平均し、列点列を
        `make_smoothing_spline` (C2 3 次平滑化スプライン, 列点残差 ~1e-3) に通す。
        戻り値はそのスプライン (`.derivative(k)` で $M',M''$ — x_reach アンカーも
        この曲線から取ることで接続 C2 が構成的に成り立つ)。`pad` は両端の微分を
        安定させるための延長窓。band 既定 = コンストラクタの axis_band_frac
        (= 測定と同じ帯)。"""
        if band is None:
            band = self._axis_band_frac
        m = (self._seg_x >= x_lo - pad) & (self._seg_x <= x_hi + pad) & (self._seg_r < band)
        xs, rs, Ms = self._seg_x[m], self._seg_r[m], self._seg_M[m]
        xcol = np.round(xs, 3)
        uniq = np.unique(xcol)
        if len(uniq) < 8:
            raise ValueError(f"axis_segment_curve: 列不足 ({len(uniq)})")
        pts = []
        for xc in uniq:
            mm = xcol == xc
            pts.append((float(np.mean(xs[mm])), float(np.mean(Ms[mm]))))
        pts = np.asarray(sorted(pts))
        from scipy.interpolate import make_smoothing_spline
        return make_smoothing_spline(pts[:, 0], pts[:, 1], lam=lam)

    def local_anchor(self, x_c: float, halfwidth: float = 0.4, degree: int = 2) -> tuple:
        r"""$x_c$ 近傍の狭窓局所フィットで $(M,M',M'')$ を返す (B7)。

        コンストラクタの $x_0$ 点アンカー (`_axis_c`) と同じ手法を任意の中心点に
        一般化したもの — $x_{\mathrm{reach}}$ でのアンカー取得に使う。

        **経緯 (2026-08-15 実測)**: 当初 $[x_0,x_{\mathrm{reach}})$ 全体を 1 本の
        曲線としてフィットしようとしたが、この近スロート域は軸方向のセル列が
        非常に疎 (実測: 幅 1.4 $r_t$ の区間にわずか ~30 列、cell モードに軸上
        DOF が無いため列ごとに近軸セルも 1 個程度) で、高次多項式は Runge 振動、
        bin 平均は列疎密の混入で非単調、等調回帰は物理的な勾配まで潰して
        $x_{\mathrm{reach}}$ で $M'\approx0$ になった。**全区間を 1 本でフィット
        するのをやめ**、$x_0$ の点アンカーと同じ**狭窓ローカルフィット**を
        $x_{\mathrm{reach}}$ でも取り、区間全体は両端の $(M,M',M'')$ だけで
        決まる 5 次 Hermite (`MachBezier.from_constraints`, free_cp なし) で
        繋ぐ。ただし $x_0$ 近傍 (throat_refine でセル密) と違い任意点近傍は
        セル列が疎な場合があるため、既定は $x_0$ の点アンカー (degree=4,
        halfwidth=0.4) より**低次・広窓** (degree=2, halfwidth=0.4) — 実測で
        degree=4/halfwidth=0.15 は $x_{\mathrm{reach}}$ 近傍のセル不足で $M'$ の
        符号が反転した (負の傾き)。degree=2/halfwidth=0.4 で正の傾きに安定
        (2026-08-15)。"""
        ax = np.abs(self._axis_x_wide - x_c) < halfwidth
        if int(ax.sum()) < degree + 4:
            raise ValueError(f"local_anchor: セル不足 ({int(ax.sum())}, x_c={x_c:.3f})")
        c = np.polyfit(self._axis_x_wide[ax] - x_c, self._axis_M_wide[ax], degree)
        h = 1e-4
        M = float(np.polyval(c, 0.0))
        Mp = float((np.polyval(c, h) - np.polyval(c, -h)) / (2 * h))
        Mpp = float((np.polyval(c, h) - 2.0 * M + np.polyval(c, -h)) / h ** 2)
        return M, Mp, Mpp

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


def axis_curve_true(run_dir, x_lo: float, x_hi: float, scale: float = 0.01,
                    n_last: int = 3, band: float = 0.35, r_deg: int = 4) -> tuple:
    r"""**真の軸値** $a_0(x)$ を各 $x$ station の偶関数フィットで復元する (B10)。

    cell 計算は軸上に DOF が無いため、近軸帯の単純平均は径方向勾配を混入させる
    (実測: 近スロートで真値より +0.02〜0.04)。軸対称性から $M$ は $r$ の偶関数
    $M(x,r)=a_0(x)+a_2(x)r^2+a_4(x)r^4+\cdots$ なので、station ごとに偶多項式を
    最小二乗し $a_0$ を採る。**末尾 `n_last` スナップ平均**を先に取る。
    戻り: (x_stations, a0)。
    """
    from ..geometry.cminus_cfd import load_field
    fld = load_field(run_dir, n_last=n_last, scale=scale, discretization="cell")
    m = (fld["x"] >= x_lo) & (fld["x"] <= x_hi) & (fld["r"] < band)
    xs, rs, Ms = fld["x"][m], fld["r"][m], fld["M"][m]
    col = np.round(xs, 3)
    out = []
    for xc in np.unique(col):
        s = col == xc
        rr, MM = rs[s], Ms[s]
        k = min(r_deg // 2 + 1, max(len(rr), 1))
        B = np.stack([rr ** (2 * j) for j in range(k)], axis=1)
        if len(rr) < k:
            out.append((float(np.mean(xs[s])), float(np.mean(MM))))
            continue
        a, *_ = np.linalg.lstsq(B, MM, rcond=None)
        out.append((float(np.mean(xs[s])), float(a[0])))
    arr = np.asarray(sorted(out))
    return arr[:, 0], arr[:, 1]


def onesided_axis_anchor(run_dir, x_c: float, side: str = "left",
                         degree: int = 3, window: float = 0.6,
                         scale: float = 0.01, n_last: int = 3,
                         band: float = 0.35) -> dict:
    r"""$x_c$ の**片側だけ**を使う局所回帰で $(M,M',M'')$ を返す (B10 正本)。

    $M_r=M(x_{reach}^-)$ 等の**左極限**が要る。$x>x_{reach}$ はこれから設計し直す
    領域なので、左右対称窓のフィット (既存 `local_anchor`) に混ぜてはいけない —
    下流の仮壁の影響が境界条件へ漏れる。

    **メッシュ値の差分は使わない**。$\xi=x-x_c$ の局所多項式
    $M=c_0+c_1\xi+c_2\xi^2+\cdots$ を上流側 (`side="left"`) の真の軸値
    $a_0(x)$ (`axis_curve_true`) へ最小二乗し、係数から解析的に
    $M_r=c_0,\;M'_r=c_1,\;M''_r=2c_2$ を取る。
    """
    lo, hi = ((x_c - window, x_c) if side == "left" else (x_c, x_c + window))
    xs, a0 = axis_curve_true(run_dir, lo - 0.05, hi + 0.05, scale=scale,
                             n_last=n_last, band=band)
    m = (xs >= lo) & (xs <= hi)
    xi, Mi = xs[m] - x_c, a0[m]
    if len(xi) < degree + 2:
        return {"ok": False, "reason": f"station 不足 ({len(xi)}, deg={degree}, win={window})"}
    c = np.polyfit(xi, Mi, degree)          # 高次→低次
    return {"ok": True, "M": float(c[-1]), "Mp": float(c[-2]), "Mpp": float(2.0 * c[-3]),
            "n_station": int(len(xi)), "side": side, "degree": degree,
            "window": window, "rms_fit": float(np.sqrt(np.mean(
                (np.polyval(c, xi) - Mi) ** 2)))}


def anchor_sensitivity(run_dir, x_c: float, scale: float = 0.01,
                       n_last: int = 3, x_c_err: float = 0.01) -> dict:
    """B10 が要求する $(M,M',M'')$ の感度一式を測る。

    次数 (2/3/4) × 上流窓幅 / 末尾スナップ数 / $x_{reach}$ トレース誤差 /
    **左右片側フィットの勾配差** (格子収束後も残るなら仮壁に本当に接続波がある)。
    """
    out = {"x_c": x_c, "left": {}, "right": {}, "snap": {}, "x_c_shift": {}}
    for deg in (2, 3, 4):
        for win in (0.3, 0.45, 0.6, 0.9):
            r = onesided_axis_anchor(run_dir, x_c, "left", deg, win, scale, n_last)
            if r.get("ok"):
                out["left"][f"deg{deg}_win{win}"] = (r["M"], r["Mp"], r["Mpp"], r["rms_fit"])
            r2 = onesided_axis_anchor(run_dir, x_c, "right", deg, win, scale, n_last)
            if r2.get("ok"):
                out["right"][f"deg{deg}_win{win}"] = (r2["M"], r2["Mp"], r2["Mpp"], r2["rms_fit"])
    for n in (1, 3, 5):
        r = onesided_axis_anchor(run_dir, x_c, "left", 3, 0.6, scale, n)
        if r.get("ok"):
            out["snap"][f"n{n}"] = (r["M"], r["Mp"], r["Mpp"])
    for dx in (-x_c_err, 0.0, x_c_err):
        r = onesided_axis_anchor(run_dir, x_c + dx, "left", 3, 0.6, scale, n_last)
        if r.get("ok"):
            out["x_c_shift"][f"{dx:+.3f}"] = (r["M"], r["Mp"], r["Mpp"])
    return out


def from_run(run_dir, scale: float, R: float, gamma: float = 1.4) -> CFDThroat:
    """run ディレクトリ (nozzle.h5 + 末尾 res_*.h5 最大 3 個の平均) から作る。

    末尾平均は中域リミットサイクル (片振幅 ~0.003) の位相ノイズ対策 — 軸 M の
    測定を末尾スナップ平均で行うのと同じ規約。"""
    rd = Path(run_dir)
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    res = [f for f in res if _steps_of_name(f) > 0]
    if not res:
        raise FileNotFoundError(f"{rd} に res_*.h5 が無い")
    return CFDThroat(rd / "nozzle.h5", res[-3:], scale, R, gamma)


def _steps_of_name(f: Path) -> int:
    return int("".join(c for c in f.stem if c.isdigit()))
