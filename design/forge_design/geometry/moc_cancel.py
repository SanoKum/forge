r"""wave-cancellation MOC (Phase W4): D 発 C⁻ データ線 → 下流壁の従属生成。

plan: plans/active/tooling-nozzle-walldriven-chain.md §4 (W4 cancellation MOC の仕様)。
CONTUR の「壁は MOC の従属出力」思想 (調査 §2.1) の実装。$M_c(x)$ を事前規定せず、
データ線 (D から出る C⁻、CFD 実場から抽出) + 軸対称のみを Cauchy データとして
場を前進マーチし、**壁を「追加の compatibility を課さない」流線として自然に得る**。

## 数理

データ線は軸対称・非回転・等エントロピー MOC の意味で「両方の Riemann 不変量
$J^+=\theta-\nu$, $J^-=\theta+\nu$ が各点で既知」な Cauchy データ (CFD 実場から
両方読み取っているため)。データ線自身が (ほぼ) C⁻ 方向に沿っていることは、
標準の interior 単位過程 (2 点 A, B から C+ (A 発) と C- (B 発) の交点を作る) を
内点の充填に適用する妨げにならない。

**役割の割当て (2026-08-15 導出)**: 軸側 (r 小) を A (C+ 担体)、壁側 (r 大) を B
(C- 担体) とする。これは `moc_kernel.KernelMOC` の楔充填・壁ステーションマーチが
実際に使っている割当て (`front[j]` = 既存側 = A) と一致し、`moc_inverse.fill()` が
縦の starting line 区間で使う固定割当て (逆) とは異なる — 後者は近軸への丸め誤差で
文書化済みの限界を持つ (moc_inverse.py 冒頭の楔欠損の記録)。データ線は C⁻ 方向に
沿う (starting line より軸傾向が強い) ため、この割当てが正しい。

## 状態 (2026-08-15): データ線抽出は検証済み・production 品質。**壁生成
(`_wall_cancel`/`march_wall_from_dataline`/`cancellation_contour`) は未検証**

**検証済み**: `extract_dataline_cfd` (実 CFD run から D 発 C⁻ + P0 診断を抽出) と
`build_field` (interior()+軸反射での内点充填) — 放射源流厳密解照合で個々の点は
機械精度近く一致 (相対誤差 ~1e-4)、格子収束も確認済み (n=41→81 で <0.1%)。

**未検証 (既知の課題)**: 壁点列を生成する専用単位過程は、文献で確認した高位の
定式化 (壁で $J^-$ を直前の壁点から引き継ぎ、$J^+$ を隣接内点から得る —
`θ=(K^-+K^+)/2, ν=(K^--K^+)/2` の標準 MLN 型公式) は正しいはずだが、**軸対称
源項の付加・フロントの更新方式 (どの点を「隣接内点」に使うか) の実装を 3 通り
試して、いずれも放射源流厳密解照合で有意な乖離が残った** (最良でも壁角が
15 反復で厳密解の 3 倍近い速さで減衰。詳細検証ログは
`plans/active/tooling-nozzle-walldriven-chain.md` W4 節を参照):

1. 定長フロントを毎回全体縮約し `front[-2]` を使い回す方式 — 流線として不連続
   (隣接壁点間で $dr/dx\not\approx\tan\theta$、厳密解とも大きく乖離)。
2. `interior()+軸反射` だけで場を作り事後的に `wall_streamline` で積分する方式 —
   数学的には健全 (円錐データ線で流線が機械精度で厳密解と一致) だが、
   **interior()+軸反射だけの構築は $r\le r_D$ 側にしか点を作らず**、
   $\theta_D>0$ の壁が実際に進む $r>r_D$ 側を全く被覆しないため場外 (NaN) に落ちる。
3. `_march_stations` (moc_kernel) 型の「壁側からはしごを軸側へ下る」構成に
   `_wall_cancel` (下記) を組み込んだ方式 — 流線としては連続に近づいたが、
   個々の壁点が厳密解と有意に乖離 (数十%) したまま。原因未特定
   (軸対称源項の符号・隣接点選択・はしごの向きのいずれかに残存バグの疑い)。

`_wall_cancel(A, W_prev)` の現在の実装 (試行 3 相当) は、
$J^+_W=J^+_A+S_p$ (通常の C+ compatibility) と
$J^-_W=J^-_{W_{prev}}+S_m$ (壁上で $J^-$ を引き継ぐ、$S_m$ は軸対称源項補正) の
連立を `_interior` と同型の predictor-corrector で解く。退化チェック
($A=W_{prev}$ で $\theta_W\to\theta_A,\nu_W\to\nu_A$) は機械精度で通るが、
実際の march では厳密解と一致しない。**この関数群は production では使わないこと
— 呼び出すと `march_wall_from_dataline`/`cancellation_contour` が動きはするが、
出力される壁座標は検証されておらず信用できない。** 次セッションでの継続時は、
一次資料 (Anderson *Modern Compressible Flow* Ch.11.7) に当たって壁点単位過程の
軸対称版を確認するか、2 点だけの孤立ケース (単純波・既知の解析壁) で
`_wall_cancel` 単体をまず検証してから march へ組み込むことを推奨する。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cminus_cfd import _interpolators, load_field, trace_cminus_from_wall
from .moc_kernel import KernelMOC, _Pt, _mu, _sin_over_r, pm_mach, pm_nu


@dataclass
class DataLine:
    """D 発 C⁻ データ線 (軸→壁 順, r* 単位)。"""

    x: np.ndarray
    r: np.ndarray
    theta: np.ndarray
    M: np.ndarray
    gamma: float
    p0_diag: dict
    x_j: float
    r_j: float
    x_axis: float

    def to_pts(self) -> list:
        nu = pm_nu(self.M, self.gamma)
        return [_Pt(float(self.x[i]), float(self.r[i]), float(self.theta[i]),
                    float(nu[i]), self.gamma) for i in range(len(self.x))]


def _p0_over_pt(P, M, gamma):
    return P * (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0))


def extract_dataline_cfd(run_dir, x_j: float, r_j: float, gamma: float = 1.4,
                         scale: float = 0.01, n_line: int = 31,
                         n_last: int = 3, wall_theta_geom: float | None = None
                         ) -> DataLine:
    r"""D=(x_j,r_j) から C⁻ を壁面状態起点で軸まで追跡し、$(x,r,\theta,M)$ + $P_0$
    診断を持つ `DataLine` (軸→壁 順、n_line 点に再標本化) を返す。

    `wall_theta_geom`: 壁接線角 (幾何値) を渡すと `trace_cminus_from_wall` が
    壁面状態出力の $\theta$ の代わりに使う (slip 壁では流れ角と厳密に一致するはずの
    量)。
    """
    fld = load_field(run_dir, n_last=n_last, scale=scale, discretization="node",
                     with_pressure=True)
    tr = trace_cminus_from_wall(fld, run_dir, x_j, r_j, scale=scale, gamma=gamma,
                                th_wall_geom=wall_theta_geom)
    if not tr.get("ok", False):
        raise RuntimeError(f"データ線抽出失敗: {tr.get('reason')}")
    path = np.asarray(tr["path"])  # [wall(D)...axis] 順、(n,2) = (x,r)
    x_axis = float(tr["x_axis"])

    # 軸端点: 極近軸は補間の凸包縁で不安定なため、ロバストな帯フィットで置き換える
    from ..feedback.cfd_anchor import axis_curve_true
    xs_a, a0 = axis_curve_true(run_dir, x_axis - 0.15, x_axis + 0.15, scale=scale,
                               n_last=n_last, band=0.35)
    if len(xs_a) == 0:
        raise RuntimeError(f"軸端点 x={x_axis:.4f} 近傍の帯フィットが空")
    M_axis = float(a0[np.argmin(np.abs(xs_a - x_axis))])

    iM, iTh = _interpolators(fld)
    # path は wall→axis 順で密 (ds=2e-3)。壁近傍過密を避けるため r で等間隔に間引く。
    r_wall_path = float(path[0, 1])
    r_grid = np.linspace(0.0, r_wall_path, n_line)[1:-1]  # 端は個別に扱う
    x_of_r = np.interp(r_grid, path[::-1, 1], path[::-1, 0])
    xs = np.concatenate([[x_axis], x_of_r, [x_j]])
    rs = np.concatenate([[0.0], r_grid, [r_j]])
    M = np.empty_like(xs)
    th = np.empty_like(xs)
    M[0], th[0] = M_axis, 0.0
    for i in range(1, len(xs) - 1):
        M[i], th[i] = float(iM(xs[i], rs[i])), float(iTh(xs[i], rs[i]))
    # 壁端点は壁面状態 (境界出力) を使う — 内部補間の凸包縁より信頼できる
    M[-1] = float(tr["wall_state"]["M_wall"])
    th[-1] = float(wall_theta_geom if wall_theta_geom is not None
                   else tr["wall_state"]["th_wall_cfd_deg"] * np.pi / 180.0)
    if not np.all(np.isfinite(M)) or not np.all(np.isfinite(th)):
        raise RuntimeError("データ線の再標本化点に非有限値 (n_line を減らすか場を確認)")

    # P0 診断 (壁直近 5% は境界補間混入を避け除外)
    Pv = np.array([_p0_at(fld, x, r) for x, r in zip(xs, rs)])
    valid = np.isfinite(Pv)
    n_skip = max(1, int(0.05 * valid.sum()))
    Pv_core = Pv[valid][n_skip:]
    Mv_core = M[valid][n_skip:]
    P0 = _p0_over_pt(Pv_core, Mv_core, gamma)
    dP0_rel = float((P0.max() - P0.min()) / np.mean(P0)) if len(P0) > 1 else 0.0
    gate = "PASS" if dP0_rel <= 1e-3 else ("WARN" if dP0_rel <= 5e-3 else "FAIL")
    p0_diag = {"dP0_rel": dP0_rel, "gate": gate, "n_pts": int(len(P0)),
              "P0_mean_Pa": float(np.mean(P0))}

    return DataLine(x=xs, r=rs, theta=th, M=M, gamma=gamma, p0_diag=p0_diag,
                    x_j=float(x_j), r_j=float(r_j), x_axis=x_axis)


def _p0_at(fld, x, r):
    """`fld` (load_field 出力) の最近傍点から局所 P [Pa] を取る簡易実装。"""
    if "P" not in fld:
        return float("nan")
    idx = np.argmin((fld["x"] - x) ** 2 + (fld["r"] - r) ** 2)
    return float(fld["P"][idx])


class _CancelKernel(KernelMOC):
    """`_interior`/`_axis` (内点充填用) + `_wall_cancel` (壁専用単位過程)。"""

    def _wall_cancel(self, A: _Pt, W_prev: _Pt) -> _Pt:
        r"""壁点 W を A (C+ 担体) と直前の壁点 W_prev から構成する。

        $J^+_W = J^+_A + S_p$ (通常の C+ compatibility)。
        $J^-_W = J^-_{W_{prev}} + S_m$ (壁上で $J^-$ を引き継ぐ — 反射なしの意味)。
        位置は C+ レイ (A 発) と壁接線レイ (W_prev 発, 傾き $\tan(\theta$ 平均$)$、
        通常の C⁻ の $-\mu$ 項なし) の交点。
        """
        g, d = self.g, self.delta
        thW, nuW = W_prev.th, A.nu + (W_prev.th - A.th)  # 初期推定
        xW = rW = None
        for _ in range(2 + self.n_corr):
            MW = pm_mach(max(nuW, 1e-8), g)
            muW = _mu(MW)
            mp = np.tan(0.5 * (A.th + thW) + 0.5 * (_mu(A.M) + muW))
            mw = np.tan(0.5 * (W_prev.th + thW))
            if abs(mw - mp) < 1e-12:
                raise RuntimeError("壁専用単位過程: C+ と壁接線が平行")
            xW = (A.r - W_prev.r + mw * W_prev.x - mp * A.x) / (mw - mp)
            rW = A.r + mp * (xW - A.x)
            Pn = _Pt(xW, max(rW, 0.0), thW, max(nuW, 0.0), g)
            fA = 0.5 * (np.sin(_mu(A.M)) * _sin_over_r(A, Pn)
                       + np.sin(muW) * _sin_over_r(Pn, A))
            # 壁レイの軸対称源項 (2026-08-15 修正: sin(0.0) は常にゼロになるバグ —
            # 壁点 W_prev/Pn 自身の Mach 角を使う。角度は「壁レイ自身の傾き」
            # 0.5*(W_prev.th+thW) を使う (C- の ∓μ オフセットは持たない)
            fW = 0.5 * (np.sin(_mu(W_prev.M)) * _sin_over_r(W_prev, Pn)
                       + np.sin(muW) * _sin_over_r(Pn, W_prev))
            Sp = -d * fA / np.cos(0.5 * (A.th + thW) + 0.5 * (_mu(A.M) + muW)) * (xW - A.x)
            Sm = d * fW / np.cos(0.5 * (W_prev.th + thW)) * (xW - W_prev.x)
            Jp = (A.th - A.nu) + Sp
            Jm = (W_prev.th + W_prev.nu) + Sm
            thW, nuW = 0.5 * (Jp + Jm), 0.5 * (Jm - Jp)
        return _Pt(xW, max(rW, 0.0), thW, max(nuW, 0.0), g)


def build_field(dataline: DataLine, n_rounds: int, gamma: float | None = None,
                n_corr: int = 2) -> list:
    """データ線から interior()+軸反射を `n_rounds` 回反復し内点を蓄積する
    (診断・格子収束確認用。壁生成には `march_wall_from_dataline` を使う)。
    """
    g = gamma if gamma is not None else dataline.gamma
    kmoc = _CancelKernel(gamma=g, delta=1.0, n_corr=n_corr)
    front = dataline.to_pts()
    if len(front) < 3:
        raise ValueError("データ線の点数が不足 (interior 単位過程に最低 3 点必要)")
    all_pts = list(front)
    for _ in range(n_rounds):
        new = []
        for i in range(len(front) - 1):
            try:
                P = kmoc._interior(front[i], front[i + 1])
            except RuntimeError:
                continue
            if not (np.isfinite(P.x) and np.isfinite(P.r)) or P.r < -1e-9:
                continue
            new.append(P)
        if len(new) < 2:
            break
        ax = kmoc._axis(new[0])
        front = [ax] + new
        all_pts.extend(front)
    return all_pts


def march_wall_from_dataline(dataline: DataLine, theta_tol_deg: float = 0.01,
                             M_target: float | None = None,
                             eps_M_rel: float = 5e-3, max_iter: int = 4000,
                             n_corr: int = 2) -> dict:
    r"""**未検証 — production では使わないこと** (モジュール docstring 「状態」節参照)。
    放射源流厳密解照合で壁座標が厳密解と有意に乖離することを確認済みで、
    原因は未特定。ここでの実装は継続調査のための足場として残してある。

    データ線 (軸→壁 順) から下流へ `_march_stations` (moc_kernel) と同型の
    はしご法でマーチする (§4 準拠)。`_wall_prescribed` (既知壁形状) の代わりに
    `_wall_cancel` (非反射・形状未知) を使う以外は完全に同じ構造。

    各反復: (1) `_wall_cancel(front[-2], front[-1])` で新しい壁点 W を作る、
    (2) W を種に、旧フロントを壁側から軸側へ interior() のはしごで下る
    (`front[j]` (j=len-2..0) を A、はしご先端を B として `interior(front[j], B)`)、
    (3) はしご終端 (軸側最先端) を `_axis()` で軸反射、(4) 新フロント =
    reversed([軸反射点, はしご各段..., W]) (長さ不変)。

    **2026-08-15、3 回の実装訂正を経て確定**: 試行 1 (前フロントを interior() で
    全体縮約し front[-2] を毎回使い回す) は θ が厳密解よりずっと速く減衰した
    (15 反復で 10°→5.16°、厳密解は 9.45°) — 原因は A 自身が縮約バイアスを
    引きずること。試行 2 (軸反射点を前フロントの近軸点 `front[1]` から作る
    「はしご」) は、本データ線の軸端が壁端よりずっと下流 ($x_{axis}>x_D$、
    C⁻ が浅い角度で降りるため) という座標の向きを見落とし、`front[1]` が
    実質データ線自身の軸端と同じ位置になり退化した。**採用方式**は
    `_march_stations` を素直に模倣し、はしごを**壁側 (front[-2]) から軸側
    (front[0]) へ**下る — 既存の検証済み構造をそのまま流用するのが最も安全。

    戻り値: `wall` (n,4 [x,r,θ,M])、`axis` (n,3 [x,M,θ])、`all_pts`、
    `terminated_by` ("theta_tol" / "max_iter")、`n_iter`。
    """
    import warnings
    warnings.warn(
        "march_wall_from_dataline は未検証 (放射源流厳密解照合で有意な乖離あり)。"
        "出力される壁座標を production で使わないこと — moc_cancel.py モジュール"
        "docstring の「状態」節を参照。", stacklevel=2)
    g = dataline.gamma
    kmoc = _CancelKernel(gamma=g, delta=1.0, n_corr=n_corr)
    front = dataline.to_pts()
    if len(front) < 3:
        raise ValueError("データ線の点数が不足 (単位過程に最低 3 点必要)")
    all_pts = list(front)
    wall = [(front[-1].x, front[-1].r, front[-1].th, front[-1].M)]
    axis = [(front[0].x, front[0].M, front[0].th)]
    theta_tol = np.deg2rad(theta_tol_deg)
    terminated_by = "max_iter"
    for it in range(max_iter):
        try:
            W = kmoc._wall_cancel(front[-2], front[-1])
        except RuntimeError as exc:
            raise RuntimeError(f"march が iter={it} で壁点構築に失敗: {exc}") from exc
        if not (np.isfinite(W.x) and np.isfinite(W.r)):
            raise RuntimeError(f"march が iter={it} で壁点が非有限")
        new = [W]
        ok = True
        for j in range(len(front) - 2, -1, -1):
            try:
                P = kmoc._interior(front[j], new[-1])
            except RuntimeError:
                ok = False
                break
            if not (np.isfinite(P.x) and np.isfinite(P.r)) or P.r < -1e-9:
                ok = False
                break
            new.append(P)
        if not ok:
            raise RuntimeError(f"march が iter={it} ではしご構築に失敗")
        ax = kmoc._axis(new[-1])
        new.append(ax)
        front = new[::-1]
        all_pts.extend(front)
        wall.append((W.x, W.r, W.th, W.M))
        axis.append((ax.x, ax.M, ax.th))
        theta_ok = abs(W.th) <= theta_tol
        M_ok = M_target is None or abs(W.M - M_target) <= eps_M_rel * M_target
        if theta_ok and M_ok:
            terminated_by = "theta_tol"
            break
    return {"wall": np.asarray(wall), "axis": np.asarray(axis),
            "all_pts": all_pts, "terminated_by": terminated_by,
            "n_iter": len(wall) - 1}


def cancellation_contour(dataline: DataLine, theta_tol_deg: float = 0.01,
                         M_target: float | None = None, eps_M_rel: float = 5e-3,
                         max_iter: int = 4000) -> np.ndarray:
    """`march_wall_from_dataline` を実行し壁点列 (n,4 [x,r,θ,M]) のみ返す薄いラッパ。"""
    res = march_wall_from_dataline(dataline, theta_tol_deg=theta_tol_deg,
                                   M_target=M_target, eps_M_rel=eps_M_rel,
                                   max_iter=max_iter)
    if res["terminated_by"] != "theta_tol":
        raise RuntimeError(
            f"cancellation march が max_iter={max_iter} まで θ_tol 未達 "
            f"(最終壁角 {np.rad2deg(res['wall'][-1, 2]):.4f}°)")
    return res["wall"]
