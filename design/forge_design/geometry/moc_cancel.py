r"""wall design rule (Phase W4 の暫定近似): D 発 C⁻ データ線 → 下流壁の生成。

plan: plans/active/tooling-nozzle-walldriven-chain.md §4 (W4 cancellation MOC の仕様)。
**このモジュールの現行実装は「軸対称 wave-cancellation MOC」ではない
(2026-08-15、外部レビュー [Codex] で指摘・確認)。CONTUR の「壁は MOC の従属
出力」思想の**近似**にとどまる — 詳細は下記「限界 (重要・未解決)」参照。**

## 何をしているか

$D$ 発 C⁻ データ線 (CFD 実場から抽出) から内点充填 (`build_field`) で内側
($r\le r_D$) の場を厳密に再構成しつつ、壁 ($r>r_D$) は「隣接内点由来の $J^+$」
と「**目標出口一様状態 $M=M_{\rm target},\theta=0$ が持つ $J^-=\nu(M_{\rm
target})$ を全壁点に一定値として直接代入**」を組み合わせて生成する
(`_wall_cancel`)。

## 限界 (重要・未解決、2026-08-15 外部レビューで確認)

**$J^-_W=\nu(M_{\rm target})$ を全壁点で一定と置くのは、出口から真に逆走した
結果ではない。** 特性線を引くには局所 $M,\theta$ (→ $\mu=\sin^{-1}(1/M)$、
$C^\mp:\,dr/dx=\tan(\theta\mp\mu)$) が要るが、$J^-$ の値だけでは特性線の
向きは分からない。真に出口から逆走するなら、出口の $C^-$ を実際に引き、
反対族 $C^+$ との交点で局所状態を更新しながら**逐次**壁点まで到達させる
必要がある。現行実装はこれをせず、出口の $J^-$ を**直接、位置によらず**
壁へ課しているだけ ($J^-_W=\nu(M_{\rm target})$)。

さらに**軸対称では厳密に誤り**である: $J^-$ は特性線に沿って一般に

$$dJ^- = S^-(M,\theta,r)\,ds$$

という幾何学的源項 (軸対称の 3D 発散効果) を持ち、平面流のように保存量
ではない。真に出口から逆走するなら

$$J^-_W = \nu(M_{\rm target}) - \int_W^{e} S^-(M,\theta,r)\,ds$$

でなければならないが、積分内の $M,\theta$ は経路上未知であり単純代入では
閉じない (循環参照— forward/backward MOC の同時マッチングか反復解法が要る)。
現行の $J^-_W=\nu(M_{\rm target})$ は、**この積分をゼロと置く平面流的な
近似**であり、$D$ の位置・$R_t$・$\theta_D$ が同じでも軸対称成分が大きい
条件では誤差が大きくなりうる。

**終了条件も自己充足的**: `march_wall_from_dataline` は壁点自身の $M_W,
\theta_W$ が目標に達したかで打ち切っており、**出口断面全体の一様性**
($M(r)\approx M_{\rm target}$, $\theta(r)\approx0$、$r$ 方向全域) は
検証していない。壁終端だけ目標に一致しても、断面全体が一様である保証はない。

**したがって現状の正しい位置づけ**: 「$J^-=\nu(M_{\rm target})$ は軸対称
cancellation MOC の解ではなく、平面 simple-wave を模した**暫定的な壁設計則**」
である。W4 を「軸対称 wave-cancellation MOC 完成」として扱ってはならない。

**真に厳密化する経路 (未実装、次段階)**:
1. 出口の終端特性線から**源項込みで** backward MOC し、$D$ 側の forward
   MOC (`build_field`) と突き合わせて matching する (出口位置・出口断面・
   壁・$\theta_D$ を同時調整する自由境界問題になる)。
2. 仮の壁 (現行の近似解でよい) で全場を解き (実 CFD)、**出口全断面の
   非一様性を目的関数として** 壁または $\theta_D$ を反復調整する
   (W5 の outer loop と統合可能)。

## 数理 (内点充填)

内点充填 (`build_field`) は、データ線が「両方の Riemann 不変量
$J^+=\theta-\nu$, $J^-=\theta+\nu$ が各点で既知」な Cauchy データ (CFD 実場
から両方読み取っているため) であることを使う。**役割の割当て**: 軸側 (r 小)
を A (C+ 担体)、壁側 (r 大) を B (C- 担体) とする。これは `moc_kernel.
KernelMOC` の楔充填・壁ステーションマーチが実際に使っている割当て
(`front[j]` = 既存側 = A) と一致し、`moc_inverse.fill()` が縦の starting
line 区間で使う固定割当て (逆) とは異なる — 後者は近軸への丸め誤差で文書化
済みの限界を持つ (moc_inverse.py 冒頭の楔欠損の記録)。この部分 (内点充填)
は放射源流厳密解照合で相対誤差 ~1e-4 と検証済み — **限界は壁生成の閉じ方
(上記) にある**。

## 実装の健全性 (数値的には正しく動く — 上記モデル近似とは別の話)

`_wall_cancel`/`march_wall_from_dataline` は、**上記の「平面近似の設計則」
としては**数値的に正しく動作する (退化チェック機械精度、流線自己整合性
~1e-8、格子収束 <0.5%):

- **退化チェック** (平面単純波 $\delta=0$、$J^-$ 厳密一定の解析場 — この場合
  近似は厳密解と一致する): 隣接内点 A が既に目標 $J^-$ と同じ値を持つとき、
  壁が A の状態を厳密に再現する (機械精度)。
- **物理的整合性** (放射源流データ線 + 任意の $M_{\rm target}$): march した壁は
  流線条件 $dr/dx\approx\tan\theta$ を満たす (自己整合、誤差 ~1e-8)、θ は
  単調に減衰、格子収束する (データ線点数 41→161 で最終壁座標の変化 <0.5%)。
- **既知の実装上の注意 (解決済み)**: (1) `front[-2]` を `_wall_cancel` の A と
  はしごの最初の段の両方に使うと march が 1 反復で固定点に収束して停止する
  バグがあった (front[-2] の二重使用) — はしごは `front[-3]` から始める。
  (2) $\theta_{\rm tol}$ をデータ線のステップ幅より厳しく指定すると、march が
  収束点 ($\theta=0$) を通り過ぎて振動域に入り最終的に破綻する (実測: 0.05°
  指定で $\theta$ が $+2°\to-1.5°\to+0.2°\to\ldots$ と振動し近軸で破綻)。
  **θ=0 交差を検出したら線形補間で即座に打ち切る安全策を実装済み**
  (`march_wall_from_dataline` 内、`theta_tol_deg` の値によらず振動域に入らない)。
  (3) 実 CFD データ (run_0043, $M_D\approx2.26\to M_{\rm target}=4$、$J^-$
  ギャップ ~20°) では、目標が局所値から大きく離れているとはしごが発散する
  ケースを確認 (iter5 で $r\approx-33$)。原因未特定 — 上記モデル近似の限界
  ($D$ 直後で一気に目標 $J^-$ へ飛ぶ非物理な第一歩) が主因の可能性が高い。
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

    # P0 診断 (壁直近 5% は境界補間混入を避け除外)。**線形補間を使う**
    # (2026-08-15 修正: 最近傍探索だと格子間隔由来の偽ノイズが乗り dP0_rel が
    # 過大評価される — 実測 0.196%→2.3% と一桁変わった。`analyze_walldriven_w3.py`
    # の p0_along_path と同じ LinearNDInterpolator 方式に統一)。
    from scipy.interpolate import LinearNDInterpolator
    if "P" in fld:
        pts = np.c_[fld["x"], fld["r"]]
        iP = LinearNDInterpolator(pts, fld["P"])
        Pv = iP(xs, rs)
    else:
        Pv = np.full_like(xs, np.nan)
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


class _CancelKernel(KernelMOC):
    """`_interior`/`_axis` (内点充填用) + `_wall_cancel` (壁専用単位過程)。"""

    def _wall_cancel(self, A: _Pt, W_prev: _Pt, nu_target: float) -> _Pt:
        r"""壁点 W を A (C+ 担体, $J^+$ の運び手) と目標出口一様状態
        ($M=M_{\rm target}, \theta=0$) が持つ反射側不変量 $\nu_{\rm target}$
        (=$\nu(M_{\rm target})$、呼び出し側で計算して渡す) から構成する
        (2026-08-15 訂正版 — plan §4 参照。旧版は $J^-$ を直前壁点から引き継ぐ
        方式だったが、C⁻ 1 本 (+軸鏡映) だけでは壁の位置が構造的に未決定と判明
        したため撤回した)。

        $J^+_W = J^+_A + S_p$ (通常の C+ compatibility)。
        $J^-_W = \nu_{\rm target}$ (**位置に依らない定数** — 出口一様状態が
        すべての壁点に共通して持つ反射側不変量。軸対称源項補正は付けない:
        古典 MLN の straightening section の簡略化であり、$D$ 直後で厳密ではない
        ことを承知の上での近似)。位置は C+ レイ (A 発) と「壁の接線方向レイ」
        (W_prev 発、傾き $\tan(\tfrac12(\theta_{prev}+\theta_W))$ — 通常の C⁻ の
        $\mp\mu$ 項を持たない、$W_{prev}\to W$ が流線であることの表現) の交点。
        """
        g, d = self.g, self.delta
        thW, nuW = W_prev.th, nu_target - W_prev.th  # 初期推定
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
            Sp = -d * fA / np.cos(0.5 * (A.th + thW) + 0.5 * (_mu(A.M) + muW)) * (xW - A.x)
            Jp = (A.th - A.nu) + Sp
            Jm = nu_target
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
            if not (np.isfinite(P.x) and np.isfinite(P.r)) or P.r < -1e-4:
                continue
            P.r = max(P.r, 0.0)
            new.append(P)
        if len(new) < 2:
            break
        ax = kmoc._axis(new[0])
        front = [ax] + new
        all_pts.extend(front)
    return all_pts


def march_wall_from_dataline(dataline: DataLine, M_target: float,
                             theta_tol_deg: float = 0.01,
                             eps_M_rel: float = 5e-3, max_iter: int = 4000,
                             n_corr: int = 2) -> dict:
    r"""データ線 (軸→壁 順) から下流へ `_march_stations` (moc_kernel) と同型の
    はしご法でマーチする (§4 準拠、2026-08-15 訂正版)。`_wall_prescribed`
    (既知壁形状) の代わりに `_wall_cancel` (目標出口一様状態 $M_{\rm target}$ が
    閉じる) を使う以外は完全に同じ構造。**`M_target` は必須** — 1 本の C⁻ (+ 軸
    鏡映) だけでは壁の位置が構造的に未決定であること (下記「訂正の経緯」) が
    判明したため、閉じるための第二の情報として要求する。

    各反復: (1) `_wall_cancel(front[-2], front[-1], nu_target)` で新しい壁点 W
    を作る (front[-2] を A、front[-1]=直前の壁点を W_prev として消費)、
    (2) W を種に、旧フロントの**残り** (front[-3] から軸側へ、front[j] を A、
    はしご先端を B として `interior(front[j], B)`) をはしごで下る (front[-2]
    は (1) で既に消費済みなので**再利用しない** — 二重使用すると
    `interior(front[-2], W)` が W とほぼ同一の点を再生成し march が 1 反復で
    固定点に収束して停止するバグがあった)、(3) はしご終端 (軸側最先端) を
    `_axis()` で軸反射、(4) 新フロント = reversed([軸反射点, はしご各段...,
    W]) (長さ不変)。

    **訂正の経緯 (2026-08-15)**: 当初「壁点で追加の compatibility は課さない」
    としていたが、C⁻ 1 本上では同族の $J^-$ が実質拘束されており (2D 平面なら
    厳密に一定)、真の 2 特性族 Cauchy データにならない — 壁の位置は依存領域の
    外側にあり構造的に未決定だった (3 通りの「非反射」試行がいずれも放射源流
    厳密解照合で乖離した根本原因、ユーザ指摘)。訂正後は $J^-_W=\nu(M_{\rm
    target})$ を壁の閉じ条件として明示的に課す (`_wall_cancel` の docstring
    参照)。平面単純波 ($\delta=0$、$J^-$ 厳密一定) での退化検証 (内点 A が
    既に目標と同じ $J^-$ を持つ場合、壁が厳密に $A$ と同じ状態を再現する) は
    機械精度で通過する。$\theta$ が 0 を交差したら振動域に入る前に線形補間で
    打ち切る (モジュール docstring の「検証状態」参照)。

    戻り値: `wall` (n,4 [x,r,θ,M])、`axis` (n,3 [x,M,θ])、`all_pts`、
    `terminated_by` ("theta_tol" / "max_iter")、`n_iter`。
    """
    g = dataline.gamma
    nu_target = float(pm_nu(M_target, g))
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
            W = kmoc._wall_cancel(front[-2], front[-1], nu_target)
        except RuntimeError as exc:
            raise RuntimeError(f"march が iter={it} で壁点構築に失敗: {exc}") from exc
        if not (np.isfinite(W.x) and np.isfinite(W.r)):
            raise RuntimeError(f"march が iter={it} で壁点が非有限")
        new = [W]
        ok = True
        # front[-2] は既に _wall_cancel の A として消費済みなので、はしごは
        # front[-3] から始める (前回のバグ: front[-2] の二重使用で
        # interior(front[-2], W) が W 自身とほぼ同一の点を再生成し march が
        # 1 反復で固定点に収束していた)。
        for j in range(len(front) - 3, -1, -1):
            try:
                P = kmoc._interior(front[j], new[-1])
            except RuntimeError:
                ok = False
                break
            if not (np.isfinite(P.x) and np.isfinite(P.r)) or P.r < -1e-4:
                ok = False
                break
            P.r = max(P.r, 0.0)
            new.append(P)
        if not ok:
            raise RuntimeError(f"march が iter={it} ではしご構築に失敗")
        ax = kmoc._axis(new[-1])
        new.append(ax)
        front = new[::-1]
        all_pts.extend(front)
        # θ=0 交差検出 (2026-08-15 追加): theta_tol_deg がステップ幅より厳しいと
        # 収束点を通り過ぎて振動域に入り最終的に破綻する (実測: 0.05° 指定で
        # ステップ幅~0.1-0.5° の march が θ=+2°→-1.5°→+0.2°... と振動し
        # 41 反復で近軸の interior() が破綻)。**W が前回の壁点と符号反転したら
        # 通り過ぎたとみなし、線形補間で θ=0 点を求めて即座に打ち切る** —
        # theta_tol_deg の指定値によらず振動域へ入らない安全策。
        W_prev_th = wall[-1][2]
        if np.sign(W.th) != np.sign(W_prev_th) and W_prev_th != 0.0:
            s = W_prev_th / (W_prev_th - W.th)
            xc = wall[-1][0] + s * (W.x - wall[-1][0])
            rc = wall[-1][1] + s * (W.r - wall[-1][1])
            Mc = wall[-1][3] + s * (W.M - wall[-1][3])
            wall.append((xc, rc, 0.0, Mc))
            axis.append((ax.x, ax.M, ax.th))
            terminated_by = "theta_tol"
            break
        wall.append((W.x, W.r, W.th, W.M))
        axis.append((ax.x, ax.M, ax.th))
        theta_ok = abs(W.th) <= theta_tol
        M_ok = abs(W.M - M_target) <= eps_M_rel * M_target
        if theta_ok and M_ok:
            terminated_by = "theta_tol"
            break
    return {"wall": np.asarray(wall), "axis": np.asarray(axis),
            "all_pts": all_pts, "terminated_by": terminated_by,
            "n_iter": len(wall) - 1}


def cancellation_contour(dataline: DataLine, M_target: float,
                         theta_tol_deg: float = 0.01, eps_M_rel: float = 5e-3,
                         max_iter: int = 4000) -> np.ndarray:
    """`march_wall_from_dataline` を実行し壁点列 (n,4 [x,r,θ,M]) のみ返す薄いラッパ。"""
    res = march_wall_from_dataline(dataline, M_target, theta_tol_deg=theta_tol_deg,
                                   eps_M_rel=eps_M_rel, max_iter=max_iter)
    if res["terminated_by"] != "theta_tol":
        raise RuntimeError(
            f"cancellation march が max_iter={max_iter} まで θ_tol 未達 "
            f"(最終壁角 {np.rad2deg(res['wall'][-1, 2]):.4f}°)")
    return res["wall"]
