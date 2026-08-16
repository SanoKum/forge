r"""軸 Mach 則 B / C (滑らかさ比較, A6 knot 則の代替候補)。

計画: plans/active/tooling-nozzle-axislaw-smoothness.md。knot 則 (`axis_law.KnotQuinticAxisLaw`,
A) は knot で $M,M',M''$ が C² だが $M'''$ が跳ぶ。ここでは:

- **B** (`MonotoneBSplineAxisLaw`): $M(x)$ を直接、次数 5 (quintic) の clamped B-spline
  (単純内部ノットで構成的 C⁴) で表現し、制御点を「端点条件 + 単調性 (制御多角形非減少)」
  の下で $J_{\rm axis}=\int(M''')^2dx$ 最小化により決める凸 QP。
- **C** (`NonnegDnuBSplineAxisLaw`): Prandtl–Meyer 勾配 $q=d\nu/dx$ を次数 4 (quartic, C³)
  の B-spline (非負制御点で構成的 $q\ge0$) で表現し、$J_\nu=\int(q'')^2dx$ 最小化。
  $\nu(x)$ は $q$ の B-spline 反導関数、$M(x)=$ `gas.mach_of_nu(nu(x))`。

両方とも `QuinticHermiteAxisLaw`/`KnotQuinticAxisLaw` と同じ `__call__`/`deriv`/`gates`
インターフェースを持ち、`inverse_design(target=law)` にそのまま渡せる。
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import BSpline, make_interp_spline
from scipy.optimize import LinearConstraint, minimize

from .moc_kernel import dnu_dM, pm_nu


# --- B-spline 共通ヘルパ ------------------------------------------------------
def _clamped_knots(x_A: float, x_E: float, n_interior: int, k: int,
                   cluster_alpha: float = 0.5) -> np.ndarray:
    r"""clamped 節点ベクトル (端 $k{+}1$ 重複)。内部節点はスロート側にクラスタ
    ($s_i=(i/(n{+}1))^\alpha$, $\alpha<1$ で $x_A$ 寄りが密)。$n_{\rm interior}=0$ なら
    単一 Bezier セグメント (次数 $k$ の多項式 1 本、B の $M''(x_E)$ 自由・単一 5 次と等価)。"""
    if n_interior < 0:
        raise ValueError("n_interior は 0 以上")
    if n_interior > 0:
        s = (np.arange(1, n_interior + 1) / (n_interior + 1)) ** cluster_alpha
        interior = x_A + s * (x_E - x_A)
    else:
        interior = np.empty(0)
    return np.concatenate([np.full(k + 1, x_A), interior, np.full(k + 1, x_E)])


def _basis_matrix(t: np.ndarray, k: int, x: np.ndarray, nu: int = 0) -> np.ndarray:
    """基底関数 (の nu 階微分) を (len(x), n_cp) 行列で返す (one-hot 制御ベクトル評価。
    このリポジトリの `wall_axismach.lsq_bspline_wall._row` と同じ流儀 — インストール済み
    scipy は `BSpline.design_matrix(nu=)` 未対応のため)。"""
    n_cp = len(t) - k - 1
    x = np.atleast_1d(np.asarray(x, dtype=float))
    B = np.empty((len(x), n_cp))
    eye = np.eye(n_cp)
    for i in range(n_cp):
        spl = BSpline(t, eye[i], k, extrapolate=False)
        B[:, i] = spl(x, nu) if nu else spl(x)
    return np.nan_to_num(B, nan=0.0)


def _gram_matrix(t: np.ndarray, k: int, nu_order: int, n_gauss: int = 6) -> np.ndarray:
    r"""3 階 (または nu_order 階) 微分基底の Gram 行列 $H_{ij}=\int B_i^{(\nu)}B_j^{(\nu)}dx$
    (各節点区間を n_gauss 点 Gauss–Legendre で厳密積分 — 区間内は次数 $k-\nu$ の多項式なので
    積の次数 $2(k-\nu)$ に対し十分な点数を取れば厳密)。"""
    knots = np.unique(t)
    gx, gw = leggauss(n_gauss)
    xs_list, ws_list = [], []
    for a, b in zip(knots[:-1], knots[1:]):
        if b - a < 1e-14:
            continue
        xm, xr = 0.5 * (a + b), 0.5 * (b - a)
        xs_list.append(xm + xr * gx)
        ws_list.append(xr * gw)
    xs = np.concatenate(xs_list)
    ws = np.concatenate(ws_list)
    B = _basis_matrix(t, k, xs, nu=nu_order)
    return B.T @ (ws[:, None] * B)


def _cumsum_transform(n: int) -> np.ndarray:
    """`c_i = z_0 + z_1 + ... + z_i` (z_0 自由, z_{i>0}>=0 が c の単調非減少と等価)
    の変換行列 T (c = T z)。"""
    T = np.zeros((n, n))
    T[:, 0] = 1.0
    for i in range(1, n):
        T[i:, i] = 1.0
    return T


def _solve_bounded_qp(H: np.ndarray, A_eq: np.ndarray, b_eq: np.ndarray,
                      lb: np.ndarray, ub: np.ndarray, T: np.ndarray | None = None,
                      ridge: float = 1e-10, A_extra: np.ndarray | None = None,
                      lb_extra: np.ndarray | None = None,
                      ub_extra: np.ndarray | None = None) -> np.ndarray:
    r"""凸 QP を**単純境界制約**に帰着させて解く: min 0.5 c^T H c s.t. A_eq c = b_eq,
    lb <= c <= ub (`T` を渡すと $c=Tz$ で変数変換し、bounds は $z$ 空間で与える —
    単調性を「制御点の隣接差分 $u_i=c_i-c_{i-1}\ge0$」に変換すると一般不等式行列より
    大幅に解きやすい、実測)。

    一般不等式 (`LinearConstraint`) より境界制約 (`Bounds`) の方が trust-constr の
    内部射影が厳密・高速で収束しやすい (実測: 一般不等式は n_interior<9 で
    非収束、境界制約化後は n_interior=3 から安定収束)。
    """
    from scipy.optimize import Bounds, lsq_linear
    n = H.shape[0]
    if T is None:
        T = np.eye(n)
    Hz = T.T @ (H + ridge * np.eye(n)) @ T
    Az = A_eq @ T
    # 初期点: 境界を厳密に満たす重み付き最小二乗 (等式は強く重み付け)
    W = np.sqrt(1.0e6)
    A_aug = np.vstack([W * Az, 1e-3 * np.eye(n)])
    b_aug = np.concatenate([W * b_eq, np.zeros(n)])
    z0 = lsq_linear(A_aug, b_aug, bounds=(lb, ub), max_iter=3000).x

    def fun(z):
        return 0.5 * float(z @ Hz @ z)

    def jac(z):
        return Hz @ z

    cons = [LinearConstraint(Az, lb=b_eq, ub=b_eq)]
    if A_extra is not None:
        cons.append(LinearConstraint(A_extra @ T,
                                     lb=lb_extra if lb_extra is not None else -np.inf,
                                     ub=ub_extra if ub_extra is not None else np.inf))
    bnds = Bounds(lb, ub)
    best = None
    for method, opts in (
        ("trust-constr", {"maxiter": 4000, "gtol": 1e-13, "xtol": 1e-15}),
        ("SLSQP", {"maxiter": 1000, "ftol": 1e-16}),
    ):
        kw = {"hess": lambda z: Hz} if method == "trust-constr" else {}
        res = minimize(fun, z0, jac=jac, method=method, constraints=cons,
                       bounds=bnds, options=opts, **kw)
        viol_eq = float(np.max(np.abs(Az @ res.x - b_eq))) if len(b_eq) else 0.0
        viol_bnd = float(max(0.0, np.max(lb - res.x), np.max(res.x - ub)))
        viol_extra = 0.0
        if A_extra is not None:
            v = A_extra @ T @ res.x
            lo = lb_extra if lb_extra is not None else np.full_like(v, -np.inf)
            hi = ub_extra if ub_extra is not None else np.full_like(v, np.inf)
            viol_extra = float(max(0.0, np.max(lo - v), np.max(v - hi)))
        viol = viol_eq + viol_bnd + viol_extra
        if viol_eq < 1e-7 and viol_bnd < 1e-9 and viol_extra < 1e-7:
            return T @ res.x
        if best is None or viol < best[0]:
            best = (viol, res.x)
    raise RuntimeError(f"QP が収束しない (最良 viol={best[0]:.3e})")


# --- B: monotone quintic B-spline (M(x) 直接) --------------------------------
class MonotoneBSplineAxisLaw:
    r"""$M(x)$ を次数 5 (quintic) の clamped B-spline (単純内部ノットで構成的 C⁴) で
    直接表現する。制御点は端点条件 (等式) + 単調性 (制御多角形非減少、不等式) の下で
    $J_{\rm axis}=\int(M''')^2dx$ を最小化して決める (凸 QP)。

    `exit_curvature`: 'hard' (既定, $M''(x_E)=0$ を等式拘束) / 'soft' (2 次ペナルティ
    `exit_curvature_weight` で目的関数に加える — 過拘束の疑いがあるときの比較用)。

    **前倒れ (front-loading) 対策の spread 拘束** (2026-08-16 実測): 端点条件 + 単調性 +
    $J_{\rm axis}$ 最小化**だけ**では、最適解が $x_E$ よりずっと手前 (実測 M6/R3/$L_c$45 で
    $x\approx20$、$L_c$ の半分未満) で $M_d$ に到達し残りを平坦域として済ませる
    ("前倒れ") — 平坦域は $M'''=0$ で目的関数を下げるため、むしろ**最適化に好まれる**。
    これは単一 quintic が起こした fold と同じ病理 (スロート近傍に曲げが集中) を
    再現する (`AxisMachCFDWall.validate` の spline リンギングで実測検出)。
    $M(x_A+f_x L_c)\le M_A+f_M(M_d-M_A)$ ($f_x$=`spread_x_frac`, $f_M$=`spread_M_frac`)
    を不等式拘束として追加し、中間点で早期到達しすぎないよう縛る (knot 則の $M_K$ と
    同じ発想を不等式として一般化 — 特定の値に固定しない)。`spread_M_frac=None` で
    無効化 (前倒れの実害を確認する比較用)。"""

    def __init__(self, x_A: float, x_E: float, M_A: float, Mp_A: float, Mpp_A: float,
                 M_d: float, n_interior: int = 15, k: int = 5,
                 cluster_alpha: float = 0.5, exit_curvature: str = "hard",
                 exit_curvature_weight: float = 1.0e4,
                 spread_x_frac: float = 0.5,
                 spread_M_frac: float | None = 0.75) -> None:
        if not np.isfinite([x_A, x_E, M_A, Mp_A, Mpp_A, M_d]).all():
            raise ValueError("MonotoneBSplineAxisLaw: 非有限の入力")
        if x_E <= x_A:
            raise ValueError(f"x_E={x_E:.6g} <= x_A={x_A:.6g}")
        if M_d <= M_A:
            raise ValueError(f"M_d={M_d:.6g} <= M_A={M_A:.6g}")
        if exit_curvature not in ("hard", "soft"):
            raise ValueError("exit_curvature は 'hard' か 'soft'")
        self.x_A, self.x_E, self.M_A, self.Mp_A, self.Mpp_A = \
            float(x_A), float(x_E), float(M_A), float(Mp_A), float(Mpp_A)
        self.M_d, self.k, self.n_interior = float(M_d), int(k), int(n_interior)
        self.exit_curvature = exit_curvature
        t = _clamped_knots(self.x_A, self.x_E, n_interior, k, cluster_alpha)
        self.t = t
        n_cp = len(t) - k - 1
        H = _gram_matrix(t, k, nu_order=3)
        # 端点等式: value/1st/2nd at x_A、value/1st at x_E (2nd は hard/soft で分岐)
        rowsA = _basis_matrix(t, k, np.array([self.x_A]), nu=0)[0], \
            _basis_matrix(t, k, np.array([self.x_A]), nu=1)[0], \
            _basis_matrix(t, k, np.array([self.x_A]), nu=2)[0]
        rowsE0 = _basis_matrix(t, k, np.array([self.x_E]), nu=0)[0]
        rowsE1 = _basis_matrix(t, k, np.array([self.x_E]), nu=1)[0]
        rowsE2 = _basis_matrix(t, k, np.array([self.x_E]), nu=2)[0]
        A_eq = [rowsA[0], rowsA[1], rowsA[2], rowsE0, rowsE1]
        b_eq = [self.M_A, self.Mp_A, self.Mpp_A, self.M_d, 0.0]
        if exit_curvature == "hard":
            A_eq.append(rowsE2); b_eq.append(0.0)
        else:
            H = H + exit_curvature_weight * np.outer(rowsE2, rowsE2)
        A_eq = np.asarray(A_eq); b_eq = np.asarray(b_eq)
        T = _cumsum_transform(n_cp)          # c = T z, z_0 自由・z_{i>0}>=0 <=> c 非減少
        lb = np.concatenate([[-np.inf], np.zeros(n_cp - 1)])
        A_extra = ub_extra = None
        self.spread_x_frac, self.spread_M_frac = spread_x_frac, spread_M_frac
        if spread_M_frac is not None:
            self.x_mid = self.x_A + spread_x_frac * (self.x_E - self.x_A)
            self.M_cap = self.M_A + spread_M_frac * (self.M_d - self.M_A)
            A_extra = _basis_matrix(t, k, np.array([self.x_mid]), nu=0)
            ub_extra = np.array([self.M_cap])
        c = _solve_bounded_qp(H, A_eq, b_eq, lb, np.full(n_cp, np.inf), T=T,
                              A_extra=A_extra, ub_extra=ub_extra)
        self._spl = BSpline(t, c, k, extrapolate=False)
        self.c = c
        self.J_axis = float(_quad_energy(self._spl, self.x_A, self.x_E, order=3))

    def __call__(self, x):
        x = np.asarray(x, dtype=float)
        return np.where(x > self.x_E, self.M_d, self._spl(np.minimum(x, self.x_E)))

    def deriv(self, x, order: int = 1):
        x = np.asarray(x, dtype=float)
        v = self._spl(np.minimum(x, self.x_E), order)
        return np.where(x > self.x_E, 0.0, v)

    def gates(self, n: int = 4001) -> dict:
        x = np.linspace(self.x_A, self.x_E, n)
        Mp, Mpp = self.deriv(x, 1), self.deriv(x, 2)
        scale = max(1.0, float(np.max(np.abs(Mp))))
        mono_ok = bool(np.all(Mp >= -1e-9 * scale))
        sg = np.sign(Mpp)
        sg = sg[np.abs(Mpp) > 1e-10 * max(1.0, float(np.max(np.abs(Mpp))))]
        n_flip = int(np.sum(np.diff(sg) != 0.0)) if len(sg) > 1 else 0
        v = [] if mono_ok else [f"M' < 0 の区間あり (min M' = {float(Mp.min()):.4g})"]
        return {"monotone_ok": mono_ok, "mpp_sign_changes": n_flip,
                "mpp_quality_ok": n_flip <= 1, "max_dMdx": float(Mp.max()),
                "min_dMdx": float(Mp.min()), "violations": v}


# --- C: 非負 dν/dx B-spline (Prandtl-Meyer 表現) -----------------------------
class NonnegDnuBSplineAxisLaw:
    r"""$q(x)=d\nu/dx \ge0$ を次数 4 (quartic, 単純ノットで C³) の B-spline (非負制御点)
    で表現し、$J_\nu=\int(q'')^2dx$ を最小化。$\nu(x)$ は $q$ の反導関数、
    $M(x)=\nu^{-1}(\nu(x))=$ `gas.mach_of_nu`。

    `exit_slope`: 'hard' (既定, $q'(x_E)=0$) / 'soft' (2 次ペナルティ)。

    **ノット配置は B と逆向きが効率的** (2026-08-16 実測, LP 実行可能性走査):
    $q(x_A)$ はアンカー由来で大きく (M6/R3 で 0.23、区間平均 0.036 の ~6.5 倍)、
    $q'(x_A)>0$ ($M''_A>0$ 由来) で**さらに立ち上がってから**長い尾部でほぼ 0 に
    張り付く必要がある。スロート側にノットを寄せる (B の `cluster_alpha<1`) と
    **尾部の解像度が足りず**逆に実行不能域が広がる (alpha=0.15 は n_interior=120
    でも不可)。均等〜出口寄り (`cluster_alpha>=1`) だと n_interior=6–8 まで
    小さくできる。既定は 1.0 (均等)。

    **前倒れ対策の spread 拘束** (2026-08-16 実測、B と同じ病理): 端点+積分+平滑化だけの
    QP は $q$ を早期に 0 へ落として長い平坦尾部で済ませる解を好む (尾部は $q''=0$ で
    $J_\nu$ をむしろ下げる) — 実際にこの形で逆 MOC に通すと単一 quintic と同じ fold が
    再現される (θ_w>μ_w がスロート近傍に再集中、`AxisMachCFDWall.validate` で検出)。
    $\nu(x_A+f_x L_c)\le\nu(M_A)+f_M(\nu(M_d)-\nu(M_A))$ を追加して中間点で早期到達を防ぐ
    (B と同じ発想、$\nu$ 空間で)。`spread_M_frac=None` で無効化。"""

    def __init__(self, x_A: float, x_E: float, M_A: float, Mp_A: float, Mpp_A: float,
                 M_d: float, gas, n_interior: int = 20, k: int = 4,
                 cluster_alpha: float = 1.0, exit_slope: str = "hard",
                 exit_slope_weight: float = 1.0e4,
                 spread_x_frac: float = 0.5,
                 spread_M_frac: float | None = 0.75) -> None:
        if not np.isfinite([x_A, x_E, M_A, Mp_A, Mpp_A, M_d]).all():
            raise ValueError("NonnegDnuBSplineAxisLaw: 非有限の入力")
        if x_E <= x_A:
            raise ValueError(f"x_E={x_E:.6g} <= x_A={x_A:.6g}")
        if M_d <= M_A:
            raise ValueError(f"M_d={M_d:.6g} <= M_A={M_A:.6g}")
        if Mp_A <= 0.0:
            raise ValueError(f"Mp_A={Mp_A:.6g} <= 0")
        if exit_slope not in ("hard", "soft"):
            raise ValueError("exit_slope は 'hard' か 'soft'")
        self.x_A, self.x_E, self.M_A, self.Mp_A, self.Mpp_A = \
            float(x_A), float(x_E), float(M_A), float(Mp_A), float(Mpp_A)
        self.M_d, self.gas, self.k = float(M_d), gas, int(k)
        nu_A = float(pm_nu(self.M_A, gas))
        nu_d = float(pm_nu(self.M_d, gas))
        nuM = float(dnu_dM(self.M_A, gas, 1))
        nuMM = float(dnu_dM(self.M_A, gas, 2))
        q_A = nuM * self.Mp_A
        qp_A = nuMM * self.Mp_A ** 2 + nuM * self.Mpp_A
        if q_A <= 0.0:
            raise ValueError(f"NonnegDnuBSplineAxisLaw: q(x_A)={q_A:.4g} <= 0 "
                             "(アンカー M'_A が小さすぎる)")
        self.nu_A, self.nu_d, self.q_A, self.qp_A = nu_A, nu_d, q_A, qp_A
        t = _clamped_knots(self.x_A, self.x_E, n_interior, k, cluster_alpha)
        self.t = t
        n_cp = len(t) - k - 1
        H = _gram_matrix(t, k, nu_order=2)
        rowsA0 = _basis_matrix(t, k, np.array([self.x_A]), nu=0)[0]
        rowsA1 = _basis_matrix(t, k, np.array([self.x_A]), nu=1)[0]
        rowsE0 = _basis_matrix(t, k, np.array([self.x_E]), nu=0)[0]
        rowsE1 = _basis_matrix(t, k, np.array([self.x_E]), nu=1)[0]
        row_int = _basis_integral_row(t, k)          # ∫ B_i dx (各基底の全区間積分)
        A_eq = [rowsA0, rowsA1, rowsE0, row_int]
        b_eq = [q_A, qp_A, 0.0, nu_d - nu_A]
        if exit_slope == "hard":
            A_eq.append(rowsE1); b_eq.append(0.0)
        else:
            H = H + exit_slope_weight * np.outer(rowsE1, rowsE1)
        A_eq = np.asarray(A_eq); b_eq = np.asarray(b_eq)
        A_extra = ub_extra = None
        self.spread_x_frac, self.spread_M_frac = spread_x_frac, spread_M_frac
        if spread_M_frac is not None:
            self.x_mid = self.x_A + spread_x_frac * (self.x_E - self.x_A)
            self.nu_cap = nu_A + spread_M_frac * (nu_d - nu_A)
            A_extra = _basis_partial_integral_row(t, k, self.x_mid)[None, :]
            ub_extra = np.array([self.nu_cap - nu_A])
        c = _solve_bounded_qp(H, A_eq, b_eq, np.zeros(n_cp), np.full(n_cp, np.inf),
                              A_extra=A_extra, ub_extra=ub_extra)
        self._qspl = BSpline(t, c, k, extrapolate=False)
        self.c = c
        self.J_nu = float(_quad_energy(self._qspl, self.x_A, self.x_E, order=2))
        F = self._qspl.antiderivative()
        self._F0 = float(F(self.x_A))
        self._F = F
        # M(x)=ν^{-1}(ν(x)): mach_of_nu の既定実装 (6000 点区分線形テーブル) を経由すると
        # 高階微分の診断が C0 キンクを増幅してノイズだらけになる (実測: M'' 符号反転 300+ 回)。
        # 物理量 M(ν) 自体は滑らかなので、既存テーブルに次数 5 補間スプラインを当て直し
        # (`_mach_of_nu_spline`, 生データそのものへの当て直し = 数値差分ではない)、
        # M(x)=Mspl(ν(x)) と Faà di Bruno の合成関数微分で M',M'',M''' を厳密に評価する
        # (q, q', q'' は q スプラインの厳密微分、Mspl', Mspl'', Mspl''' は Mspl の厳密微分)。
        self._Mspl = _mach_of_nu_spline(gas, self.nu_d)

    def nu(self, x):
        x = np.asarray(x, dtype=float)
        xc = np.clip(x, self.x_A, self.x_E)
        v = self.nu_A + (self._F(xc) - self._F0)
        return np.where(x > self.x_E, self.nu_d, v)

    def q(self, x, order: int = 0):
        x = np.asarray(x, dtype=float)
        xc = np.minimum(x, self.x_E)
        v = self._qspl(xc, order) if order else self._qspl(xc)
        v = np.nan_to_num(v, nan=0.0)
        return np.where(x > self.x_E, 0.0, v)

    def __call__(self, x):
        x = np.asarray(x, dtype=float)
        xc = np.minimum(x, self.x_E)
        v = self._Mspl(self.nu(xc))
        return np.where(x > self.x_E, self.M_d, v)

    def deriv(self, x, order: int = 1):
        x = np.asarray(x, dtype=float)
        xc = np.minimum(x, self.x_E)
        nux = self.nu(xc)
        q0 = self.q(xc)
        Mp1 = self._Mspl(nux, 1)
        if order == 1:
            out = Mp1 * q0
        elif order == 2:
            Mp2 = self._Mspl(nux, 2)
            out = Mp2 * q0 ** 2 + Mp1 * self.q(xc, 1)
        elif order == 3:
            Mp2 = self._Mspl(nux, 2)
            Mp3 = self._Mspl(nux, 3)
            q1, q2 = self.q(xc, 1), self.q(xc, 2)
            out = Mp3 * q0 ** 3 + 3.0 * Mp2 * q0 * q1 + Mp1 * q2
        else:
            raise ValueError("deriv: order は 1..3")
        return np.where(x > self.x_E, 0.0, out)

    def gates(self, n: int = 4001) -> dict:
        x = np.linspace(self.x_A, self.x_E, n)
        Mp, Mpp = self.deriv(x, 1), self.deriv(x, 2)
        scale = max(1.0, float(np.max(np.abs(Mp))))
        mono_ok = bool(np.all(Mp >= -1e-6 * scale))     # M(x) は fit spline 経由で緩めの床
        sg = np.sign(Mpp)
        sg = sg[np.abs(Mpp) > 1e-10 * max(1.0, float(np.max(np.abs(Mpp))))]
        n_flip = int(np.sum(np.diff(sg) != 0.0)) if len(sg) > 1 else 0
        v = [] if mono_ok else [f"M' < 0 の区間あり (min M' = {float(Mp.min()):.4g})"]
        return {"monotone_ok": mono_ok, "mpp_sign_changes": n_flip,
                "mpp_quality_ok": n_flip <= 1, "max_dMdx": float(Mp.max()),
                "min_dMdx": float(Mp.min()), "violations": v}


def _mach_of_nu_spline(gas, nu_max: float):
    r"""$\nu \to M$ の**滑らかな** (次数 5) 補間スプライン。既存 `gas.mach_of_nu`
    (`np.interp` 区分線形) は 6000 点ごとに $C^0$ キンクを持ち、それを合成関数微分
    (`NonnegDnuBSplineAxisLaw.deriv`) に使うと $M''$ が雑音まみれになる (実測)。
    semi-perfect はガスモデル内部の生テーブル $(M,\nu)$ に直接次数 5 スプラインを
    当て直す (`gas.nu`/`mach_of_nu` の丸めを経由しない、キャッシュ 1 個)。
    CPG (float γ、または `GasCPG` — テーブルなし) は `pm_mach` (Newton 収束、滑らか) を
    密サンプルしてスプライン化する。"""
    from .moc_kernel import pm_mach_vec
    if hasattr(gas, "_nu") and hasattr(gas, "_M"):
        spl = getattr(gas, "_M_of_nu_spline", None)
        if spl is None:
            spl = make_interp_spline(gas._nu, gas._M, k=5)
            gas._M_of_nu_spline = spl
        return spl
    nus = np.linspace(0.0, max(nu_max, 1e-3) * 1.05, 4001)
    Ms = pm_mach_vec(nus, gas)
    return make_interp_spline(nus, Ms, k=5)


def _basis_partial_integral_row(t: np.ndarray, k: int, x_upper: float) -> np.ndarray:
    r"""各基底関数の $[t_0, x_{\rm upper}]$ 部分積分 $\int_{t_0}^{x_{\rm upper}} B_i dx$
    (反導関数は scipy 標準 `BSpline.antiderivative` — 手で次数上げ公式を導出しない)。"""
    n_cp = len(t) - k - 1
    row = np.empty(n_cp)
    eye = np.eye(n_cp)
    for i in range(n_cp):
        spl = BSpline(t, eye[i], k, extrapolate=False)
        F = spl.antiderivative()
        row[i] = float(F(x_upper) - F(t[0]))
    return row


def _basis_integral_row(t: np.ndarray, k: int) -> np.ndarray:
    """各基底関数の $[t_0, t_{-1}]$ 全区間積分 (`_basis_partial_integral_row` の x_upper=t[-1])。"""
    return _basis_partial_integral_row(t, k, t[-1])


def _quad_energy(spl: BSpline, x_A: float, x_E: float, order: int, n_gauss: int = 8) -> float:
    """$\int (\mathrm{spl}^{(order)})^2\,dx$ を節点区間ごとの Gauss–Legendre で評価。"""
    knots = np.unique(spl.t)
    knots = knots[(knots >= x_A - 1e-12) & (knots <= x_E + 1e-12)]
    gx, gw = leggauss(n_gauss)
    total = 0.0
    for a, b in zip(knots[:-1], knots[1:]):
        if b - a < 1e-14:
            continue
        xm, xr = 0.5 * (a + b), 0.5 * (b - a)
        xs = xm + xr * gx
        v = spl(xs, order)
        total += xr * float(np.sum(gw * v * v))
    return total
