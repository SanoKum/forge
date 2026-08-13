"""中心線マッハ分布の Bézier 表現。

CP の x 座標を等間隔固定 (x(t)=x0+L·t) とするため M(x) は次数 n の Bernstein
多項式で、端点条件 1 個 = CP 1 個の拘束として陽に解ける
(plans/active/tooling-nozzle-design-tool.md §4.6(b))。両端 C2 拘束なら自由 CP
数 k = n-5。単調 CP → 単調カーブ (十分条件)。
"""
from __future__ import annotations

import numpy as np


def _bernstein_eval(cp: np.ndarray, t: np.ndarray) -> np.ndarray:
    """de Casteljau (ベクトル化)。"""
    t = np.asarray(t, dtype=float)
    b = np.repeat(cp[:, None], t.size, axis=1).astype(float)
    n = len(cp) - 1
    for r in range(n):
        b = b[:-1] * (1.0 - t) + b[1:] * t
    return b[0]


class MachBezier:
    """M(x) = Bernstein(cp; t), t=(x-x0)/L。"""

    def __init__(self, x0: float, x1: float, cp) -> None:
        if x1 <= x0:
            raise ValueError("x1 <= x0")
        self.x0 = float(x0)
        self.x1 = float(x1)
        self.cp = np.asarray(cp, dtype=float)
        if self.cp.size < 2:
            raise ValueError("CP は 2 個以上")

    # -- 構築 ---------------------------------------------------------------
    @classmethod
    def from_constraints(cls, x0: float, x1: float, start, free_cp, end=None) -> "MachBezier":
        """端点条件 + 自由 CP から構築する。

        start / end: (M,) | (M, dM/dx) | (M, dM/dx, d2M/dx2)。end=None は拘束なし。
        次数 n = len(start)+len(free_cp)+len(end)-1。
        """
        start = tuple(start)
        endc = tuple(end) if end is not None else ()
        free = list(free_cp)
        n = len(start) + len(free) + len(endc) - 1
        if n < 1:
            raise ValueError("次数が不足")
        L = x1 - x0
        cp = np.empty(n + 1)
        # 始端: P0..P2 (M, M', M'' の順で確定)
        cp[0] = start[0]
        if len(start) >= 2:
            cp[1] = cp[0] + start[1] * L / n
        if len(start) >= 3:
            cp[2] = start[2] * L * L / (n * (n - 1)) + 2.0 * cp[1] - cp[0]
        # 終端: Pn..Pn-2
        if len(endc) >= 1:
            cp[n] = endc[0]
        if len(endc) >= 2:
            cp[n - 1] = cp[n] - endc[1] * L / n
        if len(endc) >= 3:
            cp[n - 2] = endc[2] * L * L / (n * (n - 1)) + 2.0 * cp[n - 1] - cp[n]
        # 自由 CP を中央に充填
        i0 = len(start)
        for k, v in enumerate(free):
            cp[i0 + k] = v
        return cls(x0, x1, cp)

    @classmethod
    def fit_free_cp(cls, x0: float, x1: float, start, n_free: int, end,
                    ref_fn, n_sample: int = 400) -> "MachBezier":
        r"""端点拘束の下で **自由 CP を参照曲線へ拘束付き最小二乗射影**する (B10)。

        Bernstein 基底は CP について線形なので、端点拘束で確定する CP を右辺へ
        移せば自由 CP は 3 変数の線形最小二乗で**決定的に**決まる (手調整不要):
        $\min_{P_{free}}\sum_k\bigl(\sum_i P_iB_{i,n}(t_k)-M_{ref}(x_k)\bigr)^2$。

        用途: target の起点を $x_0\to x_{reach,CFD}$ へ移す際、旧起点用に調整
        された自由 CP 値をインデックスずらしで流用すると制御多角形が壊れる
        (実測: 旧値 2.0/2.9/3.6 を order=2 へ流用して $P_2{=}3.67\to P_3{=}2.00$
        のジグザグ)。参照曲線への射影なら旧 target の形を保ったまま新起点・新拘束
        へ載せ替えられる。`ref_fn`: 参照 $M(x)$ (呼び出し可能)。
        """
        start, endc = tuple(start), tuple(end)
        n = len(start) + int(n_free) + len(endc) - 1
        L = x1 - x0
        base = cls.from_constraints(x0, x1, start, [0.0] * int(n_free), endc)
        fixed = np.zeros(n + 1, dtype=bool)
        fixed[:len(start)] = True
        fixed[n + 1 - len(endc):] = True
        i_free = np.where(~fixed)[0]
        xs = np.linspace(x0, x1, int(n_sample))
        t = (xs - x0) / L
        # 基底行列 B[k,i] = B_{i,n}(t_k)
        from math import comb
        B = np.stack([comb(n, i) * t ** i * (1 - t) ** (n - i) for i in range(n + 1)], axis=1)
        rhs = np.asarray([float(np.atleast_1d(ref_fn(float(x)))[0]) for x in xs])
        rhs = rhs - B[:, fixed] @ base.cp[fixed]
        sol, *_ = np.linalg.lstsq(B[:, i_free], rhs, rcond=None)
        cp = base.cp.copy()
        cp[i_free] = sol
        return cls(x0, x1, cp)

    # -- 評価 ---------------------------------------------------------------
    def _t(self, x):
        return (np.asarray(x, dtype=float) - self.x0) / (self.x1 - self.x0)

    def __call__(self, x):
        return _bernstein_eval(self.cp, self._t(x))

    def deriv(self, x, order: int = 1):
        """d^k M/dx^k。導関数 Bézier の CP は差分 × n/L。"""
        cp = self.cp.copy()
        L = self.x1 - self.x0
        n = len(cp) - 1
        for _ in range(order):
            cp = np.diff(cp) * n / L
            n -= 1
            if n < 0:
                return np.zeros_like(np.asarray(x, dtype=float))
        return _bernstein_eval(cp, self._t(x))

    # -- 検査 (事前フィルタ) -------------------------------------------------
    def cp_monotone(self, tol: float = 0.0) -> bool:
        """CP 単調 (→ カーブ単調の十分条件)。"""
        return bool(np.all(np.diff(self.cp) >= -tol))

    def max_dMdx(self, nsample: int = 400) -> float:
        x = np.linspace(self.x0, self.x1, nsample)
        return float(np.max(self.deriv(x, 1)))
