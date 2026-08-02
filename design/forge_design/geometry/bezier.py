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
