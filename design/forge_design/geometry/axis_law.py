r"""軸中心 Mach 分布の 5 次 Hermite 表現 (axis-Mach チェーン)。

計画: plans/active/tooling-nozzle-axismach-chain.md §5.1。
モード F の Bézier (自由 CP 3 個 — B10-c で始端 C² と出口 C² の両立不能が確定) を
両端 6 条件をちょうど消費する 5 次 Hermite に差し替え、残る自由度を加速区間長
L_c の 1 個に集約する。

s = (x - x_A)/L_c ∈ [0,1],  M(s) = Σ a_i s^i。両端 6 条件
(M_A, M'_A, M''_A / M_d, 0, 0) の閉形式:

    a0 = M_A,  a1 = L_c M'_A,  a2 = (1/2) L_c² M''_A,  ΔM = M_d − M_A
    a3 = 10ΔM − 6a1 − 3a2
    a4 = −15ΔM + 8a1 + 3a2
    a5 = 6ΔM − 3a1 − a2

x > x_E (= x_A + L_c) では M = M_d 一定 (一様軸流への C² 接続は端条件が保証)。
x < x_A は定義域外 (呼び出し側が実測軸曲線を使う — B10 の帳簿方針)。

設計ゲート (原方針 §8, §28):
- hard: M'(x) ≥ 0 全域 (⇒ M_A ≤ M ≤ M_d、overshoot 0 が自動保証)
- 品質指標: M'' の符号反転回数 ≤ 1 (M' の一峰性)
端点条件は構成的に満足する (テストで機械精度確認)。
"""
from __future__ import annotations

import numpy as np


class QuinticHermiteAxisLaw:
    """M(x): [x_A, x_E] は 5 次 Hermite、x > x_E は M_d 一定。"""

    def __init__(self, x_A: float, L_c: float, M_A: float, Mp_A: float,
                 Mpp_A: float, M_d: float) -> None:
        if not np.isfinite([x_A, L_c, M_A, Mp_A, Mpp_A, M_d]).all():
            raise ValueError("QuinticHermiteAxisLaw: 非有限の入力")
        if L_c <= 0.0:
            raise ValueError(f"QuinticHermiteAxisLaw: L_c = {L_c:.6g} ≤ 0")
        if M_d <= M_A:
            raise ValueError(f"QuinticHermiteAxisLaw: M_d = {M_d:.6g} ≤ M_A = {M_A:.6g}")
        self.x_A, self.L_c = float(x_A), float(L_c)
        self.M_A, self.Mp_A, self.Mpp_A = float(M_A), float(Mp_A), float(Mpp_A)
        self.M_d = float(M_d)
        self.x_E = self.x_A + self.L_c
        a1 = self.L_c * self.Mp_A
        a2 = 0.5 * self.L_c ** 2 * self.Mpp_A
        dM = self.M_d - self.M_A
        self._a = np.array([
            self.M_A, a1, a2,
            10.0 * dM - 6.0 * a1 - 3.0 * a2,
            -15.0 * dM + 8.0 * a1 + 3.0 * a2,
            6.0 * dM - 3.0 * a1 - a2,
        ])

    # -- 評価 -----------------------------------------------------------------
    def _s(self, x):
        return (np.asarray(x, dtype=float) - self.x_A) / self.L_c

    def __call__(self, x):
        s = np.clip(self._s(x), None, 1.0)   # x > x_E は M_d (s=1 で厳密に M_d)
        return np.polynomial.polynomial.polyval(s, self._a)

    def deriv(self, x, order: int = 1):
        """d^k M/dx^k。x > x_E では 0。"""
        s = self._s(x)
        c = self._a.copy()
        for _ in range(order):
            c = c[1:] * np.arange(1, len(c))
        out = np.polynomial.polynomial.polyval(np.clip(s, None, 1.0), c) / self.L_c ** order
        return np.where(s > 1.0, 0.0, out)

    # -- 設計ゲート -----------------------------------------------------------
    def gates(self, n: int = 4001) -> dict:
        """設計ゲート一式。violations が空なら合格。

        M' ≥ 0 (hard) / M'' 符号反転 ≤ 1 (品質) / 最大 dM/dx (情報)。"""
        x = np.linspace(self.x_A, self.x_E, n)
        Mp = self.deriv(x, 1)
        Mpp = self.deriv(x, 2)
        scale = max(1.0, float(np.max(np.abs(Mp))))
        mono_ok = bool(np.all(Mp >= -1e-12 * scale))
        sg = np.sign(Mpp)
        sg = sg[np.abs(Mpp) > 1e-10 * max(1.0, float(np.max(np.abs(Mpp))))]
        n_flip = int(np.sum(np.diff(sg) != 0.0)) if len(sg) > 1 else 0
        v = []
        if not mono_ok:
            v.append(f"M' < 0 の区間あり (min M' = {float(Mp.min()):.4g})")
        return {
            "monotone_ok": mono_ok,
            "mpp_sign_changes": n_flip,
            "mpp_quality_ok": n_flip <= 1,
            "max_dMdx": float(Mp.max()),
            "min_dMdx": float(Mp.min()),
            "violations": v,
        }

    @classmethod
    def admissible_Lc_range(cls, x_A: float, M_A: float, Mp_A: float,
                            Mpp_A: float, M_d: float, Lc_hi: float = 200.0,
                            n_scan: int = 400, tol: float = 1e-5) -> tuple:
        r"""M'(x) ≥ 0 を満たす L_c の連結許容窓 (lo, hi) を返す。

        **重要 (2026-08-15 実装時の発見)**: 単調ゲートの構造は原方針 §10 の想定と
        **逆向き**である。

        - $L_c\to0$ では $a_1,a_2\to0$ で純粋な 0→ΔM smoothstep
          ($M'\propto s^2(1-s)^2\ge0$) に退化するため、**短い側は常に単調合格**。
          実用下限は勾配の大きさ (max dM/dx → 壁角・曲率) として MOC 後の壁 QA
          が課す。
        - $a_2=\frac12 L_c^2 M''_A$ が $L_c^2$ で成長するため、$M''_A>0$ の接続では
          **長すぎる $L_c$ が単調性を破る** (序盤に上がりすぎ → 終盤で $M'<0$)。
          実測例: B10 アンカー $(M,M',M'')=(1.946,0.290,0.190)$ の窓上限は
          $L_c\approx9.8$ で、B8 相当長 12.3 は不合格。

        したがって単純 bisection (単調性前提) ではなく log 格子スキャン + 端の
        bisection で連結窓を確定する。下端がスキャン床 (1e-2) のままなら「下限なし」
        と読む。窓が Lc_hi で切れずに続く場合は hi = Lc_hi。許容 L_c が無ければ
        ValueError。"""
        def ok(Lc: float) -> bool:
            try:
                return cls(x_A, Lc, M_A, Mp_A, Mpp_A, M_d).gates()["monotone_ok"]
            except ValueError:
                return False

        grid = np.geomspace(1e-2, Lc_hi, n_scan)
        mask = np.array([ok(L) for L in grid])
        if not mask.any():
            raise ValueError("許容 L_c が存在しない (接続条件が強すぎる)")
        i0 = int(np.argmax(mask))                      # 最初の許容点
        # 下端: grid[i0-1] (不合格) と grid[i0] の間を bisection
        if i0 == 0:
            lo = float(grid[0])
        else:
            a, b = float(grid[i0 - 1]), float(grid[i0])
            while b - a > tol * max(1.0, b):
                m = 0.5 * (a + b)
                if ok(m):
                    b = m
                else:
                    a = m
            lo = b
        # 上端: i0 以降で最初に不合格になる点
        rest = mask[i0:]
        if rest.all():
            hi = float(Lc_hi)
        else:
            j = i0 + int(np.argmin(rest))              # 最初の不合格 index
            a, b = float(grid[j - 1]), float(grid[j])
            while b - a > tol * max(1.0, b):
                m = 0.5 * (a + b)
                if ok(m):
                    a = m
                else:
                    b = m
            hi = a
        return lo, hi

    @classmethod
    def min_Lc(cls, x_A: float, M_A: float, Mp_A: float, Mpp_A: float,
               M_d: float, Lc_hi: float = 200.0) -> float:
        """許容窓の下限 (原方針 §10 の「許される最小 L_c」)。"""
        return cls.admissible_Lc_range(x_A, M_A, Mp_A, Mpp_A, M_d, Lc_hi)[0]
