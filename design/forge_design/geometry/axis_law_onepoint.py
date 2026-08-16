r"""軸 Mach 則 D: 端点アンカー + 内部補間点 1 点の C⁴ 区分 5 次則 (`OnePointC4AxisLaw`)。

計画: plans/active/tooling-nozzle-axislaw-onepoint.md。A/B/C 比較 (axislaw-smoothness) の
教訓「軸の数学的滑らかさを最小化しても MOC 後の壁・特性線網は改善しない」を受け、
**自由度を増やさず、上下流アンカーを固定したまま、内部補間点 1 点だけで膨張配分を動かす**。

- 無次元設計変数: $L_c=x_E-x_A$ (外側)、$\xi_P=(x_P-x_A)/L_c$、$\eta_P=(M_P-M_A)/(M_d-M_A)$。
  $x_P$ は候補自身の $L_c$ で無次元化する。
- 曲線: $[x_A,x_P]$ (長さ $L_1$) と $[x_P,x_E]$ (長さ $L_2$) の 2 区間、各区間を局所
  $t\in[0,1]$ の 5 次 Bernstein 基底で表す。条件 12 (A 側 3 / E 側 3 / P の C⁰..C⁴ 5 /
  $M(x_P)=M_P$ 1) = 係数 12 で一意。
- **P は Bezier 制御点ではなく、曲線が実際に通る内部補間点**。

Bernstein 係数 ($b_i$: 区間 1、$c_i$: 区間 2):
  A 側: $b_0=M_A$, $b_1=b_0+L_1M'_A/5$, $b_2=2b_1-b_0+L_1^2M''_A/20$
  E 側: $c_5=M_d$, $c_4=c_5$, $c_3=c_5$ ($M'=M''=0$)
  P:    $b_5=c_0=M_P$
  残り $(b_3,b_4,c_1,c_2)$ を P での C¹..C⁴ 連続 4 式で解く (4×4 線形、区間長のべきを含む):
    $5\Delta^1b_4/L_1=5\Delta^1c_0/L_2$, $20\Delta^2b_3/L_1^2=20\Delta^2c_0/L_2^2$,
    $60\Delta^3b_2/L_1^3=60\Delta^3c_0/L_2^3$, $120\Delta^4b_1/L_1^4=120\Delta^4c_0/L_2^4$
  ($\Delta^k$ は前進差分)。
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import comb


def _bernstein_eval(coef: np.ndarray, t, order: int = 0, L: float = 1.0):
    r"""Bernstein 係数 `coef` (次数 n=len-1) の曲線の $t\in[0,1]$ での値、または
    物理座標 $x=x_0+Lt$ に関する `order` 階微分 (差分公式を再帰適用し $L^{order}$ で割る)。"""
    c = np.asarray(coef, dtype=float)
    n = len(c) - 1
    for _ in range(order):
        if n == 0:
            return np.zeros_like(np.asarray(t, dtype=float))
        c = n * np.diff(c)
        n -= 1
    t = np.asarray(t, dtype=float)
    B = np.zeros_like(t)
    for i in range(n + 1):
        B += c[i] * comb(n, i) * t ** i * (1.0 - t) ** (n - i)
    return B / L ** order


class OnePointC4AxisLaw:
    r"""端点アンカー + 内部補間点 1 点 (ξ_P, η_P) の C⁴ 区分 5 次軸 Mach 則 (D 案)。

    インターフェースは `QuinticHermiteAxisLaw`/`KnotQuinticAxisLaw` と同じ
    (`__call__`, `deriv`, `gates`, `x_A`, `x_E`, `L_c`) に加え `x_P`, `M_P`, `xi_P`,
    `eta_P`, `L_1`, `L_2`, `ell_P` (= x_P - x_A)。"""

    def __init__(self, x_A: float, L_c: float, M_A: float, Mp_A: float, Mpp_A: float,
                 M_d: float, xi_P: float, eta_P: float) -> None:
        vals = [x_A, L_c, M_A, Mp_A, Mpp_A, M_d, xi_P, eta_P]
        if not np.isfinite(vals).all():
            raise ValueError("OnePointC4AxisLaw: 非有限の入力")
        if L_c <= 0.0:
            raise ValueError(f"OnePointC4AxisLaw: L_c = {L_c:.6g} <= 0")
        if M_d <= M_A:
            raise ValueError(f"OnePointC4AxisLaw: M_d = {M_d:.6g} <= M_A = {M_A:.6g}")
        if not (0.0 < xi_P < 1.0):
            raise ValueError(f"OnePointC4AxisLaw: xi_P = {xi_P:.6g} は (0,1) 内でなければならない")
        if not (0.0 < eta_P < 1.0):
            raise ValueError(f"OnePointC4AxisLaw: eta_P = {eta_P:.6g} は (0,1) 内でなければならない")
        self.x_A, self.L_c = float(x_A), float(L_c)
        self.M_A, self.Mp_A, self.Mpp_A, self.M_d = float(M_A), float(Mp_A), float(Mpp_A), float(M_d)
        self.xi_P, self.eta_P = float(xi_P), float(eta_P)
        self.x_E = self.x_A + self.L_c
        self.x_P = self.x_A + self.xi_P * self.L_c
        self.M_P = self.M_A + self.eta_P * (self.M_d - self.M_A)
        self.L_1 = self.x_P - self.x_A
        self.L_2 = self.x_E - self.x_P
        self.ell_P = self.L_1
        L1, L2 = self.L_1, self.L_2
        # 直接決まる係数
        b0 = self.M_A
        b1 = b0 + L1 * self.Mp_A / 5.0
        b2 = 2.0 * b1 - b0 + L1 ** 2 * self.Mpp_A / 20.0
        b5 = self.M_P
        c0 = self.M_P
        c3 = c4 = c5 = self.M_d
        # 未知 u = [b3, b4, c1, c2] の 4×4 線形系 (C1..C4 at P)
        # C1: (b5-b4)/L1 = (c1-c0)/L2
        # C2: (b5-2b4+b3)/L1^2 = (c2-2c1+c0)/L2^2
        # C3: (b5-3b4+3b3-b2)/L1^3 = (c3-3c2+3c1-c0)/L2^3
        # C4: (b5-4b4+6b3-4b2+b1)/L1^4 = (c4-4c3+6c2-4c1+c0)/L2^4
        A = np.array([
            [0.0, -1.0 / L1, -1.0 / L2, 0.0],
            [1.0 / L1 ** 2, -2.0 / L1 ** 2, 2.0 / L2 ** 2, -1.0 / L2 ** 2],
            [3.0 / L1 ** 3, -3.0 / L1 ** 3, -3.0 / L2 ** 3, 3.0 / L2 ** 3],
            [6.0 / L1 ** 4, -4.0 / L1 ** 4, 4.0 / L2 ** 4, -6.0 / L2 ** 4],
        ])
        rhs = np.array([
            -b5 / L1 - c0 / L2,
            -b5 / L1 ** 2 + c0 / L2 ** 2,
            (-b5 + b2) / L1 ** 3 + (c3 - c0) / L2 ** 3,
            (-b5 + 4.0 * b2 - b1) / L1 ** 4 + (c4 - 4.0 * c3 + c0) / L2 ** 4,
        ])
        self.cond_number = float(np.linalg.cond(A))
        b3, b4, c1, c2 = np.linalg.solve(A, rhs)
        self._b = np.array([b0, b1, b2, b3, b4, b5])
        self._c = np.array([c0, c1, c2, c3, c4, c5])
        self.constraint_residual = float(np.max(np.abs(A @ np.array([b3, b4, c1, c2]) - rhs)))

    # -- 評価 -----------------------------------------------------------------
    def _eval(self, x, order: int):
        x = np.asarray(x, dtype=float)
        scalar = x.ndim == 0
        x = np.atleast_1d(x)
        out = np.empty_like(x)
        m1 = x <= self.x_P
        m2 = ~m1
        if m1.any():
            t1 = (x[m1] - self.x_A) / self.L_1
            out[m1] = _bernstein_eval(self._b, t1, order, self.L_1)
        if m2.any():
            t2 = np.clip((x[m2] - self.x_P) / self.L_2, None, 1.0)
            out[m2] = _bernstein_eval(self._c, t2, order, self.L_2)
        if order > 0:
            out = np.where(x > self.x_E, 0.0, out)
        else:
            out = np.where(x > self.x_E, self.M_d, out)
        return out[0] if scalar else out

    def __call__(self, x):
        return self._eval(x, 0)

    def deriv(self, x, order: int = 1):
        """d^k M/dx^k。x > x_E では 0。"""
        return self._eval(x, order)

    # -- 設計ゲート -----------------------------------------------------------
    def gates(self, n: int = 4001) -> dict:
        """端点拘束残差 / P での 0〜4 階 jump / min M' (hard: >=0) / M'' 符号反転 /
        max|M'''| / ∫(M''')² dx (区間ごと Gauss、M''' は t の 2 次なので厳密)。"""
        e = 1e-7
        end_res = max(
            abs(float(self(self.x_A)) - self.M_A),
            abs(float(self.deriv(self.x_A, 1)) - self.Mp_A),
            abs(float(self.deriv(self.x_A, 2)) - self.Mpp_A),
            abs(float(self(self.x_E)) - self.M_d),
            abs(float(self.deriv(self.x_E, 1))),
            abs(float(self.deriv(self.x_E, 2))),
            abs(float(self(self.x_P)) - self.M_P))
        jumps = []
        for k in range(5):
            l = float(self._eval(self.x_P - e, k)) if k else float(self(self.x_P - e))
            r = float(self._eval(self.x_P + e, k)) if k else float(self(self.x_P + e))
            jumps.append(r - l)
        x = np.linspace(self.x_A, self.x_E, n)
        Mp = self.deriv(x, 1)
        Mpp = self.deriv(x, 2)
        Mppp = self.deriv(x, 3)
        scale = max(1.0, float(np.max(np.abs(Mp))))
        mono_ok = bool(np.all(Mp >= -1e-10 * scale))
        sg = np.sign(Mpp)
        sg = sg[np.abs(Mpp) > 1e-10 * max(1.0, float(np.max(np.abs(Mpp))))]
        n_flip = int(np.sum(np.diff(sg) != 0.0)) if len(sg) > 1 else 0
        gx, gw = leggauss(6)
        J = 0.0
        for coef, L, x0 in ((self._b, self.L_1, self.x_A), (self._c, self.L_2, self.x_P)):
            t = 0.5 * (gx + 1.0)
            v = _bernstein_eval(coef, t, 3, L)
            J += 0.5 * L * float(np.sum(gw * v * v))
        v = []
        if not mono_ok:
            v.append(f"M' < 0 の区間あり (min M' = {float(Mp.min()):.4g})")
        if end_res > 1e-8:
            v.append(f"端点拘束残差 {end_res:.2e}")
        if max(abs(j) for j in jumps) > 1e-6 * max(1.0, float(np.max(np.abs(Mppp)))):
            v.append(f"P での導関数 jump {max(abs(j) for j in jumps):.2e}")
        return {
            "monotone_ok": mono_ok, "min_dMdx": float(Mp.min()), "max_dMdx": float(Mp.max()),
            "mpp_sign_changes": n_flip, "mpp_quality_ok": n_flip <= 1,
            "max_abs_Mppp": float(np.max(np.abs(Mppp))), "J_axis": float(J),
            "endpoint_residual": end_res, "P_jumps": jumps, "P_jump_max": max(abs(j) for j in jumps),
            "cond_number": self.cond_number, "constraint_residual": self.constraint_residual,
            "xi_P": self.xi_P, "eta_P": self.eta_P, "x_P": self.x_P, "M_P": self.M_P,
            "ell_P": self.ell_P, "L_1": self.L_1, "L_2": self.L_2,
            "violations": v,
        }

    @classmethod
    def monotone_feasible(cls, x_A, L_c, M_A, Mp_A, Mpp_A, M_d, xi_P, eta_P) -> bool:
        """(ξ_P, η_P) で構成した則が全域 M'≥0 か (MOC 前の篩い)。"""
        try:
            return cls(x_A, L_c, M_A, Mp_A, Mpp_A, M_d, xi_P, eta_P).gates()["monotone_ok"]
        except (ValueError, np.linalg.LinAlgError):
            return False
