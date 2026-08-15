r"""wall-driven チェーン (①風洞) のスロート域多項式壁 (無次元, r* = 1)。

現行モード F (目標軸マッハ → 逆 MOC) と並存する新チェーンの Phase W1
(plans/active/tooling-nozzle-walldriven-chain.md)。凍結円弧を持たず、
U (上流直管端) → T (スロート) → D (最大壁角点) を低自由度多項式で設計する。
スロート曲率は幾何区分でなく条件 r''(0) = 1/R_t として入る。

- U→T: 両端 6 条件の quintic Hermite。U 側 (r_U, 0, 0) → T 側 (1, 0, 1/R_t)。
  単調収縮 ⇔ μ = L_U²/(R_t Δr) ≤ 20 (必要十分。勾配閉形式
  dr/dx = Δr/(2L_U) ξ²(ξ-1)[5μξ - 3μ - 60ξ + 60] から。数値検証 2026-08-15)。
- T→D: 壁勾配 q = r' を cubic Hermite で規定した 4 次壁。
  q(ξ) = (L_D/R_t) ξ(1-ξ)² + tanθ_D ξ²(3-2ξ),  ξ = x/L_D。
  D 半径は従属: r_D = 1 + L_D²/(12 R_t) + L_D tanθ_D / 2。
  途中変曲なし (r'' ≥ 0 全区間) ⇔ L_D ≤ 3 R_t tanθ_D (必要十分)。

r, r', r'', r''' は全て解析式で持ち、θ, κ, dκ/ds も解析評価する (数値微分しない)。
T で κ は 1/R_t に連続だが κ' は一般に不連続 — 跳び量は診断値として常時出せる
(kappa_prime_jump_at_throat)。滑らかさの許容可否は CFD で判定する方針
(調査 §2.4: 「局所を滑らかにすれば波が消える」という直感には反例がある)。

座標: x = 軸 (スロート = 0, U は x = -L_U, D は x = +L_D), r = 半径。
角度は全て radian。

エラー方針 (2026-08-15 レビュー反映): **構築時 ValueError = 入力が不正**
(非有限・非正の長さ/半径/曲率半径・角度域外・出口半径非正) — 係数生成前に検査し
NaN 座標が下流 (メッシュ生成) へ流れることを防ぐ。**validate() = 幾何は定義できるが
設計条件違反** (μ > 20、L_D > 3 R_t tanθ_D、λ 窓外、任意の max|dκ/ds| 上限超過)。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .wall import _hermite_quintic


def _check_inputs(name: str, positive: dict, finite_only: dict | None = None,
                  angles: dict | None = None) -> None:
    """係数生成前の入力検査。不正なら ValueError (validate() には回さない)。

    positive: 有限かつ > 0 を要求。finite_only: 有限のみ要求。
    angles: 有限かつ 0 < θ < π/2 (radian) を要求。
    """
    for k, v in {**positive, **(finite_only or {}), **(angles or {})}.items():
        if not np.isfinite(v):
            raise ValueError(f"{name}: {k} = {v!r} が有限でない")
    for k, v in positive.items():
        if v <= 0.0:
            raise ValueError(f"{name}: {k} = {v:.6g} ≤ 0")
    for k, v in (angles or {}).items():
        if not 0.0 < v < np.pi / 2.0:
            raise ValueError(f"{name}: {k} = {v:.6g} rad ∉ (0, π/2)")


def _poly_eval(c, xi, deriv=0):
    """ξ 基底多項式 c[j] ξ^j の deriv 階微分 (ξ 微分) を評価する。"""
    xi = np.asarray(xi, dtype=float)
    out = np.zeros_like(xi)
    for j in range(deriv, len(c)):
        fac = 1.0
        for k in range(deriv):
            fac *= j - k
        out += c[j] * fac * xi ** (j - deriv)
    return out


def _curvature(rp, rpp):
    return rpp / (1.0 + rp * rp) ** 1.5


def _dkappa_ds(rp, rpp, rppp):
    """dκ/ds = [r'''(1+r'²) - 3 r' r''²] / (1+r'²)³ (弧長微分, 解析式)。"""
    g = 1.0 + rp * rp
    return (rppp * g - 3.0 * rp * rpp * rpp) / g**3


@dataclass
class UpstreamThroatPoly:
    """U→T quintic Hermite: (r_U, 0, 0) → (r_t, 0, 1/R_t)。x ∈ [-L_U, 0]。"""

    r_U: float
    R_t: float
    L_U: float
    r_t: float = 1.0

    def __post_init__(self):
        _check_inputs("UpstreamThroatPoly",
                      positive={"R_t": self.R_t, "L_U": self.L_U, "r_t": self.r_t},
                      finite_only={"r_U": self.r_U})
        if self.r_U <= self.r_t:
            raise ValueError(f"UpstreamThroatPoly: Δr = r_U - r_t = "
                             f"{self.r_U - self.r_t:.6g} ≤ 0 (収縮が定義できない)")
        self._c = _hermite_quintic(
            -self.L_U, 0.0, self.r_U, 0.0, 0.0, self.r_t, 0.0, 1.0 / self.R_t)

    @property
    def mu(self) -> float:
        """単調収縮の無次元判定量 μ = L_U² / (R_t Δr)。単調 ⇔ μ ≤ 20。"""
        return self.L_U**2 / (self.R_t * (self.r_U - self.r_t))

    def _xi(self, x):
        return (np.asarray(x, dtype=float) + self.L_U) / self.L_U

    def r(self, x, deriv: int = 0):
        return _poly_eval(self._c, self._xi(x), deriv) / self.L_U**deriv

    def validate(self, n: int = 2001) -> list[str]:
        v = []
        if self.mu > 20.0:
            v.append(f"μ = {self.mu:.4g} > 20 (非単調収縮: L_U ≤ "
                     f"√(20 R_t Δr) = {np.sqrt(20*self.R_t*(self.r_U-self.r_t)):.4g})")
        # 数値サンプル検査 (表現差し替え時の保険。解析条件と独立に判定)
        x = np.linspace(-self.L_U, 0.0, n)
        rp = self.r(x, 1)
        scale = max(1.0, float(np.max(np.abs(rp))))
        if np.any(rp > 1e-12 * scale):
            v.append("数値検査: dr/dx > 0 の区間あり (非単調収縮)")
        if np.any(self.r(x) <= 0.0):
            v.append("数値検査: r ≤ 0 の区間あり")
        return v


@dataclass
class ThroatExpansionPoly:
    """T→D 4 次壁 (勾配 cubic Hermite)。x ∈ [0, L_D]。theta_D は radian。"""

    R_t: float
    L_D: float
    theta_D: float
    r_t: float = 1.0

    def __post_init__(self):
        _check_inputs("ThroatExpansionPoly",
                      positive={"R_t": self.R_t, "L_D": self.L_D, "r_t": self.r_t},
                      angles={"theta_D": self.theta_D})

    @property
    def m_D(self) -> float:
        return float(np.tan(self.theta_D))

    @property
    def r_D(self) -> float:
        """D 半径 (従属量): r_t + L_D²/(12 R_t) + L_D tanθ_D / 2。"""
        return self.r_t + self.L_D**2 / (12.0 * self.R_t) + 0.5 * self.L_D * self.m_D

    @property
    def L_D_max(self) -> float:
        """途中変曲を作らない上限 3 R_t tanθ_D。"""
        return 3.0 * self.R_t * self.m_D

    def r(self, x, deriv: int = 0):
        xi = np.asarray(x, dtype=float) / self.L_D
        L, R, m = self.L_D, self.R_t, self.m_D
        if deriv == 0:
            return (self.r_t
                    + L * L / R * (xi**2 / 2.0 - 2.0 * xi**3 / 3.0 + xi**4 / 4.0)
                    + L * m * (xi**3 - xi**4 / 2.0))
        if deriv == 1:  # q(ξ)
            return L / R * xi * (1.0 - xi) ** 2 + m * xi**2 * (3.0 - 2.0 * xi)
        if deriv == 2:  # (1-ξ)/R [1 + ξ(6Rm/L - 3)]
            return (1.0 - xi) / R * (1.0 + xi * (6.0 * R * m / L - 3.0))
        if deriv == 3:  # (1/L²) d²q/dξ² = (1/L²)[(L/R)(6ξ-4) + 6m(1-2ξ)]
            return ((L / R) * (6.0 * xi - 4.0) + 6.0 * m * (1.0 - 2.0 * xi)) / L**2
        raise ValueError("deriv は 0..3")

    def validate(self, n: int = 2001) -> list[str]:
        v = []
        if self.L_D > self.L_D_max:
            v.append(f"L_D = {self.L_D:.4g} > 3 R_t tanθ_D = {self.L_D_max:.4g} "
                     "(T–D 間に余分な変曲点)")
        x = np.linspace(0.0, self.L_D, n)
        rpp = self.r(x, 2)
        scale = max(1.0, float(np.max(np.abs(rpp))))
        if np.any(rpp < -1e-12 * scale):
            v.append("数値検査: r'' < 0 の区間あり (途中変曲)")
        return v


class WallDrivenThroatRegion:
    """U→T→D 合成壁: 解析 r/r'/r''/r''' → θ, κ, dκ/ds。x ∈ [-L_U, L_D]。"""

    def __init__(self, r_U: float, R_t: float, L_U: float, L_D: float,
                 theta_D: float, r_t: float = 1.0) -> None:
        self.up = UpstreamThroatPoly(r_U=r_U, R_t=R_t, L_U=L_U, r_t=r_t)
        self.dn = ThroatExpansionPoly(R_t=R_t, L_D=L_D, theta_D=theta_D, r_t=r_t)

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        return np.where(x < 0.0, self.up.r(np.minimum(x, 0.0), deriv),
                        self.dn.r(np.maximum(x, 0.0), deriv))

    def theta(self, x):
        return np.arctan(self.r(x, 1))

    def kappa(self, x):
        return _curvature(self.r(x, 1), self.r(x, 2))

    def dkappa_ds(self, x):
        return _dkappa_ds(self.r(x, 1), self.r(x, 2), self.r(x, 3))

    @property
    def kappa_prime_jump_at_throat(self) -> float:
        """T での dκ/ds の跳び (下流側 − 上流側)。κ 自体は 1/R_t に連続。"""
        up = _dkappa_ds(0.0, 1.0 / self.up.R_t, float(self.up.r(0.0, 3)))
        dn = _dkappa_ds(0.0, 1.0 / self.dn.R_t, float(self.dn.r(0.0, 3)))
        return float(dn - up)

    def validate(self, n: int = 2001, dkds_max: float | None = None) -> list[str]:
        """設計条件の検査。dkds_max を与えると max|dκ/ds| の上限も課す
        (数値上限の既定値は未定 — CFD 感度で決める。plan W3/W5)。"""
        v = (["U→T: " + m for m in self.up.validate(n)]
             + ["T→D: " + m for m in self.dn.validate(n)])
        if dkds_max is not None:
            m = float(np.max(np.abs(self.sample(n)["dkappa_ds"])))
            if m > dkds_max:
                v.append(f"max|dκ/ds| = {m:.4g} > 上限 {dkds_max:.4g}")
        return v

    def sample(self, n: int = 401) -> dict:
        """診断用サンプル (弧長比例でなく x 等間隔。source_region 付き)。"""
        x = np.linspace(-self.up.L_U, self.dn.L_D, n)
        return {
            "x": x,
            "r": self.r(x),
            "drdx": self.r(x, 1),
            "theta": self.theta(x),
            "kappa": self.kappa(x),
            "dkappa_ds": self.dkappa_ds(x),
            "source_region": np.where(x < 0.0, "UPSTREAM_HERMITE", "THROAT_TO_D"),
        }

    def diagnostics(self, n: int = 2001, dkds_max: float | None = None) -> dict:
        s = self.sample(n)
        return {
            "mu": self.up.mu,
            "L_D_over_max": self.dn.L_D / self.dn.L_D_max,
            "r_D": self.dn.r_D,
            "kappa_prime_jump_at_throat": self.kappa_prime_jump_at_throat,
            "max_abs_dkappa_ds": float(np.max(np.abs(s["dkappa_ds"]))),
            "violations": self.validate(n, dkds_max=dkds_max),
        }


@dataclass
class ContractionToConeQuintic:
    """geometry option (原案 §10): 直管 (r_in, 0, 0) → 一定傾斜点 (r_in-Δr, -tanθ_a, 0)。

    曲率符号一定 (r'' ≤ 0 全区間) ⇔ λ = L tanθ_a / Δr ∈ [5/3, 5/2]。
    λ = 2 で 5 次項が消え r(s) = r_in - 2Δr s³ + Δr s⁴ に退化する。
    ①風洞の主経路では未使用 (Korte 型入口 contour との遷移などに使う保存知見)。
    x ∈ [0, L]。
    """

    r_in: float
    dr: float
    theta_a: float
    L: float

    def __post_init__(self):
        _check_inputs("ContractionToConeQuintic",
                      positive={"r_in": self.r_in, "dr": self.dr, "L": self.L},
                      angles={"theta_a": self.theta_a})
        if self.r_in - self.dr <= 0.0:
            raise ValueError(f"ContractionToConeQuintic: 出口半径 r_in - dr = "
                             f"{self.r_in - self.dr:.6g} ≤ 0")
        self._c = _hermite_quintic(
            0.0, self.L, self.r_in, 0.0, 0.0,
            self.r_in - self.dr, -float(np.tan(self.theta_a)), 0.0)

    @property
    def lam(self) -> float:
        return self.L * float(np.tan(self.theta_a)) / self.dr

    def r(self, x, deriv: int = 0):
        xi = np.asarray(x, dtype=float) / self.L
        return _poly_eval(self._c, xi, deriv) / self.L**deriv

    def validate(self, n: int = 2001) -> list[str]:
        v = []
        if not (5.0 / 3.0 - 1e-12 <= self.lam <= 2.5 + 1e-12):
            v.append(f"λ = {self.lam:.4g} ∉ [5/3, 5/2] (収縮部で曲率符号が反転)")
        x = np.linspace(0.0, self.L, n)
        rpp = self.r(x, 2)
        scale = max(1.0, float(np.max(np.abs(rpp))))
        if np.any(rpp > 1e-12 * scale):
            v.append("数値検査: r'' > 0 の区間あり")
        return v
