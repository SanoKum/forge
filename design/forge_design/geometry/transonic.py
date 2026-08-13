"""円弧スロートの遷音速解 (Sauer 一次・軸対称)。

無次元化: 長さ = r* (スロート半径)、速度摂動 ubar=(u-a*)/a*。Sauer 内部座標 xs
(軸上ソニック点 = 0) で
    ubar = A xs + (g+1)/4 A^2 r^2
    vbar = (g+1)/2 A^2 xs r + (g+1)^2/16 A^3 r^3
    A    = sqrt(2 / ((g+1) R)),  R = 下流円弧曲率半径 / r*
壁接線条件から A と「幾何スロート位置」が決まる: 壁最小断面は xs_t = -R(g+1)^2 A^3/16
(< 0)、すなわち**軸上ソニック点は幾何スロートより下流** (古典的事実)。
本クラスの公開 API は幾何座標 x (壁最小断面 = 0) で受け、内部で xs = x + xs_t に
変換する。導出は plans §4.7(b)、一次資料 Sauer NACA TM-1147 / Hall 1962。
Kliegel–Levine 高次 (1/(R+1) 級数) は未実装 — R < 1 は精度警告 (サーベイ A2.5(A))。
"""
from __future__ import annotations

import numpy as np


class SauerThroat:
    def __init__(self, R: float, gamma: float = 1.4,
                 theta_mode: str = "linearized") -> None:
        if R < 1.0:
            import warnings

            warnings.warn(
                f"R={R:.3g} < 1: Sauer 一次解の適用域外 (Kliegel-Levine 高次未実装)。"
                "アンカー精度は帰還で吸収されるが初期形状の質は落ちる"
            )
        if theta_mode not in ("linearized", "exact_ratio"):
            raise ValueError("theta_mode は 'linearized' か 'exact_ratio'")
        self.theta_mode = theta_mode
        self.R = float(R)
        self.g = float(gamma)
        self.A = np.sqrt(2.0 / ((self.g + 1.0) * self.R))
        # 幾何スロート (壁最小断面) の Sauer 内部座標
        self.xs_throat = -self.R * (self.g + 1.0) ** 2 * self.A ** 3 / 16.0
        # 軸上ソニック点の幾何座標 (>0: スロート下流)
        self.x_sonic_axis = -self.xs_throat

    # -- 摂動場 (引数は幾何座標 x: 壁最小断面 = 0) ----------------------------
    def _xs(self, x):
        return np.asarray(x, dtype=float) + self.xs_throat

    def ubar(self, x, r):
        g, A = self.g, self.A
        return A * self._xs(x) + (g + 1.0) / 4.0 * A * A * np.asarray(r) ** 2

    def vbar(self, x, r):
        g, A = self.g, self.A
        r = np.asarray(r)
        return (g + 1.0) / 2.0 * A * A * self._xs(x) * r + (g + 1.0) ** 2 / 16.0 * A ** 3 * r ** 3

    def mach(self, x, r):
        """|q| から等エントロピー厳密関係で M (小摂動は q のみ近似)。"""
        g = self.g
        qb = np.hypot(1.0 + self.ubar(x, r), self.vbar(x, r))  # q/a*
        ab2 = 1.0 - 0.5 * (g - 1.0) * (qb * qb - 1.0)          # (a/a*)^2
        return qb / np.sqrt(np.maximum(ab2, 1e-12))

    def theta(self, x, r):
        r"""流れ角。`theta_mode` で摂動の扱いを選ぶ (2026-08-14)。

        - `"linearized"` (既定): $\theta=\arctan\bar v$。Sauer 解は**線形化壁 BC**
          $\bar v = dr_w/dx$ を満たすよう構成されているので、この定義なら
          **壁流線が壁に沿う** (幾何整合)。
        - `"exact_ratio"`: $\theta=\arctan(\bar v/(1+\bar u))$。流れ角の厳密定義だが、
          Sauer の $\bar v$ が線形化 BC しか満たさないため壁で $1+\bar u\approx1.25$ の
          分だけ角度が**過小**になる (R=2 で実測 7.39°→5.89°, −25.7%)。これは
          starting line の時点で円弧と 1.5° のキンクを作り、圧縮波の源になる。

        両者は $O(\varepsilon)$ で一致し、差は 2 次。根治は Kliegel–Levine 高次
        (未実装) — 本オプションは 1 次解の枠内で幾何整合を優先する暫定。"""
        v = self.vbar(x, r)
        if self.theta_mode == "linearized":
            return np.arctan(v)
        return np.arctan2(v, 1.0 + self.ubar(x, r))

    # -- MOC 初期値線 -------------------------------------------------------
    def qbar_of_mach(self, M: float) -> float:
        g = self.g
        return float(np.sqrt(M * M * (2.0 + (g - 1.0)) / (2.0 + (g - 1.0) * M * M)))

    def x_axis_of_mach(self, M: float) -> float:
        """軸上 (r=0) で M に達する幾何座標 x。"""
        return (self.qbar_of_mach(M) - 1.0) / self.A + self.x_sonic_axis

    def starting_line(self, M_start: float = 1.05, n: int = 25):
        """x = x0 の縦線 (軸→壁) 上の (x0, r, M, theta)。軸上マッハ = M_start。"""
        x0 = self.x_axis_of_mach(M_start)
        r_w = 1.0 + x0 * x0 / (2.0 * self.R)  # 下流円弧の壁半径 (放物近似)
        r = np.linspace(0.0, r_w, n)
        return x0, r, self.mach(x0, r), self.theta(x0, r)

    # -- 流量係数の初期見込み --------------------------------------------------
    def cd_estimate(self, n: int = 200) -> float:
        """幾何スロート断面の質量流束数値積分 (Sauer 一次の範囲)。初期見込み用。"""
        g = self.g
        r = np.linspace(0.0, 1.0, n)
        ub = self.ubar(0.0, r)
        flux = 1.0 - 0.5 * (g + 1.0) * ub * ub  # ρu/(ρ*a*) の M=1 近傍展開
        return float(np.trapezoid(flux * 2.0 * r, r))
