"""円弧スロートの遷音速解 (Sauer 一次 / Hall–Kliegel–Levine 高次・軸対称)。

無次元化: 長さ = r* (スロート半径)、速度摂動 ubar=(u-a*)/a*。Sauer 内部座標 xs
(軸上ソニック点 = 0) で
    ubar = A xs + (g+1)/4 A^2 r^2
    vbar = (g+1)/2 A^2 xs r + (g+1)^2/16 A^3 r^3
    A    = sqrt(2 / ((g+1) R)),  R = 下流円弧曲率半径 / r*
壁接線条件から A と「幾何スロート位置」が決まる: 壁最小断面は xs_t = -R(g+1)^2 A^3/16
(< 0)、すなわち**軸上ソニック点は幾何スロートより下流** (古典的事実)。
本クラスの公開 API は幾何座標 x (壁最小断面 = 0) で受け、内部で xs = x + xs_t に
変換する。導出は plans §4.7(b)、一次資料 Sauer NACA TM-1147 / Hall 1962。

`HallThroat` (2026-08-15, axis-Mach チェーン A2) は Hall (1962) の軸対称 3 次系列を
Kliegel–Levine の置換 R⁻¹ = S⁻¹+S⁻²+S⁻³ (S = R+1) 込みで実装したもの。式の出典は
CONTUR (Sivells AEDC-TR-78-63, papers/nozzle_design/ ローカル) Appendix A:
一般式 (A-1) [u]・(A-26) [v] の σ=1 (軸対称)、係数 (A-14)〜(A-25)・(A-34)〜(A-40)、
λ = sqrt((1+σ)/((γ+1)S)) (本文 Eq. 10 — 原文の "(γ-1)" は Sauer 退化・Eq. (3)/(4)
との整合から (γ+1) の誤植/判読と確定)。座標は x=0 が幾何スロート (壁最小断面、
そこで壁 v=0)・y0=r*=1 規格化、速度は音速点規格化 (u = q_x/a*)。抽出の検算
(2026-08-15): スロート壁 u(0,1) = 本文 Eq. (8) を係数恒等式で確認、v(0,1)=0 が
V42−V22+V02=0 / V63−V43+V23−V03=0 の厳密恒等式として成立、du/dx(0,1) = Eq. (9)
(GT+UP2−UP0 = −(64γ²+117γ−1026)/1152 を確認)。
"""
from __future__ import annotations

import numpy as np


class SauerThroat:
    def __init__(self, R: float, gamma: float = 1.4,
                 theta_mode: str = "exact_ratio") -> None:
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

        - `"exact_ratio"` (**既定**): $\theta=\operatorname{atan2}(\bar v, 1+\bar u)$ —
          流れ角の**物理的定義**。`mach()` が速度ベクトル $(1+\bar u,\bar v)$ から
          計算しているので、M と θ が**同一の速度ベクトルを表す**という整合が保たれる。
        - `"linearized"` (**実験オプション**): $\theta=\arctan\bar v$。Sauer 解が満たす
          線形化壁 BC $\bar v=dr_w/dx$ と揃うため壁足で円弧接線に近くなる (R=2 で
          7.64° vs 円弧 7.39°、`exact_ratio` は 5.89°) が、**M と θ が別の速度ベクトルを
          指す**ため非線形 MOC へ渡す状態としては不整合。一次資料照合と整合試験が
          済むまで既定にしない (2026-08-14 レビュー指摘)。

        両者は $O(\varepsilon)$ で一致し差は 2 次。なお**壁足 θ は `inverse_design`
        の `th_wall0` で円弧接線に上書きされる**ので、「従来は starting line が直接
        1.5° のキンクを渡していた」という説明は不正確 (撤回)。"""
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


class HallThroat:
    """Hall–Kliegel–Levine 軸対称遷音速解 (3 次系列, S = R+1)。

    API は `SauerThroat` 互換 (`mach` / `theta` / `x_axis_of_mach` /
    `starting_line` / `qbar_of_mach`) + 軸上アンカーの解析微分 `axis_anchor`。
    座標は幾何スロート原点 (x=0 = 壁最小断面)、y0 = r* = 1 規格化。
    出典・検算はモジュール docstring 参照。
    """

    def __init__(self, R: float, gamma: float = 1.4) -> None:
        if R < 0.5:
            import warnings

            warnings.warn(f"R={R:.3g} < 0.5: Hall–KL 系列でも精度低下域 "
                          "(Cuffel らの実験照合範囲外)")
        self.R = float(R)
        self.g = g = float(gamma)
        self.S = S = self.R + 1.0
        self.lam = np.sqrt(2.0 / ((g + 1.0) * S))
        # --- 係数 (CONTUR Appendix A, 軸対称 σ=1) ---
        self.GR = (15.0 - 10.0 * g) / 288.0
        self.GS = (2708.0 * g * g + 2079.0 * g + 2115.0) / 82944.0
        self.GT = (92.0 * g * g + 180.0 * g - 9.0) / 1152.0
        self.GV = (g + 1.0) / 8.0
        self.GK = (4.0 * g * g - 57.0 * g + 27.0) / 48.0
        self.U42 = (2.0 * g + 9.0) / 24.0
        self.U22 = (4.0 * g + 3.0) / 24.0
        self.U63 = (556.0 * g * g + 1737.0 * g + 3069.0) / 10368.0
        self.U43 = (388.0 * g * g + 777.0 * g + 153.0) / 2304.0
        self.U23 = (304.0 * g * g + 255.0 * g - 54.0) / 1728.0
        self.UP2 = (52.0 * g * g + 51.0 * g + 327.0) / 384.0
        self.UP0 = (52.0 * g * g + 75.0 * g - 9.0) / 192.0
        self.V42 = (g + 3.0) / 9.0
        self.V22 = (20.0 * g + 27.0) / 96.0
        self.V02 = (28.0 * g - 15.0) / 288.0
        self.V63 = (6836.0 * g * g + 23031.0 * g + 30627.0) / 82944.0
        self.V43 = (3380.0 * g * g + 7551.0 * g + 3771.0) / 13824.0
        self.V23 = (3424.0 * g * g + 4071.0 * g - 972.0) / 13824.0
        self.V03 = (7100.0 * g * g + 2151.0 * g + 2169.0) / 82944.0
        # 軸上 u(x) = c0 + c1 z + c2 z² + c3 z³, z = λx (解析微分・根探索に使う)
        self._c_ax = (
            1.0 - 1.0 / (4.0 * S) - self.GR / S ** 2 - self.GS / S ** 3,
            1.0 - 1.0 / (8.0 * S) + self.GT / S ** 2,
            0.5 * (1.0 - 2.0 * g / 3.0 - self.GV / S),
            self.GK / 3.0,
        )
        # 軸上ソニック点 (u=1) の幾何座標 (>0: スロート下流)
        self.x_sonic_axis = self._x_of_u_axis(1.0)

    # -- 速度場 (音速点規格化の全量: u = q_x/a*, v = q_r/a*) -------------------
    def u(self, x, y):
        g, S = self.g, self.S
        z = self.lam * np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        y2 = y * y
        c0, c1, c2, c3 = self._c_ax
        return (c0 + c1 * z + c2 * z * z + c3 * z ** 3
                + y2 / (2.0 * S) + (self.U42 * y2 * y2 - self.U22 * y2) / S ** 2
                + (self.U63 * y2 ** 3 - self.U43 * y2 * y2 + self.U23 * y2) / S ** 3
                + z * (y2 / S + (self.UP2 * y2 * y2 - self.UP0 * y2) / S ** 2)
                + 0.5 * z * z * y2 * (3.0 - 7.0 * g) / (4.0 * S))

    def v(self, x, y):
        g, S = self.g, self.S
        z = self.lam * np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        y2 = y * y
        return (y / (self.lam * S)) * (
            (y2 - 1.0) / (4.0 * S)
            + (self.V42 * y2 * y2 - self.V22 * y2 + self.V02) / S ** 2
            + (self.V63 * y2 ** 3 - self.V43 * y2 * y2 + self.V23 * y2 - self.V03) / S ** 3
            + z * (1.0 + ((2.0 * g + 9.0) * y2 - 2.0 * g - 1.5) / (6.0 * S)
                   + (6.0 * self.U63 * y2 * y2 - 4.0 * self.U43 * y2
                      + 2.0 * self.U23) / S ** 2)
            + 0.5 * z * z * (2.0 + (4.0 * self.UP2 * y2 - 2.0 * self.UP0) / S)
            + (z ** 3 / 3.0) * (3.0 - 13.0 * g) / 4.0)

    # -- SauerThroat 互換 API --------------------------------------------------
    def mach(self, x, r):
        g = self.g
        qb = np.hypot(self.u(x, r), self.v(x, r))
        ab2 = 1.0 - 0.5 * (g - 1.0) * (qb * qb - 1.0)
        return qb / np.sqrt(np.maximum(ab2, 1e-12))

    def theta(self, x, r):
        return np.arctan2(self.v(x, r), self.u(x, r))

    def qbar_of_mach(self, M: float) -> float:
        g = self.g
        return float(np.sqrt(M * M * (2.0 + (g - 1.0)) / (2.0 + (g - 1.0) * M * M)))

    def _x_of_u_axis(self, u_target: float) -> float:
        """軸上 u(x) = u_target の根 (Newton, スロート近傍の枝)。"""
        c0, c1, c2, c3 = self._c_ax
        z = (float(u_target) - c0) / c1
        for _ in range(60):
            f = c0 + c1 * z + c2 * z * z + c3 * z ** 3 - u_target
            df = c1 + 2.0 * c2 * z + 3.0 * c3 * z * z
            dz = f / df
            z -= dz
            if abs(dz) < 1e-15:
                break
        return z / self.lam

    def x_axis_of_mach(self, M: float) -> float:
        """軸上 (r=0) で M に達する幾何座標 x。"""
        return self._x_of_u_axis(self.qbar_of_mach(M))

    def starting_line(self, M_start: float = 1.05, n: int = 25):
        """x = x0 の縦線 (軸→壁) 上の (x0, r, M, theta)。壁は骨接放物線。"""
        x0 = self.x_axis_of_mach(M_start)
        r_w = 1.0 + x0 * x0 / (2.0 * self.R)
        r = np.linspace(0.0, r_w, n)
        return x0, r, self.mach(x0, r), self.theta(x0, r)

    # -- 軸上アンカー (M, M', M'') の解析微分 -----------------------------------
    def axis_anchor(self, x: float) -> tuple:
        """軸上 x での (M, dM/dx, d²M/dx²) を同一級数の解析微分で返す。

        M = f(u), f(u) = u·ab2^(-1/2), ab2 = 1 − (γ−1)(u²−1)/2。
        f' = ab2^(-1/2)·(1 + (γ−1)u²/(2·ab2)),
        f'' = 2g' + u·g'',  g = ab2^(-1/2) (導出は素直な連鎖律)。
        """
        g = self.g
        c0, c1, c2, c3 = self._c_ax
        z = self.lam * float(x)
        u = c0 + c1 * z + c2 * z * z + c3 * z ** 3
        up = self.lam * (c1 + 2.0 * c2 * z + 3.0 * c3 * z * z)
        upp = self.lam ** 2 * (2.0 * c2 + 6.0 * c3 * z)
        ab2 = 1.0 - 0.5 * (g - 1.0) * (u * u - 1.0)
        gg = ab2 ** -0.5
        ggp = 0.5 * (g - 1.0) * u * ab2 ** -1.5                    # dg/du
        ggpp = 0.5 * (g - 1.0) * ab2 ** -1.5 * (1.0 + 1.5 * (g - 1.0) * u * u / ab2)
        f = u * gg
        fp = gg + u * ggp
        fpp = 2.0 * ggp + u * ggpp
        return float(f), float(fp * up), float(fpp * up * up + fp * upp)

    # -- 流量係数 --------------------------------------------------------------
    def cd_series(self) -> float:
        """CONTUR Eq. (14) の解析級数 (軸対称)。"""
        g, S = self.g, self.S
        return float(1.0 - (g + 1.0) / (96.0 * S * S)
                     * (1.0 - (8.0 * g - 27.0) / (24.0 * S)
                        + (754.0 * g * g - 757.0 * g + 3615.0) / (2880.0 * S * S)))

    def cd_estimate(self, n: int = 400) -> float:
        """幾何スロート断面 (x=0) の質量流束数値積分 (級数場の厳密等エントロピー密度)。"""
        g = self.g
        r = np.linspace(0.0, 1.0, n)
        uu, vv = self.u(0.0, r), self.v(0.0, r)
        q2 = uu * uu + vv * vv
        ab2 = np.maximum(1.0 - 0.5 * (g - 1.0) * (q2 - 1.0), 1e-12)
        flux = ab2 ** (1.0 / (g - 1.0)) * uu       # ρ u_x / (ρ* a*)
        return float(np.trapezoid(flux * 2.0 * r, r))

    # -- スロート特性線 (CONTUR の throat characteristic) -----------------------
    def throat_characteristic(self, n: int = 61, ds: float = 2e-4,
                              max_steps: int = 400000) -> tuple:
        r"""**スロート壁点 $(0, r_t)$ から軸へ下る C⁻** を Hall 場の中で追跡する。

        CONTUR (AEDC-TR-78-63) §2 の "throat characteristic": 幾何スロートの壁は
        (曲率のあるスロートでは) 既に超音速なので ($R=2$, $\gamma=1.4$ で $M=1.129$)、
        そこから右進特性線 $dr/dx=\tan(\theta-\mu)$ を軸まで下ろせる。これが古典的な
        超音速 MOC の**初期値線**であり、その軸着地点が「壁の設計が軸に影響を及ぼせる
        最初の位置」になる。

        **なぜ縦線 starting line ではだめか** (2026-08-15 ユーザ指摘): $M_{start}=1.05$ の
        縦線は特性線ではないため、(a) 古典手法に対応物がなく (調査 §9-5)、(b) その壁足
        から下ろした C⁻ の軸着地は $x=1.81$ と、スロート特性線の $x=0.54$ より **1.3 $r_t$
        も下流**になる。つまり縦線構成では「スロート直後の壁が軸に効く区間」が設計対象
        から丸ごと外れていた。さらに (c) 縦線は特性線でないので単位過程の担体割当てが
        区間内で反転し、スロート直後に未計算の楔が残る (`moc_inverse.InverseMOC.fill`
        の docstring)。

        戻り値: `(x, r, M, theta)` の 4 配列、**軸→壁**の順 (`starting_line` と同じ規約)、
        各 n 点 ($r$ 等間隔にリサンプル)。
        """
        x, r = 0.0, 1.0
        path = [(x, r)]
        for _ in range(max_steps):
            M = float(self.mach(x, r))
            if M <= 1.0:
                raise RuntimeError(
                    f"throat_characteristic: スロート特性線が亜音速へ入った "
                    f"(M={M:.4f} @ x={x:.4f}, r={r:.4f})。R={self.R} が小さすぎる可能性")
            th = float(self.theta(x, r))
            slope = np.tan(th - np.arcsin(1.0 / M))
            r_n, x_n = r + ds * slope, x + ds
            if r_n <= 0.0:
                t = r / max(r - r_n, 1e-30)
                path.append((x + t * ds, 0.0))
                break
            x, r = x_n, r_n
            path.append((x, r))
        else:
            raise RuntimeError("throat_characteristic: 軸に到達しない")
        p = np.asarray(path)[::-1]              # 軸→壁
        r_q = np.linspace(0.0, 1.0, n)
        x_q = np.interp(r_q, p[:, 1], p[:, 0])
        return x_q, r_q, self.mach(x_q, r_q), self.theta(x_q, r_q)
