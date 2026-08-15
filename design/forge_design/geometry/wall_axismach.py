r"""axis-Mach チェーンの CFD ドメイン壁と壁 QA (①風洞, 無次元 r* = 1)。

計画: plans/accepted/tooling-nozzle-axismach-chain.md §5.3 +
plans/active/tooling-nozzle-axismach-throat-characteristic.md (A8)。壁構成 (上流→下流):

  入口直管 (r_U) → U→T 5次 Hermite (`UpstreamThroatPoly` 流用, r''(T)=1/R)
  → 逆 MOC 壁流線 [x_T, x_F] (端条件クランプ 5 次 B-spline = `ModeFWall` と同じ流儀)

**スロート T (x=0) から下流は全て MOC の出力**で、円弧も放物線も挟まない
(ユーザ指定 2026-08-15)。初期値線を**スロート特性線**にしたことで壁流線が
スロート壁点そのものから始まるため、旧実装にあった [T, x0] の骨接放物線区間は
不要になった (2026-08-15)。T での接続は spline 左端クランプ (r'=0, r''=1/R) と
Hermite の端点条件 (同じ 0 と 1/R) の一致で C²。

旧・縦 starting line 構成では壁流線が x0>0 から始まるため、[T, x0] を Hall 模型が
仮定する骨接放物線 r = 1 + x²/(2R) で埋めていた (`throat_start=False` で当時の
挙動を再現できる)。

`wall_qa` は逆 MOC 壁テーブルの品質指標 (単調性・最大壁角・出口角・x_F/r_F の
理論比較・曲率) を返す (原方針 §8.3, §28)。
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline, make_interp_spline

from .wall_walldriven import UpstreamThroatPoly


def area_ratio_isentropic(M_d: float, gamma: float = 1.4) -> float:
    """A_e/A_t (1D 等エントロピー・完全気体)。r_F/r_t = sqrt(この値)。"""
    g = gamma
    return float((1.0 / M_d) * ((2.0 + (g - 1.0) * M_d * M_d) / (g + 1.0))
                 ** ((g + 1.0) / (2.0 * (g - 1.0))))


class AxisMachCFDWall:
    """直管 + U→T Hermite + 逆 MOC 壁 (クランプ 5 次 B-spline)。
    旧・縦 starting line 構成では U→T Hermite と設計壁の間に骨接放物線が入る。

    `mesh2d.generate_axisym_mesh` / `paste_isentropic_ic` 互換
    (`x_in` / `x_e` / `r(x, deriv)`)。
    """

    def __init__(self, wall_pts, R: float, r_U: float = 2.5, L_U: float = 3.5,
                 L_pipe: float = 0.5) -> None:
        """wall_pts: (n,>=2) [x, r] — 逆 MOC 壁テーブル。先頭点は
        スロート (x=0, r=1) か、旧構成では x0>0 (骨接放物線上)。"""
        wall_pts = np.asarray(wall_pts, dtype=float)
        if wall_pts.ndim != 2 or len(wall_pts) < 10:
            raise ValueError("wall_pts は (n>=10, >=2) のテーブル")
        self.R = float(R)
        self.up = UpstreamThroatPoly(r_U=float(r_U), R_t=self.R, L_U=float(L_U))
        self.L_pipe = float(L_pipe)
        self.x_in = -self.up.L_U - self.L_pipe
        self.x0 = float(wall_pts[0, 0])
        if self.x0 < -1e-9:
            raise ValueError(f"設計壁始点 x0 = {self.x0:.4g} < 0 (スロート上流)")
        self.throat_start = self.x0 < 1e-6
        if self.throat_start:
            if abs(float(wall_pts[0, 1]) - 1.0) > 5e-3:
                raise ValueError(f"スロート始点なのに r = {wall_pts[0, 1]:.5f} ≠ 1")
        else:
            # [旧] 縦 starting line 構成: 始点が骨接放物線に乗っているか
            r0_par = 1.0 + self.x0 ** 2 / (2.0 * self.R)
            if abs(float(wall_pts[0, 1]) - r0_par) > 5e-3:
                raise ValueError(f"設計壁始点 r = {wall_pts[0, 1]:.5f} が放物線 "
                                 f"{r0_par:.5f} から乖離")
        self.x_e = float(wall_pts[-1, 0])
        # 左端 = スロート (r'=0, r''=1/R) / 旧構成では放物線の解析微分にクランプ
        d0 = self.x0 / self.R
        s0 = 1.0 / self.R
        _cs = CubicSpline(wall_pts[:, 0], wall_pts[:, 1])
        self._spl = make_interp_spline(
            wall_pts[:, 0], wall_pts[:, 1], k=5,
            bc_type=([(1, d0), (2, s0)],
                     [(1, float(_cs(self.x_e, 1))), (2, float(_cs(self.x_e, 2)))]))

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        xU = -self.up.L_U
        out = np.empty_like(x)
        m_pipe = x < xU
        m_up = (x >= xU) & (x < 0.0)
        m_par = (x >= 0.0) & (x < self.x0)
        m_dsg = x >= self.x0
        if deriv == 0:
            out[m_pipe] = self.up.r_U
            out[m_par] = 1.0 + x[m_par] ** 2 / (2.0 * self.R)
        elif deriv == 1:
            out[m_pipe] = 0.0
            out[m_par] = x[m_par] / self.R
        elif deriv == 2:
            out[m_pipe] = 0.0
            out[m_par] = 1.0 / self.R
        else:
            raise ValueError("deriv は 0..2")
        if m_up.any():
            out[m_up] = self.up.r(x[m_up], deriv)
        if m_dsg.any():
            out[m_dsg] = self._spl(np.minimum(x[m_dsg], self.x_e), deriv) \
                if deriv else self._spl(np.minimum(x[m_dsg], self.x_e))
        return out

    def theta(self, x):
        return np.arctan(self.r(x, 1))

    def kappa(self, x):
        rp, rpp = self.r(x, 1), self.r(x, 2)
        return rpp / (1.0 + rp * rp) ** 1.5

    def validate(self, n: int = 4000) -> list:
        msgs = ["U→T: " + m for m in self.up.validate()]
        xs = np.linspace(self.x_in, self.x_e, n)
        rv = self.r(xs)
        if np.any(rv <= 0.0):
            msgs.append("壁半径が非正")
        if abs(float(rv.min()) - 1.0) > 5e-3:
            msgs.append(f"最小半径 {float(rv.min()):.4f} != 1")
        m_dsg = xs >= self.x0
        if np.any(np.diff(rv[m_dsg]) < -1e-9):
            msgs.append("設計壁区間で半径が非単調")
        # 接合の C1/C2 (構成的に成り立つはずだが実測で保証 — ModeFWall と同じ流儀)
        h = 1e-6
        joints = [("直管/U→T", -self.up.L_U)]
        joints += ([("U→T/設計壁 (スロート)", 0.0)] if self.throat_start
                   else [("U→T/放物線", 0.0), ("放物線/設計壁", self.x0)])
        for name, xc in joints:
            dl = float(self.r(np.array([xc - h]), 1)[0])
            dr_ = float(self.r(np.array([xc + h]), 1)[0])
            if abs(dl - dr_) > 5e-3:
                msgs.append(f"{name} の接線不連続 ({dl:.4f} vs {dr_:.4f})")
            cl = float(self.r(np.array([xc - h]), 2)[0])
            cr = float(self.r(np.array([xc + h]), 2)[0])
            if abs(cl - cr) > 0.05 * max(abs(cl), abs(cr), 1.0):
                msgs.append(f"{name} の曲率不連続 ({cl:.4f} vs {cr:.4f})")
        return msgs


def wall_qa(wall, M_d: float, x_E: float, gamma: float = 1.4) -> dict:
    r"""逆 MOC 壁テーブル (n,4)[x,r,θ,M] の品質指標 (原方針 §8.3, §28)。

    - 単調性 / 最大壁角 / 出口壁角 (リップで θ→0 に戻っているか)
    - $x_F, r_F$ と理論予測 ($r_F=\sqrt{A_e/A_t}$、$x_F-x_E\approx r_F\sqrt{M_d^2-1}$)
      の比較 — terminal Mach line 近似は一様域でのみ厳密なので比率は情報値
    - 壁 M の終端値 (設計 $M_d$ との差)
    violations が空なら合格。
    """
    w = np.asarray(wall, dtype=float)
    x, r, th = w[:, 0], w[:, 1], w[:, 2]
    r_F, x_F = float(r[-1]), float(x[-1])
    r_F_pred = float(np.sqrt(area_ratio_isentropic(M_d, gamma)))
    dxEF_pred = r_F_pred * float(np.sqrt(M_d * M_d - 1.0))
    v = []
    if np.any(np.diff(r) < -1e-9):
        v.append("壁半径が非単調")
    th_exit_deg = float(np.rad2deg(th[-1]))
    if abs(th_exit_deg) > 0.2:
        v.append(f"出口壁角 {th_exit_deg:.3f}° (|θ|>0.2° — F 未到達)")
    if abs(r_F / r_F_pred - 1.0) > 0.03:
        v.append(f"r_F = {r_F:.4f} が 1D 理論 {r_F_pred:.4f} から 3% 超乖離")
    out = {
        "x_F": x_F, "r_F": r_F, "r_F_pred_1d": r_F_pred,
        "r_F_err_rel": float(r_F / r_F_pred - 1.0),
        "xF_minus_xE": x_F - float(x_E),
        "xF_minus_xE_pred": dxEF_pred,
        "theta_max_deg": float(np.rad2deg(th.max())),
        "theta_exit_deg": th_exit_deg,
        "M_wall_exit": float(w[-1, 3]),
        "violations": v,
    }
    return out
