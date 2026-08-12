"""モード F 壁 (①風洞): 上流解析区分 + 逆設計壁テーブルの複合壁。

区間: 入口直管 → Bell–Mehta 5 次収縮 (両端 C2) → 単一円弧 R (φu → スロート →
starting line 壁足 s_w0) → **逆設計壁 (点列テーブル, 3 次スプライン)** → リップ。
逆設計壁は starting line の壁足 (円弧上, 接線 = s_w0) から始まるため接続は
構成的に C1。API は NozzleWall 互換 (r / drdx / x_in / x_e / validate)。
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from .wall import _hermite_quintic, _poly_eval


class ModeFWall:
    def __init__(self, wall_pts, R: float = 2.0, r_inlet: float = 2.5,
                 L_pipe: float = 0.5, L_contract: float = 3.0,
                 phi_u: float = np.deg2rad(25.0)) -> None:
        """wall_pts: (n,>=2) [x, r] — 逆設計壁 (無次元 r*=1)。先頭点は円弧上。"""
        wall_pts = np.asarray(wall_pts, dtype=float)
        self.R, self.r_inlet = float(R), float(r_inlet)
        self.phi_u = float(phi_u)
        # 上流円弧始端と収縮
        self.x_au = -R * np.sin(phi_u)
        self.r_au = 1.0 + R * (1.0 - np.cos(phi_u))
        self.x_cs = self.x_au - L_contract
        self.x_in = self.x_cs - L_pipe
        curv_au = 1.0 / (R * np.cos(phi_u) ** 3)
        self._c_con = _hermite_quintic(self.x_cs, self.x_au,
                                       r_inlet, 0.0, 0.0,
                                       self.r_au, -np.tan(phi_u), curv_au)
        # 逆設計壁テーブル (円弧終端 = 先頭点)
        self.x_w0 = float(wall_pts[0, 0])
        self._spl = CubicSpline(wall_pts[:, 0], wall_pts[:, 1])
        self.x_e = float(wall_pts[-1, 0])

    def r(self, x):
        x = np.asarray(x, dtype=float)
        out = np.empty_like(x)
        m_pipe = x < self.x_cs
        m_con = (x >= self.x_cs) & (x < self.x_au)
        m_arc = (x >= self.x_au) & (x < self.x_w0)
        m_dsg = x >= self.x_w0
        out[m_pipe] = self.r_inlet
        out[m_con] = _poly_eval(self._c_con, self.x_cs, self.x_au, x[m_con])
        out[m_arc] = 1.0 + self.R - np.sqrt(np.maximum(self.R ** 2 - x[m_arc] ** 2, 0.0))
        out[m_dsg] = self._spl(np.minimum(x[m_dsg], self.x_e))
        return out

    def drdx(self, x):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        m_con = (x >= self.x_cs) & (x < self.x_au)
        m_arc = (x >= self.x_au) & (x < self.x_w0)
        m_dsg = x >= self.x_w0
        out[m_con] = _poly_eval(self._c_con, self.x_cs, self.x_au, x[m_con], order=1)
        out[m_arc] = x[m_arc] / np.sqrt(np.maximum(self.R ** 2 - x[m_arc] ** 2, 1e-30))
        out[m_dsg] = self._spl(np.minimum(x[m_dsg], self.x_e), 1)
        return out

    def validate(self) -> list:
        msgs = []
        xs = np.linspace(self.x_in, self.x_e, 3000)
        rv = self.r(xs)
        if np.any(rv <= 0.0):
            msgs.append("壁半径が非正")
        if abs(float(rv.min()) - 1.0) > 5e-3:
            msgs.append(f"最小半径 {rv.min():.4f} != 1")
        m_dsg = xs >= self.x_w0
        if np.any(np.diff(rv[m_dsg]) < -1e-9):
            msgs.append("設計壁区間で半径が非単調")
        # 円弧→設計壁の C1 (接線差)
        h = 1e-6
        dl = float(self.drdx(np.array([self.x_w0 - h]))[0])
        dr_ = float(self.drdx(np.array([self.x_w0 + h]))[0])
        if abs(dl - dr_) > 5e-3:
            msgs.append(f"円弧/設計壁の接線不連続 ({dl:.4f} vs {dr_:.4f})")
        return msgs
