"""モード F 壁 (①風洞): 上流解析区分 + 逆設計壁テーブルの複合壁。

区間: 入口直管 → Bell–Mehta 5 次収縮 (両端 C2) → 単一円弧 R (φu → スロート →
starting line 壁足 s_w0) → **逆設計壁 (点列テーブルを通す not-a-knot CubicSpline)**
→ リップ。円弧と設計壁の接続は **C1 が事後検査、C2 は未保証** (`blend_len>0` で
曲率ブレンドを挟めるが既定 off・撤回済み — 下記)。

**曲率ブレンドは既定 off・撤回済み (2026-08-14)**: 以下は経緯の記録。実測で
健全な設計壁を上書きして出口指標を悪化させたため既定 0.0。C2 接続は
(B1) 中心線 Cauchy データの C2 整合 と (B2) 拘束付き B-spline 表現 で保証する
方針に変更した (plan §9.2)。以下は旧説明 — 円弧終端で接線だけ合わせて (C1) 直接
スプラインに繋ぐと、曲率が円弧の $1/R$ から逆設計壁側の局所曲率へ**符号反転
込みで跳ぶ** (実測: +0.51 → −3.18, +719%)。Sivells/CONTUR の古典法は
「曲率連続の非粘性コンタ」を要件にしており (サーベイ A0)、この跳びはスロート
直後に不要な圧縮波を生む — 目標 Bézier 側の $M''$ を Sauer 解に整合させても
(逆設計は非線形 2D プロセスのため) 曲率跳びは解消しない/悪化する
(実測 +1116%) ため、**幾何側で直接ブレンド**する。手法は収縮→円弧接続と同じ
`_hermite_quintic` (両端で 値・接線・曲率の 3 条件, 計 6 条件 = 5 次多項式)。
ブレンド終端の曲率は逆設計スプラインを `blend_len` だけ内側 (局所ノイズを
避ける) で評価して採用する。
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from .wall import _hermite_quintic, _poly_eval


class ModeFWall:
    def __init__(self, wall_pts, R: float = 2.0, r_inlet: float = 2.5,
                 L_pipe: float = 0.5, L_contract: float = 3.0,
                 phi_u: float = np.deg2rad(25.0), blend_len: float = 0.0) -> None:
        """wall_pts: (n,>=2) [x, r] — 逆設計壁 (無次元 r*=1)。先頭点は円弧上。

        blend_len: 円弧→逆設計壁の曲率ブレンド長 (無次元 r* 単位)。"""
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
        # 円弧終端 (値・接線・曲率, 解析) → ブレンド終端 (スプラインを内側で評価)
        # blend_len <= 0 でブレンド無効 (C1 のみで直接スプライン接続 = 旧挙動)
        if float(blend_len) <= 0.0:
            self.x_bl = self.x_w0
            self._c_bl = None
            return
        self.x_bl = self.x_w0 + float(blend_len)
        r0 = 1.0 + R - np.sqrt(max(R ** 2 - self.x_w0 ** 2, 0.0))
        d0 = self.x_w0 / np.sqrt(max(R ** 2 - self.x_w0 ** 2, 1e-30))
        s0 = R ** 2 / max(R ** 2 - self.x_w0 ** 2, 1e-30) ** 1.5
        r1 = float(self._spl(self.x_bl))
        d1 = float(self._spl(self.x_bl, 1))
        s1 = float(self._spl(self.x_bl, 2))
        self._c_bl = _hermite_quintic(self.x_w0, self.x_bl, r0, d0, s0, r1, d1, s1)

    def r(self, x):
        x = np.asarray(x, dtype=float)
        out = np.empty_like(x)
        m_pipe = x < self.x_cs
        m_con = (x >= self.x_cs) & (x < self.x_au)
        m_arc = (x >= self.x_au) & (x < self.x_w0)
        m_bl = (x >= self.x_w0) & (x < self.x_bl)
        m_dsg = x >= self.x_bl
        out[m_pipe] = self.r_inlet
        out[m_con] = _poly_eval(self._c_con, self.x_cs, self.x_au, x[m_con])
        out[m_arc] = 1.0 + self.R - np.sqrt(np.maximum(self.R ** 2 - x[m_arc] ** 2, 0.0))
        if self._c_bl is not None and m_bl.any():
            out[m_bl] = _poly_eval(self._c_bl, self.x_w0, self.x_bl, x[m_bl])
        out[m_dsg] = self._spl(np.minimum(x[m_dsg], self.x_e))
        return out

    def drdx(self, x):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        m_con = (x >= self.x_cs) & (x < self.x_au)
        m_arc = (x >= self.x_au) & (x < self.x_w0)
        m_bl = (x >= self.x_w0) & (x < self.x_bl)
        m_dsg = x >= self.x_bl
        out[m_con] = _poly_eval(self._c_con, self.x_cs, self.x_au, x[m_con], order=1)
        out[m_arc] = x[m_arc] / np.sqrt(np.maximum(self.R ** 2 - x[m_arc] ** 2, 1e-30))
        if self._c_bl is not None and m_bl.any():
            out[m_bl] = _poly_eval(self._c_bl, self.x_w0, self.x_bl, x[m_bl], order=1)
        out[m_dsg] = self._spl(np.minimum(x[m_dsg], self.x_e), 1)
        return out

    def curvature(self, x):
        """r''(x) (診断用, r/drdx と独立の数値微分ではなく各区間の解析式)。"""
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        m_con = (x >= self.x_cs) & (x < self.x_au)
        m_arc = (x >= self.x_au) & (x < self.x_w0)
        m_bl = (x >= self.x_w0) & (x < self.x_bl)
        m_dsg = x >= self.x_bl
        out[m_con] = _poly_eval(self._c_con, self.x_cs, self.x_au, x[m_con], order=2)
        out[m_arc] = self.R ** 2 / np.maximum(self.R ** 2 - x[m_arc] ** 2, 1e-30) ** 1.5
        if self._c_bl is not None and m_bl.any():
            out[m_bl] = _poly_eval(self._c_bl, self.x_w0, self.x_bl, x[m_bl], order=2)
        out[m_dsg] = self._spl(np.minimum(x[m_dsg], self.x_e), 2)
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
        # 円弧/ブレンド/設計壁の接線・曲率連続 (構成的に一致するはずだが実測で保証)
        h = 1e-6
        joints = ([("円弧/設計壁", self.x_w0)] if self._c_bl is None
                  else [("円弧/ブレンド", self.x_w0), ("ブレンド/設計壁", self.x_bl)])
        for name, xc in joints:
            dl = float(self.drdx(np.array([xc - h]))[0])
            dr_ = float(self.drdx(np.array([xc + h]))[0])
            if abs(dl - dr_) > 5e-3:
                msgs.append(f"{name} の接線不連続 ({dl:.4f} vs {dr_:.4f})")
            if self._c_bl is None:
                continue   # ブレンド無効時は曲率跳びが既知の前提 (接線のみ検査)
            cl = float(self.curvature(np.array([xc - h]))[0])
            cr = float(self.curvature(np.array([xc + h]))[0])
            if abs(cl - cr) > 0.05 * max(abs(cl), abs(cr), 1.0):
                msgs.append(f"{name} の曲率不連続 ({cl:.4f} vs {cr:.4f})")
        return msgs
