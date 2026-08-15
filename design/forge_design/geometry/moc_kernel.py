"""kernel MOC: 円弧スロート下流の特性曲線法マーチ (軸対称・非回転・等エントロピー)。

目的: スロート円弧 (R_d, theta_a) が支配する kernel 領域の軸上マッハ分布と、
kernel 境界 x_k (円弧終端 P_a からの最終 C- が軸に当たる点) での M, M', M''
(= 目標 Bézier のアンカー) を、設計変数のみの決定的関数として与える
(親計画 §4.7(b))。

特性方程式 (導出は plans §4.7 / 標準教科書 Zucrow-Hoffman):
    C±: dr/dx = tan(θ±μ)
    C+ に沿って d(θ-ν) = -δ (sinμ sinθ / r) dσ+
    C- に沿って d(θ+ν) = +δ (sinμ sinθ / r) dσ-
δ=1 (軸対称)。マーチは「連続する C- フロント」方式: 各フロント = 壁点から軸まで
の 1 本の C-。壁点は前フロントの壁直下点からの C+ で決める。予測子-修正子
(係数平均) の 2 次精度単位プロセス。
"""
from __future__ import annotations

import numpy as np


# --- Prandtl-Meyer -----------------------------------------------------------
def pm_nu(M, g=1.4):
    M = np.asarray(M, dtype=float)
    a = np.sqrt((g + 1.0) / (g - 1.0))
    b = np.sqrt(np.maximum(M * M - 1.0, 0.0))
    return a * np.arctan(b / a) - np.arctan(b)


def pm_mach(nu, g=1.4, tol=1e-12):
    """ν → M (Newton)。"""
    nu = float(nu)
    M = 1.0 + max(nu, 1e-6)  # 初期値
    for _ in range(60):
        f = pm_nu(M, g) - nu
        df = np.sqrt(M * M - 1.0) / (1.0 + 0.5 * (g - 1.0) * M * M) / M
        dM = f / max(df, 1e-14)
        M -= dM
        if M <= 1.0:
            M = 1.0 + 1e-9
        if abs(dM) < tol:
            break
    return M


def _mu(M):
    return np.arcsin(1.0 / max(M, 1.0 + 1e-12))


class _Pt:
    __slots__ = ("x", "r", "th", "nu", "M")

    def __init__(self, x, r, th, nu, g, M=None):
        self.x = x
        self.r = r
        self.th = th
        self.nu = nu
        # M は ν から一意に決まるが、既に計算済みなら渡して Newton を省ける
        # (ベクトル充填が 100 万点規模の _Pt を作るため — 2026-08-15)
        self.M = pm_mach(nu, g) if M is None else M


def _sin_over_r(p: "_Pt", other: "_Pt"):
    """軸対称項 sinθ/r。

    軸近傍では θ の丸め/打ち切り誤差 δθ が δθ/r に増幅されるため、
    (p, other) のうち**軸から遠い方**の点で比を評価する (l'Hôpital 的に
    sinθ/r → ∂θ/∂r で、滑らかな場では両者の比はほぼ共通)。
    """
    a, b = (p, other) if p.r >= other.r else (other, p)
    if a.r > 1e-9:
        return np.sin(a.th) / a.r
    return 0.0


def trace_cminus(field_M, field_th, x_w: float, r_w: float, n: int, g: float):
    """解析場 (M(x,r), θ(x,r)) の中で C- を壁点から軸まで逆トレースする。

    初期値線を鉛直線でなく C- 特性線にするための前処理。鉛直初期線だと最初の
    C- フロントとの楔領域を長い C+ セグメントで跨ぎ O(Δ²) バイアスが解像度に
    依らず残る (放射流照合で ~1% を確認)。r を独立変数に RK2。
    戻り値: [(x, r, th, nu)] 軸→壁の順。
    """
    pts = [(x_w, r_w)]
    x, r = x_w, r_w
    dr = -r_w / n
    for _ in range(n):
        M1 = field_M(x, r)
        th1 = field_th(x, r)
        s1 = 1.0 / np.tan(th1 - _mu(M1))
        xm = x + s1 * dr * 0.5
        rm = r + dr * 0.5
        M2 = field_M(xm, max(rm, 1e-12))
        th2 = field_th(xm, max(rm, 1e-12))
        x = x + dr / np.tan(th2 - _mu(M2))
        r = r + dr
        pts.append((x, max(r, 0.0)))
    out = []
    for (x, r) in pts[::-1]:
        th = 0.0 if r < 1e-12 else float(field_th(x, r))
        out.append((x, r, th, float(pm_nu(field_M(x, r), g))))
    return out


class KernelMOC:
    """円弧 (または一般壁曲線) 下の kernel を C- フロントでマーチする。"""

    def __init__(self, gamma=1.4, delta=1.0, n_corr=2):
        self.g = gamma
        self.delta = delta
        self.n_corr = n_corr

    # -- 単位プロセス --------------------------------------------------------
    def _interior(self, A: _Pt, B: _Pt) -> _Pt:
        """A (下, C+ 担体) と B (上, C- 担体) から新点。"""
        g, d = self.g, self.delta
        thP, nuP = 0.5 * (A.th + B.th), 0.5 * (A.nu + B.nu)
        xP = rP = None
        for _ in range(1 + self.n_corr):
            MP = pm_mach(nuP, g)
            muP = _mu(MP)
            mp = np.tan(0.5 * (A.th + thP) + 0.5 * (_mu(A.M) + muP))
            mm = np.tan(0.5 * (B.th + thP) - 0.5 * (_mu(B.M) + muP))
            if abs(mm - mp) < 1e-12:
                raise RuntimeError("特性線が平行 (M~1?)")
            xP = (A.r - B.r + mm * B.x - mp * A.x) / (mm - mp)
            rP = A.r + mp * (xP - A.x)
            Pn = _Pt(xP, max(rP, 0.0), thP, max(nuP, 0.0), g)
            fA = 0.5 * (np.sin(_mu(A.M)) * _sin_over_r(A, B) + np.sin(muP) * _sin_over_r(Pn, A))
            fB = 0.5 * (np.sin(_mu(B.M)) * _sin_over_r(B, A) + np.sin(muP) * _sin_over_r(Pn, B))
            Sp = -d * fA / np.cos(0.5 * (A.th + thP) + 0.5 * (_mu(A.M) + muP)) * (xP - A.x)
            Sm = +d * fB / np.cos(0.5 * (B.th + thP) - 0.5 * (_mu(B.M) + muP)) * (xP - B.x)
            Jp = (A.th - A.nu) + Sp   # θP - νP
            Jm = (B.th + B.nu) + Sm   # θP + νP
            thP, nuP = 0.5 * (Jp + Jm), 0.5 * (Jm - Jp)
        return _Pt(xP, rP, thP, nuP, self.g)

    def _axis(self, B: _Pt) -> _Pt:
        """B から C- で軸 (r=0, θ=0) へ。"""
        g, d = self.g, self.delta
        nuP = B.th + B.nu
        xP = B.x
        for _ in range(1 + self.n_corr):
            MP = pm_mach(nuP, g)
            muP = _mu(MP)
            mm = np.tan(0.5 * B.th - 0.5 * (_mu(B.M) + muP))
            xP = B.x - B.r / mm
            Pn = _Pt(xP, 0.0, 0.0, nuP, g)
            fB = 0.5 * (np.sin(_mu(B.M)) * _sin_over_r(B, B) + np.sin(muP) * _sin_over_r(Pn, B))
            Sm = +d * fB / np.cos(0.5 * B.th - 0.5 * (_mu(B.M) + muP)) * (xP - B.x)
            nuP = (B.th + B.nu) + Sm
        return _Pt(xP, 0.0, 0.0, nuP, self.g)

    def _wall_prescribed(self, s: float, wall_xy, wall_th, front) -> _Pt:
        """規定の壁ステーション s の壁点。C+ の足を前フロント折れ線上に補間で求める。

        壁ステップを規定側で制御する (前フロントの点間隔任せだと近壁解像度で
        ステップが決まり kernel 全域で数フロントしか取れず破綻する)。
        """
        g, d = self.g, self.delta
        xw, rw = wall_xy(s)
        thw = wall_th(s)
        nuW = max(front[-1].nu + (thw - front[-1].th), 1e-8)  # 2D 単純波の初期推定
        A, i_foot = front[-2], len(front) - 2
        for _ in range(2 + self.n_corr):
            muW = _mu(pm_mach(nuW, g))
            m = np.tan(0.5 * (A.th + thw) + 0.5 * (_mu(A.M) + muW))
            # W から後方 C+ 直線 r = rw + m (x - xw) と front 折れ線の交点
            A, i_foot = self._foot_on_front(front, xw, rw, m)
            Pn = _Pt(xw, rw, thw, nuW, g)
            fA = 0.5 * (np.sin(_mu(A.M)) * _sin_over_r(A, Pn) + np.sin(muW) * _sin_over_r(Pn, A))
            Sp = -d * fA / np.cos(0.5 * (A.th + thw) + 0.5 * (_mu(A.M) + muW)) * (xw - A.x)
            nuW = max(thw - (A.th - A.nu) - Sp, 1e-8)
        return _Pt(xw, rw, thw, nuW, g), i_foot

    def _foot_on_front(self, front, xw, rw, m):
        """直線 r = rw + m (x - xw) と front 折れ線 (軸→壁) の交点を線形補間。

        戻り値 (点, セグメント下端 index)。interior ペアリングは index 以下の
        front 点のみ使う (足より上の点の C+ は新 C- でなく壁に当たるため)。
        """
        for i in range(len(front) - 2, -1, -1):
            P0, P1 = front[i], front[i + 1]
            dx, drr = P1.x - P0.x, P1.r - P0.r
            den = drr - m * dx
            if abs(den) < 1e-14:
                continue
            t = (rw + m * (P0.x - xw) - P0.r) / den
            if -1e-9 <= t <= 1.0 + 1e-9:
                t = min(max(t, 0.0), 1.0)
                return _Pt(P0.x + t * dx, P0.r + t * drr,
                           P0.th + t * (P1.th - P0.th), P0.nu + t * (P1.nu - P0.nu), self.g), i
        # 交点なし → 壁直下点で代用 (ステップが極端に小さい場合)
        return front[-2], len(front) - 2

    def _wall(self, A: _Pt, wall_xy, wall_th, s0: float) -> tuple:
        """A から C+ で壁曲線 (param s) へ。戻り (点, s)。"""
        g, d = self.g, self.delta
        s = s0
        thP = wall_th(s)
        nuP = A.nu
        for _ in range(2 + self.n_corr):
            MP = pm_mach(max(nuP, 1e-8), g)
            muP = _mu(MP)
            mp = np.tan(0.5 * (A.th + thP) + 0.5 * (_mu(A.M) + muP))
            # C+ 直線と壁曲線の交点 (s を Newton)
            for _ in range(30):
                xw, rw = wall_xy(s)
                f = rw - (A.r + mp * (xw - A.x))
                h = 1e-7
                xw2, rw2 = wall_xy(s + h)
                df = (rw2 - (A.r + mp * (xw2 - A.x)) - f) / h
                if abs(df) < 1e-14:
                    break
                ds = f / df
                s -= ds
                if abs(ds) < 1e-13:
                    break
            xw, rw = wall_xy(s)
            thP = wall_th(s)
            Pn = _Pt(xw, rw, thP, max(nuP, 1e-8), g)
            fA = 0.5 * (np.sin(_mu(A.M)) * _sin_over_r(A, Pn) + np.sin(muP) * _sin_over_r(Pn, A))
            Sp = -d * fA / np.cos(0.5 * (A.th + thP) + 0.5 * (_mu(A.M) + muP)) * (xw - A.x)
            nuP = thP - (A.th - A.nu) - Sp
        return _Pt(xw, rw, thP, nuP, self.g), s

    # -- kernel マーチ (円弧壁) ----------------------------------------------
    def march_arc(self, throat, theta_a: float, n_start: int = 41, n_wall: int = 80):
        """throat: SauerThroat。円弧終端角 theta_a まで規定壁ステーションでマーチ。

        手順: (1) Sauer 有効域内の鉛直 starting line (軸上 M=1.05) を初期値線に、
        **三角形充填** (interior 列 + 軸反射) で初期壁点からの C- (楔境界) を
        細かい折れ線として構築 (Sauer 場を有効域外へ外挿しない)。
        (2) その C- をフロントに壁ステーションマーチ (θa まで)。
        戻り値 dict: axis_x, axis_M, x_k, anchor (M, M', M''), n_fronts。
        """
        g = self.g
        Rd = throat.R

        def wall_xy(s):
            return Rd * np.sin(s), 1.0 + Rd * (1.0 - np.cos(s))

        def wall_th(s):
            return s

        x0, rr, MM, tt = throat.starting_line(M_start=1.05, n=n_start)
        if x0 >= 0.8 * Rd * np.sin(theta_a):
            raise ValueError(
                f"遷音速パッチ (x0={x0:.3f}) が円弧終端に迫る: R={Rd} は Sauer 一次の"
                "適用域外 (Kliegel-Levine 高次が必要 — 既知の制約)")
        s_w0 = float(np.arcsin(min(x0 / Rd, 1.0)))
        front = [_Pt(x0, rr[i], tt[i], float(pm_nu(MM[i], g)), g) for i in range(n_start)]
        front[-1].th = wall_th(s_w0)
        axis_pts = [(front[0].x, front[0].M)]

        # (1) 楔充填: 初期値線の各点 k から出る C- を下から順にフロント化する
        #     (front_k = init[k] から軸までの C-)。最後のフロント = 初期壁点からの
        #     C- = 楔境界で、以後の壁ステーションマーチの初期フロントになる。
        wedge = [front[0]]
        for k in range(1, n_start):
            new = [front[k]]
            for j in range(len(wedge) - 1, -1, -1):
                try:
                    P = self._interior(wedge[j], new[-1])
                except RuntimeError:
                    continue
                if P.r < 0.0 or P.x < min(wedge[j].x, new[-1].x) - 1e-9:
                    continue
                new.append(P)
            ax = self._axis(new[-1])
            new.append(ax)
            wedge = new[::-1]
            axis_pts.append((ax.x, ax.M))
        front = wedge

                # (2) 壁ステーションマーチ
        stations = np.linspace(s_w0, theta_a, n_wall + 1)[1:]
        front = self._march_stations(front, stations, wall_xy, wall_th, axis_pts)
        axis_x = np.array([p[0] for p in axis_pts])
        axis_M = np.array([p[1] for p in axis_pts])
        x_k = float(axis_x[-1])
        n_fit = max(6, int(0.3 * len(axis_x)))
        c = np.polyfit(axis_x[-n_fit:] - x_k, axis_M[-n_fit:], 3)
        anchor = (float(np.polyval(c, 0.0)),
                  float(np.polyval(np.polyder(c), 0.0)),
                  float(np.polyval(np.polyder(c, 2), 0.0)))
        return {"axis_x": axis_x, "axis_M": axis_M, "x_k": x_k, "anchor": anchor,
                "n_fronts": len(axis_x)}

    def _march_stations(self, front, stations, wall_xy, wall_th, axis_pts):
        """規定壁ステーション列で C- フロントを順に構築する共通ループ。"""
        for s_k in stations:
            W, i_foot = self._wall_prescribed(float(s_k), wall_xy, wall_th, front)
            new = [W]
            for j in range(i_foot, -1, -1):
                try:
                    P = self._interior(front[j], new[-1])
                except RuntimeError:
                    continue
                if P.r < 0.0 or P.x < min(front[j].x, new[-1].x) - 1e-9:
                    continue
                new.append(P)
            ax = self._axis(new[-1])
            new.append(ax)
            front = new[::-1]
            axis_pts.append((ax.x, ax.M))
        return front

    # -- 検証用: 一般壁 (円錐等) --------------------------------------------
    def march_wall(self, front_pts, wall_xy, wall_th, stations):
        """既製フロントから一般壁下を規定ステーションでマーチ (テスト用)。"""
        g = self.g
        front = [_Pt(*p, g) for p in front_pts]  # (x, r, th, nu)
        axis_pts = [(front[0].x, front[0].M)]
        self._march_stations(front, np.asarray(stations, dtype=float), wall_xy, wall_th, axis_pts)
        return np.array(axis_pts)


# --- ベクトル化ユーティリティ (fill の O(n²) スカラーループ解消, 2026-08-15) -----
def pm_mach_vec(nu, g=1.4, tol=1e-13, iters=60):
    """ν → M の配列版 Newton (`pm_mach` と同じ反復・同じ初期値)。"""
    nu = np.asarray(nu, dtype=float)
    M = 1.0 + np.maximum(nu, 1e-6)
    for _ in range(iters):
        f = pm_nu(M, g) - nu
        df = np.sqrt(np.maximum(M * M - 1.0, 0.0)) / (1.0 + 0.5 * (g - 1.0) * M * M) / M
        dM = f / np.maximum(df, 1e-14)
        M = M - dM
        M = np.where(M <= 1.0, 1.0 + 1e-9, M)
        if np.all(np.abs(dM) < tol):
            break
    return M


def mu_vec(M):
    return np.arcsin(1.0 / np.maximum(M, 1.0 + 1e-12))


def _sin_over_r_vec(r_p, th_p, r_o, th_o):
    """`_sin_over_r` の配列版 (軸から遠い方の点で sinθ/r を評価)。"""
    far = r_p >= r_o
    ra = np.where(far, r_p, r_o)
    tha = np.where(far, th_p, th_o)
    return np.where(ra > 1e-9, np.sin(tha) / np.where(ra > 1e-9, ra, 1.0), 0.0)


def interior_vec(Ax, Ar, Ath, Anu, AM, Bx, Br, Bth, Bnu, BM,
                 g: float, delta: float, n_corr: int):
    """`KernelMOC._interior` の**配列版** (A=C⁺担体, B=C⁻担体 の対を一括処理)。

    スカラー版と同一の予測子-修正子・同一の係数平均・同一の軸対称源項評価。
    戻り: (xP, rP, thP, nuP, ok) — ok=False は特性線が平行/非有限で棄却すべき対。
    """
    thP = 0.5 * (Ath + Bth)
    nuP = 0.5 * (Anu + Bnu)
    muA, muB = mu_vec(AM), mu_vec(BM)
    xP = np.zeros_like(thP)
    rP = np.zeros_like(thP)
    ok = np.ones(thP.shape, dtype=bool)
    for _ in range(1 + n_corr):
        MP = pm_mach_vec(nuP, g)
        muP = mu_vec(MP)
        ang_p = 0.5 * (Ath + thP) + 0.5 * (muA + muP)
        ang_m = 0.5 * (Bth + thP) - 0.5 * (muB + muP)
        mp, mm = np.tan(ang_p), np.tan(ang_m)
        den = mm - mp
        good = np.abs(den) >= 1e-12
        den = np.where(good, den, 1.0)
        xP = (Ar - Br + mm * Bx - mp * Ax) / den
        rP = Ar + mp * (xP - Ax)
        rPc = np.maximum(rP, 0.0)
        fA = 0.5 * (np.sin(muA) * _sin_over_r_vec(Ar, Ath, Br, Bth)
                    + np.sin(muP) * _sin_over_r_vec(rPc, thP, Ar, Ath))
        fB = 0.5 * (np.sin(muB) * _sin_over_r_vec(Br, Bth, Ar, Ath)
                    + np.sin(muP) * _sin_over_r_vec(rPc, thP, Br, Bth))
        Sp = -delta * fA / np.cos(ang_p) * (xP - Ax)
        Sm = +delta * fB / np.cos(ang_m) * (xP - Bx)
        Jp = (Ath - Anu) + Sp
        Jm = (Bth + Bnu) + Sm
        thP, nuP = 0.5 * (Jp + Jm), 0.5 * (Jm - Jp)
        ok &= good
    ok &= np.isfinite(xP) & np.isfinite(rP) & np.isfinite(thP) & np.isfinite(nuP)
    return xP, rP, thP, nuP, ok
