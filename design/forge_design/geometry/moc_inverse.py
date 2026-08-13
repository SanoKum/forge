"""逆 MOC マーチ (モード F 逆設計): 軸目標 M(x) → 場充填 → 壁流線抽出。

親計画 §4.7 v1 の実体。数理 (2026-08-13 の設計議論で確定):

- 軸は条件 2 個 (M 指定 + θ=0) を持つ Cauchy データ線。初期データ =
  starting line (遷音速パッチ出力, 軸→壁) + 軸目標点列 (下流へ延長) の L 字。
- **三角充填**: データ列の隣接ペア (右=C⁻ 担体, 左=C⁺ 担体) から interior
  単位プロセスで 1 段上の点を作る。レベルごとに 1 点ずつ減り、決定領域 =
  「左端 C⁺ (starting line 壁点から) と右端 C⁻ (最下流軸点へ降りる) に挟まれた
  三角形」を過不足なく埋める。壁点の依存域は軸の両側に足を持つ (C⁺ 足=上流,
  C⁻ 足=下流) ため、**軸目標は壁端の C⁻ 足まで下流延長が必要** (一様出口では
  M_d 一定でタダ)。
- **壁抽出**: 場を線形補間し、starting line 壁端から流線 dr/dx = tanθ を積分。
  質量流量の整合 (∫ρu dA vs スロート) は診断として併記。

単位プロセス (interior/軸対称源項/予測子修正子) は moc_kernel.KernelMOC を継承。
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import LinearNDInterpolator

from .moc_kernel import KernelMOC, _Pt, pm_nu, pm_mach


def _mass_flux_density(M, g):
    """等エントロピー流の ρ|V|/(ρ0 a0) (よどみ量規格化)。"""
    M = np.asarray(M, dtype=float)
    t = 1.0 + 0.5 * (g - 1.0) * M * M
    return M * t ** (-0.5 * (g + 1.0) / (g - 1.0))


class InverseMOC(KernelMOC):
    """軸 Cauchy データからの三角充填と壁流線抽出。"""

    # -- 場充填 ---------------------------------------------------------------
    def fill(self, init_front) -> list:
        """init_front: _Pt 列。**最下流軸点 → 上流軸点 → starting line を壁へ**
        の順 (L 字に沿って単調)。戻り値: 全計算点 (init 含む)。

        単位プロセスの役割 (A=C⁺ 担体 / B=C⁻ 担体) は前線の向きに依存する
        (軸区間では上流側が C⁺ 担体、縦の starting line 区間では下の点が C⁺
        担体)。固定割当てだと縦線区間で幾何的に逆転し、交点が後方に出て棄却
        され**スロート直後の楔領域が無計算のまま残る** (2026-08-14 に実測で
        発見 — 壁流線の付け根の曲率ゴミ/スロート直後圧縮波の真因)。ここでは
        ペアごとに両割当てを試し、**前方交点** (両点より下流) を与える側を
        採る向き非依存の充填にする。"""
        pts = list(init_front)
        front = list(init_front)
        while len(front) >= 2:
            new = []
            for i in range(len(front) - 1):
                P = None
                for A, B in ((front[i + 1], front[i]), (front[i], front[i + 1])):
                    try:
                        cand = self._interior(A, B)
                    except RuntimeError:
                        continue
                    if (np.isfinite(cand.x) and np.isfinite(cand.r)
                            and cand.r > -1e-12
                            and cand.x >= min(A.x, B.x) - 1e-9):
                        P = cand
                        break
                if P is not None:
                    new.append(P)
            front = new
            pts.extend(new)
        return pts

    # -- 壁流線 ---------------------------------------------------------------
    def wall_streamline(self, pts, x_start: float, r_start: float,
                        dx: float = 0.02) -> np.ndarray:
        """(x_start, r_start) から dr/dx = tanθ を RK2 積分。場外に出たら終了。

        戻り値: (n,4) [x, r, theta, M]。"""
        xy = np.array([[p.x, p.r] for p in pts])
        th = np.array([p.th for p in pts])
        nu = np.array([p.nu for p in pts])
        itp_th = LinearNDInterpolator(xy, th)
        itp_nu = LinearNDInterpolator(xy, nu)
        out = []
        x, r = float(x_start), float(r_start)
        t0 = float(itp_th(x, r))
        n0 = float(itp_nu(x, r))
        while np.isfinite(t0):
            out.append((x, r, t0, pm_mach(max(n0, 1e-9), self.g)))
            xm, rm = x + 0.5 * dx, r + 0.5 * dx * np.tan(t0)
            tm = float(itp_th(xm, rm))
            if not np.isfinite(tm):
                break
            x, r = x + dx, r + dx * np.tan(tm)
            t0, n0 = float(itp_th(x, r)), float(itp_nu(x, r))
        return np.array(out)

    # -- 質量流量診断 -----------------------------------------------------------
    def massflux_across(self, pts, x_plane: float, r_top: float, n: int = 120) -> float:
        """x=x_plane の縦断面 (0..r_top) の ∫ρu·2πr dr (よどみ量規格化)。"""
        xy = np.array([[p.x, p.r] for p in pts])
        th = np.array([p.th for p in pts])
        nu = np.array([p.nu for p in pts])
        itp_th = LinearNDInterpolator(xy, th)
        itp_nu = LinearNDInterpolator(xy, nu)
        rr = np.linspace(1e-6, r_top, n)
        tt = itp_th(np.full(n, x_plane), rr)
        vv = itp_nu(np.full(n, x_plane), rr)
        ok = np.isfinite(tt) & np.isfinite(vv)
        M = np.array([pm_mach(max(v, 1e-9), self.g) for v in vv[ok]])
        f = _mass_flux_density(M, self.g) * np.cos(tt[ok])
        return float(np.trapezoid(f * 2.0 * np.pi * rr[ok], rr[ok]))


def _flux_along(line, g: float) -> np.ndarray:
    """折れ線 (軸→上) を横切る累積質量流束 (よどみ量規格化, 2πr 込み)。

    セグメント法線 n dℓ = (dr, −dx) (下流向き) で
    Δflux = ρV (cosθ·dr − sinθ·dx) · 2π r_mid。戻り値は各点までの累積 (len(line))。"""
    acc = [0.0]
    for P0, P1 in zip(line[:-1], line[1:]):
        Mm = 0.5 * (P0.M + P1.M)
        thm = 0.5 * (P0.th + P1.th)
        rm = 0.5 * (P0.r + P1.r)
        f = float(_mass_flux_density(Mm, g))
        acc.append(acc[-1] + f * (np.cos(thm) * (P1.r - P0.r)
                                  - np.sin(thm) * (P1.x - P0.x)) * 2.0 * np.pi * rm)
    return np.asarray(acc)


class _CPlusMarch(InverseMOC):
    """C⁺ 線マーチ (Foelsch/CONTUR 型): 各軸点から C⁺ 線を立ち上げ、
    **壁点 = 累積質量流束が ṁ* に達する点** (C⁻ 担体不要の閉包) で切る。

    レベル三角充填は前線が縮む構造のため「壁自身の依存領域」(starting line
    頂点の C⁻ より上の帯) を埋められず、壁流線の付け根曲率が壊れる
    (2026-08-14 実測: 円弧 0.51 に対し 0.03)。壁点を流束閉包で直接決める
    本方式が古典の正解で、付け根は自然に円弧の続きになる。"""

    @staticmethod
    def _mu2(p: _Pt) -> float:
        return np.arcsin(1.0 / max(p.M, 1.0 + 1e-12))

    def _foot_index(self, line, xp, rp, m):
        """点 (xp,rp) を通る傾き m の直線と折れ線 line の交点 (セグメント index)。"""
        for i in range(len(line) - 1):
            P0, P1 = line[i], line[i + 1]
            dx, dr = P1.x - P0.x, P1.r - P0.r
            den = dr - m * dx
            if abs(den) < 1e-14:
                continue
            t = (rp + m * (P0.x - xp) - P0.r) / den
            if -1e-9 <= t <= 1.0 + 1e-9:
                return i
        return 0

    def _seg_flux(self, A: _Pt, P: _Pt) -> float:
        Mm, thm, rm = 0.5 * (A.M + P.M), 0.5 * (A.th + P.th), 0.5 * (A.r + P.r)
        return (float(_mass_flux_density(Mm, self.g))
                * (np.cos(thm) * (P.r - A.r) - np.sin(thm) * (P.x - A.x))
                * 2.0 * np.pi * rm)

    def march(self, init_line, axis_x, target, mdot_star: float):
        """init_line: starting line (軸→壁足, 壁足 θ 補正済み)。axis_x: 軸点列 (昇順)。

        各軸点から C⁺ 線を「はしご」で立ち上げる: 新線の次点 = interior(直前点,
        前線の次の点)。前線の点のうち軸点の C⁻ 足より上のものだけが担体になる
        (足より下の C⁻ は A の上流の軸に落ちる)。累積流束が ṁ* に達した
        セグメントで壁点を補間し、前線トップまで使っても不足する場合は
        **C⁺ 方向へ凍結状態で延長して流束閉包** (壁間の薄い帯の古典的処理 —
        不足分は前線間隔 O(Δ) で小さい)。戻り: wall (n,4), 全点。"""
        g = self.g
        prev = list(init_line)
        wall = [(prev[-1].x, prev[-1].r, prev[-1].th, prev[-1].M)]
        pts = list(prev)
        for xk in axis_x:
            A0 = _Pt(float(xk), 0.0, 0.0, float(pm_nu(float(target(xk)), g)), g)
            i0 = self._foot_index(prev, A0.x, A0.r, np.tan(A0.th - self._mu2(A0)))
            line = [A0]
            flux = 0.0
            top = None
            for B in prev[i0 + 1:]:
                try:
                    P = self._interior(line[-1], B)
                except RuntimeError:
                    break  # 担体を飛ばすと はしごの順序が壊れる — 線を打ち切り延長閉包へ
                if (not (np.isfinite(P.x) and np.isfinite(P.r))
                        or P.r <= line[-1].r or P.x < line[-1].x - 0.5):
                    break
                df = self._seg_flux(line[-1], P)
                if flux + df >= mdot_star:
                    s = (mdot_star - flux) / max(df, 1e-30)
                    A = line[-1]
                    top = _Pt(A.x + s * (P.x - A.x), A.r + s * (P.r - A.r),
                              A.th + s * (P.th - A.th), A.nu + s * (P.nu - A.nu), g)
                    line.append(top)
                    break
                flux += df
                line.append(P)
            if top is None:
                # 前線トップまで使って ṁ* 未達 — C⁺ 方向へ凍結状態で延長して閉包。
                # 単位 dr あたり flux = c·r (c = f(M)(cosθ − sinθ/mp)·2π) の線形則から
                # r_top を陽に解く。
                A = line[-1]
                mp = np.tan(A.th + self._mu2(A))
                c = (float(_mass_flux_density(A.M, g))
                     * (np.cos(A.th) - np.sin(A.th) / mp) * 2.0 * np.pi)
                rem = mdot_star - flux
                r_top = float(np.sqrt(max(A.r ** 2 + 2.0 * rem / max(c, 1e-30), 0.0)))
                top = _Pt(A.x + (r_top - A.r) / mp, r_top, A.th, A.nu, g)
                line.append(top)
            wall.append((top.x, top.r, top.th, top.M))
            pts.extend(line[1:])
            prev = line
        return np.asarray(wall), pts


def inverse_design(throat, target, x_axis_end: float, n_axis: int = 260,
                   n_start: int = 41, gamma: float = 1.4, dx_wall: float = 0.02,
                   th_wall0: float | None = None):
    """starting line (throat: SauerThroat) + 軸目標 target(x) から壁を逆設計する。

    target: 呼び出し可能 M(x) — starting line の軸点 x0 で場に C1 整合していること
    (モード F: MachBezier.from_constraints の start に Sauer 軸微分を渡す)。
    x_axis_end: 軸目標の下流端 (壁端の C⁻ 足より下流まで — 一様出口なら
    M_d 一定を伸ばすだけ)。
    戻り値 dict: wall (n,4 [x,r,θ,M]), pts, mdot_start, mdot_exit。
    """
    # 注: _CPlusMarch (壁点=流束閉包の古典法) は WIP — 壁帯の被覆は正しいが
    # 壁閉包の凍結延長が粗く壁角が崩れる。現行は向き修正済み三角充填 + 壁流線
    # (既知の限界: 壁足の曲率が円弧 1/R でなく ~0 に寝る — 接合は ModeFWall の
    # 曲率ブレンドで C2 化する)。
    inv = InverseMOC(gamma=gamma, delta=1.0)
    x0, rr, MM, tt = throat.starting_line(M_start=1.05, n=n_start)
    g = gamma
    ax = np.linspace(x0, x_axis_end, n_axis)[1:][::-1]
    init = [_Pt(float(x), 0.0, 0.0, float(pm_nu(float(target(x)), g)), g) for x in ax]
    init += [_Pt(float(x0), float(rr[i]), float(tt[i]), float(pm_nu(float(MM[i]), g)), g)
             for i in range(n_start)]
    if th_wall0 is not None:
        # 壁足の θ を厳密壁接線で上書き (Sauer 線形化 vbar は過小 — kernel と同じ補正)
        init[-1].th = float(th_wall0)
    mdot_star = float(_flux_along(init[len(ax):], g)[-1])
    pts = inv.fill(init)
    wall = inv.wall_streamline(pts, x0, float(rr[-1]), dx=dx_wall)
    # 一様出口設計ではリップ (θ が最大値を経て ~0 に戻る点) で壁が完成する。
    # そこで切る (以降は一様菱形領域の水平流線で物理壁ではない)
    if len(wall) > 10:
        i_pk = int(np.argmax(wall[:, 2]))
        after = wall[i_pk:, 2] <= np.deg2rad(0.05)
        if np.any(after):
            wall = wall[: i_pk + int(np.argmax(after)) + 1]
    # 質量流量診断: 壁 (=流管であるべき曲線) までの断面積分を ṁ* と比較
    mdot1 = (inv.massflux_across(pts, wall[-1, 0] - 0.3,
                                 float(np.interp(wall[-1, 0] - 0.3, wall[:, 0], wall[:, 1])))
             if len(wall) > 10 else float("nan"))
    return {"wall": wall, "pts": pts, "x0": float(x0),
            "mdot_start": mdot_star, "mdot_exit": mdot1}
