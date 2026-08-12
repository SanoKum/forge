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
        の順 (L 字に沿って単調)。戻り値: 全計算点 (init 含む)。"""
        pts = list(init_front)
        front = list(init_front)
        while len(front) >= 2:
            new = []
            for i in range(len(front) - 1):
                B, A = front[i], front[i + 1]  # B=右(下流)側 C⁻ 担体, A=左側 C⁺ 担体
                try:
                    P = self._interior(A, B)
                except RuntimeError:
                    continue
                if not (np.isfinite(P.x) and np.isfinite(P.r)) or P.r < -1e-12:
                    continue
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


def inverse_design(throat, target, x_axis_end: float, n_axis: int = 260,
                   n_start: int = 41, gamma: float = 1.4, dx_wall: float = 0.02):
    """starting line (throat: SauerThroat) + 軸目標 target(x) から壁を逆設計する。

    target: 呼び出し可能 M(x) — starting line の軸点 x0 で場に C1 整合していること
    (モード F: MachBezier.from_constraints の start に Sauer 軸微分を渡す)。
    x_axis_end: 軸目標の下流端 (壁端の C⁻ 足より下流まで — 一様出口なら
    M_d 一定を伸ばすだけ)。
    戻り値 dict: wall (n,4 [x,r,θ,M]), pts, mdot_start, mdot_exit。
    """
    inv = InverseMOC(gamma=gamma, delta=1.0)
    x0, rr, MM, tt = throat.starting_line(M_start=1.05, n=n_start)
    g = gamma
    # L 字初期データ: 最下流軸点 → x0 直後の軸点 → starting line (軸側は重複させない)
    ax = np.linspace(x0, x_axis_end, n_axis)[1:][::-1]
    init = [_Pt(float(x), 0.0, 0.0, float(pm_nu(float(target(x)), g)), g) for x in ax]
    init += [_Pt(float(x0), float(rr[i]), float(tt[i]), float(pm_nu(float(MM[i]), g)), g)
             for i in range(n_start)]
    pts = inv.fill(init)
    wall = inv.wall_streamline(pts, x0, float(rr[-1]), dx=dx_wall)
    # 一様出口設計ではリップ (θ が最大値を経て ~0 に戻る点) で壁が完成する。
    # そこで切る (以降は一様菱形領域の水平流線で物理壁ではない)
    if len(wall) > 10:
        i_pk = int(np.argmax(wall[:, 2]))
        after = wall[i_pk:, 2] <= np.deg2rad(0.05)
        if np.any(after):
            wall = wall[: i_pk + int(np.argmax(after)) + 1]
    # 質量流量診断: 壁 (=流管であるべき曲線) までを積分して比較する
    mdot0 = inv.massflux_across(pts, x0, float(rr[-1]))
    if len(wall) > 10:
        xm = wall[-1, 0] - 0.3
        rm = float(np.interp(xm, wall[:, 0], wall[:, 1]))
        mdot1 = inv.massflux_across(pts, xm, rm)
    else:
        mdot1 = float("nan")
    return {"wall": wall, "pts": pts, "x0": float(x0),
            "mdot_start": mdot0, "mdot_exit": mdot1}
