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

from .moc_kernel import (KernelMOC, _Pt, interior_vec, pm_mach,
                         pm_mach_vec, pm_nu)


def field_interpolator(pts):
    r"""全計算点の Delaunay を **1 回だけ**作り、$(\theta, \nu)$ をまとめて返す補間器。

    充填をベクトル化した後 (2026-08-15) は**三角形分割が支配コスト**になった
    (実測 n_axis=1000 で 1 個 8.4 s、n_axis=2000 で ~35 s)。従来は θ 用と ν 用に
    別々の `LinearNDInterpolator` を作り、さらに壁流線・質量流束診断・終端特性線が
    それぞれ作り直していたため同じ分割を 6 回構築していた。値を (n,2) にまとめ、
    生成した補間器を呼び出し側で使い回すことで 1 回に減らす (数値は完全に同一 —
    同じ点集合の同じ線形補間)。"""
    xy = np.array([[p.x, p.r] for p in pts])
    val = np.array([[p.th, p.nu] for p in pts])
    return LinearNDInterpolator(xy, val)


def _mass_flux_density(M, g):
    """等エントロピー流の ρ|V|/(ρ0 a0) (よどみ量規格化)。ガスモデルなら ρV/(ρ*a*)
    (規格化定数が違うが、本モジュールでは**比**しか使わないので等価)。"""
    if hasattr(g, "mass_flux_density"):
        return np.asarray(g.mass_flux_density(M), dtype=float)
    M = np.asarray(M, dtype=float)
    t = 1.0 + 0.5 * (g - 1.0) * M * M
    return M * t ** (-0.5 * (g + 1.0) / (g - 1.0))


class InverseMOC(KernelMOC):
    """軸 Cauchy データからの三角充填と壁流線抽出。"""

    # -- 場充填 ---------------------------------------------------------------
    def fill(self, init_front) -> list:
        """init_front: _Pt 列。**最下流軸点 → 上流軸点 → starting line を壁へ**
        の順 (L 字に沿って単調)。戻り値: 全計算点 (init 含む)。

        **既知の限界と、棄却された「修正」の記録 (2026-08-14)**: 単位プロセスの
        役割 (A=C⁺ 担体 / B=C⁻ 担体) は本来前線の向きに依存し、縦の starting
        line 区間では下の点が C⁺ 担体であるべき。この固定割当てでは縦線区間で
        逆転するため交点が棄却され、**スロート直後の楔領域が無計算のまま残る**
        (壁流線はそこを凸包補間で跨ぐ)。そこで「両割当てを試し前方交点を採る」
        向き非依存の充填を試したが、**実設計は悪化した** — 楔が埋まる代わりに
        誤ったデータで埋まり、壁流線がそれを拾うため:

        | 設計 | 出口 ε_M | コア M (目標 4.0) |
        | --- | --- | --- |
        | 本実装 (楔は空・凸包補間) | 0.173% | 3.9995 |
        | 向き非依存 fill | 0.932% | 4.0353 |

        よって**本実装 (楔を空のまま残す) を維持**する。なお初期値線をスロート
        特性線にした (A8) 時点で縦線区間が無くなり、この楔自体が構造的に消えた。
        壁を古典的に閉じる経路は `cplus_flux_wall` (A10) を使う。"""
        arr = self.fill_arrays(init_front)
        g = self.g
        return [_Pt(float(a[0]), float(a[1]), float(a[2]), float(a[3]), g, float(a[4]))
                for a in arr]

    def fill_arrays(self, init_front) -> np.ndarray:
        r"""`fill` の**フロント一括ベクトル版**。戻り: (n,5) [x, r, th, nu, M]。

        レベルごとに隣接ペアを**まとめて**単位過程に通す (計算内容はスカラー版と
        同一で、ループの順序だけが変わる)。スカラー版は 1 点ずつ Python で回すため
        点数 $\sim n^2/2$ に比例した Python オーバヘッドが支配していた
        (実測: $n_{axis}$=2000 で 369 秒)。ベクトル版は**レベル数 $n$ 回**の
        numpy 呼び出しで済む (同 2.4 秒, 150 倍)。numba/C++ を持ち込まずに済むのは、
        1 レベル内のペアが互いに独立だから (レベル間の依存だけが逐次)。
        """
        x = np.array([p.x for p in init_front], dtype=float)
        r = np.array([p.r for p in init_front], dtype=float)
        th = np.array([p.th for p in init_front], dtype=float)
        nu = np.array([p.nu for p in init_front], dtype=float)
        M = np.array([p.M for p in init_front], dtype=float)
        out = [np.column_stack([x, r, th, nu, M])]
        while len(x) >= 2:
            # B=右(下流)側 C⁻ 担体 = front[:-1], A=左側 C⁺ 担体 = front[1:]
            xP, rP, thP, nuP, ok = interior_vec(
                x[1:], r[1:], th[1:], nu[1:], M[1:],
                x[:-1], r[:-1], th[:-1], nu[:-1], M[:-1],
                self.g, self.delta, self.n_corr)
            ok &= rP >= -1e-12
            if not np.any(ok):
                break
            x, r, th, nu = xP[ok], rP[ok], thP[ok], nuP[ok]
            M = pm_mach_vec(nu, self.g)
            out.append(np.column_stack([x, r, th, nu, M]))
        return np.vstack(out)

    def fill_levels(self, init_front) -> np.ndarray:
        r"""`fill_arrays` の**レベル構造を保ったまま**返す版。戻り: (n_lev, n_pt, 5)。

        `[k, i]` = レベル $k$ の位置 $i$ の点 $[x, r, \theta, \nu, M]$、
        欠損 (棄却された対・前線の縮み) は NaN。**添字を詰めない**のが要点で、
        これにより特性線が添字だけで読み出せる (`cplus_lines`)。

        $L_k[i]$ は $L_{k-1}[i]$ (C⁻ 担体) と $L_{k-1}[i+1]$ (C⁺ 担体) から作るので:

        - **C⁻ 線** = 添字固定の列 $L_k[i],\ k=0,1,\dots$
        - **C⁺ 線** = 反対角線 $L_k[m-k],\ k=0,1,\dots,m$ (起点 $L_0[m]$)

        つまり三角充填の網と特性線網は同じもので、走査順が違うだけ。
        """
        x = np.array([p.x for p in init_front], dtype=float)
        r = np.array([p.r for p in init_front], dtype=float)
        th = np.array([p.th for p in init_front], dtype=float)
        nu = np.array([p.nu for p in init_front], dtype=float)
        M = np.array([p.M for p in init_front], dtype=float)
        n = len(x)
        out = np.full((n, n, 5), np.nan)
        out[0] = np.column_stack([x, r, th, nu, M])
        live = np.ones(n, dtype=bool)          # レベル k で有効な添字
        for k in range(1, n):
            # 対 (i, i+1) がともに有効なときだけ新点を作る (添字は i を継承)
            pair = live[:-1] & live[1:]
            if not np.any(pair):
                break
            xP, rP, thP, nuP, ok = interior_vec(
                x[1:], r[1:], th[1:], nu[1:], M[1:],
                x[:-1], r[:-1], th[:-1], nu[:-1], M[:-1],
                self.g, self.delta, self.n_corr)
            ok &= pair & (rP >= -1e-12)
            if not np.any(ok):
                break
            MP = pm_mach_vec(np.where(ok, nuP, 0.0), self.g)
            x = np.where(ok, xP, np.nan)
            r = np.where(ok, rP, np.nan)
            th = np.where(ok, thP, np.nan)
            nu = np.where(ok, nuP, np.nan)
            M = np.where(ok, MP, np.nan)
            out[k, :len(x)] = np.column_stack([x, r, th, nu, M])
            live = ok
        return out

    # -- 壁流線 ---------------------------------------------------------------
    def wall_streamline(self, pts, x_start: float, r_start: float,
                        dx: float = 0.02, itp=None) -> np.ndarray:
        """(x_start, r_start) から dr/dx = tanθ を RK2 積分。場外に出たら終了。

        `itp`: `field_interpolator(pts)` を使い回す場合に渡す (省略時は自前で構築)。
        戻り値: (n,4) [x, r, theta, M]。"""
        itp = field_interpolator(pts) if itp is None else itp
        out = []
        x, r = float(x_start), float(r_start)
        t0, n0 = (float(v) for v in itp(x, r))
        while np.isfinite(t0):
            out.append((x, r, t0, pm_mach(max(n0, 1e-9), self.g)))
            xm, rm = x + 0.5 * dx, r + 0.5 * dx * np.tan(t0)
            tm = float(itp(xm, rm)[0])
            if not np.isfinite(tm):
                break
            x, r = x + dx, r + dx * np.tan(tm)
            t0, n0 = (float(v) for v in itp(x, r))
        return np.array(out)

    # -- 質量流量診断 -----------------------------------------------------------
    def massflux_across(self, pts, x_plane: float, r_top: float, n: int = 120,
                        itp=None) -> float:
        """x=x_plane の縦断面 (0..r_top) の ∫ρu·2πr dr (よどみ量規格化)。"""
        itp = field_interpolator(pts) if itp is None else itp
        rr = np.linspace(1e-6, r_top, n)
        val = itp(np.full(n, x_plane), rr)
        tt, vv = val[:, 0], val[:, 1]
        ok = np.isfinite(tt) & np.isfinite(vv)
        M = np.array([pm_mach(max(v, 1e-9), self.g) for v in vv[ok]])
        f = _mass_flux_density(M, self.g) * np.cos(tt[ok])
        return float(np.trapezoid(f * 2.0 * np.pi * rr[ok], rr[ok]))


def flux_closure_wall(itp, mdot_star: float, xs, r_top, gamma: float = 1.4,
                      n_r: int = 400) -> np.ndarray:
    r"""**壁 = 累積質量流束が $\dot m^*$ に達する半径** (各 $x$ 断面, ベクトル化)。

    計画: plans/archived/tooling-nozzle-moc-flux-closure-wall.md (A9)。

    $$\int_0^{r_w}\rho(\mathbf u\cdot\mathbf n)\,dA
      = \int_0^{r_w} \rho V\cos\theta\,2\pi r\,dr = \dot m^*$$

    **⚠ 既定では使わない — CFD で棄却済み (2026-08-16)**。定常流では流線面と質量流束
    一定面は同一の対象なので「保存則を構成的に満たす閉包の方が良い壁が出る」と期待して
    実装したが、**実測は逆だった**。同一場 ($n_{\rm axis}$=2000) で `wall_mode` だけを
    振った CFD A/B (case/41 run_0050 [流線] vs run_0052 [本関数]):

    | 指標 | 流線 | 流束閉包 |
    | --- | --- | --- |
    | ‖ΔM‖∞ [% M_d] | 0.353 | 0.507 |
    | overshoot | +0.068% | +0.217% |
    | 出口 ε_M rms | 0.020% | 0.057% |

    厳密解を持つ放射源流でも流線 6.1e-4 に対し本閉包は +1.43e-3 (2.3 倍悪い、しかも
    $x$ によらずほぼ一定の系統バイアス)。実装は正しい (一様流で $r_w\equiv1$ を 1.7e-7)。

    **理由**: MOC 場は特性線適合関係 (運動量+等エントロピー) から作られており、
    離散レベルでは質量保存を満たさない。断面積分で壁を決めてもその不整合は消えず、
    **壁半径の誤差から壁形状の誤差へ移るだけ**。しかも断面積分は各断面の半径方向
    プロファイル全体 (軸対称ソース項 $1/r$ が効く近軸域を含む) の精度に依存するのに対し、
    流線は $\theta$ を 1 本の曲線上で積分するだけなので誤差が小さい。

    **紛らわしい指標に注意**: 本閉包を使うと `mdot_ratio_moc` は構成的に 1 になる
    (循環指標で品質の証拠にならない) し、$r_F$ も $C_D$ 整合値に近づく。どちらも
    改善に見えるが CFD は悪化する。**判定は必ず CFD の軸 M / 出口一様性で行うこと**。

    **適用範囲** (使う場合): 断面がまるごと計算領域に入る $x$ (初期値線の軸着地より
    下流) のみ。上流側は呼び出し側で流線壁と合成する。

    `itp`: `field_interpolator` の補間器。`xs`: 断面位置。`r_top`: 各断面の探索上限
    (スカラーまたは `xs` と同形)。戻り: (n,) の $r_w$。
    """
    xs = np.asarray(xs, dtype=float)
    r_top = np.broadcast_to(np.asarray(r_top, dtype=float), xs.shape)
    s = np.linspace(0.0, 1.0, n_r)                       # 断面ごとの正規化半径
    rr = r_top[:, None] * s[None, :]                     # (n_x, n_r)
    val = itp(np.repeat(xs, n_r), rr.ravel())            # 補間器呼び出しは 1 回だけ
    th = val[:, 0].reshape(rr.shape)
    nu = val[:, 1].reshape(rr.shape)
    M = pm_mach_vec(np.where(np.isfinite(nu), nu, 0.0), gamma)
    f = _mass_flux_density(M, gamma) * np.cos(th) * 2.0 * np.pi * rr
    good = np.isfinite(th) & np.isfinite(nu)
    # 場の外 (凸包外) は NaN。軸から連続する有限区間だけを積分に使う
    good &= np.cumprod(good, axis=1).astype(bool)
    f = np.where(good, f, 0.0)
    dr = np.diff(rr, axis=1)
    cum = np.concatenate([np.zeros((len(xs), 1)),
                          np.cumsum(0.5 * (f[:, 1:] + f[:, :-1]) * dr, axis=1)], axis=1)
    cum = np.where(good, cum, np.nan)
    out = np.empty(len(xs))
    for i in range(len(xs)):
        c, r_i, ok = cum[i], rr[i], good[i]
        n_ok = int(ok.sum())
        if n_ok < 3:
            out[i] = np.nan
            continue
        c, r_i = c[:n_ok], r_i[:n_ok]
        if c[-1] >= mdot_star:
            out[i] = float(np.interp(mdot_star, c, r_i))
        else:
            # 場の縁まで使っても不足 — 縁の状態を凍結して外挿 (薄い帯の古典的処理)。
            # dṁ/dr = k·r (k = ρV cosθ·2π 一定) から r_w を陽に解く。
            k = f[i, n_ok - 1] / max(r_i[-1], 1e-30)
            out[i] = float(np.sqrt(max(r_i[-1] ** 2
                                       + 2.0 * (mdot_star - c[-1]) / max(k, 1e-30), 0.0)))
    return out


def cplus_lines(levels: np.ndarray, m: int) -> np.ndarray:
    r"""レベル配列から**起点 $L_0[m]$ の C⁺ 線**を取り出す (`fill_levels` の反対角線)。

    $L_k[i]$ は $A=L_{k-1}[i+1]$ を C⁺ 担体として作られるので、$L_0[m]$ から出る
    C⁺ 線は $L_k[m-k]$ ($k=0,\dots,m$)。**計算は一切せず添字を読むだけ**。
    戻り: (L,5)、NaN の手前まで。"""
    n = levels.shape[1]
    d = np.diagonal(levels[:, ::-1], offset=n - 1 - m, axis1=0, axis2=1).T
    ok = np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])
    cut = int(np.argmin(ok)) if not ok.all() else len(ok)
    return d[:cut]


def cplus_flux_wall(levels: np.ndarray, init_cum: np.ndarray, mdot_star: float,
                    gamma: float = 1.4) -> np.ndarray:
    r"""**壁点 = C⁺ 線上で累積質量流束が $\dot m^*$ に達する点** (古典的な逆設計閉包)。

    計画: plans/active/tooling-nozzle-moc-wall-unit-process.md。

    定常流では壁は流線であり、流線は質量流束一定面なので、「軸から数えて
    $\dot m^*$ 分の流れが通った高さ」がその位置の壁になる。各 C⁺ 線に沿って

    $$\Delta\dot m = \rho V(\cos\theta\,\Delta r - \sin\theta\,\Delta x)\cdot 2\pi\bar r$$

    を足し上げ、$\dot m^*$ を跨いだ線分を線形補間して切る。

    **`flux_closure_wall` (A9, 棄却済) との決定的な違いは積分する線と使う値**:
    あちらは鉛直断面 × Delaunay 補間場だったので半径方向プロファイル全体の
    補間誤差を拾った。本関数は **C⁺ 線 × 網の点そのもの**で、補間を一切通さない。
    壁が網を作る過程から出てくるので、Delaunay も流線 ODE 積分も不要になる
    (`wall_streamline` は補間場の中で $dr/dx=\tan\theta$ を積分するため、
    誤差が下流へ累積し、それを抑えるために網を極端に細かくする必要があった)。

    `init_cum[m]`: 起点 $L_0[m]$ までに**すでに軸から流れている**流束。初期値線上の
    点なら軸からその点までの累積、軸上の点なら 0。
    戻り: (n,4) [x, r, theta, M] を x 昇順で。"""
    n = levels.shape[1]
    out = []
    for m in range(n - 1, -1, -1):
        line = cplus_lines(levels, m)
        if len(line) < 2:
            continue
        x, r, th, nu, M = (line[:, i] for i in range(5))
        Mm = 0.5 * (M[1:] + M[:-1])
        thm = 0.5 * (th[1:] + th[:-1])
        rm = 0.5 * (r[1:] + r[:-1])
        d = (_mass_flux_density(Mm, gamma)
             * (np.cos(thm) * np.diff(r) - np.sin(thm) * np.diff(x))
             * 2.0 * np.pi * rm)
        cum = init_cum[m] + np.concatenate([[0.0], np.cumsum(d)])
        if cum[0] >= mdot_star:                 # 起点が既に壁 (初期値線の壁足)
            out.append((x[0], r[0], th[0], M[0]))
            continue
        j = int(np.argmax(cum >= mdot_star))
        if cum[j] < mdot_star:                  # 線が尽きた = 壁に届かない
            continue
        s = (mdot_star - cum[j - 1]) / max(cum[j] - cum[j - 1], 1e-30)
        out.append((x[j - 1] + s * (x[j] - x[j - 1]),
                    r[j - 1] + s * (r[j] - r[j - 1]),
                    th[j - 1] + s * (th[j] - th[j - 1]),
                    M[j - 1] + s * (M[j] - M[j - 1])))
    w = np.asarray(out, dtype=float)
    return w[np.argsort(w[:, 0])] if len(w) else w


def _smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


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



def terminal_exit(pts, wall, x_E: float, M_d: float, gamma: float = 1.4,
                  ds: float = 2e-3, n_max: int = 400000, itp=None) -> dict:
    r"""**終端特性線による物理出口 F の決定** (2026-08-15, ユーザ指摘)。

    軸上で $M=M_d$ に達する点 $E=(x_E,0)$ から出る C⁺ 特性線
    ($dr/dx=\tan(\theta+\mu)$) を MOC 場の中で追跡し、**壁流線との交点**を物理出口
    $F$ とする。この特性線は一様域の上流境界であり、$F$ より下流の断面は
    $(M_d,\theta=0)$ で埋まる — これが「一様出口」の定義そのもの。

    **なぜ壁角しきい値ではだめか**: 従来は壁流線の $\theta$ が 0.05° を切った点で
    切っていた。壁角は $F$ に向けて**漸近的に**ゼロへ近づくので、しきい値は $F$ の
    数 $r_t$ 上流で発火する (実測 M_d=4, R=2 で 2.5$r_t$ 手前)。半径差は
    $10^{-3}r_t$ と無視できるが、**その分だけ一様コアが出口面に届かない**
    (実測: コア半径が出口半径の 80%)。長さの定義として不正確。

    場の外へ出た区間は一様出口状態 $(\theta=0, M=M_d)$ で凍結して直線延長する
    (終端特性線の下流側は定義上一様域なので近似ではない)。戻り値に凍結ステップ数を
    含めるので、追跡の大半が凍結なら「場が短すぎる」と分かる。
    """
    itp = field_interpolator(pts) if itp is None else itp
    wx, wr = wall[:, 0], wall[:, 1]
    mu_d = float(np.arcsin(1.0 / max(M_d, 1.0 + 1e-12)))
    x, r, n_frozen = float(x_E), 0.0, 0
    path = [(x, r)]
    for _ in range(n_max):
        th, nu = (float(v) for v in itp(x, r))
        if np.isfinite(th) and np.isfinite(nu):
            M = pm_mach(max(nu, 1e-9), gamma)
            slope = np.tan(th + np.arcsin(1.0 / max(M, 1.0 + 1e-12)))
        else:
            slope = np.tan(mu_d)
            n_frozen += 1
        x_n, r_n = x + ds, r + ds * slope
        rw, rw_n = float(np.interp(x, wx, wr)), float(np.interp(x_n, wx, wr))
        if r_n >= rw_n:                       # 壁を越えた — 線形補間で交点
            t = (rw - r) / max((r_n - r) - (rw_n - rw), 1e-30)
            x_F = x + t * ds
            return {"x_F": float(x_F), "r_F": float(np.interp(x_F, wx, wr)),
                    "ok": True, "n_frozen": n_frozen, "path": np.asarray(path)}
        x, r = x_n, r_n
        path.append((x, r))
        if x > wx[-1]:
            break
    return {"x_F": float("nan"), "r_F": float("nan"), "ok": False,
            "n_frozen": n_frozen, "path": np.asarray(path)}


def _design_cplus(inv, init, n_ax, ax, mdot_star, g, x0, x_wall0,
                  exit_mode, x_E, M_d, start_line) -> dict:
    r"""`wall_mode='cplus'` の設計本体 — **補間構造を一切作らない**経路。

    壁は各 C⁺ 線上の流束閉包 (`cplus_flux_wall`)、物理出口 $F$ は $E=(x_E,0)$ を
    起点とする**網の C⁺ 線**と壁の交点。どちらも網の点だけで決まるので Delaunay も
    流線 ODE も不要 (設計コストの 82% を占めていた三角形分割が丸ごと消える)。

    診断は `mdot_ratio_moc` を返さない — 本閉包では構成的に 1 になる**循環指標**に
    なるため (A9 の教訓)。代わりに**流線整合残差** $\max|dr/dx-\tan\theta_{\rm net}|$
    を返す: 壁は流線でもあるはずなので、独立な 2 つの閉包の食い違いを測っている。
    """
    lev = inv.fill_levels(init)
    cum0 = np.zeros(len(init))
    cum0[n_ax:] = _flux_along(init[n_ax:], g)
    wall = cplus_flux_wall(lev, cum0, mdot_star, g)
    if len(wall) < 10:
        raise RuntimeError(f"cplus: 壁点が {len(wall)} 個しか取れない (場が不足)")
    wall_full = wall.copy()
    exit_info: dict = {"mode": exit_mode}
    if exit_mode == "characteristic":
        if x_E is None or M_d is None:
            raise ValueError("exit_mode='characteristic' には x_E と M_d が必要")
        m_E = int(np.argmin(np.abs(np.asarray(ax) - float(x_E))))
        term = cplus_lines(lev, m_E)
        if len(term) < 3:
            raise RuntimeError("cplus: 終端 C⁺ が短すぎる")
        # 終端 C⁺ と壁の交点 (どちらも網由来の折れ線)
        rw_on_term = np.interp(term[:, 0], wall[:, 0], wall[:, 1],
                               left=np.nan, right=np.nan)
        d = term[:, 1] - rw_on_term
        j = int(np.argmax(d >= 0.0)) if np.any(d >= 0.0) else -1
        if j <= 0:
            raise RuntimeError("終端特性線が壁に到達しない (軸 target の延長不足)")
        s = -d[j - 1] / max(d[j] - d[j - 1], 1e-30)
        x_F = float(term[j - 1, 0] + s * (term[j, 0] - term[j - 1, 0]))
        exit_info.update({"x_F": x_F,
                          "r_F": float(np.interp(x_F, wall[:, 0], wall[:, 1])),
                          "x_E_index": m_E, "n_frozen": 0})
        wall = wall[wall[:, 0] <= x_F]
    elif exit_mode == "lip":
        i_pk = int(np.argmax(wall[:, 2]))
        after = wall[i_pk:, 2] <= np.deg2rad(0.05)
        if np.any(after):
            wall = wall[: i_pk + int(np.argmax(after)) + 1]
        exit_info.update({"x_F": float(wall[-1, 0]), "r_F": float(wall[-1, 1])})
    else:
        raise ValueError("exit_mode は 'characteristic' か 'lip'")
    # 独立診断: 壁は流線でもあるはず (循環しない整合チェック)
    m = wall[:, 0] > wall[0, 0] + 0.5
    res = (np.gradient(wall[m, 1], wall[m, 0]) - np.tan(wall[m, 2])) if m.sum() > 5 \
        else np.array([np.nan])
    return {"wall": wall, "wall_full": wall_full, "pts": None, "x0": float(x0),
            "start_line": start_line, "x_wall0": float(x_wall0),
            "wall_mode": {"mode": "cplus", "n_wall": int(len(wall)),
                          "streamline_residual_max": float(np.max(np.abs(res))),
                          "streamline_residual_rms": float(np.sqrt(np.mean(res ** 2)))},
            "levels": lev, "exit": exit_info,
            "mdot_start": mdot_star, "mdot_exit": float("nan")}


def inverse_design(throat, target, x_axis_end: float, n_axis: int = 260,
                   n_start: int = 41, gamma: float = 1.4, dx_wall: float = 0.02,
                   th_wall0: float | None = None, M_start: float = 1.05,
                   exit_mode: str = "lip", x_E: float | None = None,
                   M_d: float | None = None, start_line: str = "vertical",
                   wall_mode: str = "streamline", blend_width: float = 1.0):
    """starting line (throat: SauerThroat) + 軸目標 target(x) から壁を逆設計する。

    target: 呼び出し可能 M(x) — starting line の軸点 x0 で場に C1 整合していること
    (モード F: MachBezier.from_constraints の start に Sauer 軸微分を渡す)。
    x_axis_end: 軸目標の下流端 (壁端の C⁻ 足より下流まで — 一様出口なら
    M_d 一定を伸ばすだけ)。
    戻り値 dict: wall (n,4 [x,r,θ,M]), pts, mdot_start, mdot_exit。
    """
    # wall_mode: 'cplus' = 壁点を C⁺ 線上の流束閉包で決める古典法 (A10, `_design_cplus`。
    # Delaunay も流線 ODE も通らない) / 'streamline' = 三角充填 + 補間場の流線積分 (旧既定。
    # 向き非依存 fill への「修正」は出口指標を悪化させ撤回 — plan §9.1)。
    # streamline の既知の限界: 壁足の曲率が円弧 1/R でなく ~0 に寝る。
    inv = InverseMOC(gamma=gamma, delta=1.0)
    g = gamma
    if start_line == "throat_char":
        # **スロート特性線を初期値線にする** (CONTUR 流, 2026-08-15 ユーザ指摘)。
        # 線自体が C⁻ なので担体割当てが区間内で反転せず、縦線構成で残っていた
        # スロート直後の未計算楔が構造的に消える。壁流線はスロート壁点から始まる。
        xs_c, rr, MM, tt = throat.throat_characteristic(n=n_start)
        x0 = float(xs_c[0])            # 軸着地 = 設計が軸に効き始める位置
        x_line = np.asarray(xs_c, dtype=float)
    elif start_line == "vertical":
        x0, rr, MM, tt = throat.starting_line(M_start=M_start, n=n_start)
        x_line = np.full(n_start, float(x0))
    else:
        raise ValueError("start_line は 'throat_char' か 'vertical'")
    ax = np.linspace(x0, x_axis_end, n_axis)[1:][::-1]
    init = [_Pt(float(x), 0.0, 0.0, float(pm_nu(float(target(x)), g)), g) for x in ax]
    init += [_Pt(float(x_line[i]), float(rr[i]), float(tt[i]),
                 float(pm_nu(float(MM[i]), g)), g) for i in range(n_start)]
    if th_wall0 is not None and start_line == "vertical":
        # 壁足の θ を厳密壁接線で上書き (Sauer 線形化 vbar は過小 — kernel と同じ補正)
        # throat_char では壁足 = 幾何スロートで θ=0 が厳密に成り立つので不要。
        init[-1].th = float(th_wall0)
    mdot_star = float(_flux_along(init[len(ax):], g)[-1])
    if wall_mode == "cplus":
        return _design_cplus(inv, init, len(ax), ax, mdot_star, g, x0,
                             float(x_line[-1]), exit_mode, x_E, M_d, start_line)
    pts = inv.fill(init)
    # Delaunay は充填ベクトル化後の支配コスト → 1 個だけ作って全診断で使い回す
    itp = field_interpolator(pts)
    wall = inv.wall_streamline(pts, float(x_line[-1]), float(rr[-1]), dx=dx_wall,
                               itp=itp)
    wall_diag: dict = {"mode": wall_mode}
    if wall_mode == "flux":
        # **壁を流束閉包で置き換える** (A9)。断面が計算領域に入る x0 より下流だけ
        # 置換し、[x0, x0+blend_width] で流線壁と smoothstep 合成する
        # (r(x_wall0)=r_t を厳密に保つため — 合成帯では両者の差はまだ小さい)。
        x_lo = float(x0) + 0.05
        m = wall[:, 0] >= x_lo
        if int(m.sum()) < 10:
            raise RuntimeError("wall_mode='flux': 置換区間が短すぎる")
        xs_f = wall[m, 0]
        r_f = flux_closure_wall(itp, mdot_star, xs_f, 1.08 * wall[m, 1], gamma=g)
        if not np.all(np.isfinite(r_f)):
            raise RuntimeError("wall_mode='flux': 流束閉包が解けない断面がある")
        r_s = wall[m, 1].copy()
        w = _smoothstep((xs_f - x_lo) / max(blend_width, 1e-12))
        wall = wall.copy()
        wall[m, 1] = (1.0 - w) * r_s + w * r_f
        v = itp(wall[:, 0], wall[:, 1])          # θ, M は新しい壁上で取り直す
        okv = np.isfinite(v[:, 0]) & np.isfinite(v[:, 1])
        wall[okv, 2] = v[okv, 0]
        wall[okv, 3] = pm_mach_vec(v[okv, 1], g)
        wall_diag.update({
            "x_lo": x_lo, "blend_width": float(blend_width),
            # 合成帯で両者がどれだけ違うか (小さいほど合成の恣意性が小さい)
            "d_blend_max": float(np.max(np.abs((r_f - r_s)[w < 1.0]))) if np.any(w < 1.0) else 0.0,
            # 置換で壁がどれだけ動いたか (= 流線積分が溜め込んでいた誤差)
            "dr_max": float(np.max(np.abs(r_f - r_s))),
            "dr_exit": float((r_f - r_s)[-1])})
    elif wall_mode != "streamline":
        raise ValueError("wall_mode は 'streamline' か 'flux'")
    wall_full = wall.copy()
    exit_info: dict = {"mode": exit_mode}
    if exit_mode == "characteristic":
        # **物理出口 = 終端 C⁺ (E 発) と壁流線の交点** (正しい定義, `terminal_exit`)
        if x_E is None or M_d is None:
            raise ValueError("exit_mode='characteristic' には x_E と M_d が必要")
        te = terminal_exit(pts, wall, float(x_E), float(M_d), gamma=gamma, itp=itp)
        if not te["ok"]:
            raise RuntimeError("終端特性線が壁に到達しない (軸 target の延長不足の疑い"
                               f" — x_axis_end={x_axis_end:.3g}, wall_end={wall[-1,0]:.3g})")
        i_cut = int(np.searchsorted(wall[:, 0], te["x_F"]))
        wall = wall[:max(i_cut, 10)]
        exit_info.update({k: te[k] for k in ("x_F", "r_F", "n_frozen")})
        exit_info["term_path"] = te["path"]
    elif exit_mode == "lip":
        # [旧] リップ (θ が最大値を経て 0.05° を切る点) で切る。壁角は F へ漸近的に
        # 落ちるため、この閾値は F の数 r_t 上流で発火する (`terminal_exit` docstring)。
        if len(wall) > 10:
            i_pk = int(np.argmax(wall[:, 2]))
            after = wall[i_pk:, 2] <= np.deg2rad(0.05)
            if np.any(after):
                wall = wall[: i_pk + int(np.argmax(after)) + 1]
        exit_info.update({"x_F": float(wall[-1, 0]), "r_F": float(wall[-1, 1])})
    else:
        raise ValueError("exit_mode は 'characteristic' か 'lip'")
    # 質量流量診断: 壁 (=流管であるべき曲線) までの断面積分を ṁ* と比較
    mdot1 = (inv.massflux_across(pts, wall[-1, 0] - 0.3,
                                 float(np.interp(wall[-1, 0] - 0.3, wall[:, 0], wall[:, 1])),
                                 itp=itp)
             if len(wall) > 10 else float("nan"))
    return {"wall": wall, "wall_full": wall_full, "pts": pts, "x0": float(x0),
            "start_line": start_line, "x_wall0": float(x_line[-1]),
            "wall_mode": wall_diag,
            "exit": exit_info, "mdot_start": mdot_star, "mdot_exit": mdot1}
