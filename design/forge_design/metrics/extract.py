"""res_*.h5 からの目的関数抽出 (③ベル最小セット)。

規約 (サーベイ B4.5 / case/29 postprocess_thrust.py):
- 出口面の状態は境界ダンプでなく res h5 を BCONDS/<outlet>/iCells で引く。
- 軸対称の PLANES/surfArea は 2D 辺長 (dr) — dA = 2πr·dr を手で構成する。
- P は roe から再構成 (CPG)。
- 全指標は最小化形でなく生値も併記 (η は最大化量なので optimizer へは -η)。
"""
from __future__ import annotations

import h5py
import numpy as np

from ..evaluate.ic import invert_area_ratio


def thrust_metrics(mesh_h5, res_h5, gamma, Pt, p_ambient, r_throat,
                   outlet_physid=2) -> dict:
    """推力・C_F・η (= C_F / C_F,ideal 同面積比)。

    出口面の運動量法 F = ∫(ρu² + (P - Pa)) dA。壁摩擦・壁圧の影響は出口状態に
    既に織り込まれているため**別加算しない** (加算すると二重計上 — サーベイ
    B4.5 の「+壁摩擦」はこの抽出法では不採用。子 plan §9 参照)。
    """
    g = gamma
    with h5py.File(mesh_h5, "r") as nz:
        ip = nz[f"/BCONDS/{outlet_physid}/iPlanes"][:]
        ic = nz[f"/BCONDS/{outlet_physid}/iCells"][:]
        pc = nz["/PLANES/centCoords"][:].reshape(-1, 3)
        rf = pc[ip, 1]
        dr = nz["/PLANES/surfArea"][:][ip]  # 軸対称: 2D 辺長 (dr)
        dA = 2.0 * np.pi * rf * dr
    with h5py.File(res_h5, "r") as f:
        ro = f["/VALUE/ro"][:][ic]
        u = f["/VALUE/roUx"][:][ic] / ro
        v = f["/VALUE/roUy"][:][ic] / ro
        roe = f["/VALUE/roe"][:][ic]
        P = (g - 1.0) * (roe - 0.5 * ro * (u * u + v * v))
    F = float(np.sum((ro * u * u + (P - p_ambient)) * dA))
    At = np.pi * r_throat ** 2
    Ae = float(np.sum(dA))
    CF = F / (Pt * At)
    CF_id = cf_ideal(Ae / At, g, p_ambient / Pt)
    return {
        "thrust_N": F,
        "CF": CF,
        "CF_ideal": CF_id,
        "eta_cf": CF / CF_id,
        "eps_measured": Ae / At,
        "mdot_kg_s": float(np.sum(ro * u * dA)),
    }


def cf_ideal(eps, gamma, pa_over_pt) -> float:
    """1D 等エントロピーの理想推力係数 (同面積比・同背圧)。"""
    g = gamma
    Me = float(invert_area_ratio(np.array([eps]), np.array([True]), g)[0])
    pe_pt = (1.0 + 0.5 * (g - 1.0) * Me * Me) ** (-g / (g - 1.0))
    gterm = np.sqrt(
        (2.0 * g * g / (g - 1.0))
        * (2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))
        * (1.0 - pe_pt ** ((g - 1.0) / g))
    )
    return float(gterm + (pe_pt - pa_over_pt) * eps)


def core_radius_traced(mesh_h5, res_h5, x_d, x_plane, r_start_frac=0.02,
                       ds_frac=2e-4, gamma=1.4) -> float:
    r"""**有効菱形の上流境界**を実場でトレースして求める (テストコア半径)。

    軸が $M_d$ に達する点 $(x_d,0)$ から出る C⁺ ($dr/dx=\tan(\theta+\mu)$) を
    **実際の CFD 場**の中で前進積分し、`x_plane` での半径を返す。この線より下が
    一様域 (有効菱形の内側)。

    直線近似 $(x-x_d)\tan\mu_d$ は特性線が**一様域でしか直線でない**ことを無視する
    ため近似にすぎない (本ケースでは 0.824 vs 実トレース 0.830 で 0.7% 差)。
    軸上ちょうどは cell 中心が無く補間の凸包外なので、`r_start_frac`×r_max から
    開始する (node 化すれば 0 から開始できる)。
    """
    with h5py.File(mesh_h5, "r") as nz:
        cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
    with h5py.File(res_h5, "r") as f:
        Ux, Uy, son = f["/VALUE/Ux"][:], f["/VALUE/Uy"][:], f["/VALUE/sonic"][:]
    from scipy.interpolate import LinearNDInterpolator
    x, r = cc[:, 0], cc[:, 1]
    itp_M = LinearNDInterpolator(np.c_[x, r], np.hypot(Ux, Uy) / np.maximum(son, 1e-9))
    itp_th = LinearNDInterpolator(np.c_[x, r], np.arctan2(Uy, Ux))
    r_max = float(r.max())
    xx, rr = float(x_d), max(r_start_frac * r_max, float(r.min()) * 1.05)
    ds = ds_frac * float(x.max() - x.min())
    for _ in range(200000):
        if xx >= x_plane or rr >= r_max:
            break
        Mv, tv = float(itp_M(xx, rr)), float(itp_th(xx, rr))
        if not (np.isfinite(Mv) and np.isfinite(tv)) or Mv <= 1.0:
            break
        mu = np.arcsin(1.0 / max(Mv, 1.0 + 1e-12))
        rr += ds * np.tan(tv + mu)
        xx += ds
    return float(min(rr, r_max))


def core_radius_straight(x_plane: float, x_d: float, M_design: float,
                         r_wall: float) -> float:
    r"""有効菱形上流境界の**直線近似** $(x-x_d)\tan\mu_d$ (壁でクリップ)。

    一様域外では特性線が曲がるので近似。トレース版との差は本ケースで 0.7%。"""
    mu = np.arcsin(1.0 / max(float(M_design), 1.0 + 1e-12))
    return float(min(max((float(x_plane) - float(x_d)) * np.tan(mu), 0.0), float(r_wall)))


def exit_uniformity(mesh_h5, res_h5, M_design, x_d=None, outlet_physid=2,
                    core_mode: str = "traced", core_radius=None,
                    gamma=1.4, n_grid: int | None = None) -> dict:
    r"""出口一様性 $\varepsilon_M$ / $\varepsilon_\theta$ (サーベイ B4.5 定義)。

    **幾何は出口境界面の実 connectivity から取る** (`BCONDS/<id>/iPlanes` +
    `PLANES/surfArea`)。軸対称では面中心が半径中点なので
    $2\pi r_c\,dr = \pi(r_{out}^2-r_{in}^2)$ が**厳密**に成り立ち、環状面積は
    近似なしで得られる (セル中心から仮想境界を作る旧実装は近似だった)。
    壁半径も最外面の $r_{out}$ = **真の壁半径** (最外セル中心ではない)。

    - $\varepsilon_M$ = コア上の**質量流束重み RMS 偏差**
      $\sqrt{\langle (M-M_d)^2\rangle_w}/M_d$、$w_i=\rho_i u_{x,i} A_i$。
    - $\varepsilon_\theta$ = コア上の $\max|\theta|$ [deg] (独立の目的)。

    `core_mode`: `"traced"` (既定, 実場で C⁺ をトレース) / `"straight"` (直線近似)
    / `"full"` (コア制限なし) / `core_radius` 明示。
    `n_grid` を指定すると $r/r_{wall}$ 固定格子へ補間して積分し直す
    (サーベイの別法。厳密面積との一致を確認する交差検証用)。
    """
    with h5py.File(mesh_h5, "r") as nz:
        ip = nz[f"/BCONDS/{outlet_physid}/iPlanes"][:]
        ic = nz[f"/BCONDS/{outlet_physid}/iCells"][:]
        pc = nz["/PLANES/centCoords"][:].reshape(-1, 3)
        rf, xf = pc[ip, 1], pc[ip, 0]
        dr = nz["/PLANES/surfArea"][:][ip]
    with h5py.File(res_h5, "r") as f:
        Ux, Uy = f["/VALUE/Ux"][:][ic], f["/VALUE/Uy"][:][ic]
        son, ro = f["/VALUE/sonic"][:][ic], f["/VALUE/ro"][:][ic]
    o = np.argsort(rf)
    rf, dr, Ux, Uy, son, ro = rf[o], dr[o], Ux[o], Uy[o], son[o], ro[o]
    A = 2.0 * np.pi * rf * dr              # = π(r_out²−r_in²) (厳密)
    r_wall = float(rf[-1] + 0.5 * dr[-1])  # 真の壁半径 (最外面の外縁)
    x_plane = float(np.mean(xf))
    M = np.hypot(Ux, Uy) / np.maximum(son, 1e-9)
    th_deg = np.degrees(np.arctan2(Uy, Ux))
    Md = float(M_design)

    if core_radius is not None:
        rc, how = float(core_radius), "explicit"
    elif core_mode == "full" or x_d is None:
        rc, how = r_wall, "full(no_core_limit)"
    elif core_mode == "straight":
        rc, how = core_radius_straight(x_plane, x_d, Md, r_wall), "rhombus_straight"
    else:
        rc, how = core_radius_traced(mesh_h5, res_h5, x_d, x_plane, gamma=gamma), "rhombus_traced"
    core = rf <= rc
    if core.sum() < 4:
        raise ValueError(f"コア内の出口面が不足 ({int(core.sum())}, r_core={rc:.4g})")
    w = ro[core] * Ux[core] * A[core]
    eM = float(np.sqrt(np.average((M[core] - Md) ** 2, weights=w)) / Md)
    out = {"x_plane": x_plane, "r_core": rc, "core_def": how, "r_wall": r_wall,
           "eps_M_rms": eM,
           "eps_M_max": float(np.max(np.abs(M[core] - Md)) / Md),
           "eps_theta_max_deg": float(np.max(np.abs(th_deg[core]))),
           "M_core_massflux_avg": float(np.average(M[core], weights=w)),
           "n_core_faces": int(core.sum()), "n_exit_faces": int(len(rf))}
    if n_grid:   # 交差検証: 固定 r 格子へ補間して同じ量を積分し直す
        rg = np.linspace(0.0, rc, int(n_grid))
        Mg = np.interp(rg, rf, M)
        wg = np.interp(rg, rf, ro * Ux) * 2.0 * np.pi * rg
        out["eps_M_rms_grid"] = float(np.sqrt(np.average((Mg - Md) ** 2,
                                                         weights=np.maximum(wg, 0))) / Md)
    return out


def axis_mach(mesh_h5, res_h5, axis_band) -> tuple:
    r"""軸上/軸近傍の (x, M) を x 昇順で返す。**node では真の軸値**を返す。

    **cell/node で座標の意味が変わる罠** (2026-08-17 実測):
    node モードでは `/CELLS/centCoords` が**双対 CV 重心**に書き換わり VALUE と
    同じ長さになる。したがって cell 用の実装がエラーも出さずに「動いてしまう」が、
    軸上ノード (`MESH/COORD` で $r=0$) の双対 CV 重心は $r>0$ (実測 中央値
    0.022 $r_t$、最大 0.029) なので、`cy < axis_band` の帯選択には**非軸ノードが
    混入する** (実測: band=0.09 で 567 DOF 中 246 個が非軸、$M$ が最大 +0.001
    ずれる)。node は軸上に DOF を持つのが利点なので、それを使う。

    - node (`/MESH/COORD` の長さ = VALUE 長): $r=0$ の**軸ノードのみ**を返す。
      `axis_band` は無視される (真の軸値なので帯平均が不要)。
    - cell: 従来どおり `/CELLS/centCoords` の $r<$`axis_band` 帯。cell には軸上に
      DOF が無いため帯平均が避けられない (真の軸値が要る用途では
      `feedback.cfd_anchor.axis_curve_true` の偶関数フィットを使う)。
    """
    with h5py.File(mesh_h5, "r") as nz:
        cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
        node_c = (nz["/MESH/COORD"][:].reshape(-1, 3) if "/MESH/COORD" in nz else None)
    with h5py.File(res_h5, "r") as f:
        Ux = f["/VALUE/Ux"][:]
        Uy = f["/VALUE/Uy"][:]
        son = f["/VALUE/sonic"][:]
    is_node = node_c is not None and len(node_c) == len(Ux)
    if is_node:
        m = node_c[:, 1] < 1e-12 * max(float(np.max(node_c[:, 1])), 1.0)
        x = node_c[m, 0]
    else:
        m = cc[:, 1] < axis_band
        x = cc[m, 0]
    M = np.hypot(Ux[m], Uy[m]) / np.maximum(son[m], 1e-9)
    o = np.argsort(x)
    return x[o], M[o]
