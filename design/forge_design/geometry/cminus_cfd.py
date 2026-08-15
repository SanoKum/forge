r"""固定円弧終端 $J$ から出る **CFD C⁻ 特性線**の追跡と $x_{\mathrm{reach,CFD}}$ 抽出
(plan §9.2 B10 の正本)。

$x_{\mathrm{reach,CFD}}$ = 「自由壁が軸へ影響できる最初の位置」。可変壁の始点
$J=(x_j,r_j)$ (固定円弧の終端) から出る C⁻ が軸へ着地する $x$ で定義する。
設計 target と自由設計変数はここから始める。

**逆 MOC 場の `cminus_map.build_cminus_map()` は新正本では使わない** — 設計場は
仮定 (完全気体・等エントロピー・Sauer/CFD アンカー) の産物であり、「実際に壁が
軸へ届く位置」は実流で決めるべきだから。

`feedback/euler_loop.build_cminus_map_cfd()` は帰還マップ用の**診断実装**で、
最新 1 スナップ・壁の 0.985 倍から出発・軸近傍を直線外挿、と正本化前提の検証を
経ていない。本モジュールはそれらを設定可能にし、**感度を測って記録する**。

積分は $\dfrac{dr}{dx}=\tan(\theta-\mu)$、$\mu=\arcsin(1/M)$ (下流向き左走行波)。
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np


def _steps_of(f: Path) -> int:
    return int("".join(c for c in f.stem if c.isdigit()))


def load_field(run_dir, n_last: int = 3, scale: float = 0.01,
               discretization: str = "auto", with_pressure: bool = False) -> dict:
    """run から **末尾 n_last スナップ平均**の (x, r, M, θ) を返す (r* 単位)。

    末尾平均は中域リミットサイクル (片振幅 ~0.003) の位相ノイズ対策 — 軸 M 測定
    や CFD アンカー抽出と同じ規約。`discretization`: "cell" は CELLS/centCoords、
    "node" は MESH/COORD を座標に使う ("auto" は node を優先)。
    `with_pressure=True` で `/VALUE/P` も末尾平均して "P" キーに追加する
    (W4 の $P_0$ 診断用 — `moc_cancel.py`)。
    """
    rd = Path(run_dir)
    with h5py.File(rd / "nozzle.h5") as nz:
        if discretization == "cell":
            cc = nz["/CELLS/centCoords"][:].reshape(-1, 3)
        elif discretization == "node":
            cc = nz["/MESH/COORD"][:].reshape(-1, 3)
        else:
            cc = (nz["/MESH/COORD"][:].reshape(-1, 3) if "/MESH/COORD" in nz
                  else nz["/CELLS/centCoords"][:].reshape(-1, 3))
    res = sorted(rd.glob("res_[0-9]*.h5"), key=_steps_of)
    res = [f for f in res if _steps_of(f) > 0][-n_last:]
    if not res:
        raise FileNotFoundError(f"{rd} に res_*.h5 (step>0) が無い")
    Ms, ths, Ps = [], [], []
    for f in res:
        with h5py.File(f) as h:
            Ux, Uy, son = h["/VALUE/Ux"][:], h["/VALUE/Uy"][:], h["/VALUE/sonic"][:]
            if with_pressure:
                Ps.append(h["/VALUE/P"][:])
        if len(Ux) != len(cc):
            raise ValueError(f"座標 {len(cc)} と値 {len(Ux)} の長さ不一致 "
                             f"(discretization 指定を確認)")
        Ms.append(np.hypot(Ux, Uy) / np.maximum(son, 1e-9))
        ths.append(np.arctan2(Uy, Ux))
    out = {"x": cc[:, 0] / scale, "r": cc[:, 1] / scale,
           "M": np.mean(Ms, axis=0), "th": np.mean(ths, axis=0),
           "snaps": [f.name for f in res], "n_snap": len(res)}
    if with_pressure:
        out["P"] = np.mean(Ps, axis=0)
    return out


def _interpolators(fld: dict, x_max: float = 30.0):
    from scipy.interpolate import LinearNDInterpolator
    m = fld["x"] < x_max
    pts = np.c_[fld["x"][m], fld["r"][m]]
    return (LinearNDInterpolator(pts, fld["M"][m]),
            LinearNDInterpolator(pts, fld["th"][m]))


def trace_cminus(fld: dict, x_j: float, r_j: float, start_offset: float = 0.01,
                 ds: float = 2e-3, r_stop: float | None = None,
                 max_steps: int = 400000) -> dict:
    r"""$J$ の少し内側から C⁻ を軸方向へ前進積分する。

    `start_offset`: 出発点を壁から内側へずらす量 [r*]。壁ちょうどは補間の凸包
    境界で値が退化しうるため必要だが、**結果がこれに依存しては困る**ので
    呼び出し側 (`extract_x_reach`) が複数値で感度を取る。
    `r_stop`: ここまで積分して以降は最終勾配で $r=0$ へ外挿する (cell モードは
    軸上に DOF が無いため。node なら 0 近くまで積分できる)。既定は場の r 最小値。
    戻り: 軌跡 `path` (n,2)、`x_axis` (軸着地点)、`extrap`(外挿で稼いだ x 量)。
    """
    iM, iTh = _interpolators(fld)
    if r_stop is None:
        r_stop = float(np.min(fld["r"][fld["r"] > 0])) * 1.5 if np.any(fld["r"] > 0) else 0.0
    xx, rr = float(x_j), float(r_j) - float(start_offset)
    path = [(xx, rr)]
    slope = None
    for _ in range(max_steps):
        if rr <= r_stop:
            break
        Mv, tv = float(iM(xx, rr)), float(iTh(xx, rr))
        if not (np.isfinite(Mv) and np.isfinite(tv)) or Mv <= 1.0:
            return {"ok": False, "reason": f"M={Mv} th={tv} @({xx:.4f},{rr:.4f})",
                    "path": np.asarray(path), "x_axis": float("nan"), "extrap": float("nan")}
        mu = np.arcsin(1.0 / max(Mv, 1.0 + 1e-12))
        slope = np.tan(tv - mu)          # dr/dx (負: 下流へ進むと軸へ近づく)
        if abs(slope) < 1e-12:
            return {"ok": False, "reason": "slope~0", "path": np.asarray(path),
                    "x_axis": float("nan"), "extrap": float("nan")}
        rr += ds * slope
        xx += ds
        path.append((xx, rr))
    x_axis = xx - rr / slope if slope else float("nan")   # 最終勾配で r=0 へ外挿
    return {"ok": True, "path": np.asarray(path), "x_axis": float(x_axis),
            "extrap": float(x_axis - xx), "r_stop": float(r_stop)}


def wall_state_at(run_dir, x_j: float, scale: float = 0.01, gamma: float = 1.4,
                  n_last: int = 3) -> dict:
    r"""$x_j$ における**壁面境界出力**から $(M,\theta)$ を取る (末尾平均)。

    `res_wall_3_*.h5` (BCONDS/3 の面出力) を使う。slip 壁なので流れは壁接線
    方向 — $\theta$ は幾何 (壁接線) と一致するはずで、両者の差は診断値になる。
    場の内部補間と違い**凸包端の退化やセル中心オフセットの影響を受けない**。
    """
    rd = Path(run_dir)
    fs = sorted(rd.glob("res_wall_3_[0-9]*.h5"),
                key=lambda f: int(f.stem.split("_")[-1]))
    fs = [f for f in fs if int(f.stem.split("_")[-1]) > 0][-n_last:]
    if not fs:
        raise FileNotFoundError(f"{rd} に res_wall_3_*.h5 が無い")
    Ms, ths, xs = [], [], None
    for f in fs:
        with h5py.File(f) as h:
            cc = h["MESH/COORD"][:].reshape(-1, 3)
            Ux, Uy = h["VALUE/Ux"][:], h["VALUE/Uy"][:]
            ro, Ps = h["VALUE/ro"][:], h["VALUE/Ps"][:]
        xf = 0.5 * (cc[:-1, 0] + cc[1:, 0]) / scale if len(cc) == len(Ux) + 1 else cc[:len(Ux), 0] / scale
        a = np.sqrt(gamma * np.maximum(Ps, 1e-9) / np.maximum(ro, 1e-12))
        Ms.append(np.hypot(Ux, Uy) / a)
        ths.append(np.arctan2(Uy, Ux))
        xs = xf
    o = np.argsort(xs)
    xs = xs[o]
    M = np.mean(Ms, axis=0)[o]
    th = np.mean(ths, axis=0)[o]
    return {"M": float(np.interp(x_j, xs, M)), "th": float(np.interp(x_j, xs, th)),
            "n_snap": len(fs)}


def trace_cminus_from_wall(fld: dict, run_dir, x_j: float, r_j: float,
                           scale: float = 0.01, gamma: float = 1.4,
                           th_wall_geom: float | None = None,
                           first_step: float = 0.01, ds: float = 2e-3,
                           n_last: int = 3, **kw) -> dict:
    r"""**壁面状態を起点**に C⁻ を追跡する (B10 正本経路)。

    場の内部補間だけで壁近傍から出発すると、出発オフセットに強く依存する
    (実測 0.005→0.04 r* で $x_{reach}$ が 1.697→1.597 と 0.10 r* 動く) —
    内側から出発すると **$J$ より上流の壁点由来の別の C⁻ に乗る**ため。
    そこで最初の一歩だけは壁面状態 ($M_w$ は境界出力、$\theta_w$ は壁接線 —
    slip 壁なので流れ角と一致) が与える特性線勾配で内側へ入り、以降は場の
    補間で積分する。こうすると `first_step` 依存は 1 次で消える。
    """
    ws = wall_state_at(run_dir, x_j, scale=scale, gamma=gamma, n_last=n_last)
    th_w = float(th_wall_geom) if th_wall_geom is not None else ws["th"]
    mu_w = np.arcsin(1.0 / max(ws["M"], 1.0 + 1e-12))
    slope_w = np.tan(th_w - mu_w)
    x1 = x_j + first_step
    r1 = r_j + first_step * slope_w
    sub = trace_cminus(fld, x1, r1, start_offset=0.0, ds=ds, **kw)
    sub["wall_state"] = {"M_wall": ws["M"], "th_wall_cfd_deg": float(np.rad2deg(ws["th"])),
                         "th_wall_geom_deg": (float(np.rad2deg(th_wall_geom))
                                              if th_wall_geom is not None else None),
                         "first_step": first_step}
    if len(sub.get("path", [])):
        sub["path"] = np.vstack([[x_j, r_j], sub["path"]])
    return sub


def extract_x_reach(run_dir, x_j: float, r_j: float, scale: float = 0.01,
                    n_last: int = 3, discretization: str = "auto") -> dict:
    r"""$x_{\mathrm{reach,CFD}}$ を抽出し、**必須の感度を全部測って返す** (B10)。

    感度項目 (plan §9.2 B10): 出発オフセット / 積分刻み / 軸近傍外挿量 /
    スナップ平均数。呼び出し側は返り値の `spread` を見て「基準値として使える
    ばらつきか」を判断すること (単一値を無検証で確定しない)。
    """
    fld = load_field(run_dir, n_last=n_last, scale=scale, discretization=discretization)
    base = trace_cminus(fld, x_j, r_j)
    out = {"run_dir": str(run_dir), "x_j": float(x_j), "r_j": float(r_j),
           "n_snap": fld["n_snap"], "snaps": fld["snaps"],
           "x_reach": base["x_axis"], "ok": base["ok"],
           "extrap_len": base.get("extrap"), "r_stop": base.get("r_stop"),
           "path": base["path"]}
    sens = {}
    sens["start_offset"] = {f"{o}": trace_cminus(fld, x_j, r_j, start_offset=o)["x_axis"]
                            for o in (0.005, 0.01, 0.02, 0.04)}
    sens["ds"] = {f"{s}": trace_cminus(fld, x_j, r_j, ds=s)["x_axis"]
                  for s in (5e-3, 2e-3, 1e-3, 5e-4)}
    sens["n_snap"] = {}
    for n in (1, 3, 5):
        try:
            f2 = load_field(run_dir, n_last=n, scale=scale, discretization=discretization)
            sens["n_snap"][f"{n}"] = trace_cminus(f2, x_j, r_j)["x_axis"]
        except Exception as exc:
            sens["n_snap"][f"{n}"] = f"ERR {type(exc).__name__}"
    out["sensitivity"] = sens
    vals = [v for grp in sens.values() for v in grp.values() if isinstance(v, float) and np.isfinite(v)]
    out["spread"] = float(max(vals) - min(vals)) if vals else float("nan")
    return out
