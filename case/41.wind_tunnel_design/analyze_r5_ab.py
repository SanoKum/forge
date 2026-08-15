#!/usr/bin/env python3
"""W0 判別実験 (R=2 vs R=5) の計測一式 (walldriven plan W0 / 調査 §8)。

対照 (R=2, node): run_0034_node_base — 公表値 こぶ振幅 0.0180@1.44 / 谷 −0.0231。
被験 (R=5, node): run_0039_r5_cold。

計測項目:
  1. 意図非依存うねり指標: node 軸 DOF の M(x) (末尾3スナップ平均) の
     多項式トレンド (3次/5次) 残差の山と谷。**同一コードを両 run に適用**。
  2. x_reach,CFD (壁面状態起点の C⁻ 追跡, 感度込み) と、そこでの (M, M', M'')
     (軸 M の左片側局所多項式フィット, 次数・窓の感度併記)。
  3. 壁 dM/dx の J 両側線形近似 (res_wall_3 境界出力)。
  4. mdot 比 (スロート断面の ∫ρu·2πr dr / 等エントロピー 1D)。
  5. 出口 ε_M / |θ|max (core traced)。
使い方: design/.venv-opt/bin/python analyze_r5_ab.py
"""
import json
import sys
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE.parents[1] / "design"))

from forge_design.geometry.cminus_cfd import (  # noqa: E402
    extract_x_reach, load_field, wall_state_at)
from forge_design.metrics.extract import exit_uniformity  # noqa: E402

GAMMA, CP, SCALE = 1.4, 1004.5, 0.01
RUNS = {
    "R2": {"dir": CASE / "run_0034_node_base", "R": 2.0,
           "x_j": 0.24654255216465334, "r_j": 1.0152539784717696, "x_d": 14.0},
    "R5": {"dir": CASE / "run_0039_r5_cold", "R": 5.0,
           "x_j": 0.22325764104830179, "r_j": 1.0049843974286452, "x_d": 14.0},
}


def axis_M_node(run_dir, n_last=3):
    """node 軸 DOF (r=0) の M(x) 末尾スナップ平均。"""
    fld = load_field(run_dir, n_last=n_last, scale=SCALE, discretization="node")
    m = fld["r"] < 1e-8
    o = np.argsort(fld["x"][m])
    return fld["x"][m][o], fld["M"][m][o]


def waviness(x, M, x_lo, x_hi, deg):
    w = (x >= x_lo) & (x <= x_hi)
    c = np.polyfit(x[w], M[w], deg)
    r = M[w] - np.polyval(c, x[w])
    i, j = int(np.argmax(r)), int(np.argmin(r))
    return {"deg": deg, "window": [x_lo, x_hi],
            "amp": float(r[i]), "x_amp": float(x[w][i]),
            "trough": float(r[j]), "x_trough": float(x[w][j]),
            "n_station": int(w.sum())}


def local_anchor_left(x, M, x_c, degree=3, window=0.6):
    m = (x >= x_c - window) & (x <= x_c)
    if m.sum() < degree + 2:
        return None
    c = np.polyfit(x[m] - x_c, M[m], degree)
    return {"M": float(c[-1]), "Mp": float(c[-2]), "Mpp": float(2.0 * c[-3]),
            "degree": degree, "window": window, "n": int(m.sum())}


def wall_dMdx_jump(run_dir, x_j, half=0.15, n_last=3):
    """res_wall_3 の壁 M(x) を J 両側 half 幅で線形近似し勾配比を返す。"""
    import h5py
    rd = Path(run_dir)
    fs = sorted(rd.glob("res_wall_3_[0-9]*.h5"),
                key=lambda f: int(f.stem.split("_")[-1]))
    fs = [f for f in fs if int(f.stem.split("_")[-1]) > 0][-n_last:]
    Ms, xs = [], None
    for f in fs:
        with h5py.File(f) as h:
            cc = h["MESH/COORD"][:].reshape(-1, 3)
            Ux, Uy = h["VALUE/Ux"][:], h["VALUE/Uy"][:]
            ro, Ps = h["VALUE/ro"][:], h["VALUE/Ps"][:]
        xf = (0.5 * (cc[:-1, 0] + cc[1:, 0]) / SCALE
              if len(cc) == len(Ux) + 1 else cc[:len(Ux), 0] / SCALE)
        a = np.sqrt(GAMMA * np.maximum(Ps, 1e-9) / np.maximum(ro, 1e-12))
        Ms.append(np.hypot(Ux, Uy) / a)
        xs = xf
    o = np.argsort(xs)
    xs, M = xs[o], np.mean(Ms, axis=0)[o]
    up = (xs >= x_j - half) & (xs < x_j)
    dn = (xs > x_j) & (xs <= x_j + half)
    su = float(np.polyfit(xs[up], M[up], 1)[0])
    sd = float(np.polyfit(xs[dn], M[dn], 1)[0])
    return {"slope_up": su, "slope_dn": sd, "ratio": sd / su, "half": half,
            "n_up": int(up.sum()), "n_dn": int(dn.sum())}


def mdot_ratio_cfd(run_dir, n_last=3):
    """最終スナップのスロート最寄り列で ∫ρUx 2πr dr / 等エントロピー 1D。"""
    import h5py
    rd = Path(run_dir)
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    with h5py.File(rd / "nozzle.h5") as nz:
        cc = nz["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(res) as f:
        ro, Ux = f["/VALUE/ro"][:], f["/VALUE/Ux"][:]
    x, r = cc[:, 0] / SCALE, cc[:, 1] / SCALE
    cols = np.round(x, 6)
    xc = np.unique(cols)
    x0 = xc[np.argmin(np.abs(xc))]          # スロート最寄りの列
    s = cols == x0
    o = np.argsort(r[s])
    rr, f_ru = r[s][o] * SCALE, (ro[s] * Ux[s])[o]
    mdot = float(np.trapezoid(f_ru * 2.0 * np.pi * rr, rr))
    # 等エントロピー 1D (Pt=1e6, Tt=800, A*=π r_t²)
    Pt, Tt = 1.0e6, 800.0
    Rg = CP * (GAMMA - 1.0) / GAMMA
    rost = Pt / (Rg * Tt) * (2.0 / (GAMMA + 1.0)) ** (1.0 / (GAMMA - 1.0))
    ast = np.sqrt(GAMMA * Rg * Tt * 2.0 / (GAMMA + 1.0))
    mdot_1d = rost * ast * np.pi * SCALE ** 2
    return {"x_col": float(x0), "mdot": mdot, "mdot_1d": float(mdot_1d),
            "ratio": mdot / mdot_1d}


def main():
    out = {}
    for key, cfg in RUNS.items():
        rd = cfg["dir"]
        r = {"run": rd.name, "R": cfg["R"]}
        x, M = axis_M_node(rd)
        np.savetxt(rd / "axis_M_node.csv", np.c_[x, M], delimiter=",",
                   header="x_rt,M", comments="")
        # 1. うねり (基準窓 0.55-2.5 + 次数感度)
        r["waviness"] = {f"deg{d}": waviness(x, M, 0.55, 2.5, d) for d in (3, 5)}
        # 2. x_reach + アンカー
        xr = extract_x_reach(rd, cfg["x_j"], cfg["r_j"], scale=SCALE,
                             discretization="node")
        xr.pop("path", None)
        r["x_reach"] = {k: xr[k] for k in ("x_reach", "ok", "extrap_len", "r_stop", "sens")
                        if k in xr}
        xrv = xr["x_reach"]
        if xrv is not None:
            r["anchor_at_x_reach"] = {
                f"deg{d}_win{wn}": local_anchor_left(x, M, xrv, d, wn)
                for d in (2, 3) for wn in (0.4, 0.6)}
        # 3. 壁 dM/dx 跳び
        r["wall_dMdx_J"] = wall_dMdx_jump(rd, cfg["x_j"])
        # 4. mdot
        r["mdot"] = mdot_ratio_cfd(rd)
        # 5. 出口一様性
        res = sorted(rd.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        try:
            r["exit"] = exit_uniformity(rd / "nozzle.h5", res, 4.0, x_d=cfg["x_d"],
                                        gamma=GAMMA)
        except Exception as exc:
            r["exit"] = {"error": f"{type(exc).__name__}: {exc}"}
        out[key] = r
    (CASE / "r5_ab_metrics.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
