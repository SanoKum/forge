#!/usr/bin/env python3
"""W3 検証: run_0040_walldriven_w3 (poly R_t=5) の受け入れ計測。

1. D 発 C⁻ (W4 のデータ線) を壁面状態起点でトレース → 軸着地 x_land が
   ドメイン内 (x_e より margin 上流) にあるか = L_ext_factor の実測検証。
2. データ線上の P0 変動 (等エントロピー診断) → W4 ゲート初期値
   (≤0.1% PASS / 0.1–0.5% WARN / >0.5% FAIL) の予備判定。
3. 軸 M のうねり (同一 3 次トレンド, 窓 0.55–2.5) を run_0039 (円弧 R=5) と比較
   — 「同じ R でも接合特異性が J (超音速域) から T (遷音速点) へ移ると波源が
   弱いか」(walldriven plan §2.2 仮説) の初回データ点。
使い方: design/.venv-opt/bin/python analyze_walldriven_w3.py
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

CASE = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE.parents[1] / "design"))

from forge_design.geometry.cminus_cfd import (  # noqa: E402
    load_field, trace_cminus_from_wall)

GAMMA, SCALE = 1.4, 0.01
RUN = CASE / "run_0042_walldriven_w3c"
ARC = CASE / "run_0039_r5_cold"


def axis_M_node(run_dir, n_last=3):
    fld = load_field(run_dir, n_last=n_last, scale=SCALE, discretization="node")
    m = fld["r"] < 1e-8
    o = np.argsort(fld["x"][m])
    return fld["x"][m][o], fld["M"][m][o]


def waviness(x, M, x_lo, x_hi, deg=3):
    w = (x >= x_lo) & (x <= x_hi)
    c = np.polyfit(x[w], M[w], deg)
    r = M[w] - np.polyval(c, x[w])
    i, j = int(np.argmax(r)), int(np.argmin(r))
    return {"amp": float(r[i]), "x_amp": float(x[w][i]),
            "trough": float(r[j]), "x_trough": float(x[w][j]),
            "deg": deg, "window": [x_lo, x_hi]}


def p0_along_path(run_dir, path, n_last=3):
    """データ線 (x,r 点列, r* 単位) 上の P0/Pt 変動。node 座標 + 末尾平均 P, M。"""
    from scipy.interpolate import LinearNDInterpolator
    rd = Path(run_dir)
    with h5py.File(rd / "nozzle.h5") as nz:
        cc = nz["/MESH/COORD"][:].reshape(-1, 3)
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    res = [f for f in res if int("".join(c for c in f.stem if c.isdigit())) > 0][-n_last:]
    Ps, Ms = [], []
    for f in res:
        with h5py.File(f) as h:
            P = h["/VALUE/P"][:]
            Ux, Uy, son = h["/VALUE/Ux"][:], h["/VALUE/Uy"][:], h["/VALUE/sonic"][:]
        Ps.append(P)
        Ms.append(np.hypot(Ux, Uy) / np.maximum(son, 1e-9))
    P = np.mean(Ps, axis=0)
    M = np.mean(Ms, axis=0)
    pts = np.c_[cc[:, 0] / SCALE, cc[:, 1] / SCALE]
    fP = LinearNDInterpolator(pts, P)
    fM = LinearNDInterpolator(pts, M)
    Pp, Mp = fP(path[:, 0], path[:, 1]), fM(path[:, 0], path[:, 1])
    ok = np.isfinite(Pp) & np.isfinite(Mp)
    P0 = Pp[ok] * (1.0 + 0.5 * (GAMMA - 1.0) * Mp[ok] ** 2) ** (GAMMA / (GAMMA - 1.0))
    rel = (P0.max() - P0.min()) / np.mean(P0)
    verdict = "PASS" if rel <= 1e-3 else ("WARN" if rel <= 5e-3 else "FAIL")
    return {"n_pts": int(ok.sum()), "P0_mean_Pa": float(np.mean(P0)),
            "dP0_rel": float(rel), "gate": verdict,
            "P0_over_Pt_mean": float(np.mean(P0) / 1.0e6)}


def main():
    info = json.loads((RUN / "prepare_info.json").read_text())
    x_D, r_D = info["x_D"], info["r_D"]
    th_D = np.deg2rad(info["theta_D_deg"])
    out = {"run": RUN.name}

    # 1. D 発 C⁻ (壁面状態起点)
    fld = load_field(RUN, n_last=3, scale=SCALE, discretization="node")
    tr = trace_cminus_from_wall(fld, RUN, x_D, r_D, scale=SCALE, gamma=GAMMA,
                                th_wall_geom=th_D)
    path = np.asarray(tr["path"])
    out["dataline"] = {"ok": bool(tr["ok"]), "x_land": float(tr["x_axis"]),
                       "x_e": info["x_e"],
                       "margin_to_exit": float(info["x_e"] - tr["x_axis"]),
                       "wall_state": tr["wall_state"]}
    np.savetxt(RUN / "dataline_cminus_from_D.csv", path, delimiter=",",
               header="x_rt,r_rt", comments="")

    # 2. P0 診断 (壁直近 5% は境界補間の混入を避けて除外)
    out["p0_diag"] = p0_along_path(RUN, path[int(0.05 * len(path)):])

    # 3. うねり: poly R_t=5 vs arc R=5 (同一コード・同一窓)
    x, M = axis_M_node(RUN)
    np.savetxt(RUN / "axis_M_node.csv", np.c_[x, M], delimiter=",",
               header="x_rt,M", comments="")
    xa, Ma = axis_M_node(ARC)
    out["waviness"] = {
        "poly_Rt5_win055_25": waviness(x, M, 0.55, 2.5),
        "arc_R5_win055_25": waviness(xa, Ma, 0.55, 2.5),
        "poly_Rt5_unreach": waviness(x, M, 0.55, float(tr["x_axis"])),
    }
    (RUN / "w3_metrics.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
