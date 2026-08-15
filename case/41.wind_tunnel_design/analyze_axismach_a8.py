#!/usr/bin/env python3
"""A8 A/B: 初期値線 (縦線 vs スロート特性線) の同条件比較。

plan: plans/active/tooling-nozzle-axismach-throat-characteristic.md §6

対照 run_0049_axismach_a8v (start_line: vertical) / 本命 run_0050_axismach_a8c
(start_line: throat_char)。両 run は Hall アンカー pass0・n_axis 2000・
exit_mode characteristic・Lc_mode max で、差分は初期値線と ni (壁長比例) のみ。

計測:
  1. 設計側 (prepare_info.json): x_A / x_E / x_F / r_F 誤差 / MOC mdot リーク / θ_max
  2. CFD 軸 M vs target (metrics.json): ΔM_max/Md, ΔM rms, overshoot
  3. 意図非依存うねり (analyze_r5_ab と同一コード): 3 次トレンド残差の山谷 (窓 0.55–2.5)
  4. **到達不能域**: 設計壁始点発 C⁻ の CFD 実測軸着地 x_reach,CFD と、
     軸 target 開始点 x_A の差 = 「設計が指定しているのに壁が到達できない帯」
  5. 出口一様性 (metrics.json の exit_uniformity)
使い方: design/.venv-opt/bin/python analyze_axismach_a8.py
"""
import json
import sys
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE.parents[1] / "design"))

from forge_design.geometry.cminus_cfd import extract_x_reach, load_field  # noqa: E402

SCALE = 0.01
RUNS = {"vertical": CASE / "run_0049_axismach_a8v",
        "throat_char": CASE / "run_0050_axismach_a8c"}
WAVE_WINDOW = (0.55, 2.5)   # W0/A5 と同一 (円弧こぶ帯)


def axis_M_node(run_dir, n_last=3):
    fld = load_field(run_dir, n_last=n_last, scale=SCALE, discretization="node")
    m = fld["r"] < 1e-8
    o = np.argsort(fld["x"][m])
    return fld["x"][m][o], fld["M"][m][o]


def waviness(x, M, x_lo, x_hi, deg=3):
    w = (x >= x_lo) & (x <= x_hi)
    r = M[w] - np.polyval(np.polyfit(x[w], M[w], deg), x[w])
    i, j = int(np.argmax(r)), int(np.argmin(r))
    return {"deg": deg, "window": [x_lo, x_hi], "n_station": int(w.sum()),
            "amp": float(r[i]), "x_amp": float(x[w][i]),
            "trough": float(r[j]), "x_trough": float(x[w][j]),
            "p2p": float(r[i] - r[j])}


def main():
    out = {}
    for tag, rd in RUNS.items():
        info = json.loads((rd / "prepare_info.json").read_text())
        met = json.loads((rd / "metrics.json").read_text())
        wall = np.loadtxt(rd / "wall_design.csv", delimiter=",", skiprows=1)
        x_w0, r_w0 = float(wall[0, 0]) / SCALE, float(wall[0, 1]) / SCALE
        xa, Ma = axis_M_node(rd)
        reach = extract_x_reach(rd, x_w0, r_w0, scale=SCALE, discretization="node")
        sens = [v for d in reach["sensitivity"].values() for v in d.values()
                if isinstance(v, float) and np.isfinite(v)]
        out[tag] = {
            "run": rd.name,
            "design": {"x_A": info["x_A"], "x_E": info["x_E"],
                       "x_F": info["qa"]["x_F"], "r_F_err_rel": info["qa"]["r_F_err_rel"],
                       "theta_max_deg": info["qa"]["theta_max_deg"],
                       "mdot_ratio_moc": info["mdot_ratio_moc"],
                       "wall_start": [x_w0, r_w0]},
            "cfd": {"dM_max_rel_Md": met["dM_max_rel_Md"], "dM_rms": met["dM_rms"],
                    "overshoot_rel": met["overshoot_rel"],
                    "exit": {k: met["exit_uniformity"].get(k)
                             for k in ("eps_M_rms", "eps_M_max",
                                       "eps_theta_max_deg", "r_core", "r_wall",
                                       "n_core_faces", "n_exit_faces")}},
            "waviness": waviness(xa, Ma, *WAVE_WINDOW),
            "reach": {"x_reach_cfd": reach["x_reach"],
                      "spread": reach["spread"],
                      "unreachable_band": float(reach["x_reach"] - info["x_A"])},
        }
    (CASE / "axismach_a8_ab.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    print("\n--- 要約 ---")
    hdr = f"{'指標':28s} {'vertical':>14s} {'throat_char':>14s}"
    print(hdr)
    rows = [("設計壁始点 x_w0", lambda d: d["design"]["wall_start"][0]),
            ("軸 target 開始 x_A", lambda d: d["design"]["x_A"]),
            ("x_reach,CFD", lambda d: d["reach"]["x_reach_cfd"]),
            ("到達不能帯 (x_reach−x_A)", lambda d: d["reach"]["unreachable_band"]),
            ("x_E", lambda d: d["design"]["x_E"]),
            ("x_F", lambda d: d["design"]["x_F"]),
            ("r_F 誤差 [%]", lambda d: 100 * d["design"]["r_F_err_rel"]),
            ("MOC mdot 比", lambda d: d["design"]["mdot_ratio_moc"]),
            ("ΔM_max [% Md]", lambda d: 100 * d["cfd"]["dM_max_rel_Md"]),
            ("ΔM rms", lambda d: d["cfd"]["dM_rms"]),
            ("うねり amp", lambda d: d["waviness"]["amp"]),
            ("うねり p2p", lambda d: d["waviness"]["p2p"]),
            ("出口 ε_M rms [%]", lambda d: 100 * (d["cfd"]["exit"]["eps_M_rms"] or np.nan)),
            ("出口 |θ|max [deg]", lambda d: d["cfd"]["exit"]["eps_theta_max_deg"]),
            ("出口コア面/全面", lambda d: d["cfd"]["exit"]["n_core_faces"]
             / max(d["cfd"]["exit"]["n_exit_faces"], 1)),
            ]
    for name, f in rows:
        print(f"{name:28s} {f(out['vertical']):14.5f} {f(out['throat_char']):14.5f}")


if __name__ == "__main__":
    main()
