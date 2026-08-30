#!/usr/bin/env python3
"""Lc_mode: from_length (dv.L_total → L_c 逆算) の統合 test。

検証対象 (plans/accepted/tooling-nozzle-axismach-length-dv.md):
  1. dv.L_total 指定で design_chain が |x_F − L_total| ≤ tol の壁を返す
     (L_c は許容窓内・x_E = x_A + L_c・設計パス数が少数)
  2. 往復一致: 解かれた L_c を Lc_mode: explicit に渡すと同一の x_F になる
  3. 窓内で実現不能な L_total は明示エラー
  4. knot 則でも同じ逆算が動く (窓が開放でも可)

実行時間 ~30 s (逆 MOC cplus n500 を数パス)。
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.evaluate.runner_axismach import design_chain  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


BASE = {
    "name": "test_length_dv",
    "type": "wind_tunnel_axisym_axismach",
    "gas": {"gamma": 1.4, "cp": 1004.5},
    "spec": {"Pt": 1.0e6, "Tt": 300.0, "M_design": 4.0, "r_throat": 0.05},
    "dv": {},
    "geometry": {"R": 2.0, "start_line": "throat_char", "wall_mode": "cplus",
                 "n_axis_inv": 500, "n_start": 41},
    "mesh": {},
    "evaluate": {},
}


def make_problem(tmp, name, dv, geom_extra):
    d = {**BASE, "name": name, "dv": dv,
         "geometry": {**BASE["geometry"], **geom_extra}}
    path = Path(tmp) / f"{name}.yaml"
    path.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
    return load_problem(path)


with tempfile.TemporaryDirectory() as tmp:
    # --- 1. from_length: L_total = 20 (M4/R2, x_F 可到達域 ~[13, 23]) ---
    L_target = 20.0
    p = make_problem(tmp, "fl20", {"L_total": {"value": L_target}},
                     {"Lc_mode": "from_length"})
    d = design_chain(p)
    xF = float(d["exit"]["x_F"])
    lo, hi = d["Lc_window"]
    print(f"info from_length: L_c={d['L_c']:.4f} 窓=({lo:.3f},{hi:.3f}) "
          f"x_F={xF:.5f} solve={d['Lc_solve']}")
    check(f"|x_F − L_total| = {abs(xF - L_target):.2e} ≤ 0.05 (n500 ノイズ床)",
          abs(xF - L_target) <= 0.05)
    check("L_c が許容窓内", lo <= d["L_c"] <= hi)
    check(f"x_E = x_A + L_c ({d['x_E']:.4f})",
          abs(d["x_E"] - (d["x_A"] + d["L_c"])) < 1e-12)
    check(f"設計パス数 {d['Lc_solve']['n_design_evals']} ≤ 6",
          d["Lc_solve"]["n_design_evals"] <= 6)

    # --- 2. 往復一致: explicit に解かれた L_c を渡すと同一 x_F ---
    p2 = make_problem(tmp, "rt", {"L_c": {"value": d["L_c"]}},
                      {"Lc_mode": "explicit"})
    d2 = design_chain(p2)
    check(f"往復一致 |Δx_F| = {abs(float(d2['exit']['x_F']) - xF):.2e} < 1e-9",
          abs(float(d2["exit"]["x_F"]) - xF) < 1e-9)

    # --- 3. 実現不能な L_total (窓上限超え) は明示エラー ---
    p3 = make_problem(tmp, "inf30", {"L_total": {"value": 30.0}},
                      {"Lc_mode": "from_length"})
    try:
        design_chain(p3)
        check("L_total=30 で ValueError", False)
    except ValueError as e:
        check(f"L_total=30 で ValueError ({str(e)[:40]}…)", "実現不能" in str(e))

    # --- 4. knot 則 (窓が開放, 生産域 M6/R3) でも逆算が動く ---
    d4dict = {**BASE, "name": "knot_m6",
              "spec": {**BASE["spec"], "M_design": 6.0},
              "dv": {"L_total": {"value": 100.0}},
              "geometry": {**BASE["geometry"], "R": 3.0,
                           "Lc_mode": "from_length", "axis_law": "knot",
                           "axis_dx0": 0.03}}
    p4_path = Path(tmp) / "knot_m6.yaml"
    p4_path.write_text(yaml.safe_dump(d4dict, sort_keys=False, allow_unicode=True))
    d4 = design_chain(load_problem(p4_path))
    xF4 = float(d4["exit"]["x_F"])
    print(f"info knot: L_c={d4['L_c']:.4f} M_K={d4['M_knot']:.3f} "
          f"x_F={xF4:.5f} solve={d4['Lc_solve']}")
    check(f"knot M6: |x_F − 100| = {abs(xF4 - 100.0):.2e} ≤ 0.05 (n500 ノイズ床)",
          abs(xF4 - 100.0) <= 0.05)

print("FAILED" if FAIL else "ALL OK", f"(fail = {FAIL})")
sys.exit(1 if FAIL else 0)
