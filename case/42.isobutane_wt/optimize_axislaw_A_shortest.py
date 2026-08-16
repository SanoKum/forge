"""軸 Mach 則 A の最短ロバスト長さを探索する。

目的は壁角の単純最小化ではなく、hard gate・単峰・mu-theta余裕・特性線交差角余裕を
制約として満たす候補の閉包後端 x_F を最小化すること。M6/R3では x_F-x_E がほぼ固定なので、
探索上は L_c の最小化と同値になる。

解像度は全探索 n_axis=600、最短境界だけ n_axis=1200。2400点は使わない。

再生成:
  design/.venv-opt/bin/python case/42.isobutane_wt/optimize_axislaw_A_shortest.py
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


CASE = Path(__file__).resolve().parent
_SWEEP_PATH = CASE / "sweep_axislaw_A_MK.py"
_SPEC = importlib.util.spec_from_file_location("sweep_axislaw_A_MK_for_shortest", _SWEEP_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"感度モジュールを読み込めない: {_SWEEP_PATH}")
_SWEEP = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _SWEEP
_SPEC.loader.exec_module(_SWEEP)

MARGIN_FLOOR_DEG = 1.0
TOPO_FLOOR = 0.02


def _grid(lo: float, hi: float, step: float) -> list[float]:
    """浮動小数の端点漏れを避けた閉区間格子。"""
    n = int(round((hi - lo) / step))
    return [round(lo + i * step, 10) for i in range(n + 1)]


def _valid_for_topology(rec: dict) -> bool:
    return _SWEEP._valid(rec) and _SWEEP._topo(rec) >= TOPO_FLOOR


def _robust(rec: dict) -> bool:
    return _valid_for_topology(rec) and _SWEEP._margin(rec) >= MARGIN_FLOOR_DEG


def _best_margin(records: list[dict], L_c: float) -> dict | None:
    candidates = [r for r in records if r.get("L_c") == L_c and _valid_for_topology(r)]
    return max(candidates, key=_SWEEP._margin) if candidates else None


def _first_feasible_lc(records: list[dict], lcs: list[float]) -> tuple[float, dict]:
    for L_c in lcs:
        best = _best_margin(records, L_c)
        if best is not None and _SWEEP._margin(best) >= MARGIN_FLOOR_DEG:
            return L_c, best
    raise RuntimeError(f"探索範囲 {lcs[0]}..{lcs[-1]} にロバスト候補がない")


def _run_stage(name: str, lcs: list[float], mks: list[float]) -> list[dict]:
    print(f"--- {name}: n_axis=600, {len(lcs)} x {len(mks)} ---", flush=True)
    records = _SWEEP.evaluate_many([(L_c, M_K, 600) for L_c in lcs for M_K in mks], 600)
    first_lc, best = _first_feasible_lc(records, lcs)
    print(
        f"{name}: first feasible L_c={first_lc:.3f}, M_K={best['M_K']:.3f}, "
        f"margin={_SWEEP._margin(best):.6f}, topo={_SWEEP._topo(best):.6f}",
        flush=True,
    )
    for rec in records:
        rec["search_stage"] = name
    return records


def _compact_summary(rec: dict) -> dict:
    return {
        "L_c": rec["L_c"],
        "M_K": rec["M_K"],
        "n_axis": rec["n_axis"],
        "x_F": rec["x_F"],
        "x_E": rec["x_E"],
        "xF_minus_xE": rec["xF_minus_xE"],
        "theta_max_deg": rec["theta_max_deg"],
        "min_mu_minus_theta_deg": _SWEEP._margin(rec),
        "min_sin_angle_interior": _SWEEP._topo(rec),
        "hard_gate_pass": rec["hard_gate_pass"],
        "mpp_quality_ok": rec["axis_gates"]["mpp_quality_ok"],
        "interior_flip_count": rec["topology"]["n_orientation_flip_interior"],
        "crossing_count": rec["topology"]["n_crossings_sampled"],
    }


def main() -> None:
    all_records: list[dict] = []

    coarse_lcs = _grid(35.0, 40.0, 0.5)
    coarse_mks = _grid(2.2, 3.2, 0.1)
    coarse = _run_stage("coarse", coarse_lcs, coarse_mks)
    all_records.extend(coarse)
    coarse_lc, coarse_best = _first_feasible_lc(coarse, coarse_lcs)

    fine_lcs = _grid(coarse_lc - 0.5, coarse_lc, 0.05)
    fine_mks = _grid(float(coarse_best["M_K"]) - 0.30,
                     float(coarse_best["M_K"]) + 0.30, 0.02)
    fine = _run_stage("fine", fine_lcs, fine_mks)
    all_records.extend(fine)
    fine_lc, fine_best = _first_feasible_lc(fine, fine_lcs)

    micro_lcs = _grid(fine_lc - 0.05, fine_lc + 0.10, 0.01)
    micro_mks = _grid(float(fine_best["M_K"]) - 0.15,
                      float(fine_best["M_K"]) + 0.15, 0.01)
    micro = _run_stage("micro", micro_lcs, micro_mks)
    all_records.extend(micro)
    micro_lc, _ = _first_feasible_lc(micro, micro_lcs)

    # 600点で得た最短境界から長い側へ、各長さの最大margin候補だけを1200点で先に確認する。
    probe_sources = [_best_margin(micro, L_c) for L_c in micro_lcs if L_c >= micro_lc]
    probe_sources = [r for r in probe_sources if r is not None]
    print(f"--- boundary probe: n_axis=1200, {len(probe_sources)} candidates ---", flush=True)
    probes = _SWEEP.evaluate_many(
        [(float(r["L_c"]), float(r["M_K"]), 1200) for r in probe_sources], 1200)
    for rec in probes:
        rec["search_stage"] = "boundary_probe"
    all_records.extend(probes)
    passing_probes = sorted((r for r in probes if _robust(r)), key=lambda r: r["L_c"])
    if not passing_probes:
        raise RuntimeError("1200点で合格する境界候補がない。micro探索範囲を長い側へ拡張する必要がある")
    verified_lc = float(passing_probes[0]["L_c"])

    # 最短長さが確定した後、その長さのmicro M_K全点を1200点で再評価し、二次目的を決める。
    print(f"--- final M_K at L_c={verified_lc:.3f}: n_axis=1200 ---", flush=True)
    final_1200 = _SWEEP.evaluate_many(
        [(verified_lc, M_K, 1200) for M_K in micro_mks], 1200)
    for rec in final_1200:
        rec["search_stage"] = "final_MK"
    all_records.extend(final_1200)
    final_candidates = [r for r in final_1200 if _robust(r)]
    if not final_candidates:
        raise RuntimeError("確認済み最短長さでロバストなM_Kがない")
    # 同一 L_c では x_F の1e-8程度の数値差を目的差とみなさず、壁角を二次目的にする。
    winner = min(final_candidates, key=lambda r: (r["theta_max_deg"], -_SWEEP._margin(r)))

    json_path = CASE / "optimize_axislaw_A_shortest.json"
    csv_path = CASE / "optimize_axislaw_A_shortest_summary.csv"
    result = {
        "objective": "hard gates and engineering margins, then min x_F; at equal L_c min theta_max",
        "constraints": {
            "min_mu_minus_theta_deg": MARGIN_FLOOR_DEG,
            "min_sin_angle_interior": TOPO_FLOOR,
            "single_peak_Mpp": True,
            "hard_gate_pass": True,
        },
        "resolution_policy": {"search": 600, "final": 1200, "n_axis_2400_used": False},
        "grids": {
            "coarse": {"L_c": coarse_lcs, "M_K": coarse_mks},
            "fine": {"L_c": fine_lcs, "M_K": fine_mks},
            "micro": {"L_c": micro_lcs, "M_K": micro_mks},
        },
        "winner": _compact_summary(winner),
        "passing_MK_at_winner_Lc": [_compact_summary(r) for r in final_candidates],
        "records": all_records,
    }
    json_path.write_text(json.dumps(result, indent=1, allow_nan=True) + "\n")

    header = [
        "stage", "L_c", "M_K", "n_axis", "robust", "hard_gate", "single_peak",
        "x_F", "x_E", "xF_minus_xE", "theta_max_deg", "min_mu_minus_theta_deg",
        "min_sin_angle_interior", "flip_interior", "crossings",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(header)
        for rec in all_records:
            if "error" in rec:
                continue
            writer.writerow([
                rec["search_stage"], rec["L_c"], rec["M_K"], rec["n_axis"], _robust(rec),
                rec["hard_gate_pass"], rec["axis_gates"]["mpp_quality_ok"], rec["x_F"],
                rec["x_E"], rec["xF_minus_xE"], rec["theta_max_deg"], _SWEEP._margin(rec),
                _SWEEP._topo(rec), rec["topology"]["n_orientation_flip_interior"],
                rec["topology"]["n_crossings_sampled"],
            ])

    print("WINNER", json.dumps(_compact_summary(winner), indent=1), flush=True)
    print(f"saved {json_path}", flush=True)
    print(f"saved {csv_path}", flush=True)


if __name__ == "__main__":
    main()
