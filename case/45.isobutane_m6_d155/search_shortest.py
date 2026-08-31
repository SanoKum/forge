"""M6 / 出口径 1.55 m (物理壁) の最短ロバスト軸則探索 (R∈{2,3} × knot 則)。

条件: イソブタン燃焼ガス φ=0.9 (case/42 と同組成), Pt 5.5 MPa, Tt 1600 K,
semi-perfect (NASA-9)。判定基準は case/42 の最短ロバスト study
(plans/accepted/tooling-nozzle-shortest-robust-axis.md) を踏襲:

  hard gate (壁QA + CFDWall validate + 壁単調 + flip/交差/堆積 0)
  + 単峰 (M'' 反転 ≤1) + min(μ_w−θ_w) ≥ 1° + min sin∠(C+,C−) ≥ 0.02

の合格集合で x_F (閉包後端 = スロート→物理出口長, r_t 単位) を最小化する。
探索は n_axis=600、最短境界のみ n_axis=1200 で確認 (600 の sin 角は楽観側)。

再生成:
  design/.venv-opt/bin/python case/45.isobutane_m6_d155/search_shortest.py
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

DESIGN = Path("/home/sano/work/forge/design")
sys.path.insert(0, str(DESIGN))

from forge_design.gas import GasSemiPerfect                                    # noqa: E402
from forge_design.geometry.axis_law import KnotQuinticAxisLaw                  # noqa: E402
from forge_design.geometry.moc_diagnostics import (                            # noqa: E402
    characteristic_topology_diagnostics, wall_margin_diagnostics)
from forge_design.geometry.moc_inverse import inverse_design                   # noqa: E402
from forge_design.geometry.transonic import HallThroat                         # noqa: E402
from forge_design.geometry.wall_axismach import (AxisMachCFDWall,              # noqa: E402
                                                 area_ratio_isentropic, wall_qa)

CASE = Path(__file__).resolve().parent
Y = {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.022, "N2": 0.7245}
TT, MD = 1600.0, 6.0
L_U, L_PIPE = 12.0, 0.5
AXIS_DX0 = 0.03
MARGIN_FLOOR_DEG = 1.0
TOPO_FLOOR = 0.02

gas = GasSemiPerfect(Y, Tt=TT)
RF = float(np.sqrt(area_ratio_isentropic(MD, gas)))
# r_inlet: 入口配管の物理半径 0.5 m / r_t。r_t は δ* 込み確定前なので探索では
# 仮値 0.0771 m を使う (r_inlet は壁上流にしか効かず x_F 最小化に無関係)。
R_T_EST = 0.0771
R_U = 0.5 / R_T_EST

_ANCHOR = {}
for _R in (2.0, 3.0):
    _ht = HallThroat(R=_R, gamma=gas.gamma_star)
    _xs = _ht.throat_characteristic(n=41)[0]
    _xA = float(_xs[0])
    _ANCHOR[_R] = (_ht, _xA, *_ht.axis_anchor(_xA))


def evaluate(R: float, L_c: float, M_K: float, n_axis: int) -> dict:
    """knot 則 → 逆 MOC → 全ゲート診断 (compare_axislaw_D.evaluate と同判定)。"""
    ht, X_A, M_A, MP_A, MPP_A = _ANCHOR[R]
    base = {"R": float(R), "L_c": float(L_c), "M_K": float(M_K), "n_axis": int(n_axis)}
    try:
        law = KnotQuinticAxisLaw(X_A, L_c, M_A, MP_A, MPP_A, MD, M_K)
        gates = law.gates()
        if gates["violations"]:
            return base | {"error": "axis gate: " + "; ".join(gates["violations"])}
        x_E = law.x_E
        x_end = x_E + 2.3 * RF * np.sqrt(MD * MD - 1.0)
        res = inverse_design(ht, lambda x: float(law(max(x, X_A))), x_axis_end=x_end,
                             n_axis=n_axis, n_start=41, gamma=gas, M_start=1.05,
                             exit_mode="characteristic", x_E=x_E, M_d=MD,
                             start_line="throat_char", wall_mode="cplus",
                             axis_dx0=AXIS_DX0)
        w = res["wall"]
        qa = wall_qa(w, MD, x_E, gas, R=R)
        try:
            cw = AxisMachCFDWall(w, R=R, r_U=R_U, L_U=L_U, L_pipe=L_PIPE)
            cfd_msgs = cw.validate()
        except Exception as e:                                   # noqa: BLE001
            cfd_msgs = [f"CFDWall 構築失敗: {e}"]
        topo = characteristic_topology_diagnostics(
            res["levels"], n_axis - 1,
            n_crossing_sample=(20 if n_axis <= 600 else 40),
            seam_i_pad=90, seam_k_max=110, seam_x_max=X_A + 1.0, wall=w[:, :2])
        margin = wall_margin_diagnostics(w)
        wall_mono = bool(np.all(np.diff(w[:, 1]) >= -1e-9))
        hard = (not qa["violations"]) and (not cfd_msgs) and wall_mono \
            and topo["n_orientation_flip_interior"] == 0 \
            and topo["n_crossings_sampled"] == 0 \
            and topo["pileup_pairs_interior"] == 0
        return base | {
            "x_A": X_A, "x_E": float(x_E), "x_F": float(qa["x_F"]),
            "r_F": float(qa["r_F"]), "theta_max_deg": float(qa["theta_max_deg"]),
            "min_mu_minus_theta_deg": float(margin["min_mu_minus_theta_deg"]),
            "min_sin_angle_interior": float(topo.get("min_sin_angle_interior",
                                                     float("nan"))),
            "interior_flip": int(topo["n_orientation_flip_interior"]),
            "crossings": int(topo["n_crossings_sampled"]),
            "pileups": int(topo["pileup_pairs_interior"]),
            "mpp_quality_ok": bool(gates["mpp_quality_ok"]),
            "hard_gate_pass": bool(hard),
            "wall_qa_violations": qa["violations"], "cfd_wall_violations": cfd_msgs,
        }
    except Exception as exc:                                     # noqa: BLE001
        return base | {"error": f"{type(exc).__name__}: {exc}"}


def robust(rec: dict) -> bool:
    if "error" in rec or not rec["hard_gate_pass"] or not rec["mpp_quality_ok"]:
        return False
    ts = rec["min_sin_angle_interior"]
    return (rec["min_mu_minus_theta_deg"] >= MARGIN_FLOOR_DEG
            and np.isfinite(ts) and ts >= TOPO_FLOOR)


def _job(args):
    return evaluate(*args)


def run_grid(jobs, n_axis):
    workers = min({600: 6, 1200: 4}.get(n_axis, 2), os.cpu_count() or 1, len(jobs))
    if workers <= 1:
        return [_job(j) for j in jobs]
    with ProcessPoolExecutor(max_workers=workers,
                             mp_context=mp.get_context("fork")) as pool:
        return list(pool.map(_job, jobs, chunksize=1))


def grid(lo, hi, step):
    n = int(round((hi - lo) / step))
    return [round(lo + i * step, 10) for i in range(n + 1)]


def stage(R, lcs, mks, n_axis, name):
    t0 = time.time()
    recs = run_grid([(R, L, M, n_axis) for L in lcs for M in mks], n_axis)
    ok = [r for r in recs if robust(r)]
    print(f"R={R} {name}: {len(recs)} pts ({time.time()-t0:.1f} s), "
          f"robust {len(ok)}", flush=True)
    return recs, ok


def shortest(recs):
    ok = [r for r in recs if robust(r)]
    return min(ok, key=lambda r: r["x_F"]) if ok else None


def main():
    all_recs = []
    winners = {}
    for R in (2.0, 3.0):
        # coarse: L_c 36–48 (R3 の既往最短 39.05 @Tt1550 を跨ぐ)
        c, _ = stage(R, grid(36.0, 48.0, 1.0), grid(2.2, 3.2, 0.1), 600, "coarse")
        all_recs += c
        best = shortest(c)
        if best is None:
            print(f"R={R}: coarse でロバスト候補なし", flush=True)
            continue
        print(f"  coarse best: L_c={best['L_c']} M_K={best['M_K']} "
              f"x_F={best['x_F']:.3f}", flush=True)
        # fine: 最短 L_c の下 1.0 を 0.1 刻み、M_K ±0.3 を 0.05 刻み
        f, _ = stage(R, grid(best["L_c"] - 1.0, best["L_c"], 0.1),
                     grid(best["M_K"] - 0.3, best["M_K"] + 0.3, 0.05), 600, "fine")
        all_recs += f
        bf = shortest(f) or best
        # 確認: 最短境界の候補 (L_c 昇順・各 L_c の margin 最大) を n1200 で再評価し、
        # 最初に robust が立つ L_c を採る (600 の sin 角楽観の補正)
        cand = sorted([r for r in f if robust(r)], key=lambda r: (r["L_c"], -r["min_mu_minus_theta_deg"]))
        seen, probe = set(), []
        for r in cand:
            if r["L_c"] not in seen:
                seen.add(r["L_c"])
                probe.append(r)
            if len(probe) >= 6:
                break
        v = run_grid([(R, r["L_c"], r["M_K"], 1200) for r in probe], 1200)
        all_recs += v
        vw = shortest(v)
        winners[R] = vw or bf
        tag = "n1200" if vw else "n600(1200未確認)"
        w = winners[R]
        print(f"R={R} shortest [{tag}]: L_c={w['L_c']} M_K={w['M_K']} "
              f"x_F={w['x_F']:.4f} margin={w['min_mu_minus_theta_deg']:.3f}° "
              f"topo={w['min_sin_angle_interior']:.4f} "
              f"theta_max={w['theta_max_deg']:.2f}°", flush=True)

    out = {"conditions": {"Y": Y, "Tt": TT, "Md": MD, "gamma_star": gas.gamma_star,
                          "rF_rt": RF, "L_U": L_U, "margin_floor_deg": MARGIN_FLOOR_DEG,
                          "topo_floor": TOPO_FLOOR},
           "winners": winners,
           "records": [{k: v for k, v in r.items()
                        if k not in ("wall_qa_violations", "cfd_wall_violations")
                        or v} for r in all_recs]}
    (CASE / "search_shortest.json").write_text(json.dumps(out, indent=1))
    print("saved:", CASE / "search_shortest.json", flush=True)


if __name__ == "__main__":
    main()
