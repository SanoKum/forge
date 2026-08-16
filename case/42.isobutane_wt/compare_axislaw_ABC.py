"""軸 Mach 則 A (knot) / B (monotone B-spline) / C (非負 dν/dx B-spline) の比較。

計画: plans/accepted/tooling-nozzle-axislaw-smoothness.md。基準条件: case/42
problem_ib_m6_R3_Lc45_MK2.5.yaml 相当 (R=3, L_c=45 r_t, M_d=6, semi-perfect Tt=1550K,
Hall throat anchor, start_line=throat_char, wall_mode=cplus, exit_mode=characteristic,
axis_dx0=0.03)。n_axis in (600, 1200, 2400) の 3 解像度で軸則・特性線網トポロジ・
壁 fairness を比較する。

再生成: design/.venv-opt/bin/python case/42.isobutane_wt/compare_axislaw_ABC.py
出力: case/42.isobutane_wt/compare_axislaw_ABC.json、
      /tmp (Artifact 用) の PNG 一式 (実行時に表示するパスへ)。
"""
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

DESIGN = Path("/home/sano/work/forge/design")
sys.path.insert(0, str(DESIGN))

from forge_design.gas import GasSemiPerfect                                    # noqa: E402
from forge_design.geometry.axis_law import KnotQuinticAxisLaw                  # noqa: E402
from forge_design.geometry.axis_law_bspline import (MonotoneBSplineAxisLaw,    # noqa: E402
                                                     NonnegDnuBSplineAxisLaw)
from forge_design.geometry.moc_diagnostics import (                            # noqa: E402
    axis_smoothness_diagnostics, characteristic_topology_diagnostics,
    wall_curvature_diagnostics, wall_margin_diagnostics)
from forge_design.geometry.moc_inverse import inverse_design                   # noqa: E402
from forge_design.geometry.transonic import HallThroat                        # noqa: E402
from forge_design.geometry.wall_axismach import (AxisMachCFDWall,             # noqa: E402
                                                 area_ratio_isentropic, wall_qa)

OUT = Path("/tmp/claude-1000/-home-sano-work-forge/b78d653c-d741-4dd8-bcfc-5a359da19ae8/scratchpad")
CASE = Path("/home/sano/work/forge/case/42.isobutane_wt")

# --- 基準条件 (case/42 problem_ib_m6_R3_Lc45_MK2.5.yaml と同一) ------------------
Y = {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.022, "N2": 0.7245}
TT, MD, R, LC, MK = 1550.0, 6.0, 3.0, 45.0, 2.5
R_U, L_U, L_PIPE = 9.303, 12.0, 0.5
RESOLUTIONS = (600, 1200, 2400)
AXIS_DX0 = 0.03


def build_laws(x_A, M_A, Mp_A, Mpp_A, gas):
    x_E = x_A + LC
    return {
        "A_knot": KnotQuinticAxisLaw(x_A, LC, M_A, Mp_A, Mpp_A, MD, MK),
        "B_bspline_M": MonotoneBSplineAxisLaw(x_A, x_E, M_A, Mp_A, Mpp_A, MD,
                                              n_interior=15),
        "C_bspline_dnu": NonnegDnuBSplineAxisLaw(x_A, x_E, M_A, Mp_A, Mpp_A, MD, gas,
                                                 n_interior=20),
    }


def breakpoints_of(tag, law):
    if tag == "A_knot":
        return [law.x_K]
    return None


def run_one(tag, law, ht, gas, x_A, x_E, n_axis):
    rF = float(np.sqrt(area_ratio_isentropic(MD, gas)))
    x_end = x_E + 2.3 * rF * np.sqrt(MD * MD - 1.0)
    t0 = time.time()
    res = inverse_design(ht, lambda x: float(law(max(x, x_A))), x_axis_end=x_end,
                         n_axis=n_axis, n_start=41, gamma=gas, M_start=1.05,
                         exit_mode="characteristic", x_E=x_E, M_d=MD,
                         start_line="throat_char", wall_mode="cplus",
                         axis_dx0=AXIS_DX0)
    t_moc = time.time() - t0
    w = res["wall"]
    qa = wall_qa(w, MD, x_E, gas, R=R)
    try:
        cfd_wall = AxisMachCFDWall(w, R=R, r_U=R_U, L_U=L_U, L_pipe=L_PIPE)
        cfd_msgs = cfd_wall.validate()
        wc = wall_curvature_diagnostics(cfd_wall, w[0, 0], w[-1, 0])
    except Exception as e:  # noqa: BLE001
        cfd_msgs = [f"CFDWall 構築失敗: {e}"]
        wc = {}
    ax_diag = axis_smoothness_diagnostics(law, x_A, x_E, breakpoints=breakpoints_of(tag, law))
    n_ax = n_axis - 1
    # seam_x_max: スロート特性線が初期値線として使われる throat_char 構成では、
    # x A_A 近傍 (実測 x<1.0 r_t) の三角形が「軸目標配列 / スロート特性線配列」の
    # 添字境界と重なり、健全な設計 (A) でも符号反転が局在する (plan §4 参照)。
    # 物理座標での seam 除外を index 除外と併用し、解像度に依らず同じ近傍を除く。
    topo = characteristic_topology_diagnostics(res["levels"], n_ax, n_crossing_sample=40,
                                               seam_i_pad=90, seam_k_max=110,
                                               seam_x_max=x_A + 1.0)
    margin = wall_margin_diagnostics(w)
    hard_gate = (not qa["violations"]) and (not cfd_msgs) \
        and topo["n_orientation_flip_interior"] == 0 \
        and topo["n_crossings_sampled"] == 0 and topo["pileup_pairs"] == 0 \
        and bool(np.all(np.diff(w[:, 1]) >= -1e-9))
    return {
        "tag": tag, "n_axis": n_axis, "t_moc_s": t_moc,
        "x_A": x_A, "x_E": x_E, "x_F": qa["x_F"],
        "wall_qa_violations": qa["violations"],
        "theta_max_deg": qa["theta_max_deg"], "kappa0_R_ratio": qa.get("kappa0_R_ratio"),
        "cfd_wall_violations": cfd_msgs,
        "axis": ax_diag, "topology": topo, "wall_curvature": wc, "margin": margin,
        "hard_gate_pass": bool(hard_gate),
        "wall_table": w.tolist() if n_axis == 1200 else None,   # 図用 (代表解像度のみ保存)
    }


def main():
    gas = GasSemiPerfect(Y, Tt=TT)
    ht = HallThroat(R=R, gamma=gas.gamma_star)
    xs_c, rr, MM, tt = ht.throat_characteristic(n=41)
    x_A = float(xs_c[0])
    M_A, Mp_A, Mpp_A = ht.axis_anchor(x_A)
    x_E = x_A + LC

    results = {}
    for n_axis in RESOLUTIONS:
        laws = build_laws(x_A, M_A, Mp_A, Mpp_A, gas)
        for tag, law in laws.items():
            key = f"{tag}@{n_axis}"
            print(f"=== {key} ===", flush=True)
            try:
                rec = run_one(tag, law, ht, gas, x_A, x_E, n_axis)
                results[key] = rec
                print(f"  hard_gate={rec['hard_gate_pass']}  theta_max={rec['theta_max_deg']:.2f}"
                     f"  flip_int={rec['topology']['n_orientation_flip_interior']}"
                     f"  max|Mppp|={rec['axis']['max_abs_Mppp']:.4f}"
                     f"  max|dk/ds|={rec['wall_curvature'].get('max_abs_dkappa_ds', float('nan')):.3f}"
                     f"  wall_qa={rec['wall_qa_violations']}  cfd_wall={rec['cfd_wall_violations']}",
                     flush=True)
            except Exception as e:  # noqa: BLE001
                traceback.print_exc()
                results[key] = {"tag": tag, "n_axis": n_axis, "error": str(e)}

    (CASE / "compare_axislaw_ABC.json").write_text(json.dumps(results, indent=1))
    print("saved", CASE / "compare_axislaw_ABC.json")

    # --- anchor/law メタ情報 (図・レポート用に別途保存) ---
    meta = {"x_A": x_A, "M_A": M_A, "Mp_A": Mp_A, "Mpp_A": Mpp_A, "x_E": x_E,
           "Tt": TT, "Md": MD, "R": R, "Lc": LC, "M_knot": MK,
           "r_U": R_U, "L_U": L_U, "L_pipe": L_PIPE}
    (CASE / "compare_axislaw_ABC_meta.json").write_text(json.dumps(meta, indent=1))
    print("DONE")


if __name__ == "__main__":
    main()
