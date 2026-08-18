"""va3 L_c8 形状 (新条件で設計したノズル) を **一個前の入口条件 (va2: Pt 1.137 MPa / Tt 1058 K / 旧組成)** で回す
off-design バッチ (run_0095–0097) + 新条件の平衡凝縮 (run_0094)。

形状・メッシュは `problem_va3_M4.19_Lc8_dry.yaml` の設計チェーンから作り (= run_0091 と同一壁)、
ガス (擬似種 DB)・BC (Pt/Tt/組成)・IC (等エントロピー) だけを `problem_va3old_Lc8_*.yaml` から差し替える。
`collect` は va3 の設計目標軸 M に対して評価するので、旧条件側の ‖ΔM‖∞ は「off-design のずれ」を意味する。

cd design && ./.venv-opt/bin/python ../case/44.vitiated_air_wt/run_va3old_batch.py
"""
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.evaluate import runner_axismach as R
from forge_design.evaluate.ic import paste_isentropic_ic
from forge_design.probdef import load_problem

CASE = Path("/home/sano/work/forge/case/44.vitiated_air_wt")
GEO_PROB = CASE / "problem_va3_M4.19_Lc8_dry.yaml"   # 形状 (壁・メッシュ) の供給元 = run_0091 と同一


def run_one(run_dir: Path, bc_prob: Path, geo_prob: Path = GEO_PROB):
    """geo_prob の形状で prepare し、bc_prob のガス/BC/IC に差し替えてから段階起動。"""
    p_geo = load_problem(geo_prob)
    info = R.prepare(geo_prob, run_dir)                      # 形状・メッシュ・(新条件の) config/IC
    if bc_prob != geo_prob:
        p_bc = load_problem(bc_prob)
        d = R.design_chain(p_geo)                            # 壁 (IC の面積比に使う)
        scale = float(p_geo.spec["r_throat"])
        n, out_int = info["nStepOuter"], int(p_geo.evaluate.get("outStepInterval", 4000))
        # ガス (擬似種 DB) + 凝縮設定 + BC を旧条件に差し替え
        (run_dir / "solverConfig.yaml").write_text(
            R._apply_gas_to_config(R._config_euler_node(p_bc, n, out_int, 4.0, 1), p_bc, run_dir))
        (run_dir / "bcondConfig.yaml").write_text(
            R._bcond_with_species(R._bcond(p_bc, euler=True), R._tp_species_Y(p_bc)))
        # IC も旧ガス・旧 Pt/Tt で貼り直す (面積比は va3 壁のまま)
        exit_ic = paste_isentropic_ic(run_dir / "nozzle.h5", d["wall"], scale,
                                      float(p_bc.spec["Pt"]), float(p_bc.spec["Tt"]),
                                      p_bc.gamma, p_bc.cp, gas=p_bc.gas_model,
                                      h_ref_T=float(p_bc.evaluate.get("thermo_href_temp", 298.15)),
                                      species_Y=R._tp_species_Y(p_bc))
        info["bc_problem"] = str(bc_prob)
        info["bc_gas"] = {"Pt": float(p_bc.spec["Pt"]), "Tt": float(p_bc.spec["Tt"]),
                          "Y": dict(p_bc.raw["gas"]["species"])}
        info["ic_exit"] = {k: float(v) for k, v in exit_ic.items()} if isinstance(exit_ic, dict) else None
        (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
        pp = p_bc
    else:
        pp = p_geo
    rc = R.run_staged(run_dir, cfl_main=pp.evaluate.get("cfl_main"),
                      mid_stage=bool(pp.evaluate.get("mid_stage", pp.is_semiperfect)))
    met = R.collect(geo_prob, run_dir)                       # 目標は常に va3 設計則
    return rc, met


JOBS = [
    ("run_0094_va3_M4.19_Lc8_eq",  CASE / "problem_va3_M4.19_Lc8_eq.yaml", GEO_PROB),
    ("run_0095_va3old_Lc8_dry",    CASE / "problem_va3old_Lc8_dry.yaml",   GEO_PROB),
    ("run_0096_va3old_Lc8_noneq",  CASE / "problem_va3old_Lc8_noneq.yaml", GEO_PROB),
    ("run_0097_va3old_Lc8_eq",     CASE / "problem_va3old_Lc8_eq.yaml",    GEO_PROB),
]

if __name__ == "__main__":
    for name, bc, geo in JOBS:
        rd = CASE / name
        if (rd / "metrics.json").exists():
            print("skip", name, flush=True)
            continue
        try:
            rc, met = run_one(rd, bc, geo)
            print(f"{name} rc {rc} dM {100*met['dM_max_rel_Md']:.3f}% M_exit {met['M_axis_exit']:.4f}", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"{name} FAILED: {e}", flush=True)
    print("VA3OLD BATCH DONE")
