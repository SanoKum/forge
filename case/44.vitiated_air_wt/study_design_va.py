"""M4.19 加熱空気風洞の設計のみスイープ (design_chain、CFD 無し) → study_design_va.json。

R∈{1.5,2,3} × L_c∈{7,8,9,10,11} @L_U=6、L_U∈{4,9} @R2 (L_c 7/8/10)。
基底 YAML problem_va_R2_LU6_Lc8.yaml の geometry/dv を上書きして design_chain を呼ぶ。
実行: cd design && ./.venv-opt/bin/python ../case/44.vitiated_air_wt/study_design_va.py
"""
import copy, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.probdef import load_problem, Problem
from forge_design.evaluate.runner_axismach import design_chain

CASE = Path("/home/sano/work/forge/case/44.vitiated_air_wt")
base = load_problem(CASE / "problem_va_R2_LU6_Lc8.yaml")
PTS = [(R, 6.0, Lc) for R in (1.5, 2.0, 3.0) for Lc in (7.0, 8.0, 9.0, 10.0, 11.0)]
PTS += [(2.0, LU, Lc) for LU in (4.0, 9.0) for Lc in (7.0, 8.0, 10.0)]
out = []
for R, LU, Lc in PTS:
    raw = copy.deepcopy(base.raw)
    raw["geometry"]["R"] = R; raw["geometry"]["L_U"] = LU
    raw["dv"]["L_c"]["value"] = Lc; raw["dv"]["L_c"]["max"] = 20.0
    p = Problem(name=f"va_R{R:g}_LU{LU:g}_Lc{Lc:g}", type=base.type, gamma=base.gamma, cp=base.cp,
                spec=raw["spec"], dv=raw["dv"], geometry=raw["geometry"], mesh=raw["mesh"],
                evaluate=raw["evaluate"], raw=raw)
    rt = float(p.spec["r_throat"]); rin = float(raw["geometry"]["r_inlet"]); Lp = float(raw["geometry"]["L_pipe"])
    rec = {"tag": f"R{R:g}_LU{LU:g}_Lc{Lc:g}", "R": R, "L_U": LU, "L_c": Lc,
           "mu_up": LU * LU / (R * (rin - 1.0))}
    try:
        d = design_chain(p)
        qa = d["qa"]
        rec.update(ok=True, x_A=d["x_A"], x_E=d["x_E"], x_F=qa["x_F"], theta_max=qa["theta_max_deg"],
                   kappa0R=qa["kappa0_R_ratio"], rF_err=qa["r_F_err_rel"], theta_exit=qa["theta_exit_deg"],
                   Lc_window=list(d["Lc_window"]), L_div_m=qa["x_F"] * rt,
                   L_total_m=(qa["x_F"] + LU + Lp) * rt, gates=d["gates"])
        print(f"{rec['tag']:16s} OK  th_max {qa['theta_max_deg']:6.2f}  k0R {qa['kappa0_R_ratio']:.3f}  "
              f"x_F {qa['x_F']:6.2f}  L_tot {rec['L_total_m']:.3f} m  window ({d['Lc_window'][0]:.2f},{d['Lc_window'][1]:.2f})", flush=True)
    except Exception as e:  # noqa: BLE001
        rec.update(ok=False, reason=str(e)[:200])
        print(f"{rec['tag']:16s} REJECT {str(e)[:150]}", flush=True)
    out.append(rec)
(CASE / "study_design_va.json").write_text(json.dumps(out, indent=1, default=float))
print("DESIGN SWEEP DONE", len(out))
