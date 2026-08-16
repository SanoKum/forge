"""M4.2 イソブタン風洞 R/L スタディ (TP 正本, run_0031–0038) の投入ドライバ (再生成用)。

各点: 設計 (semi-perfect MOC) → メッシュ → 段階起動 node Euler TP CFD → metrics。
実行: cd design && ./.venv-opt/bin/python ../case/42.isobutane_wt/run_study_m42_tp.py
出力: case/42.isobutane_wt/run_00NN_ib_<tag>_tp/ と study_cfd_m42_tp.json。
run 番号 (seq) は既存と衝突しないよう README の run 一覧を見て振り直すこと。
"""
import json, subprocess, sys
from pathlib import Path
CASE = Path("/home/sano/work/forge/case/42.isobutane_wt"); DESIGN = Path("/home/sano/work/forge/design")
PTS = [("R1.5_LU6_Lcmax", 31), ("R2_LU6_Lcmax", 32), ("R3_LU6_Lcmax", 33),
       ("R1.5_LU6_Lc8", 34), ("R2_LU6_Lc8", 35), ("R3_LU6_Lc8", 36),
       ("R2_LU4_Lcmax", 37), ("R2_LU9_Lcmax", 38)]
out = []
for tag, seq in PTS:
    rd = CASE / f"run_{seq:04d}_ib_{tag}_tp"
    print(f"=== {tag} -> {rd.name} ===", flush=True)
    r = subprocess.run([str(DESIGN/".venv-opt/bin/python"), "-m", "forge_design.evaluate.runner_axismach",
                        str(CASE/f"problem_ib_{tag}_tp.yaml"), str(rd)], cwd=DESIGN, capture_output=True, text=True)
    if r.returncode != 0 or not (rd/"metrics.json").exists():
        print("  FAIL:", r.stderr[-300:], flush=True); out.append({"tag": tag, "run": rd.name, "ok": False}); continue
    m = json.loads((rd/"metrics.json").read_text()); info = json.loads((rd/"prepare_info.json").read_text())
    rec = {"tag": tag, "run": rd.name, "ok": True, "R": info["R"], "L_c": info["L_c"], "x_E": info["x_E"],
           "x_F": info["qa"]["x_F"], "theta_max": info["qa"]["theta_max_deg"],
           "dM": m["dM_max_rel_Md"], "dM_rms": m["dM_rms"], "over": m["overshoot_rel"],
           "eps_M": m["exit_uniformity"].get("eps_M_rms"), "eps_th": m["exit_uniformity"].get("eps_theta_max_deg"),
           "M_exit_axis": m["M_axis_exit"], "core": f"{m['exit_uniformity'].get('n_core_faces')}/{m['exit_uniformity'].get('n_exit_faces')}"}
    out.append(rec)
    print(f"  dM={100*rec['dM']:.3f}%Md eps_M={100*rec['eps_M']:.3f}% eps_th={rec['eps_th']:.4f}° x_F={rec['x_F']:.2f} M_exit={rec['M_exit_axis']:.4f}", flush=True)
(CASE/"study_cfd_m42_tp.json").write_text(json.dumps(out, indent=1))
print("CFD STUDY TP DONE", flush=True)
