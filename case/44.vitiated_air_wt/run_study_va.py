"""M4.19 加熱空気風洞 R/L_U/L_c スタディ (node Euler + TP semi-perfect, run_0001–0015) の投入ドライバ。

各点: 設計 (semi-perfect MOC) → メッシュ → 段階起動 node Euler TP CFD (soft→mid→本段 cfl2, 12000 step) → metrics。
実行: cd design && ./.venv-opt/bin/python ../case/44.vitiated_air_wt/run_study_va.py
出力: case/44.vitiated_air_wt/run_00NN_va_<tag>/ と study_cfd_va.json。
run 番号 (seq) は既存と衝突しないよう README の run 一覧を見て振ること。
"""
import json, subprocess, sys
from pathlib import Path
CASE = Path("/home/sano/work/forge/case/44.vitiated_air_wt"); DESIGN = Path("/home/sano/work/forge/design")
PTS = [("R1.5_LU6_Lc7", 1), ("R1.5_LU6_Lc8", 2), ("R1.5_LU6_Lc10", 3),
       ("R2_LU6_Lc7", 4), ("R2_LU6_Lc8", 5), ("R2_LU6_Lc10", 6),
       ("R3_LU6_Lc7", 7), ("R3_LU6_Lc8", 8), ("R3_LU6_Lc10", 9),
       ("R2_LU4_Lc8", 10), ("R2_LU9_Lc8", 11),
       ("R1.5_LU6_Lc9", 12), ("R2_LU6_Lc9", 13), ("R3_LU6_Lc9", 14), ("R3_LU6_Lc11", 15)]
if len(sys.argv) > 1:
    sel = set(sys.argv[1:]); PTS = [t for t in PTS if t[0] in sel or str(t[1]) in sel]
outp = CASE / "study_cfd_va.json"
out = json.loads(outp.read_text()) if outp.exists() else []
out = [r for r in out if r["tag"] not in {t for t, _ in PTS}]
for tag, seq in PTS:
    rd = CASE / f"run_{seq:04d}_va_{tag}"
    print(f"=== {tag} -> {rd.name} ===", flush=True)
    r = subprocess.run([str(DESIGN/".venv-opt/bin/python"), "-m", "forge_design.evaluate.runner_axismach",
                        str(CASE/f"problem_va_{tag}.yaml"), str(rd)], cwd=DESIGN, capture_output=True, text=True)
    (rd / "runner_stdout.log").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr) if rd.exists() else None
    if r.returncode != 0 or not (rd/"metrics.json").exists():
        print("  FAIL:", r.stderr[-300:], flush=True); out.append({"tag": tag, "run": rd.name, "ok": False}); continue
    m = json.loads((rd/"metrics.json").read_text()); info = json.loads((rd/"prepare_info.json").read_text())
    rec = {"tag": tag, "run": rd.name, "ok": True, "R": info["R"], "L_U": info["L_U"], "L_c": info["L_c"], "x_E": info["x_E"],
           "x_F": info["qa"]["x_F"], "theta_max": info["qa"]["theta_max_deg"], "kappa0R": info["qa"]["kappa0_R_ratio"],
           "L_total_m": (info["qa"]["x_F"] + info["L_U"] + 0.5) * info["scale_m"],
           "dM": m["dM_max_rel_Md"], "dM_rms": m["dM_rms"], "over": m["overshoot_rel"],
           "eps_M": m["exit_uniformity"].get("eps_M_rms"), "eps_th": m["exit_uniformity"].get("eps_theta_max_deg"),
           "M_exit_axis": m["M_axis_exit"], "core": f"{m['exit_uniformity'].get('n_core_faces')}/{m['exit_uniformity'].get('n_exit_faces')}"}
    out.append(rec)
    print(f"  dM={100*rec['dM']:.3f}%Md eps_M={100*rec['eps_M']:.3f}% eps_th={rec['eps_th']:.4f}° x_F={rec['x_F']:.2f} M_exit={rec['M_exit_axis']:.4f}", flush=True)
    outp.write_text(json.dumps(out, indent=1))
print("CFD STUDY VA DONE", flush=True)
