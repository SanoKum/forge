"""M5 イソブタン風洞 R×L_U×L_c スタディ (③, design plan: design-isobutane-wt-m5-sweep) の投入ドライバ。

27 点 = R {1.5, 2, 3} × L_U {4, 6, 9} × L_c {14, 17, 20} (knot MK2.5, dry split_h2o TP)。
problem YAML はセンター (problem_ib_m5_R2_LU6_Lc17.yaml) を基底に R/L_U/L_c と name を置換生成。
センター点はインテーク smoke の run_0072 をそのまま流用 (再実行しない)。
実行: cd design && ./.venv-opt/bin/python ../case/42.isobutane_wt/run_study_m5_tp.py
出力: run_0073–0098 と study_cfd_m5_tp.json。
"""
import itertools, json, re, subprocess
from pathlib import Path

CASE = Path("/home/sano/work/forge/case/42.isobutane_wt"); DESIGN = Path("/home/sano/work/forge/design")
BASE = (CASE / "problem_ib_m5_R2_LU6_Lc17.yaml").read_text()

pts, seq = [], 73
for R, LU, Lc in itertools.product([1.5, 2, 3], [4, 6, 9], [14, 17, 20]):
    tag = f"R{R:g}_LU{LU}_Lc{Lc}"
    if (R, LU, Lc) == (2, 6, 17):
        pts.append((tag, 72, None))  # センター = smoke run_0072 流用
        continue
    pts.append((tag, seq, (R, LU, Lc))); seq += 1

out = []
for tag, n, combo in pts:
    rd = CASE / f"run_{n:04d}_ib_m5_{tag}_tp"
    if combo is not None:
        rd = CASE / f"run_{n:04d}_ib_m5_{tag}_tp"
        prob = CASE / f"problem_ib_m5_{tag}.yaml"
        R, LU, Lc = combo
        txt = BASE
        txt = re.sub(r"name: .*", f"name: isobutane_wt_m5_{tag}_tp", txt, count=1)
        txt = re.sub(r"  R: .*", f"  R: {R}", txt, count=1)
        txt = re.sub(r"  L_U: .*", f"  L_U: {LU}.0", txt, count=1)
        txt = re.sub(r"L_c: \{value: [\d.]+", f"L_c: {{value: {Lc}.0", txt, count=1)
        prob.write_text(txt)
    else:
        rd = CASE / "run_0072_ib_m5_R2_LU6_Lc17_tp"
        prob = CASE / "problem_ib_m5_R2_LU6_Lc17.yaml"
    print(f"=== {tag} -> {rd.name} ===", flush=True)
    if not (rd / "metrics.json").exists():
        r = subprocess.run([str(DESIGN/".venv-opt/bin/python"), "-m", "forge_design.evaluate.runner_axismach",
                            str(prob), str(rd)], cwd=DESIGN, capture_output=True, text=True)
        if r.returncode != 0 or not (rd/"metrics.json").exists():
            print("  FAIL:", r.stderr[-300:], flush=True)
            out.append({"tag": tag, "run": rd.name, "ok": False, "err": r.stderr[-300:]}); continue
    m = json.loads((rd/"metrics.json").read_text()); info = json.loads((rd/"prepare_info.json").read_text())
    rt = 0.080505
    rec = {"tag": tag, "run": rd.name, "ok": True, "R": info["R"], "L_U": info.get("L_U"), "L_c": info["L_c"],
           "x_F": info["qa"]["x_F"], "len_sup_m": info["qa"]["x_F"]*rt, "theta_max": info["qa"]["theta_max_deg"],
           "dM": m["dM_max_rel_Md"], "dM_rms": m["dM_rms"], "over": m["overshoot_rel"],
           "eps_M": m["exit_uniformity"].get("eps_M_rms"), "eps_th": m["exit_uniformity"].get("eps_theta_max_deg"),
           "M_exit_axis": m["M_axis_exit"], "core": f"{m['exit_uniformity'].get('n_core_faces')}/{m['exit_uniformity'].get('n_exit_faces')}"}
    out.append(rec)
    print(f"  dM={100*rec['dM']:.3f}%Md eps_M={100*rec['eps_M']:.3f}% eps_th={rec['eps_th']:.4f}° "
          f"len={rec['len_sup_m']:.3f}m M_exit={rec['M_exit_axis']:.4f}", flush=True)
    (CASE/"study_cfd_m5_tp.json").write_text(json.dumps(out, indent=1))
print("CFD STUDY M5 TP DONE", flush=True)
