"""新条件 (Pt 1.137 MPa / Tt 1058 K) × M_d 6 点 × {dry, noneq, eq} の投入ドライバ (run_0031–0048)。
cd design && ./.venv-opt/bin/python ../case/44.vitiated_air_wt/run_va2_batch.py"""
import subprocess, sys
from pathlib import Path
CASE=Path("/home/sano/work/forge/case/44.vitiated_air_wt"); DESIGN=Path("/home/sano/work/forge/design")
seq=31
for Md in ("4.19","4.3","4.4","4.5","4.75","5.0"):
    for kind in ("dry","noneq","eq"):
        rd=CASE/f"run_{seq:04d}_va2_M{Md}_{kind}"; seq+=1
        if (rd/"metrics.json").exists(): print("skip", rd.name); continue
        r=subprocess.run([str(DESIGN/".venv-opt/bin/python"),"-m","forge_design.evaluate.runner_axismach",str(CASE/f"problem_va2_M{Md}_{kind}.yaml"),str(rd)],cwd=DESIGN,capture_output=True,text=True)
        (CASE/f"{rd.name}.log").write_text(r.stdout+"\n---stderr---\n"+r.stderr)
        print(rd.name, "rc", r.returncode, flush=True)
print("VA2 BATCH DONE")
