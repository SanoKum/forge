"""va2lp の凝縮 12 run (noneq/eq) を 最終 binary (SLAU 二相面温度修正 + CEA 由来 L) で再計算 (run_0079–0090, 正本)。dry は g=0 でビット不変のため再計算しない。"""
import subprocess
from pathlib import Path
CASE=Path("/home/sano/work/forge/case/44.vitiated_air_wt"); DESIGN=Path("/home/sano/work/forge/design")
seq=79
for Md in ("4.19","4.3","4.4","4.5","4.75","5.0"):
    for kind in ("noneq","eq"):
        rd=CASE/f"run_{seq:04d}_va2c_M{Md}_{kind}"; seq+=1
        r=subprocess.run([str(DESIGN/".venv-opt/bin/python"),"-m","forge_design.evaluate.runner_axismach",str(CASE/f"problem_va2lp_M{Md}_{kind}.yaml"),str(rd)],cwd=DESIGN,capture_output=True,text=True)
        (CASE/"_logs"/f"{rd.name}.log").write_text(r.stdout+"\n---stderr---\n"+r.stderr); print(rd.name,"rc",r.returncode,flush=True)
print("VA2C DONE")
