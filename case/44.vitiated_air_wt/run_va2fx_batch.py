"""va2lp の凝縮 12 run (noneq/eq) を SLAU 二相面温度修正後 binary で再計算 (run_0067–0078)。dry は g=0 でビット不変のため再計算しない。"""
import subprocess
from pathlib import Path
CASE=Path("/home/sano/work/forge/case/44.vitiated_air_wt"); DESIGN=Path("/home/sano/work/forge/design")
seq=67
for Md in ("4.19","4.3","4.4","4.5","4.75","5.0"):
    for kind in ("noneq","eq"):
        rd=CASE/f"run_{seq:04d}_va2fx_M{Md}_{kind}"; seq+=1
        r=subprocess.run([str(DESIGN/".venv-opt/bin/python"),"-m","forge_design.evaluate.runner_axismach",str(CASE/f"problem_va2lp_M{Md}_{kind}.yaml"),str(rd)],cwd=DESIGN,capture_output=True,text=True)
        (CASE/"_logs"/f"{rd.name}.log").write_text(r.stdout+"\n---stderr---\n"+r.stderr); print(rd.name,"rc",r.returncode,flush=True)
print("VA2FX DONE")
