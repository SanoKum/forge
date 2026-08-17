import subprocess
from pathlib import Path
CASE=Path("/home/sano/work/forge/case/44.vitiated_air_wt"); DESIGN=Path("/home/sano/work/forge/design")
for seq,Md in ((33,"4.19"),(36,"4.3"),(39,"4.4"),(42,"4.5"),(45,"4.75"),(48,"5.0")):
    rd=CASE/f"run_{seq:04d}_va2_M{Md}_eq"
    r=subprocess.run([str(DESIGN/".venv-opt/bin/python"),"-m","forge_design.evaluate.runner_axismach",str(CASE/f"problem_va2_M{Md}_eq.yaml"),str(rd)],cwd=DESIGN,capture_output=True,text=True)
    (CASE/"_logs"/f"{rd.name}.log").write_text(r.stdout+"\n---stderr---\n"+r.stderr); print(rd.name,"rc",r.returncode,flush=True)
print("EQ REDO DONE")
