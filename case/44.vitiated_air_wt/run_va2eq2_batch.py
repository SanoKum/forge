"""va2 ノズル (旧条件で設計した M_d 別ベスト形状) × 旧入口条件 × 平衡凝縮 EOS 拘束形 (condEquilibrium 2) の再計算 (run_0101–0106)。
M4.19 は参照 (新ノズル va3 での旧条件 eq2 は run_0098)。
cd design && ./.venv-opt/bin/python ../case/44.vitiated_air_wt/run_va2eq2_batch.py"""
import subprocess
from pathlib import Path
CASE=Path("/home/sano/work/forge/case/44.vitiated_air_wt"); DESIGN=Path("/home/sano/work/forge/design")
seq=101
for Md in ("4.3","4.4","4.5","4.75","5.0","4.19"):
    rd=CASE/f"run_{seq:04d}_va2eq2_M{Md}_eq2"; seq+=1
    if (rd/"metrics.json").exists(): print("skip", rd.name, flush=True); continue
    r=subprocess.run([str(DESIGN/".venv-opt/bin/python"),"-m","forge_design.evaluate.runner_axismach",str(CASE/f"problem_va2eq2_M{Md}_eq2.yaml"),str(rd)],cwd=DESIGN,capture_output=True,text=True)
    (CASE/f"{rd.name}.log").write_text(r.stdout+"\n---stderr---\n"+r.stderr)
    print(rd.name, "rc", r.returncode, flush=True)
print("VA2EQ2 BATCH DONE")
