#!/usr/bin/env python3
"""同一メッシュ継続 run の生成: SRC_RUN の最終 res を interp_field で新 run の mesh h5 に貼り、nStepOuter=N で本段のみ実行。
usage: python3 continue_run.py SRC_RUN DST_RUN N"""
import sys, shutil, re, subprocess
from pathlib import Path
sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.evaluate.runner import run_forge, FORGE_TOOLS, _ENV
CASE = Path(__file__).resolve().parent
src, dst, n = CASE / sys.argv[1], CASE / sys.argv[2], int(sys.argv[3])
dst.mkdir(exist_ok=False)
for f in ["solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml", "probe.yaml", "nozzle_user_2d.h5"]:
    shutil.copy(src / f, dst / f)
res = sorted(src.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]
subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"), str(res), str(dst / "nozzle_user_2d.h5")], env=_ENV, check=True, capture_output=True, text=True)
c = (dst / "solverConfig.yaml").read_text()
c = re.sub(r"nStepOuter: \d+", f"nStepOuter: {n}", c); c = re.sub(r"outStepInterval: \d+", f"outStepInterval: {n//4}", c)
(dst / "solverConfig.yaml").write_text(c)
(dst / "CONTINUED_FROM.txt").write_text(f"{src.name}/{res.name}\n")
print("main rc", run_forge(dst))
print((dst / "CONVERGENCE_VERDICT.txt").read_text().strip().splitlines()[-1])
