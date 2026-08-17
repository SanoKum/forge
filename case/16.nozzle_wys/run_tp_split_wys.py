#!/usr/bin/env python3
"""Wyslouzil fig3 2D SST + H2O 凝縮 (Kw+HK) を **TP split (MIXDRY=N2 擬似種 + H2O 独立種)** で回す
検証 run の生成 (計画 plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md)。

基準: run_0050_fig3_2d_sst_kwhk (CPG thermalMethod 0, species [N2,H2O], Kw+HK)。同メッシュ・同 BC で
  - thermalMethod 2, species [MIXDRY, H2O] (speciesDBFile), thermoHrefTemp 298.15
  - IC = run_0048 (dry SST 発達場, CPG) を cpg_field_to_tp で TP 整合の roe に変換 (+roY)
  - condensation は run_0050 と同じ (condModel 1, condGasSpecies 1, Kantrowitz 1, HK)
usage: python3 run_tp_split_wys.py RUN_DIR [--dry]   (RUN_DIR は新規、例 run_0170_fig3_2d_sst_kwhk_tpsplit)
"""
import subprocess, sys, shutil, re
from pathlib import Path
sys.path.insert(0, "/home/sano/work/forge/design")
import yaml
from forge_design.gas.semiperfect import mixture_pseudo_species_split
from forge_design.evaluate.ic import cpg_field_to_tp

CASE = Path(__file__).resolve().parent
REF = CASE / "run_0050_fig3_2d_sst_kwhk"
DRY = CASE / "run_0048_fig3_2d_sst_dry"
TOOLS = Path("/home/sano/work/forge/solver_density_cuda/tools")
run_dir = CASE / sys.argv[1]
dry_only = "--dry" in sys.argv
run_dir.mkdir(exist_ok=False)

# 1. species DB (N2 だけの擬似種 = N2 と同係数、H2O は内蔵と同係数を明示)
db, Ys, order = mixture_pseudo_species_split({"N2": 0.98905, "H2O": 0.01095}, keep=("H2O",))
(run_dir / "species_db.yaml").write_text(yaml.safe_dump(db, sort_keys=False))
Y = [Ys[k] for k in order]

# 2. メッシュ + IC: run_0048 のメッシュ h5 に run_0048 最終 res を焼き込み → TP 変換
mesh_src = DRY / "nozzle_fig3_2d.h5"
res_src = sorted(DRY.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]
tmp = run_dir / "nozzle_fig3_2d_cpg.h5"
shutil.copy(mesh_src, tmp)
subprocess.run([sys.executable, str(CASE / "build_restart.py"), str(res_src), str(mesh_src), str(tmp), "--copyturb"], check=True)
info = cpg_field_to_tp(str(tmp), str(run_dir / "nozzle_fig3_2d.h5"), order, Y, cp_cpg=1039.0, gamma_cpg=1.4,
                       h_ref_T=298.15, db=db)
tmp.unlink()
print("IC:", info)

# 3. config: run_0050 をベースに physProp を TP split に、initial は無視されるが残す
cfg = (REF / "solverConfig.yaml").read_text()
cfg = cfg.replace("thermalMethod: 0", "thermalMethod: 2")
cfg = cfg.replace('species: ["N2", "H2O"]', 'species: ["MIXDRY", "H2O"]\n  speciesDBFile: "species_db.yaml"\n  thermoHrefTemp: 298.15')
if dry_only:
    cfg = re.sub(r"condensation: 1 ", "condensation: 0 ", cfg)
(run_dir / "solverConfig.yaml").write_text(cfg)
shutil.copy(REF / "bcondConfig.yaml", run_dir / "bcondConfig.yaml")
shutil.copy(REF / "probe.yaml", run_dir / "probe.yaml")
print("prepared", run_dir)
