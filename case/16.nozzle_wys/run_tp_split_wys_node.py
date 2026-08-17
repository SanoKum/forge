#!/usr/bin/env python3
"""Wyslouzil fig3 2D SST + H2O 凝縮 (Kw+HK) を **node (median-dual) + TP split** で回す run の生成
(計画 plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md、ユーザ指示「全部 node で」2026-08-17)。

- メッシュ: 平面 2D `mesh/nozzle_fig3_2d_planar.msh` (押し出しなし。node は 2 ノード spanwise 不可) を
  **node config を書いてから** convertGmshToForge (median-dual + no-slip 壁の wall_dist)。
- IC: run_0048 (cell, dry SST 発達場, CPG) を interp_field で cross-mesh 移植 → cpg_field_to_tp で TP 整合の roe に。
- config: node + SST 低 Re (nodeWallDirichlet 1, wallTreatmentSST 0), thermalMethod 2 [MIXDRY,H2O],
  thermoHrefTemp 298.15, 凝縮は run_0050 と同じ (condModel 1, condGasSpecies 1, Kantrowitz 1, HK)。
usage: python3 run_tp_split_wys_node.py RUN_DIR [--dry] [--cfl C] [--first-order] [--nsteps N]
"""
import subprocess, sys, shutil, os
import numpy as np
from pathlib import Path
sys.path.insert(0, "/home/sano/work/forge/design")
import yaml
from forge_design.gas.semiperfect import mixture_pseudo_species_split
from forge_design.evaluate.ic import cpg_field_to_tp, zero_wall_velocity_ic

CASE = Path(__file__).resolve().parent
REF = CASE / "run_0050_fig3_2d_sst_kwhk"
DRY = CASE / "run_0048_fig3_2d_sst_dry"
FORGE = Path("/home/sano/work/forge/solver_density_cuda")
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
args = sys.argv[1:]
run_dir = CASE / args[0]
dry_only = "--dry" in args
first_order = "--first-order" in args
cfl = float(args[args.index("--cfl") + 1]) if "--cfl" in args else 1.0
cpg = "--cpg" in args
laminar = "--laminar" in args     # 切り分け: SST を切る
slip = "--slip" in args           # 切り分け: 壁を slip (非粘性壁) に
nodir = "--nodirichlet" in args   # 切り分け: nodeWallDirichlet 0
ic_mode = args[args.index("--ic") + 1] if "--ic" in args else "isen"   # isen (1D 等エントロピー, 既定) / xmesh (cell 場の最近傍移植)
staged = "--staged" in args
mesh_name = args[args.index("--mesh") + 1] if "--mesh" in args else "planar"   # planar (壁クラスタ, 粘性用) / planar_inv (一様, 非粘性用)
cell_mode = "--cell" in args      # 対照: cell 離散化 (同一平面メッシュ)
wallfunc = "--wallfunc" in args    # node SST 壁関数 (wallTreatmentSST 1) — 緩い壁クラスタ mesh 用      # soft(1次 cfl0.5) → mid(cfl1) → 本段 の段階起動 (runner_axismach.run_staged 流用)          # 切り分け用: CPG [N2,H2O] (run_0050 と同じ熱力学) の node 版
nsteps = int(args[args.index("--nsteps") + 1]) if "--nsteps" in args else 15000
run_dir.mkdir(exist_ok=False)

db, Ys, order = mixture_pseudo_species_split({"N2": 0.98905, "H2O": 0.01095}, keep=("H2O",))
(run_dir / "species_db.yaml").write_text(yaml.safe_dump(db, sort_keys=False))
Y = [Ys[k] for k in order]

cond = "" if dry_only else """condensation:
  condensation: 1
  nCondSpecies: 1
  condModel: 1
  condGasSpecies: 1
  condKantrowitz: 1
  condGrowthModel: 0
"""
cfg = f"""mesh:
  meshFormat: "hdf5"
  discretization: "{'cell' if cell_mode else 'node'}"
  isAxisymmetric: 0
  nodeWallDirichlet: {0 if nodir else 1}
  meshFileName: "nozzle_fig3_2d_node.h5"
  valueFileName: "nozzle_fig3_2d_node.h5"
gpu: 1
solver: "SLAU"
physProp:
  isCompressible: 1
  thermalMethod: {0 if cpg else 2}
  viscMethod: {0 if slip else 1}
  ro: 1.2
  visc: 0.0000175
  thermCond: 0.0241
  cp: 1039.0
  gamma: 1.4
{'  species: ["N2", "H2O"]' if cpg else '  species: ["MIXDRY", "H2O"]' + chr(10) + '  speciesDBFile: "species_db.yaml"' + chr(10) + '  thermoHrefTemp: 298.15'}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 0.00001, cfl: {cfl}, cfl_pseudo: {2*cfl}, implicitRelax: 1.0, blockDPLUR: 1,
           dt_min: 0.00000001, dt_max: 1.0, detectNaN: 1}}
  outStepInterval: 3000
  outStepStart: 0
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: {0 if first_order else 1}, limiter: {0 if first_order else 2}}}
turbulence: {{model: "{'none' if laminar else 'sst'}", scalarDiffusion: 1, dilatationCorrection: 2, wallTreatmentSST: {1 if wallfunc else 0}, kInf: 1.0, omegaInf: 1000.0}}
initial: "nozzle_wys"
{cond}"""
(run_dir / "solverConfig.yaml").write_text(cfg)
shutil.copy(REF / "bcondConfig.yaml", run_dir / "bcondConfig.yaml")
shutil.copy(REF / "probe.yaml", run_dir / "probe.yaml")
# frontback は平面 2D では存在しない → bcond から落とす
b = (run_dir / "bcondConfig.yaml").read_text()
if "frontback:" in b:
    b = b.split("frontback:")[0]
if slip:
    b = b.replace("kind: wall ", "kind: slip ")
# 出口の逆流基準 Pt を Ps に落とす (runner_axismach と同じ)。node は壁コーナーノード (u=0 Dirichlet) が
# 出口面で逆流判定になり、Pt=59 kPa のよどみガスを注入 → 衝撃列が遡上して unstart した (run_0179/0180:
# 15000 step で場が P≈Pt に戻る)。cell の run_0050 では顕在化しなかった node 固有の経路。
import re as _re
b = _re.sub(r"(outlet:.*?floats:\s*\n(?:\s+.*\n)*?)\s+Pt: 59070\.0", lambda m: m.group(1) + "    Pt: 2000.0", b, flags=_re.S)
(run_dir / "bcondConfig.yaml").write_text(b)

# node 変換 (config を先に書いてから)
shutil.copy(CASE / "mesh" / f"nozzle_fig3_2d_{mesh_name}.msh", run_dir / "nozzle_fig3_2d_planar.msh")
subprocess.run([str(FORGE / "build" / "convertGmshToForge"), "nozzle_fig3_2d_planar.msh", "nozzle_fig3_2d_node.h5"],
               cwd=run_dir, env=ENV, check=True, capture_output=True, text=True)
if ic_mode == "isen":
    # 1D 等エントロピー IC (壁 y(x) から A/A*)。cell 場の最近傍移植は壁ノードで隣接 2 セル (壁セル u≈0,T≈Tt と
    # 第 2 セル u≈400,T≈205K) が交互に乗り T が 200↔286K でジグザグし、node が壁から発散した (run_0172–0177)。
    import h5py
    from forge_design.evaluate.ic import paste_isentropic_ic
    from forge_design.gas.semiperfect import GasSemiPerfect
    with h5py.File(run_dir / "nozzle_fig3_2d_node.h5", "r") as f:
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
    # 上壁 y(x): x ビンごとの y 最大 (平面 2D 全高メッシュ、スロート x=0, 半高 0.25 cm)
    xb = np.linspace(cc[:, 0].min(), cc[:, 0].max(), 400)
    yb = np.array([cc[(cc[:, 0] >= a) & (cc[:, 0] < b), 1].max() if np.any((cc[:, 0] >= a) & (cc[:, 0] < b)) else np.nan
                   for a, b in zip(xb[:-1], xb[1:])])
    xc = 0.5 * (xb[:-1] + xb[1:]); ok = np.isfinite(yb)

    class _Wall:
        x_throat, r_throat = 0.0, 0.0025
        def r(self, x, deriv=0):
            return np.interp(x, xc[ok], yb[ok])
    gas_ic = None if cpg else GasSemiPerfect({"N2": 0.98905, "H2O": 0.01095}, Tt=286.65)
    paste_isentropic_ic(str(run_dir / "nozzle_fig3_2d_node.h5"), _Wall(), 1.0, 59070.0, 286.65, 1.4, 1039.0,
                        k_init=1.0, omega_init=1000.0, gas=gas_ic, h_ref_T=(None if cpg else 298.15),
                        species_Y=Y)
    info = "1D isentropic IC (%s)" % ("CPG" if cpg else "TP semi-perfect, datum 298.15")
    print("IC:", info)
    print("prepared", run_dir)
    if staged:
        from forge_design.evaluate.runner_axismach import run_staged
        rc = run_staged(run_dir, cfl_main=cfl, mid_stage=True, mesh_h5="nozzle_fig3_2d_node.h5")
        print("run_staged rc", rc)
    sys.exit(0)
# cross-mesh IC: run_0048 最終 res (cell) → node h5 (nearest neighbor)
res_src = sorted(DRY.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]
subprocess.run([sys.executable, str(FORGE / "tools" / "interp_field.py"), str(res_src),
                str(run_dir / "nozzle_fig3_2d_node.h5")], env=ENV, check=True)
if cpg:
    import h5py
    with h5py.File(run_dir / "nozzle_fig3_2d_node.h5", "r+") as f:      # roY だけ貼る (CPG は roe そのまま)
        ro = f["/VALUE/ro"][:]
        for i, y in enumerate(Y):
            f.create_dataset(f"/VALUE/roY{i}", data=(ro * y).astype(np.float32))
    info = "CPG (roe unchanged, roY pasted)"
else:
    tmp = run_dir / "_cpg.h5"
    shutil.move(run_dir / "nozzle_fig3_2d_node.h5", tmp)
    info = cpg_field_to_tp(str(tmp), str(run_dir / "nozzle_fig3_2d_node.h5"), order, Y, cp_cpg=1039.0,
                           gamma_cpg=1.4, h_ref_T=298.15, db=db)
    tmp.unlink()
nw = zero_wall_velocity_ic(str(run_dir / "nozzle_fig3_2d_node.h5"))
print("IC:", info, "| wall nodes zeroed:", nw)
print("prepared", run_dir)
