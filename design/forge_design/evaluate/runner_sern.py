"""⑤ SERN 評価 runner (S3): 逆設計 → 2 バンド構造メッシュ → forge 平面 2D (Euler/SST) → 力係数。

plan: plans/active/tooling-nozzle-sern-chain.md §4.2, §4.7。問題型 `sern_2d`。
1 run = 1 作動点。作動点は spec.external で与える (多作動点束ねは S6 で driver 側)。
段階起動: soft (1 次 + cfl 0.5, 3000 step) → 本段 (2 次 + cfl_main)。IC は領域別一様
(中間線より上 = 燃焼器出口状態、下 = 外部流)。
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

from ..geometry.moc_sern import PlanarMOC, SernKernelSpec, wall_forces
from ..geometry.rao_planar import ideal_gross_thrust
from ..meshing.mesh_sern import PHYS_SERN, SernMeshParams, generate_sern_mesh, write_msh41_named
from ..metrics.sern_forces import force_history, steadiness
from ..probdef import Problem, dv_value, load_problem

FORGE_TOOLS = Path("/home/sano/work/forge/solver_density_cuda/tools")
FORGE_BUILD = Path("/home/sano/work/forge/solver_density_cuda/build")
_ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")
MESH = "sern.h5"


def _dv(p: Problem, name, default=None) -> float:
    v = dv_value(p, name, default)
    return float(v["value"] if isinstance(v, dict) else v)


def select_operating_point(p: Problem, op: str | None) -> dict:
    """spec.operating_points[] から名前で 1 点を選び spec.external (と inflow の上書き) に反映する。
    op=None で operating_points が無ければ spec.external をそのまま使う。戻り値 = 選んだ点 (重み込み)。"""
    ops = p.spec.get("operating_points")
    if not ops:
        if op not in (None, "", "default"):
            raise ValueError("spec.operating_points が無いのに --op が指定された")
        return {"name": "default", "weight": 1.0, "external": dict(p.spec["external"])}
    names = [o["name"] for o in ops]
    if op in (None, "", "default"):
        op = names[0]
    if op not in names:
        raise ValueError(f"作動点 '{op}' が無い (候補: {names})")
    o = ops[names.index(op)]
    p.spec["external"] = dict(o["external"])
    if "inflow" in o:
        p.spec["inflow"] = {**p.spec["inflow"], **o["inflow"]}
    return {"name": op, "weight": float(o.get("weight", 1.0)), "external": dict(o["external"])}


def gas_states(p: Problem) -> dict:
    g, cp = p.gamma, p.cp
    R = cp * (g - 1.0) / g
    fi, ex = p.spec["inflow"], p.spec["external"]
    if fi.get("mode", "supersonic") != "supersonic":
        raise ValueError("inflow.mode は現状 supersonic のみ (sonic_throat は S1 接続が未実装)")

    def st(M, P, T):
        ro = P / (R * T)
        u = M * np.sqrt(g * R * T)
        k = 1.5 * (0.01 * u) ** 2
        omega = ro * k / (1.8e-5 * 10.0)
        return {"M": float(M), "P": float(P), "T": float(T), "ro": float(ro), "u": float(u),
                "k": float(k), "omega": float(omega)}
    return {"exhaust": st(fi["M_in"], fi["p_in"], fi["T_in"]),
            "ext": st(ex["M_inf"], ex["p_inf"], ex["T_inf"]), "R": R}


def design_from_problem(p: Problem):
    geo = p.geometry
    st = gas_states(p)
    p_ext_ratio = st["ext"]["P"] / st["exhaust"]["P"]
    spec = SernKernelSpec(M_in=float(p.spec["inflow"]["M_in"]),
                          theta_r0=np.deg2rad(_dv(p, "theta_r0_deg")), theta_c0=np.deg2rad(_dv(p, "theta_c0_deg")),
                          L_cowl=_dv(p, "L_cowl"), gamma=p.gamma, p_ext_over_p_in=p_ext_ratio,
                          x_max=float(geo.get("x_max_kernel", 10.0)), nj=int(geo.get("nj_moc", 301)),
                          dx=float(geo.get("dx_moc", 2e-3)))
    k = PlanarMOC(spec).march()
    if str(geo.get("mode", "keypoint")) == "straight":
        d = k.straight_design(_dv(p, "L_ramp"))
    else:
        d = k.design_ramp(M_c=_dv(p, "M_c"), f=_dv(p, "f"))
    xr, yr = p.spec.get("moment_ref", [0.0, 0.0])
    fr = wall_forces(d, spec.M_in, p.gamma, pa_over_pin=p_ext_ratio, x_ref=float(xr), y_ref=float(yr))
    theta_b = float(k.TH[-1, 0])   # 自由境界の終端角 (せん断層に格子線を沿わせる)
    return k, d, fr, theta_b


def _solver_config(p: Problem, nsteps: int, out_int: int, cfl: float, p_ref: float) -> str:
    disc = p.mesh.get("discretization", "cell")
    model = p.evaluate.get("model", "euler")
    node_keys = ", nodeWallDirichlet: 1" if (disc == "node" and model != "euler") else ""
    if model == "euler":
        phys = f"physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 0, ro: 1.2, visc: 0.0, thermCond: 0.0, cp: {p.cp}, gamma: {p.gamma}}}"
        turb = 'turbulence: {model: "none"}'
    else:
        phys = (f"physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 1, ro: 1.2, visc: 1.8e-5, thermCond: 0.0257, "
                f"thermCondMethod: 1, prandtlLam: 0.72, cp: {p.cp}, gamma: {p.gamma}}}")
        turb = 'turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1, wallTreatmentSST: 1}'
    return f"""mesh: {{meshFormat: "hdf5", discretization: "{disc}", isAxisymmetric: 0{node_keys}, meshFileName: "{MESH}", valueFileName: "{MESH}"}}
gpu: 1
solver: "SLAU"
{phys}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl},
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2, pRef: {p_ref}}}
{turb}
initial: "uniform_p101325_u10"
"""


def _bcond_config(p: Problem, st: dict) -> str:
    model = p.evaluate.get("model", "euler")
    wall_kind = "slip" if model == "euler" else "wall"
    ex, en = st["exhaust"], st["ext"]

    def inlet(name, pid, s):
        return (f"{name}: {{physID: {pid}, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , "
                f"floats: {{ro: {s['ro']:.6g}, Ux: {s['u']:.6g}, Uy: 0.0, Uz: 0.0, Ps: {s['P']:.6g}, k: {s['k']:.6g}, omega: {s['omega']:.6g}}}}}\n")

    def outlet(name, pid):
        return (f"{name}: {{physID: {pid}, kind: outlet_statPress, outputHDFflg: 0, ints: , "
                f"floats: {{Ps: {en['P']:.6g}, Pt: {en['P']:.6g}, Tt: {en['T']:.6g}}}}}\n")

    def wall(name, pid):
        return f"{name}: {{physID: {pid}, kind: {wall_kind}, outputHDFflg: 1, ints: , floats: }}\n"
    P = PHYS_SERN
    return (inlet("inlet_nozzle", P["inlet_nozzle"], ex) + inlet("inlet_ext", P["inlet_ext"], en)
            + outlet("outlet", P["outlet"]) + wall("ramp", P["ramp"]) + wall("cowl_in", P["cowl_in"])
            + wall("cowl_out", P["cowl_out"]) + outlet("bottom", P["bottom"]) + outlet("top_out", P["top_out"]))


def apply_wall_offset(design, wall_offset: dict, H: float):
    """壁を法線方向 (流体と反対側) に dn(x) [m] だけ動かした SernDesign を返す (無次元化して適用)。
    ramp: 上壁なので +n = 上、cowl: 下壁 (内面が上向き) なので +n = 下。表は (x_m, dn_m) の 2 列。"""
    import copy
    d = copy.deepcopy(design)
    for name, arr in (("ramp", d.ramp_xy), ("cowl", d.cowl_xy)):
        tbl = wall_offset.get(name)
        if tbl is None:
            continue
        tbl = np.asarray(tbl, dtype=float)
        dn = np.interp(arr[:, 0], tbl[:, 0] / H, tbl[:, 1] / H, left=tbl[0, 1] / H, right=tbl[-1, 1] / H)
        t = np.gradient(arr, axis=0); t /= np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-30)[:, None]
        nrm = np.column_stack([-t[:, 1], t[:, 0]]) if name == "ramp" else np.column_stack([t[:, 1], -t[:, 0]])
        arr += dn[:, None] * nrm
    return d


def paste_region_ic(h5path, y_mid, scale: float, st: dict, gamma: float) -> None:
    ex, en = st["exhaust"], st["ext"]
    with h5py.File(h5path, "r+") as f:
        cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
        xn, yn = cc[:, 0] / scale, cc[:, 1] / scale
        upper = yn > y_mid(xn)
        ro = np.where(upper, ex["ro"], en["ro"]); u = np.where(upper, ex["u"], en["u"]); P = np.where(upper, ex["P"], en["P"])
        v = f["/VALUE"]
        v["ro"][:] = ro; v["roUx"][:] = ro * u; v["roUy"][:] = 0.0; v["roUz"][:] = 0.0
        v["roe"][:] = P / (gamma - 1.0) + 0.5 * ro * u * u
        if "roK" in v:
            v["roK"][:] = ro * np.where(upper, ex["k"], en["k"]); v["roOmega"][:] = ro * np.where(upper, ex["omega"], en["omega"])


def prepare(problem_path, run_dir, nsteps=None, op: str | None = None, wall_offset=None) -> dict:
    """op: 作動点名 (spec.operating_points)。wall_offset: {"ramp": (x_m, dn_m), "cowl": (x_m, dn_m)} の
    法線オフセット表 [m] (S5 δ* 一発補正。壁を流体と反対側へ dn だけ動かす)。"""
    p = load_problem(problem_path)
    if p.type != "sern_2d":
        raise ValueError("runner_sern は sern_2d 専用")
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=False)
    opinfo = select_operating_point(p, op)
    st = gas_states(p)
    kern, design, fr_moc, theta_b = design_from_problem(p)
    H = float(p.spec["H_m"])
    m = p.mesh
    mp = SernMeshParams(ni_up=int(m.get("ni_up", 16)), ni_noz=int(m.get("ni_noz", 120)), ni_plume=int(m.get("ni_plume", 220)),
                        nj_top=int(m.get("nj_top", 101)), nj_bot=int(m.get("nj_bot", 61)), L_up=float(m.get("L_up", 0.5)),
                        x_out_extra=float(m.get("x_out_extra", 2.0)), bot_depth=float(m.get("bot_depth", 3.0)),
                        first_wall_frac=float(m.get("first_wall_frac", 2e-3)),
                        cowl_thickness=float(m.get("cowl_thickness", 2e-3 if m.get("discretization", "cell") == "node" else 0.0)),
                        interface_angle=float(m.get("interface_angle_rad", theta_b)),
                        top_ext_angle=float(np.deg2rad(m.get("top_ext_angle_deg", np.rad2deg(design.info["theta_e"])))),
                        scale=H)
    if wall_offset:
        design = apply_wall_offset(design, wall_offset, H)
    coords, quads, bedges, minfo, y_mid = generate_sern_mesh(design, mp)
    write_msh41_named(run_dir / "sern.msh", coords, quads, bedges, PHYS_SERN)
    np.savetxt(run_dir / "ramp_contour.csv", design.ramp_xy * H, delimiter=",", header="x_m,y_m", comments="")
    np.savetxt(run_dir / "cowl_contour.csv", design.cowl_xy * H, delimiter=",", header="x_m,y_m", comments="")
    n = int(nsteps or p.evaluate.get("nStepOuter", 6000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    cfl = float(p.evaluate.get("cfl_main", 4.0))
    cfg = _solver_config(p, n, out_int, cfl, st["ext"]["P"])
    (run_dir / "bcondConfig.yaml").write_text(_bcond_config(p, st))
    (run_dir / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    disc = p.mesh.get("discretization", "cell")
    # 品質ゲートは primal (cell) 変換で
    (run_dir / "solverConfig.yaml").write_text(cfg.replace(f'discretization: "{disc}"', 'discretization: "cell"').replace(", nodeWallDirichlet: 1", ""))
    subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "sern.msh", "sern_qc.h5"], cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"), "sern_qc.h5", "--mode", "2d"], cwd=run_dir, env=_ENV, capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(q.stdout + q.stderr)
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    if disc == "cell":
        (run_dir / "sern_qc.h5").rename(run_dir / MESH)
    else:
        (run_dir / "sern_qc.h5").unlink()
        (run_dir / "solverConfig.yaml").write_text(cfg)
        subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "sern.msh", MESH], cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    for f in run_dir.glob("sern_qc.xmf"):
        f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg)
    paste_region_ic(run_dir / MESH, y_mid, H, st, p.gamma)
    ex = st["exhaust"]
    F_ideal_nd, M_e_id = ideal_gross_thrust(ex["M"], st["ext"]["P"] / ex["P"], p.gamma)
    info = {"problem": str(problem_path), "run_dir": str(run_dir), "nsteps": n, "H_m": H, "states": st,
            "operating_point": opinfo, "wall_offset": bool(wall_offset),
            "design": {"key_point": list(design.key_point), "foot_a": list(design.foot_a), "lip_e": list(design.lip_e),
                       "L_ramp": design.L_ramp, "mass_fraction_check": design.mass_fraction_check,
                       "theta_e_deg": float(np.rad2deg(design.info["theta_e"])), "theta_b_deg": float(np.rad2deg(theta_b)),
                       "p_te_over_p_in": kern.p_te_over_p_in, "warnings": design.info["warnings"]},
            "moc_forces": fr_moc, "F_ideal_N_per_m": F_ideal_nd * ex["P"] * H, "M_e_ideal": M_e_id,
            "mesh": minfo, "discretization": disc, "model": p.evaluate.get("model", "euler")}
    (run_dir / "solverConfig_main.yaml").write_text(cfg)
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def restart_by_index(res_h5, mesh_h5) -> None:
    """同一メッシュの stage 間移植: VALUE を index でコピーする (座標最近傍の `interp_field.py` は使わない)。
    理由 (2026-09-04, case/46 run_0009): スリットカウルの上下壁ノードは座標が一致し、最近傍補間が双子を同じ元
    ノードに写す → 排気側の壁ノードが外部流の圧力を持ち 2 次で発散した (interp_field の全 134 station で誤写像を確認)。"""
    with h5py.File(res_h5, "r") as src, h5py.File(mesh_h5, "r+") as dst:
        n = len(dst["VALUE/ro"])
        for k in src["VALUE"]:
            if k in dst["VALUE"] and len(src["VALUE"][k]) == n:
                dst["VALUE"][k][:] = src["VALUE"][k][:]


def run_forge(run_dir) -> int:
    r = subprocess.run([str(FORGE_TOOLS / "run_case.sh"), str(Path(run_dir).resolve())], env=_ENV, capture_output=True, text=True)
    (Path(run_dir) / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode


def run_staged(run_dir, stages: str = "full", soft_steps: int = 3000) -> int:
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig_main.yaml").read_text() if (run_dir / "solverConfig_main.yaml").exists() \
        else (run_dir / "solverConfig.yaml").read_text()
    if stages == "main":  # soft 段の res_*.h5 が残っている状態から本段だけ回す
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        restart_by_index(res[-1], run_dir / MESH)
        for f in run_dir.glob("res_*"):
            f.unlink()
        (run_dir / "solverConfig.yaml").write_text(cfg_main)
        return run_forge(run_dir)
    if stages == "none":
        return run_forge(run_dir)
    soft = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", "cfl: 0.5, cfl_pseudo: 0.5", cfg_main)
    soft = soft.replace("convMethod: 1", "convMethod: 0")
    soft = re.sub(r"nStepOuter: \d+", f"nStepOuter: {soft_steps}", soft)
    soft = re.sub(r"outStepInterval: \d+", f"outStepInterval: {soft_steps}", soft)
    (run_dir / "solverConfig.yaml").write_text(soft)
    rc = run_forge(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if rc != 0 or not res:
        raise RuntimeError("soft 段が失敗 (res_nan_*.h5 / forge_run.log を見る)")
    restart_by_index(res[-1], run_dir / MESH)
    for f in run_dir.glob("res_*"):
        f.unlink()
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    return run_forge(run_dir)


def collect(problem_path, run_dir) -> dict:
    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    st = info["states"]; ex, en = st["exhaust"], st["ext"]; H = info["H_m"]
    xr, yr = p.spec.get("moment_ref", [0.0, 0.0])
    hist = force_history(run_dir, p_a=en["P"], F_ideal=info["F_ideal_N_per_m"], H=H, x_ref=float(xr) * H, y_ref=float(yr) * H,
                         mdot_u_in=ex["ro"] * ex["u"] ** 2 * H, p_in=ex["P"],
                         twall_on_fluid=(info.get("discretization", "cell") == "cell"))
    verdict = (run_dir / "CONVERGENCE_VERDICT.txt").read_text().strip().splitlines()[-2:] if (run_dir / "CONVERGENCE_VERDICT.txt").exists() else []
    out = {"convergence_verdict": verdict, "n_snapshots": len(hist), "history": hist,
           "operating_point": info.get("operating_point"), "L_ramp": info["design"]["L_ramp"]}
    if hist:
        last = hist[-1]
        out.update({k: last[k] for k in ("step", "C_T", "C_T_wall", "C_L", "C_M", "T_wall", "L", "M_noseup")})
        for k in ("C_T_with_shear", "C_T_friction", "sep_frac_ramp", "sep_x_min_ramp"):
            if k in last:
                out[k] = last[k]
        out["steadiness"] = {k: steadiness([h[k] for h in hist]) for k in ("C_T", "C_L", "C_M")}
        out["moc_forces"] = info["moc_forces"]
        out["cfd_vs_moc"] = {k: (last[k] - info["moc_forces"][k]) for k in ("C_T", "C_L", "C_M")}
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="forge_design ⑤ SERN 評価 (1 作動点)")
    ap.add_argument("problem"); ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--stages", default="full", choices=["full", "none", "main"])
    ap.add_argument("--op", default=None, help="作動点名 (spec.operating_points)")
    ap.add_argument("--wall-offset", default=None, help="δ* オフセット JSON ({ramp: [[x_m, dn_m],...], cowl: [...]})")
    a = ap.parse_args(argv)
    wo = json.loads(Path(a.wall_offset).read_text()) if a.wall_offset else None
    info = prepare(a.problem, a.run_dir, a.steps, op=a.op, wall_offset=wo)
    print(json.dumps({k: info[k] for k in ("design", "moc_forces", "mesh")}, indent=1))
    if a.prepare_only:
        return 0
    rc = run_staged(a.run_dir, a.stages)
    out = collect(a.problem, a.run_dir)
    print(json.dumps({k: v for k, v in out.items() if k != "history"}, indent=1))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
