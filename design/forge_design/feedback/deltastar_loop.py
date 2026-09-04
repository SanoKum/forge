r"""排除厚さの固定点反復ドライバ (固定 Euler 基準・コア整合抽出)。

計画: plans/active/tooling-nozzle-deltastar-core-matched-euler.md §4.5–4.6。

1 pass = 前 pass の NS run から $\delta_r(x)$ を抽出 (`metrics.deltastar.deltastar_from_core_matched_euler`)
→ 緩和して次の入力 $\delta_{in}^{k+1} = (1-\omega)\delta_{in}^k + \omega\,\delta_{use}^k$ (hard 不合格の断面は前回値保持)
→ 半径方向オフセットの物理壁で `prepare_ns` → 段階起動 NS → `collect` (質量流量帳簿込み)
→ 新 run からも抽出して固定点差を記録 (`fixed_point.json`)。

使い方 (リポジトリルートで):
  design/.venv-opt/bin/python -m forge_design.feedback.deltastar_loop \
      --problem case/45.isobutane_m6_d155/problem_d155_ns.yaml \
      --euler-ref case/45.isobutane_m6_d155/run_0001_euler_shortest_dry \
      --prev case/45.isobutane_m6_d155/run_0004_ns_v3 \
      --run-dir case/45.isobutane_m6_d155/run_0019_ns_cm_pass1 [--omega 0.5] [--prepare-only]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..metrics.deltastar import deltastar_from_core_matched_euler, massflow_ratio


def extract_and_merge(prev_run, euler_run, omega: float = 0.5, **kw) -> dict:
    """前 pass の NS run から抽出し、次 pass の入力 δ_r(x) を作る。
    出力 (prev_run 内): delta_r_equiv.csv / delta_r_equiv_diag.json / delta_r_next.csv。"""
    prev_run = Path(prev_run)
    d = deltastar_from_core_matched_euler(prev_run, euler_run, out_dir=prev_run, **kw)
    d_in = d["delta_in"]
    d_use = d["delta_r_use"]
    held = ~np.isfinite(d_use)
    d_next = np.where(held, d_in, (1.0 - omega) * d_in + omega * np.where(held, 0.0, d_use))
    np.savetxt(prev_run / "delta_r_next.csv", np.c_[d["x"], d_next, d_in, d_use, held.astype(int)],
               delimiter=",", comments="",
               header=f"x_rt,delta_r,delta_in_prev,delta_r_use,held (omega={omega})")
    fin = np.isfinite(d_use) & (d_in > 1e-4)
    ratio = d_use[fin] / d_in[fin]
    summary = {"omega": omega, "n_stations": int(len(d["x"])), "n_ok": int(d["ok"].sum()),
               "n_hard_ok": int(d["hard_ok"].sum()), "n_held": int(held.sum()),
               "use_over_in_median": float(np.median(ratio)) if fin.any() else None,
               "use_over_in_p10_p90": [float(np.percentile(ratio, 10)), float(np.percentile(ratio, 90))] if fin.any() else None,
               "max_abs_change": float(np.nanmax(np.abs(d_use - d_in))) if fin.any() else None,
               "massflow": d["massflow"],
               "delta_r_throat_use": float(np.interp(0.0, d["x"], np.where(np.isfinite(d_use), d_use, d_in))),
               "delta_r_throat_in": float(np.interp(0.0, d["x"], d_in))}
    (prev_run / "delta_r_extract_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def run_pass(problem, euler_ref, prev_run, run_dir, omega: float = 0.5, ic_from=None,
             prepare_only: bool = False, nsteps=None) -> dict:
    from ..evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    from ..evaluate.runner import FORGE_TOOLS
    prev_run = Path(prev_run); run_dir = Path(run_dir)
    summ = extract_and_merge(prev_run, euler_ref, omega=omega)
    print("extract(prev):", json.dumps({k: v for k, v in summ.items() if k != "massflow"}), flush=True)
    print("massflow(prev):", json.dumps(summ["massflow"]), flush=True)
    info = prepare_ns(problem, run_dir, nsteps=nsteps, ic_from=ic_from or prev_run,
                      delta_r_csv=prev_run / "delta_r_next.csv", offset="radial",
                      euler_ref=euler_ref, omega=omega, prev_run=prev_run)
    print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "mesh", "nStepOuter", "cfl_main")}, indent=1), flush=True)
    print(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
    if prepare_only:
        return info
    rc = run_staged_ns(run_dir)
    print("forge rc", rc, flush=True)
    m = collect(problem, run_dir)
    print(json.dumps(m, indent=1), flush=True)
    subprocess.run([sys.executable, str(FORGE_TOOLS / "check_convergence.py"), str(run_dir)], check=False)
    # 固定点差: 新 run から抽出し、入力 (= 新 run の壁 − Euler 壁) と比較
    fp = extract_and_merge(run_dir, euler_ref, omega=omega)
    (run_dir / "fixed_point.json").write_text(json.dumps(fp, indent=1))
    print("fixed_point(new):", json.dumps({k: v for k, v in fp.items() if k != "massflow"}), flush=True)
    print("massflow(new):", json.dumps(fp["massflow"]), flush=True)
    return m


def run_pass0_integral(problem, euler_ref, run_dir, ic_from, initializer=None, prepare_only: bool = False,
                       nsteps=None, omega: float = 0.5) -> dict:
    """pass 0: 積分法 (CONTUR) 初期壁で NS を立て、終了後に抽出して次 pass 用 delta_r_next.csv まで作る。
    initializer=None なら problem YAML の `deltastar_initializer` (無ければ断熱 contur)。"""
    from ..evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    from ..evaluate.runner import FORGE_TOOLS
    run_dir = Path(run_dir)
    init = initializer if initializer is not None else {"model": "contur", "thermal_bc": {"mode": "adiabatic"}}
    info = prepare_ns(problem, run_dir, nsteps=nsteps, ic_from=ic_from, initializer=init,
                      euler_ref=euler_ref, omega=omega)
    print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "initializer", "mesh", "nStepOuter", "cfl_main")},
                     indent=1, default=str), flush=True)
    print(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
    if prepare_only:
        return info
    rc = run_staged_ns(run_dir)
    print("forge rc", rc, flush=True)
    m = collect(problem, run_dir)
    print(json.dumps(m, indent=1), flush=True)
    subprocess.run([sys.executable, str(FORGE_TOOLS / "check_convergence.py"), str(run_dir)], check=False)
    fp = extract_and_merge(run_dir, euler_ref, omega=omega)
    (run_dir / "fixed_point.json").write_text(json.dumps(fp, indent=1))
    print("extract(new):", json.dumps({k: v for k, v in fp.items() if k != "massflow"}), flush=True)
    print("massflow(new):", json.dumps(fp["massflow"]), flush=True)
    return m


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="δ_r 固定点反復 (固定 Euler 基準)")
    ap.add_argument("--problem", required=True)
    ap.add_argument("--euler-ref", required=True)
    ap.add_argument("--prev", default=None, help="前 pass の NS run (抽出元・IC 既定)。--init-integral では不要")
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--omega", type=float, default=0.5)
    ap.add_argument("--ic-from", default=None)
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--extract-only", action="store_true", help="抽出と delta_r_next.csv だけ作る")
    ap.add_argument("--init-integral", action="store_true",
                    help="pass 0: 積分法 (CONTUR) 初期壁で NS を立てる (YAML の deltastar_initializer を使う)")
    ap.add_argument("--init-thermal", default=None, help="--init-integral の熱境界条件 JSON (例 '{\"mode\":\"adiabatic\"}')")
    a = ap.parse_args(argv)
    if a.extract_only:
        s = extract_and_merge(a.prev, a.euler_ref, omega=a.omega)
        print(json.dumps(s, indent=1))
        return 0
    if a.init_integral:
        init = None
        if a.init_thermal:
            init = {"model": "contur", "thermal_bc": json.loads(a.init_thermal)}
        run_pass0_integral(a.problem, a.euler_ref, a.run_dir, a.ic_from, initializer=init,
                           prepare_only=a.prepare_only, nsteps=a.steps, omega=a.omega)
        return 0
    if not a.prev:
        ap.error("--prev が必要 (--init-integral でなければ)")
    run_pass(a.problem, a.euler_ref, a.prev, a.run_dir, omega=a.omega, ic_from=a.ic_from,
             prepare_only=a.prepare_only, nsteps=a.steps)
    return 0


if __name__ == "__main__":
    sys.exit(main())
