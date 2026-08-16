"""軸 Mach 則 A (`KnotQuinticAxisLaw`) の knot Mach 感度を逆 MOC 後まで評価する。

基準条件は `compare_axislaw_D.py` と同じ M6 / R=3。現行 A は `x_K` を直接指定せず、
`M_K` と `L_c` から `L_1=x_K-x_A` を従属決定するため、`M_K` 掃引が現行設計族における
knot 位置感度になる。

手順:
  1. L_c=30,35,40,45,50,60、M_K=1.5..4.0 (0.1 刻み) を n_axis=600 で全評価。
  2. 各 L_c について、基準 M_K=2.5、最小 theta_max、最大 mu-theta、
     margin>=1 deg かつ topology margin>=0.02 で最小 theta_max の候補を n_axis=1200 で再評価。
  3. 1200 点の代表群から基準と robust 最小壁角を選び n_axis=2400 で再評価。
  4. 全結果 JSON、表形式 CSV、感度図 PNG を保存する。

再生成:
  design/.venv-opt/bin/python case/42.isobutane_wt/sweep_axislaw_A_MK.py
"""
from __future__ import annotations

import csv
import importlib.util
import json
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CASE = Path(__file__).resolve().parent
_CMP_PATH = CASE / "compare_axislaw_D.py"
_SPEC = importlib.util.spec_from_file_location("compare_axislaw_D_for_MK_sweep", _CMP_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"比較モジュールを読み込めない: {_CMP_PATH}")
_CMP = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_CMP)

LCS = (30.0, 35.0, 40.0, 45.0, 50.0, 60.0)
MKS = np.round(np.arange(1.5, 4.0 + 1e-9, 0.1), 10)
BASELINE_MK = 2.5
TOPO_FLOOR = 0.02
MARGIN_FLOOR_DEG = 1.0
DEFAULT_WORKERS = {600: 6, 1200: 4, 2400: 2}


def _compact(rec: dict) -> dict:
    """巨大な壁点列を除き、比較に必要な診断だけを残す。"""
    rec.pop("wall_table", None)
    return rec


def evaluate(L_c: float, M_K: float, n_axis: int) -> dict:
    """A を構築して逆 MOC と全診断を実行する。"""
    try:
        law = _CMP.KnotQuinticAxisLaw(
            _CMP.X_A, L_c, _CMP.M_A, _CMP.MP_A, _CMP.MPP_A, _CMP.MD, M_K)
        gates = law.gates()
        rec = _CMP.evaluate(
            law, f"A:Lc{L_c:g}:MK{M_K:.1f}:n{n_axis}", n_axis,
            extra={
                "M_K": float(M_K), "x_K": float(law.x_K),
                "ell_K": float(law.x_K - _CMP.X_A),
                "xi_K": float((law.x_K - _CMP.X_A) / L_c),
                "L_1": float(law.L_1), "L_2": float(law.L_2),
                "s_K": float(law.s_K), "axis_gates": gates,
            })
        return _compact(rec)
    except Exception as exc:  # noqa: BLE001 -- 掃引は失敗点も結果として保存する
        return {"L_c": float(L_c), "M_K": float(M_K), "n_axis": int(n_axis),
                "error": f"{type(exc).__name__}: {exc}"}


def _evaluate_job(job: tuple[float, float, int]) -> dict:
    return evaluate(*job)


def evaluate_many(jobs: list[tuple[float, float, int]], n_axis: int) -> list[dict]:
    """独立候補を process 並列評価する。高解像度ほど一候補のメモリが増えるため並列数を絞る。"""
    cpu = os.cpu_count() or 1
    workers = min(DEFAULT_WORKERS[n_axis], cpu, len(jobs))
    if workers <= 1:
        return [_evaluate_job(job) for job in jobs]
    # Linux/WSL の設計環境では fork で構築済み gas/Hall table を copy-on-write 共有する。
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("fork")) as pool:
        return list(pool.map(_evaluate_job, jobs, chunksize=1))


def _valid(rec: dict) -> bool:
    return "error" not in rec and rec["hard_gate_pass"] \
        and rec["axis_gates"]["mpp_quality_ok"]


def _topo(rec: dict) -> float:
    return float(rec["topology"].get("min_sin_angle_interior", float("nan")))


def _margin(rec: dict) -> float:
    return float(rec["margin"]["min_mu_minus_theta_deg"])


def select_finalists(records: list[dict]) -> dict[str, dict]:
    """各 L_c で意味の異なる代表候補を選ぶ。単一の恣意的スコアには潰さない。"""
    valid = [r for r in records if _valid(r)]
    if not valid:
        return {}
    out = {
        "min_theta": min(valid, key=lambda r: r["theta_max_deg"]),
        "max_margin": max(valid, key=_margin),
    }
    baseline = [r for r in valid if abs(r["M_K"] - BASELINE_MK) < 1e-9]
    if baseline:
        out["baseline"] = baseline[0]
    robust = [r for r in valid if _margin(r) >= MARGIN_FLOOR_DEG and _topo(r) >= TOPO_FLOOR]
    if robust:
        out["robust_min_theta"] = min(robust, key=lambda r: r["theta_max_deg"])
    return out


def csv_row(rec: dict, role: str = "sweep") -> list:
    if "error" in rec:
        return [role, rec["L_c"], rec["M_K"], rec["n_axis"], "ERROR", rec["error"]]
    ax = rec["axis"]
    return [
        role, rec["L_c"], rec["M_K"], rec["n_axis"], rec["hard_gate_pass"], "",
        rec["x_K"], rec["ell_K"], rec["xi_K"], rec["s_K"],
        rec["axis_gates"]["mpp_sign_changes"], rec["axis_gates"]["mpp_quality_ok"],
        rec["theta_max_deg"], _margin(rec),
        rec["topology"]["n_orientation_flip_interior"], _topo(rec),
        rec["wall_curvature"].get("downstream_max_abs_dkappa_ds"),
        rec["wall_curvature"].get("downstream_J_dkappa_ds2"),
        ax.get("max_breakpoint_jump"), ax["max_abs_Mppp"], ax["J_axis"],
    ]


def plot_sensitivity(records: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    cmap = plt.get_cmap("viridis")
    for j, L_c in enumerate(LCS):
        rr = sorted((r for r in records if r.get("L_c") == L_c and _valid(r)),
                    key=lambda r: r["M_K"])
        if not rr:
            continue
        mk = np.array([r["M_K"] for r in rr])
        color = cmap(j / (len(LCS) - 1))
        label = f"L_c={L_c:g}"
        axes[0, 0].plot(mk, [r["theta_max_deg"] for r in rr], color=color, label=label)
        axes[0, 1].plot(mk, [_margin(r) for r in rr], color=color)
        axes[1, 0].plot(mk, [_topo(r) for r in rr], color=color)
        axes[1, 1].plot(mk, [r["xi_K"] for r in rr], color=color)
    axes[0, 0].set_ylabel(r"$\theta_{max}$ [deg]")
    axes[0, 1].set_ylabel(r"$\min(\mu_w-\theta_w)$ [deg]")
    axes[1, 0].set_ylabel("min sin(characteristic angle)")
    axes[1, 1].set_ylabel(r"$(x_K-x_A)/L_c$")
    axes[0, 1].axhline(MARGIN_FLOOR_DEG, color="0.5", ls="--", lw=1)
    axes[1, 0].axhline(TOPO_FLOOR, color="0.5", ls="--", lw=1)
    for ax in axes.flat:
        ax.axvline(BASELINE_MK, color="0.65", ls=":", lw=1)
        ax.grid(True, alpha=0.25)
    axes[1, 0].set_xlabel(r"$M_K$")
    axes[1, 1].set_xlabel(r"$M_K$")
    axes[0, 0].legend(ncol=2, fontsize=9)
    fig.suptitle("Axis law A: knot Mach sensitivity (n_axis=600 sweep)")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    all_records: list[dict] = []
    records_600 = []
    print("--- full sweep n_axis=600 ---", flush=True)
    records_600 = evaluate_many(
        [(L_c, float(M_K), 600) for L_c in LCS for M_K in MKS], 600)
    for L_c in LCS:
        valid = [r for r in records_600 if r.get("L_c") == L_c and _valid(r)]
        print(f"L_c={L_c:g}: valid {len(valid)}/{len(MKS)}", flush=True)
    all_records.extend(records_600)

    representatives_1200: list[dict] = []
    selected_600: dict[str, dict[str, dict]] = {}
    jobs_1200 = []
    roles_1200 = {}
    for L_c in LCS:
        rr = [r for r in records_600 if r.get("L_c") == L_c]
        selected = select_finalists(rr)
        selected_600[f"Lc{L_c:g}"] = selected
        seen = set()
        for role, source in selected.items():
            M_K = float(source["M_K"])
            if M_K in seen:
                continue
            seen.add(M_K)
            print(f"representative n=1200 L_c={L_c:g} M_K={M_K:.1f} ({role})", flush=True)
            jobs_1200.append((L_c, M_K, 1200))
            roles_1200[(L_c, M_K)] = [k for k, v in selected.items()
                                      if abs(float(v["M_K"]) - M_K) < 1e-9]
    representatives_1200 = evaluate_many(jobs_1200, 1200)
    for rec in representatives_1200:
        rec["selection_roles"] = roles_1200[(float(rec["L_c"]), float(rec["M_K"]))]
    all_records.extend(representatives_1200)

    selected_1200: dict[str, dict[str, dict]] = {}
    final_records_2400: list[dict] = []
    jobs_2400 = []
    roles_2400 = {}
    for L_c in LCS:
        rr = [r for r in representatives_1200 if r.get("L_c") == L_c]
        selected = select_finalists(rr)
        # 最終解像度は生産基準と、ゲート内で壁角最小の候補に限定する。
        keep = {k: v for k, v in selected.items() if k in ("baseline", "robust_min_theta")}
        selected_1200[f"Lc{L_c:g}"] = keep
        seen = set()
        for role, source in keep.items():
            M_K = float(source["M_K"])
            if M_K in seen:
                continue
            seen.add(M_K)
            print(f"final n=2400 L_c={L_c:g} M_K={M_K:.1f} ({role})", flush=True)
            jobs_2400.append((L_c, M_K, 2400))
            roles_2400[(L_c, M_K)] = [k for k, v in keep.items()
                                      if abs(float(v["M_K"]) - M_K) < 1e-9]
    final_records_2400 = evaluate_many(jobs_2400, 2400)
    for rec in final_records_2400:
        rec["selection_roles"] = roles_2400[(float(rec["L_c"]), float(rec["M_K"]))]
    all_records.extend(final_records_2400)

    summary = {
        "conditions": {
            "R": _CMP.R, "M_d": _CMP.MD, "Tt": _CMP.TT,
            "L_c_values": list(LCS), "M_K_min": float(MKS.min()),
            "M_K_max": float(MKS.max()), "M_K_step": 0.1,
            "baseline_M_K": BASELINE_MK, "topology_floor": TOPO_FLOOR,
            "margin_floor_deg": MARGIN_FLOOR_DEG,
        },
        "representatives_from_600": {
            lc: {role: {k: rec[k] for k in ("M_K", "x_K", "xi_K", "theta_max_deg")}
                 | {"min_mu_minus_theta_deg": _margin(rec), "topology_margin": _topo(rec)}
                 for role, rec in roles.items()}
            for lc, roles in selected_600.items()
        },
        "finalists_from_1200": {
            lc: {role: {k: rec[k] for k in ("M_K", "x_K", "xi_K", "theta_max_deg")}
                 | {"min_mu_minus_theta_deg": _margin(rec), "topology_margin": _topo(rec)}
                 for role, rec in roles.items()}
            for lc, roles in selected_1200.items()
        },
        "records": all_records,
    }
    json_path = CASE / "sweep_axislaw_A_MK.json"
    csv_path = CASE / "sweep_axislaw_A_MK_summary.csv"
    png_path = CASE / "sweep_axislaw_A_MK.png"
    json_path.write_text(json.dumps(summary, indent=1, allow_nan=True) + "\n")

    header = [
        "role", "L_c", "M_K", "n_axis", "hard_gate", "error", "x_K", "ell_K", "xi_K", "s_K",
        "Mpp_sign_changes", "single_peak", "theta_max_deg", "min_mu_minus_theta_deg",
        "flip_interior", "topology_margin", "downstream_max_dkds", "downstream_J_dkds2",
        "Mppp_jump", "max_Mppp", "J_axis",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(header)
        for rec in all_records:
            role = "+".join(rec.get("selection_roles", ["sweep"]))
            writer.writerow(csv_row(rec, role))
    plot_sensitivity(records_600, png_path)
    print(f"saved {json_path}", flush=True)
    print(f"saved {csv_path}", flush=True)
    print(f"saved {png_path}", flush=True)


if __name__ == "__main__":
    main()
