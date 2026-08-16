"""軸 Mach 則 D (端点アンカー + 内部補間点 1 点の C⁴ 区分 5 次) の探索・continuation・A 比較。

計画: plans/accepted/tooling-nozzle-axislaw-onepoint.md。基準条件: case/42 M6
(R=3, M_d=6, semi-perfect Tt=1550K, Hall anchor, start_line=throat_char, wall_mode=cplus,
exit_mode=characteristic, axis_dx0=0.03)。

手順:
  Phase 1  L_c=45, (ξ_P, η_P) 格子 → M'≥0 篩い → n_axis=600 で MOC → hard gate + 指標
  Phase 2  L_c continuation 60→50→45→40→35→30 (前 L_c の最良 (ξ,η) 近傍を局所再探索)
  Phase 3  各 L_c の最良 D と A (knot, M_K=2.5) を n_axis 1200/2400 で再評価
出力: compare_axislaw_D.json (全候補・全診断)、compare_axislaw_D_summary.csv、PNG (scratchpad)。

再生成: design/.venv-opt/bin/python case/42.isobutane_wt/compare_axislaw_D.py [--quick]
"""
import csv
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

DESIGN = Path("/home/sano/work/forge/design")
sys.path.insert(0, str(DESIGN))

from forge_design.gas import GasSemiPerfect                                     # noqa: E402
from forge_design.geometry.axis_law import KnotQuinticAxisLaw                   # noqa: E402
from forge_design.geometry.axis_law_onepoint import OnePointC4AxisLaw           # noqa: E402
from forge_design.geometry.moc_diagnostics import (                             # noqa: E402
    axis_smoothness_diagnostics, characteristic_topology_diagnostics,
    wall_curvature_diagnostics, wall_margin_diagnostics)
from forge_design.geometry.moc_inverse import inverse_design                    # noqa: E402
from forge_design.geometry.transonic import HallThroat                         # noqa: E402
from forge_design.geometry.wall_axismach import (AxisMachCFDWall,              # noqa: E402
                                                 area_ratio_isentropic, wall_qa)

CASE = Path("/home/sano/work/forge/case/42.isobutane_wt")
Y = {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.022, "N2": 0.7245}
TT, MD, R, MK = 1550.0, 6.0, 3.0, 2.5
R_U, L_U, L_PIPE = 9.303, 12.0, 0.5
AXIS_DX0 = 0.03
QUICK = "--quick" in sys.argv

gas = GasSemiPerfect(Y, Tt=TT)
ht = HallThroat(R=R, gamma=gas.gamma_star)
xs_c, rr, MM, tt = ht.throat_characteristic(n=41)
X_A = float(xs_c[0])
M_A, MP_A, MPP_A = ht.axis_anchor(X_A)
RF = float(np.sqrt(area_ratio_isentropic(MD, gas)))


def evaluate(law, tag, n_axis, extra=None):
    """逆 MOC + 全診断。戻り dict (JSON 可)。"""
    x_E = law.x_E
    x_end = x_E + 2.3 * RF * np.sqrt(MD * MD - 1.0)
    t0 = time.time()
    res = inverse_design(ht, lambda x: float(law(max(x, X_A))), x_axis_end=x_end,
                         n_axis=n_axis, n_start=41, gamma=gas, M_start=1.05,
                         exit_mode="characteristic", x_E=x_E, M_d=MD,
                         start_line="throat_char", wall_mode="cplus", axis_dx0=AXIS_DX0)
    t_moc = time.time() - t0
    w = res["wall"]
    qa = wall_qa(w, MD, x_E, gas, R=R)
    try:
        cw = AxisMachCFDWall(w, R=R, r_U=R_U, L_U=L_U, L_pipe=L_PIPE)
        cfd_msgs = cw.validate()
        wc = wall_curvature_diagnostics(cw, w[0, 0], w[-1, 0], x_split=0.25)
    except Exception as e:  # noqa: BLE001
        cfd_msgs = [f"CFDWall 構築失敗: {e}"]
        wc = {}
    bps = [law.x_K] if hasattr(law, "x_K") else ([law.x_P] if hasattr(law, "x_P") else None)
    ax = axis_smoothness_diagnostics(law, X_A, x_E, breakpoints=bps)
    topo = characteristic_topology_diagnostics(res["levels"], n_axis - 1,
                                               n_crossing_sample=(20 if n_axis <= 600 else 40),
                                               seam_i_pad=90, seam_k_max=110,
                                               seam_x_max=X_A + 1.0, wall=w[:, :2])
    margin = wall_margin_diagnostics(w)
    wall_mono = bool(np.all(np.diff(w[:, 1]) >= -1e-9))
    hard = (not qa["violations"]) and (not cfd_msgs) and wall_mono \
        and topo["n_orientation_flip_interior"] == 0 and topo["n_crossings_sampled"] == 0 \
        and topo["pileup_pairs_interior"] == 0
    rec = {"tag": tag, "n_axis": n_axis, "t_moc_s": t_moc, "L_c": float(law.L_c),
           "x_A": X_A, "x_E": float(x_E), "x_F": qa["x_F"], "r_F": qa["r_F"],
           "xF_minus_xE": qa["xF_minus_xE"], "xF_minus_xE_pred": qa["xF_minus_xE_pred"],
           "theta_max_deg": qa["theta_max_deg"], "theta_exit_deg": qa["theta_exit_deg"],
           "kappa0_R_ratio": qa.get("kappa0_R_ratio"), "r_F_err_rel": qa["r_F_err_rel"],
           "wall_qa_violations": qa["violations"], "cfd_wall_violations": cfd_msgs,
           "wall_monotone": wall_mono, "axis": ax, "topology": topo,
           "wall_curvature": wc, "margin": margin, "hard_gate_pass": bool(hard)}
    if extra:
        rec.update(extra)
    if n_axis >= 1200:
        rec["wall_table"] = w.tolist()
    return rec


def d_extra(law):
    return {"xi_P": law.xi_P, "eta_P": law.eta_P, "x_P": law.x_P, "M_P": law.M_P,
            "ell_P": law.ell_P, "ell_P_over_R": law.ell_P / R, "L_1": law.L_1, "L_2": law.L_2,
            "gates": {k: v for k, v in law.gates().items() if k != "P_jumps"}}


TOPO_MARGIN_FLOOR = 0.02   # min sin(角) がこれ未満なら「退化寸前」として順位を落とす


def rank_key(rec):
    """順位付け (計画 §4.3): hard gate → トポロジ余裕 → μ−θ 余裕 → 下流 fairness → θ_max → M3。

    トポロジ余裕 (内部三角形の最小 sin 角) は**しきい値**として使う: 実測 (M6/R3/L_c45,
    n=600) では隣接候補間で 0.03↔0.06 と ±50% 揺れ、最小位置も x=11↔28 と飛ぶ (網の自然な
    最小角の場所が変わるだけ) ため、連続量として細かく順位付けすると雑音で決まる。
    flip=0 かつ床 (0.02) 以上を「健全」とし、その中で μ−θ 余裕 → 下流 fairness → θ_max で並べる。"""
    if "error" in rec:
        return (1, 1, 0, 0, 0, 0)
    topo = rec["topology"]; wc = rec["wall_curvature"]
    tm = topo.get("min_sin_angle_interior", float("nan"))
    tm = -1.0 if not np.isfinite(tm) else tm
    return (0 if rec["hard_gate_pass"] else 1,
            0 if tm >= TOPO_MARGIN_FLOOR else 1,
            -round(rec["margin"]["min_mu_minus_theta_deg"], 1),
            round(wc.get("downstream_max_abs_dkappa_ds", 1e9), 2),
            round(rec["theta_max_deg"], 1),
            rec["axis"]["max_abs_Mppp"])


def make_D(L_c, xi, eta):
    return OnePointC4AxisLaw(X_A, L_c, M_A, MP_A, MPP_A, MD, xi, eta)


def make_A(L_c):
    return KnotQuinticAxisLaw(X_A, L_c, M_A, MP_A, MPP_A, MD, MK)


def feasible_band(L_c, xi, eta_lo=0.05, eta_hi=0.95, step=0.0025):
    """ξ 固定で M'≥0 となる η の連結帯 [η_lo, η_hi] (無ければ None)。帯は L_c=45 で幅
    0.01–0.07 と極めて細いので、粗い η 格子では見落とす (実測: 9×14 格子で 8/126 しか
    当たらない) — 帯を先に見つけてからその中に候補を置く。"""
    etas = np.arange(eta_lo, eta_hi + 1e-9, step)
    ok = np.array([OnePointC4AxisLaw.monotone_feasible(X_A, L_c, M_A, MP_A, MPP_A, MD, xi, e)
                   for e in etas])
    if not ok.any():
        return None
    return float(etas[ok].min()), float(etas[ok].max())


def scan(L_c, xis, etas, n_axis, label, n_in_band=3):
    """各 ξ で単調帯を求め、帯の下端/中央/上端 (n_in_band 点) を候補にする (`etas` は
    互換のため受けるが帯探索の範囲 [min,max] にだけ使う)。"""
    out = []
    n_feas = 0
    bands = {}
    eta_lo, eta_hi = float(np.min(etas)), float(np.max(etas))
    for xi in xis:
        band = feasible_band(L_c, xi, eta_lo, eta_hi)
        bands[float(xi)] = band
        if band is None:
            out.append({"tag": f"D:{label}", "L_c": L_c, "xi_P": float(xi), "eta_P": None,
                        "n_axis": n_axis, "monotone": False, "hard_gate_pass": False,
                        "error": "M'<0 (帯なし)"})
            continue
        for eta in np.linspace(band[0], band[1], n_in_band):
            eta = float(round(eta, 5))
            n_feas += 1
            law = make_D(L_c, xi, eta)
            try:
                rec = evaluate(law, f"D:{label}", n_axis, extra=d_extra(law))
                rec["monotone"] = True
                out.append(rec)
                print(f"  L_c={L_c:5.1f} xi={xi:.3f} eta={eta:.3f} M_P={law.M_P:.2f} "
                      f"gate={rec['hard_gate_pass']} θmax={rec['theta_max_deg']:.2f} "
                      f"μ-θ={rec['margin']['min_mu_minus_theta_deg']:.2f} "
                      f"flip_int={rec['topology']['n_orientation_flip_interior']} "
                      f"tm={rec['topology'].get('min_sin_angle_interior', float('nan')):.3f} "
                      f"dk_ds(dn)={rec['wall_curvature'].get('downstream_max_abs_dkappa_ds', float('nan')):.3f}",
                      flush=True)
            except Exception as e:  # noqa: BLE001
                out.append({"tag": f"D:{label}", "L_c": L_c, "xi_P": xi, "eta_P": eta,
                            "n_axis": n_axis, "monotone": True, "hard_gate_pass": False,
                            "error": str(e)[:200]})
                print(f"  L_c={L_c:5.1f} xi={xi:.3f} eta={eta:.3f} ERROR {str(e)[:80]}", flush=True)
    print(f"  [{label}] monotone-feasible candidates {n_feas} (帯あり ξ: {sum(1 for b in bands.values() if b)}/{len(xis)})", flush=True)
    for xi_, b in bands.items():
        if b:
            print(f"     ξ={xi_:.3f} ℓ={xi_*L_c:5.2f}: η∈[{b[0]:.3f},{b[1]:.3f}] M_P∈[{M_A+b[0]*(MD-M_A):.2f},{M_A+b[1]*(MD-M_A):.2f}]", flush=True)
    out.append({"tag": f"D:{label}:bands", "L_c": L_c, "bands": {str(k): v for k, v in bands.items()}, "error": "meta"})
    return out


def best_of(recs):
    ok = [r for r in recs if "error" not in r and r.get("hard_gate_pass")]
    if not ok:
        return None
    return sorted(ok, key=rank_key)[0]


def main():
    all_records = []
    summary = {}
    # ---------------- Phase 1: L_c = 45 格子 ----------------
    print("=== Phase 1: L_c=45 grid (n_axis=600) ===", flush=True)
    xis = np.round(np.arange(0.02, 0.301, 0.01 if not QUICK else 0.04), 4)
    etas = np.array([0.05, 0.95])       # 帯探索の範囲 (実際の η は帯の中に置く)
    p1 = scan(45.0, xis, etas, 600, "P1_Lc45")
    all_records += p1
    # 基準候補 (A の knot 相当)
    xi0, eta0 = 3.74 / 45.0, (3.0 - M_A) / (MD - M_A)
    print(f"  基準候補 xi={xi0:.4f} eta={eta0:.4f}", flush=True)
    if OnePointC4AxisLaw.monotone_feasible(X_A, 45.0, M_A, MP_A, MPP_A, MD, xi0, eta0):
        law0 = make_D(45.0, xi0, eta0)
        r0 = evaluate(law0, "D:P1_ref", 600, extra=d_extra(law0)); r0["monotone"] = True
        all_records.append(r0)
    b1 = best_of(p1)
    summary["phase1_best"] = None if b1 is None else {k: b1[k] for k in ("xi_P", "eta_P", "M_P", "ell_P", "theta_max_deg", "hard_gate_pass")}
    print("  Phase1 best:", summary["phase1_best"], flush=True)

    # ---------------- Phase 2: L_c continuation ----------------
    print("=== Phase 2: L_c continuation (n_axis=600) ===", flush=True)
    Lcs = [60.0, 50.0, 45.0, 40.0, 35.0, 30.0] if not QUICK else [60.0, 45.0, 30.0]
    prev = None
    cont = {}
    for L_c in Lcs:
        etas_l = np.array([0.05, 0.95])
        if prev is None:
            xis_l = np.round(np.arange(0.02, 0.301, 0.02), 4)
        else:
            # 前 L_c の最良 ξ 近傍を細かく (±0.04, 0.01 刻み) + 粗い下地 (0.02–0.30, 0.04 刻み)。
            # η は各 ξ の単調帯から取る (帯が細いので前 L_c の η をそのまま使わない)。
            xis_l = np.round(np.clip(np.arange(prev["xi_P"] - 0.04, prev["xi_P"] + 0.0401, 0.01), 0.02, 0.4), 4)
            xis_l = np.unique(np.concatenate([xis_l, np.round(np.arange(0.02, 0.301, 0.04), 4)]))
        recs = scan(L_c, xis_l, etas_l, 600, f"P2_Lc{L_c:g}")
        all_records += recs
        # A 基準 (同 L_c)
        try:
            lawA = make_A(L_c)
            rA = evaluate(lawA, f"A:Lc{L_c:g}", 600, extra={"M_knot": MK, "x_K": lawA.x_K, "ell_K": lawA.x_K - X_A})
            all_records.append(rA)
        except Exception as e:  # noqa: BLE001
            rA = {"tag": f"A:Lc{L_c:g}", "L_c": L_c, "n_axis": 600, "error": str(e)}
            all_records.append(rA)
        b = best_of(recs)
        cont[str(L_c)] = {
            "D_best": None if b is None else {k: b[k] for k in ("xi_P", "eta_P", "M_P", "ell_P", "ell_P_over_R",
                                                                "theta_max_deg", "hard_gate_pass")}
                      | ({"mu_theta": b["margin"]["min_mu_minus_theta_deg"],
                          "topo_margin": b["topology"].get("min_sin_angle_interior"),
                          "dn_dkds": b["wall_curvature"].get("downstream_max_abs_dkappa_ds"),
                          "Mppp": b["axis"]["max_abs_Mppp"]} if b else {}),
            "A": None if "error" in rA else {"theta_max_deg": rA["theta_max_deg"],
                                             "hard_gate_pass": rA["hard_gate_pass"],
                                             "mu_theta": rA["margin"]["min_mu_minus_theta_deg"],
                                             "topo_margin": rA["topology"].get("min_sin_angle_interior"),
                                             "dn_dkds": rA["wall_curvature"].get("downstream_max_abs_dkappa_ds"),
                                             "Mppp": rA["axis"]["max_abs_Mppp"]},
            "n_feasible_gate": int(sum(1 for r in recs if "error" not in r and r.get("hard_gate_pass"))),
        }
        print(f"  L_c={L_c}: D_best={cont[str(L_c)]['D_best']}  A={cont[str(L_c)]['A']}", flush=True)
        if b is not None:
            prev = b
    summary["continuation"] = cont

    # ---------------- Phase 3: 最終候補を 1200/2400 で再評価 ----------------
    print("=== Phase 3: final re-evaluation (n_axis 1200/2400) ===", flush=True)
    final = {}
    for L_c in Lcs:
        c = cont[str(L_c)]
        for n_axis in (1200, 2400):
            if c["D_best"] is not None:
                law = make_D(L_c, c["D_best"]["xi_P"], c["D_best"]["eta_P"])
                try:
                    r = evaluate(law, f"D:final_Lc{L_c:g}", n_axis, extra=d_extra(law))
                except Exception as e:  # noqa: BLE001
                    r = {"tag": f"D:final_Lc{L_c:g}", "L_c": L_c, "n_axis": n_axis, "error": str(e)}
                all_records.append(r)
                final[f"D@Lc{L_c:g}@{n_axis}"] = r
            try:
                lawA = make_A(L_c)
                rA = evaluate(lawA, f"A:final_Lc{L_c:g}", n_axis, extra={"M_knot": MK, "x_K": lawA.x_K})
            except Exception as e:  # noqa: BLE001
                rA = {"tag": f"A:final_Lc{L_c:g}", "L_c": L_c, "n_axis": n_axis, "error": str(e)}
            all_records.append(rA)
            final[f"A@Lc{L_c:g}@{n_axis}"] = rA
            for k in (f"D@Lc{L_c:g}@{n_axis}", f"A@Lc{L_c:g}@{n_axis}"):
                r = final.get(k)
                if r and "error" not in r:
                    print(f"  {k}: gate={r['hard_gate_pass']} θmax={r['theta_max_deg']:.2f} "
                          f"μ-θ={r['margin']['min_mu_minus_theta_deg']:.2f} "
                          f"flip_int={r['topology']['n_orientation_flip_interior']} "
                          f"tm={r['topology'].get('min_sin_angle_interior', float('nan')):.3f} "
                          f"dk_ds(dn)={r['wall_curvature'].get('downstream_max_abs_dkappa_ds', float('nan')):.3f} "
                          f"M'''={r['axis']['max_abs_Mppp']:.3f}", flush=True)
    summary["final"] = {k: ({kk: v[kk] for kk in ("hard_gate_pass", "theta_max_deg")} | {
        "mu_theta": v["margin"]["min_mu_minus_theta_deg"],
        "flip_int": v["topology"]["n_orientation_flip_interior"],
        "topo_margin": v["topology"].get("min_sin_angle_interior"),
        "dn_dkds": v["wall_curvature"].get("downstream_max_abs_dkappa_ds"),
        "dn_J": v["wall_curvature"].get("downstream_J_dkappa_ds2"),
        "Mppp": v["axis"]["max_abs_Mppp"], "J_axis": v["axis"]["J_axis"]}
        | ({kk: v[kk] for kk in ("xi_P", "eta_P", "ell_P", "ell_P_over_R")} if "xi_P" in v else {}))
        if "error" not in v else {"error": v["error"]} for k, v in final.items()}

    meta = {"x_A": X_A, "M_A": M_A, "Mp_A": MP_A, "Mpp_A": MPP_A, "Tt": TT, "Md": MD, "R": R,
            "M_knot_A": MK, "r_U": R_U, "L_U": L_U, "L_pipe": L_PIPE, "axis_dx0": AXIS_DX0,
            "grid_phase1": {"xi": xis.tolist(), "eta_band_search": etas.tolist(), "n_in_band": 3}, "Lc_continuation": Lcs}
    (CASE / "compare_axislaw_D.json").write_text(json.dumps(
        {"meta": meta, "summary": summary, "records": all_records}, indent=1))
    # CSV サマリ (final のみ)
    with open(CASE / "compare_axislaw_D_summary.csv", "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["key", "law", "L_c", "n_axis", "xi_P", "eta_P", "ell_P", "ell_P_over_R",
                       "hard_gate", "theta_max_deg", "min_mu_minus_theta_deg", "flip_interior",
                       "topo_margin", "downstream_max_dkds", "downstream_J_dkds2",
                       "max_Mppp", "J_axis"])
        for k, v in final.items():
            if "error" in v:
                wcsv.writerow([k, k.split("@")[0], v.get("L_c"), v.get("n_axis"), "", "", "", "", "ERR", v["error"]])
                continue
            wcsv.writerow([k, k.split("@")[0], v["L_c"], v["n_axis"],
                           v.get("xi_P", ""), v.get("eta_P", ""), v.get("ell_P", ""), v.get("ell_P_over_R", ""),
                           v["hard_gate_pass"], f"{v['theta_max_deg']:.3f}",
                           f"{v['margin']['min_mu_minus_theta_deg']:.3f}",
                           v["topology"]["n_orientation_flip_interior"],
                           f"{v['topology'].get('min_sin_angle_interior', float('nan')):.4f}",
                           f"{v['wall_curvature'].get('downstream_max_abs_dkappa_ds', float('nan')):.4f}",
                           f"{v['wall_curvature'].get('downstream_J_dkappa_ds2', float('nan')):.4f}",
                           f"{v['axis']['max_abs_Mppp']:.4f}", f"{v['axis']['J_axis']:.5f}"])
    print("saved", CASE / "compare_axislaw_D.json", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
