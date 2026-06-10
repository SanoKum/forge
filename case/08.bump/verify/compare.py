#!/usr/bin/env python3
"""08.bump 検証の合否判定。
(1) y≈0.3 の数値解 (P/Pt, Mach) が baseline と一致するか
(2) 収束度合い (残差 floor・1e-4 到達ステップ・発散有無) が baseline より悪化していないか
を各 run について判定し、PASS/FAIL を表示する。差分プロットも出力する。

usage: compare.py --workdir <dir>  (各 run_<case>_<method> と baseline_*_y03.csv を参照)
"""
import argparse, csv, math, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
# 数値解の許容相対差 (baseline floor が完全ゼロでないため 0.5% を閾値とする)
TOL_PROFILE = 0.005      # P/Pt, Mach の最大相対差
# 収束悪化の許容: floor が baseline の何倍まで許すか / 1e-4 到達ステップが何倍まで許すか
TOL_FLOOR_FACTOR = 3.0
TOL_STEP_FACTOR  = 3.0

RUNS = [  # (case, method, runDir, mesh, Pt, lastRes)
    ("loM", "exp", "run_loM_exp", "bump.h5",      120193.0),
    ("loM", "imp", "run_loM_imp", "bump.h5",      120193.0),
    ("hiM", "exp", "run_hiM_exp", "bump_4pct.h5", 100000.0),
    ("hiM", "imp", "run_hiM_imp", "bump_4pct.h5", 100000.0),
]

def load_csv(path):
    rows = list(csv.DictReader(open(path)))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}

def conv_baseline():
    out = {}
    for r in csv.DictReader(open(os.path.join(HERE, "baseline_convergence.csv"))):
        out[(r["case"], r["method"])] = r
    return out

def residual_metrics(run):
    rows = [r for r in csv.DictReader(open(os.path.join(run, "residual_history.csv")))
            if r["phase"] == "outer_begin"]
    final = None; tstep = None; nan = False
    for r in rows:
        s = int(r["step"]); v = float(r["rms_ro"])
        if math.isfinite(v):
            final = v
            if tstep is None and v < 1.0e-4: tstep = s
        else:
            nan = True
    return final, tstep, nan

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    a = ap.parse_args()
    convb = conv_baseline()
    overall_ok = True
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    print(f"{'case/meth':<10}{'profile':<26}{'convergence':<40}{'verdict':<8}")
    print("-" * 84)
    n_checked = 0
    for i, (case, meth, rd, mesh, pt) in enumerate(RUNS):
        run = os.path.join(a.workdir, rd)
        if not os.path.exists(os.path.join(run, "residual_history.csv")):
            continue  # このモードでは未実行 (quick は陰解法のみ)
        n_checked += 1
        ok = True; notes = []
        # --- profile ---
        latest = os.path.join(run, "profile_y03.csv")
        base = os.path.join(HERE, f"baseline_{case}_{meth}_y03.csv")
        dP = dM = float("nan")
        if os.path.exists(latest) and os.path.exists(base):
            L = load_csv(latest); B = load_csv(base)
            n = min(len(B["P_over_Pt"]), len(L["P_over_Pt"]))
            dP = float(np.max(np.abs(L["P_over_Pt"][:n] - B["P_over_Pt"][:n]) /
                              np.maximum(np.abs(B["P_over_Pt"][:n]), 1e-9)))
            dM = float(np.max(np.abs(L["Mach"][:n] - B["Mach"][:n]) /
                              np.maximum(np.abs(B["Mach"][:n]), 1e-9)))
            if dP > TOL_PROFILE or dM > TOL_PROFILE:
                ok = False; notes.append("profile drift")
            ax = axes.flat[i]
            ax.plot(B["x_over_L"], B["P_over_Pt"], "-", color="k", lw=1.4, label="baseline P/Pt")
            ax.plot(L["x_over_L"], L["P_over_Pt"], "--", color="tab:red", lw=1.2, label="latest P/Pt")
            ax2 = ax.twinx()
            ax2.plot(B["x_over_L"], B["Mach"], "-", color="tab:blue", lw=1.0, alpha=0.6)
            ax2.plot(L["x_over_L"], L["Mach"], ":", color="tab:green", lw=1.2)
            ax.set_title(f"{case} {meth}: ΔP={dP*100:.3f}% ΔM={dM*100:.3f}%")
            ax.set_xlabel("x/L"); ax.set_ylabel("P/Pt"); ax2.set_ylabel("Mach"); ax.legend(fontsize=7)
        else:
            ok = False; notes.append("missing profile csv")
        prof_str = f"ΔP={dP*100:.3f}% ΔM={dM*100:.3f}% (tol {TOL_PROFILE*100:.1f}%)"
        # --- convergence ---
        final, tstep, nan = residual_metrics(run)
        b = convb.get((case, meth), {})
        bfloor = float(b.get("final_rms_ro", "nan"))
        bstep = b.get("steps_to_1e-4", "None")
        bstep_i = int(bstep) if bstep not in ("None", "", None) else None
        conv_notes = []
        if nan or final is None:
            ok = False; conv_notes.append("DIVERGED")
        else:
            if final > bfloor * TOL_FLOOR_FACTOR:
                ok = False; conv_notes.append("floor worse")
            if bstep_i is not None and tstep is not None and tstep > bstep_i * TOL_STEP_FACTOR + 5:
                ok = False; conv_notes.append("slower")
            if bstep_i is not None and tstep is None:
                ok = False; conv_notes.append("never<1e-4")
        conv_str = (f"floor={final:.2e}/base{bfloor:.1e} "
                    f"s<1e-4={tstep}/base{bstep} {';'.join(conv_notes) or 'ok'}")
        verdict = "PASS" if ok else "FAIL"
        overall_ok = overall_ok and ok
        print(f"{case+'/'+meth:<10}{prof_str:<26}{conv_str:<40}{verdict:<8}"
              + ("  " + ",".join(notes) if notes else ""))

    fig.suptitle("08.bump verification: y~0.3 profile  baseline(solid) vs latest(dashed)")
    fig.tight_layout()
    out = os.path.join(a.workdir, "verify_profiles.png")
    fig.savefig(out, dpi=110); plt.close(fig)
    print("-" * 84)
    print("diff plot:", out)
    print(f"OVERALL ({n_checked} run checked):", "PASS ✅" if overall_ok else "FAIL ❌")
    return 0 if overall_ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
