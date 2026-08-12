"""前線ポリッシュ: チャンク継続 + η ドリフトゲート (子 plan §5 ステップ 5b)。

残差プラトー (STALLED) は既知のリミットサイクルで「未収束」とは限らないが、
報告量そのものの頭打ちで確認するのが規律 (準定常確認の η 直接版)。各対象 run の
最終 res から同一メッシュ継続 (+nStepOuter) を回し、|Δη| < tol (既定 3e-4 ≈
0.03%) に頭打ちするまで最大 max_chunks 回続ける。継続 run は `<tag>_pK`。

CLI: design/.venv-opt/bin/python -m forge_design.opt.polish <campaign_dir>
     (campaign の pareto.json の全点をポリッシュし polish.jsonl /
      pareto_polished.json を書く)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ..evaluate import runner
from .ehvi import nondominated_mask


def _last_res(d: Path):
    rs = sorted(Path(d).glob("res_[0-9]*.h5"),
                key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    return rs[-1] if rs else None


def polish_run(prob_yaml, run_dir, eta0: float, tol: float = 3e-4,
               max_chunks: int = 2) -> dict:
    """run_dir の解を η が頭打ちするまでチャンク継続する。"""
    run_dir = Path(run_dir)
    prev_res = _last_res(run_dir)
    eta_prev = float(eta0)
    hist = [eta_prev]
    settled = False
    for k in range(1, max_chunks + 1):
        nd = run_dir.parent / f"{run_dir.name}_p{k}"
        if not nd.exists():  # 再開可能
            runner.prepare(prob_yaml, nd, warm_from=prev_res)
            runner.run_forge(nd)
        out = runner.collect(prob_yaml, nd)
        eta = float(out["eta_cf"])
        hist.append(eta)
        drift, eta_prev, prev_res = abs(eta - eta_prev), eta, _last_res(nd)
        if drift < tol:
            settled = True
            break
    return {"eta_final": eta_prev, "eta_hist": hist, "settled": settled,
            "chunks": len(hist) - 1}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="パレート前線の η 頭打ちポリッシュ")
    ap.add_argument("campaign_dir")
    ap.add_argument("--tol", type=float, default=3e-4)
    ap.add_argument("--max-chunks", type=int, default=2)
    a = ap.parse_args(argv)
    d = Path(a.campaign_dir)
    rows = [json.loads(l) for l in (d / "ledger.jsonl").read_text().splitlines() if l]
    front = json.load(open(d / "pareto.json"))["pareto"]
    out_rows = []
    for q in front:
        row = next(r for r in rows if r["status"] == "PASS"
                   and np.allclose(r["x"], q["x"], atol=1e-12))
        res = polish_run(d / f"{row['tag']}.yaml", row["run_dir"], row["eta_cf"],
                         tol=a.tol, max_chunks=a.max_chunks)
        rec = {"tag": row["tag"], "x": row["x"], "L_over_rt": q["L_over_rt"],
               "eta_raw": row["eta_cf"], **res}
        out_rows.append(rec)
        print(f"[{row['tag']}] η {row['eta_cf']:.5f} -> {res['eta_final']:.5f} "
              f"(chunks={res['chunks']}, settled={res['settled']})", flush=True)
    (d / "polish.jsonl").write_text(
        "\n".join(json.dumps(r, ensure_ascii=False) for r in out_rows) + "\n")
    F = np.array([[-r["eta_final"], r["L_over_rt"]] for r in out_rows])
    mask = nondominated_mask(F)
    polished = [{"x": out_rows[i]["x"], "eta_cf": out_rows[i]["eta_final"],
                 "L_over_rt": out_rows[i]["L_over_rt"],
                 "settled": out_rows[i]["settled"]}
                for i in np.where(mask)[0]]
    polished.sort(key=lambda r: r["L_over_rt"])
    (d / "pareto_polished.json").write_text(
        json.dumps({"tol": a.tol, "pareto": polished}, indent=1, ensure_ascii=False))
    n_set = sum(r["settled"] for r in out_rows)
    print(f"settled {n_set}/{len(out_rows)}, polished front {len(polished)} 点", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
