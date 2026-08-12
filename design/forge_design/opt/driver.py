"""③ベル TOP 幾何 dv のバッチ MOO ドライバ (Phase 2 §5 ステップ 4–5)。

DOE → forge 評価 (VERDICT ゲート) → KRG (f1=−η のみ; f2=L/rt は dv そのもの
なので厳密値) → NSGA-II + EHVI infill → 反復。使い方:

  design/.venv-opt/bin/python -m forge_design.opt.driver \
      case/40.nozzle_design_tool/problem_bell_top.yaml \
      case/40.nozzle_design_tool/run_NNNN_moo_<tag> \
      --seed-res <収束 res_*.h5> [--n-doe 24 --n-iter 13 --batch 2]

設計 (子 plan §4):
- **ゲート**: メッシュ品質 FAIL / 発散 (DIVERGED・非有限メトリクス) / quasisteady
  非 STEADY は評価失敗としてサロゲートに入れない。残差の「NOT CONVERGED
  (stalled)」プラトーはノイズ床挙動として既知 (run_0071) なので、DIVERGED
  でなく ALL STEADY なら採用する (実装定義。子 plan 変更ログ参照)。さらに
  **物理ゲート**: ṁ/ṁ_1D (チョーク流量比 = 実効 Cd) が [0.94, 1.005]、
  η∈(0.5, 1.02) — cross-geometry warm start は発散せず**不始動流という偽
  アトラクタ** (machmax~1.4, CF~7) に落ちることがあり (run_0073 smoke で実測)、
  これは quasisteady だけでは弾けないため。
- **warm seed**: 正規化 dv 距離が最近傍の PASS 評価の最終 res から interp
  warm start。初期 DOE は --seed-res (基準 run_0072 等) から。
- **段階起動を常用**: cross-geometry warm start を陰解法 cfl4 に直投入すると
  発散 (run_0019) または不始動化 (run_0073 smoke) するため、**全評価を
  soft 段 (convMethod 0 + cfl 0.5, 3000 step) → 本設定 (12000 step)** の
  2 段で回す (run_0019/0024/0034 の実証レシピの常用化)。
- **決定性**: LHS/pymoo/KRG とも seed 固定。台帳は ledger.jsonl (逐次追記)。
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import yaml

from ..evaluate import runner
from ..geometry.wall import NozzleWall
from ..probdef import load_problem
from .doe import lhs
from .ehvi import hypervolume2d, nondominated_mask
from .moo import propose_infill
from .surrogate import KrigingSet

DV_ORDER = ("theta_a_deg", "theta_e_deg", "L_over_rt")


class _F1KrgF2Exact:
    """f1=−η は KRG、f2=L/rt は dv そのもの (厳密・σ≈0) の predict アダプタ。"""

    def __init__(self, krg_f1: KrigingSet, l_index: int) -> None:
        self.krg = krg_f1
        self.l_index = l_index

    def predict(self, X):
        X = np.atleast_2d(np.asarray(X, dtype=float))
        mu1, sg1 = self.krg.predict(X)
        L = X[:, self.l_index]
        mu = np.column_stack([mu1[:, 0], L])
        sg = np.column_stack([sg1[:, 0], np.full(L.shape, 1e-9)])
        return mu, sg


class Campaign:
    def __init__(self, base_yaml, campaign_dir, seed_res, ref=(-0.90, 13.0),
                 seed: int = 0) -> None:
        self.base_yaml = Path(base_yaml)
        self.dir = Path(campaign_dir)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.seed_res = Path(seed_res).resolve()
        self.ref = tuple(ref)
        self.seed = int(seed)
        self.base_raw = yaml.safe_load(self.base_yaml.read_text())
        self.bounds = np.array([[float(self.base_raw["dv"][k]["min"]),
                                 float(self.base_raw["dv"][k]["max"])] for k in DV_ORDER])
        self.ledger = self.dir / "ledger.jsonl"
        self.rows: list = []
        if self.ledger.exists():  # 再開: 既存台帳を読み戻す (評価済みはスキップ)
            self.rows = [json.loads(l) for l in self.ledger.read_text().splitlines() if l]
        # 1D チョーク流量 (Cd=1): 物理ゲートの分母
        g, cp = float(self.base_raw["gas"]["gamma"]), float(self.base_raw["gas"]["cp"])
        R = cp * (g - 1.0) / g
        s = self.base_raw["spec"]
        At = np.pi * float(s["r_throat"]) ** 2
        self.mdot_1d = (float(s["Pt"]) * At * np.sqrt(g / (R * float(s["Tt"])))
                        * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0))))

    # -- 幾何整合の事前フィルタ (CFD なし) ------------------------------------
    def feasible(self, x) -> bool:
        try:
            self._wall(x)
            return True
        except ValueError:
            return False

    def _wall(self, x) -> NozzleWall:
        p = load_problem(self._write_problem(x, self.dir / "_feas_probe.yaml"))
        wall, _ = runner.build_wall(p)
        return wall

    def _write_problem(self, x, path: Path) -> Path:
        raw = json.loads(json.dumps(self.base_raw))  # deep copy
        for k, v in zip(DV_ORDER, np.asarray(x, dtype=float)):
            raw["dv"][k]["value"] = float(v)
        path.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
        return path

    # -- warm seed: 正規化 dv 距離の最近傍 PASS -------------------------------
    def _warm_for(self, x) -> Path:
        lo, span = self.bounds[:, 0], self.bounds[:, 1] - self.bounds[:, 0]
        best, best_d = self.seed_res, None
        xn = (np.asarray(x, dtype=float) - lo) / span
        for r in self.rows:
            if r["status"] != "PASS":
                continue
            dn = float(np.linalg.norm((np.array(r["x"]) - lo) / span - xn))
            if best_d is None or dn < best_d:
                best, best_d = Path(r["run_dir"]) / r["res_file"], dn
        return best

    # -- 1 点評価 (ゲート込み) --------------------------------------------------
    def evaluate(self, x, tag: str) -> dict:
        x = [float(v) for v in np.asarray(x, dtype=float)]
        t0 = time.time()
        run_dir = self.dir / tag
        prob = self._write_problem(x, self.dir / f"{tag}.yaml")
        row = {"tag": tag, "x": x, "run_dir": str(run_dir), "status": "FAIL",
               "warm_from": "", "note": ""}
        try:
            warm = self._warm_for(x)
            row["warm_from"] = str(warm)
            # soft 段 (常用): convMethod 0 + cfl 0.5 で 3000 step 流れを再確立
            stage = self.dir / (tag + "_stgA")
            n_stage = 3000
            runner.prepare(prob, stage, nsteps=n_stage, warm_from=warm)
            cfg = (stage / "solverConfig.yaml").read_text()
            cfg = cfg.replace("cfl: 4.0, cfl_pseudo: 4.0", "cfl: 0.5, cfl_pseudo: 0.5")
            cfg = cfg.replace("convMethod: 1", "convMethod: 0")
            # 収束場を必ずダンプさせる (outStepInterval > n_stage だと res_0 しか
            # 出ず、本段が seed 直投入相当になって発散する — smoke で実測したバグ)
            import re
            cfg = re.sub(r"outStepInterval: \d+", f"outStepInterval: {n_stage}", cfg)
            (stage / "solverConfig.yaml").write_text(cfg)
            runner.run_forge(stage)
            sres = sorted(stage.glob("res_[0-9]*.h5"),
                          key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
            if (not sres
                    or int("".join(c for c in sres[-1].stem if c.isdigit())) < n_stage
                    or self._diverged(stage, {"eta_cf": 0.9})):
                raise RuntimeError("soft 段が発散または最終場未出力")
            # 本段 (12000 step, 全域 2 次 + cfl4)
            runner.prepare(prob, run_dir, warm_from=sres[-1])
            runner.run_forge(run_dir)
            out = runner.collect(prob, run_dir)
            verdict = (run_dir / "CONVERGENCE_VERDICT.txt").read_text()
            qs = self._quasisteady(run_dir)
            eta = out.get("eta_cf")
            mrel = (out.get("mdot_kg_s") / self.mdot_1d
                    if out.get("mdot_kg_s") is not None else None)
            row.update(eta_cf=eta, mdot=out.get("mdot_kg_s"),
                       mdot_ratio_1d=mrel,
                       thrust_N=out.get("thrust_N"), res_file=out.get("res_file", ""),
                       conv="DIVERGED" if "DIVERGED" in verdict else
                            ("PASS" if "PASS" in verdict else "STALLED"),
                       quasisteady=qs)
            n_main = int(self.base_raw["evaluate"].get("nStepOuter", 12000))
            res_step = int("".join(c for c in Path(out.get("res_file", "res_-1"))
                                   .stem if c.isdigit()) or -1)
            row["res_step"] = res_step
            ok = (row["conv"] != "DIVERGED" and qs == "ALL STEADY"
                  and res_step >= n_main                    # 最終場が最後まで到達
                  and eta is not None and np.isfinite(eta)
                  and 0.5 < eta < 1.02                      # 物理ゲート (η)
                  and mrel is not None and 0.94 < mrel < 1.005)  # 物理ゲート (Cd)
            row["status"] = "PASS" if ok else "FAIL"
        except Exception as e:  # メッシュ FAIL・幾何 NG・forge 異常終了
            row["note"] += f"{type(e).__name__}: {e}"[:300]
        row["wall_s"] = round(time.time() - t0, 1)
        self.rows.append(row)
        with open(self.ledger, "a") as f:
            f.write(json.dumps(row, ensure_ascii=False) + "\n")
        print(f"[{tag}] {row['status']} x={np.round(x,3)} eta={row.get('eta_cf')}"
              f" ({row['wall_s']}s) {row['note']}", flush=True)
        return row

    def _diverged(self, run_dir: Path, out: dict) -> bool:
        v = (run_dir / "CONVERGENCE_VERDICT.txt")
        txt = v.read_text() if v.exists() else ""
        eta = out.get("eta_cf")
        return "DIVERGED" in txt or eta is None or not np.isfinite(eta)

    def _quasisteady(self, run_dir: Path) -> str:
        r = subprocess.run([sys.executable,
                            str(runner.FORGE_TOOLS / "check_quasisteady.py"), str(run_dir)],
                           capture_output=True, text=True, env=runner._ENV)
        txt = r.stdout + r.stderr
        (run_dir / "QUASISTEADY_VERDICT.txt").write_text(txt)
        overall = [l for l in txt.splitlines() if "OVERALL:" in l]
        if not overall:
            return "UNKNOWN"
        line = overall[-1]
        if "NOT" in line:  # "NOT ALL STEADY" は "ALL STEADY" を部分文字列に含む
            return "NOT ALL STEADY"
        return "ALL STEADY" if "ALL STEADY" in line else line.split("OVERALL:")[-1].strip()

    # -- 目的ベクトル (最小化): f1=−η, f2=L/rt ---------------------------------
    def _XF(self) -> tuple:
        X = np.array([r["x"] for r in self.rows if r["status"] == "PASS"])
        F = np.array([[-r["eta_cf"], r["x"][DV_ORDER.index("L_over_rt")]]
                      for r in self.rows if r["status"] == "PASS"])
        return X, F

    # -- メインループ -----------------------------------------------------------
    def run(self, n_doe: int, n_iter: int, batch: int) -> None:
        done = {r["tag"] for r in self.rows}
        cand = lhs(self.bounds, 3 * n_doe, seed=self.seed)
        doe = [x for x in cand if self.feasible(x)][:n_doe]
        print(f"DOE {len(doe)} 点 (feasible {len(doe)}/{3*n_doe} 生成中から先頭採用)",
              flush=True)
        for i, x in enumerate(doe):
            tag = f"doe_{i:03d}"
            if tag not in done:
                self.evaluate(x, tag)
        for it in range(n_iter):
            X, F = self._XF()
            if len(X) < 4:
                raise RuntimeError("PASS 評価が少なすぎる (<4) — ゲート/レシピを確認")
            hv = hypervolume2d(F, self.ref)
            print(f"--- iter {it}: PASS {len(X)} 点, HV={hv:.5f} ---", flush=True)
            krg = KrigingSet(self.bounds).fit(X, F[:, :1])
            surro = _F1KrgF2Exact(krg, DV_ORDER.index("L_over_rt"))
            xs = propose_infill(surro, X, F, self.ref, self.bounds, n_infill=batch,
                                seed=self.seed + 1000 + it, feasible=self.feasible)
            for j, x in enumerate(xs):
                tag = f"inf_{it:02d}_{j}"
                if tag not in done:
                    self.evaluate(x, tag)
        X, F = self._XF()
        mask = nondominated_mask(F)
        pareto = [{"x": X[i].tolist(), "eta_cf": float(-F[i, 0]),
                   "L_over_rt": float(F[i, 1])} for i in np.where(mask)[0]]
        pareto.sort(key=lambda r: r["L_over_rt"])
        summary = {"n_eval": len(self.rows),
                   "n_pass": int(len(X)), "hv": hypervolume2d(F, self.ref),
                   "ref": self.ref, "pareto": pareto}
        (self.dir / "pareto.json").write_text(
            json.dumps(summary, indent=1, ensure_ascii=False))
        print(json.dumps(summary, indent=1, ensure_ascii=False), flush=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="③ベル TOP dv の MOO キャンペーン")
    ap.add_argument("problem")
    ap.add_argument("campaign_dir")
    ap.add_argument("--seed-res", required=True, help="初期 warm start の収束 res_*.h5")
    ap.add_argument("--n-doe", type=int, default=24)
    ap.add_argument("--n-iter", type=int, default=13)
    ap.add_argument("--batch", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ref", type=float, nargs=2, default=(-0.90, 13.0))
    a = ap.parse_args(argv)
    Campaign(a.problem, a.campaign_dir, a.seed_res, ref=tuple(a.ref),
             seed=a.seed).run(a.n_doe, a.n_iter, a.batch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
