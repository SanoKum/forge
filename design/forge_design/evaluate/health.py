"""CFD 健全性ゲート (Phase 3 A3) — AGENTS の収束・準定常ルールの実体化。

**既存ルールを弱めない**のが設計方針:

- `check_convergence.py` の VERDICT を**必須**として記録する。`NOT CONVERGED` を
  独自基準で「収束相当」に読み替えることはしない (2026-08-14 のレビュー指摘)。
- 加えて `check_quasisteady.py` (報告量の定常性)、残差の**絶対レベル**、
  残差の**時系列トレンド**を併用する。warm restart 継続では「低下桁数」が小さく
  出るのは自然だが、それは判定を緩める理由にはならない — 事実として両方を残し、
  未収束なら**未収束と明記**する。
- **正式評価経路では PASS のみを採る** (`gate()`): `check_convergence` が PASS で
  なければ壁もマップも更新しない。`UNKNOWN`・ツール失敗も不合格。これは
  `methods/design/overview.md` の「PASS した場のみ metrics に渡す」と整合させる
  ための実体化であり、「未収束と記録しただけで帰還を続ける」ことは**しない**
  (2026-08-14 レビュー指摘)。
- 探索的に未収束のまま回したい場合は `mode="exploratory"` を明示し、**正式経路とは
  別の run ディレクトリ**に分離する (結果を正式指標として報告しない)。
- 準定常は既定量 (shock/machmax/pmax) に加え、**報告量そのもの** ($\varepsilon_M$/
  $\varepsilon_\theta$) の時系列を `eps_series_steady()` で別途判定する。
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np

from .runner import FORGE_TOOLS, _ENV


def _verdict_of(txt: str, keys) -> str:
    for k in keys:
        if k in txt:
            return k
    return "UNKNOWN"


def residual_summary(run_dir) -> dict:
    """residual_history.csv から各列の最終絶対値・最小到達値・末尾トレンドを返す。"""
    f = Path(run_dir) / "residual_history.csv"
    if not f.exists():
        return {"error": "residual_history.csv なし"}
    rows = [l.split(",") for l in f.read_text().splitlines() if l.strip()]
    hdr = [h.strip() for h in rows[0]]
    # `phase` など非数値列があるので rms_* 列だけを個別にパースする
    out = {}
    for j, name in enumerate(hdr):
        if not name.startswith("rms_") or name.startswith("rms_dq_"):
            continue
        vals = []
        for r in rows[1:]:
            if j >= len(r):
                continue
            try:
                vals.append(float(r[j]))
            except ValueError:
                vals.append(np.nan)
        c = np.asarray(vals, dtype=float)
        if not np.any(np.isfinite(c)) or np.allclose(c, 0.0):
            continue
        n = len(c)
        tail = c[max(0, n - n // 5):]
        out[name] = {
            "final": float(c[-1]),
            "min": float(np.nanmin(c)),
            "nonfinite": bool(np.any(~np.isfinite(c))),
            # 末尾 20% の傾き (log10/step): 負なら低下継続、正なら上昇
            "tail_slope_log10": (float(np.polyfit(np.arange(len(tail)),
                                                  np.log10(np.maximum(tail, 1e-300)), 1)[0])
                                 if len(tail) > 4 else None),
        }
    return out


def eps_series_steady(run_dir, M_design, x_d=None, rel_tol: float = 2e-2,
                      n_last: int = 3) -> dict:
    """報告量そのもの ($\varepsilon_M$/$\varepsilon_\theta$) の時系列定常性。

    `check_quasisteady.py` の既定量 (shock/machmax/pmax) は報告量ではないので、
    出力スナップショット列から $\varepsilon$ を直接計算して末尾 `n_last` 個の
    相対変動が `rel_tol` 未満かを見る (2026-08-14 レビュー指摘)。"""
    from ..metrics.extract import exit_uniformity

    run_dir = Path(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if len(res) < n_last:
        return {"verdict": "UNKNOWN", "reason": f"スナップショット不足 ({len(res)})"}
    eM, eT = [], []
    for f in res[-n_last:]:
        try:
            e = exit_uniformity(run_dir / "nozzle.h5", f, M_design, x_d=x_d)
        except Exception as exc:
            return {"verdict": "UNKNOWN", "reason": f"{type(exc).__name__}: {exc}"}
        eM.append(e["eps_M_rms"])
        eT.append(e["eps_theta_max_deg"])
    def _spread(v):
        v = np.asarray(v, dtype=float)
        return float((v.max() - v.min()) / max(abs(v.mean()), 1e-30))
    sM, sT = _spread(eM), _spread(eT)
    ok = sM < rel_tol and sT < rel_tol
    return {"verdict": "STEADY" if ok else "DRIFTING",
            "eps_M_series": eM, "eps_theta_series": eT,
            "eps_M_spread": sM, "eps_theta_spread": sT, "rel_tol": rel_tol}


def gate(run_dir, M_design=None, x_d=None, mode: str = "formal") -> dict:
    """正式評価ゲート。`mode="formal"` では **PASS 以外を不合格**にする。

    戻り値 `ok=True` のときのみ、帰還ループは壁・特性線マップを更新してよい。
    """
    h = health(run_dir)
    reasons = []
    if h["convergence_verdict"] != "PASS":
        reasons.append(f"convergence={h['convergence_verdict']}")
    if h["quasisteady_verdict"] not in ("ALL STEADY",):
        reasons.append(f"quasisteady={h['quasisteady_verdict']}")
    if h["nonfinite_residual"]:
        reasons.append("residual に非有限値")
    if isinstance(h["residuals"], dict) and h["residuals"].get("error"):
        reasons.append(f"residual 読取失敗: {h['residuals']['error']}")
    if M_design is not None:
        es = eps_series_steady(run_dir, M_design, x_d=x_d)
        h["eps_steady"] = es
        if es["verdict"] != "STEADY":
            reasons.append(f"eps_series={es['verdict']}")
    h["gate_mode"] = mode
    h["gate_reasons"] = reasons
    h["ok"] = (len(reasons) == 0) if mode == "formal" else h["usable"]
    return h


def health(run_dir, quantity: str | None = None) -> dict:
    """収束 (必須) + 準定常 + 残差絶対値/トレンドをまとめて返す。

    戻り値の `usable` は「**発散していない**」だけを意味する (収束の主張ではない)。
    収束を主張したい場合は `convergence_verdict` をそのまま提示すること。
    """
    run_dir = Path(run_dir)
    conv = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_convergence.py"),
                           str(run_dir)], capture_output=True, text=True, env=_ENV)
    ctxt = conv.stdout + conv.stderr
    (run_dir / "CONVERGENCE_VERDICT.txt").write_text(ctxt)
    cver = _verdict_of(ctxt, ("DIVERGED", "NOT CONVERGED", "PASS"))

    cmd = [sys.executable, str(FORGE_TOOLS / "check_quasisteady.py"), str(run_dir)]
    if quantity:
        cmd += ["--quantity", quantity]
    qs = subprocess.run(cmd, capture_output=True, text=True, env=_ENV)
    qtxt = qs.stdout + qs.stderr
    (run_dir / "QUASISTEADY_VERDICT.txt").write_text(qtxt)
    overall = [l for l in qtxt.splitlines() if "OVERALL:" in l]
    if not overall:
        qver = "UNKNOWN"
    elif "NOT" in overall[-1]:
        qver = "NOT ALL STEADY"      # "ALL STEADY" を部分文字列に含むので先に判定
    else:
        qver = "ALL STEADY" if "ALL STEADY" in overall[-1] else overall[-1].split("OVERALL:")[-1].strip()

    res = residual_summary(run_dir)
    nonfinite = any(v.get("nonfinite") for v in res.values() if isinstance(v, dict))
    return {"convergence_verdict": cver, "quasisteady_verdict": qver,
            "residuals": res, "nonfinite_residual": nonfinite,
            # `usable` は「発散していない」だけ。正式判定は gate() を使うこと
            "usable": bool(cver not in ("DIVERGED", "UNKNOWN") and not nonfinite)}
