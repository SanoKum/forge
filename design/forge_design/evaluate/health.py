"""CFD 健全性ゲート (Phase 3 A3) — AGENTS の収束・準定常ルールの実体化。

**既存ルールを弱めない**のが設計方針:

- `check_convergence.py` の VERDICT を**必須**として記録する。`NOT CONVERGED` を
  独自基準で「収束相当」に読み替えることはしない (2026-08-14 のレビュー指摘)。
- 加えて `check_quasisteady.py` (報告量の定常性)、残差の**絶対レベル**、
  残差の**時系列トレンド**を併用する。warm restart 継続では「低下桁数」が小さく
  出るのは自然だが、それは判定を緩める理由にはならない — 事実として両方を残し、
  未収束なら**未収束と明記**する。
- 帰還ループはこの結果を**記録**し、DIVERGED / NaN のときだけパスを止める。
  それ以外は「未収束のまま進んだ」という事実を履歴に残して継続する
  (収束を主張する報告には VERDICT をそのまま貼れる形にしておく)。
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
            "usable": bool(cver != "DIVERGED" and not nonfinite)}
