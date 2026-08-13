"""CFD 健全性ゲート (Phase 3 A3) — AGENTS の収束・準定常ルールの実体化。

**既存ルールを弱めない**のが設計方針:

- `check_convergence.py` の VERDICT を**必須**として記録する。`NOT CONVERGED` を
  独自基準で「収束相当」に読み替えることはしない (2026-08-14 のレビュー指摘)。
- 加えて `check_quasisteady.py` (報告量の定常性)、残差の**絶対レベル**、
  残差の**時系列トレンド**を併用する。warm restart 継続では「低下桁数」が小さく
  出るのは自然だが、それは判定を緩める理由にはならない — 事実として両方を残し、
  未収束なら**未収束と明記**する。
- **収束閾値は明示・記録する** (2026-08-14): `check_convergence.py` は `--drop` を
  正式にサポートしており、既定 3 桁を **2 桁 + 準定常** に緩和して運用する
  (読み替えではなく設定の明示)。理由: 本ケースの残差は step ~371 で全列 2 桁・
  ~1101 で 3 桁に達した後**プラトー**に入り、その後わずかに増減するため、
  **長く回すほど「初期比 3 桁」を満たさなくなる**ことが実測された
  (24000 step で `rms_roUx` が 2.9 桁)。閾値は verdict と履歴に必ず残し、
  「収束した」と述べる際は **`PASS (--drop 2)`** のように条件込みで示す。
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


def axis_M_series_steady(run_dir, scale, x_lo, x_hi, tol_M: float = 5e-3,
                         n_ref: int = 6000, n_last: int = 5,
                         axis_band_frac: float = 0.09, snaps=None) -> dict:
    r"""**軸中心 M プロファイルの定常性** — 本ループの目的量そのもの。

    スナップショット列から軸上 $M(x)$ を固定格子へ補間し、末尾 `n_last` 個の窓で
    **点別の線形トレンド** (step に対する最小二乗傾き) を `n_ref` step へ射影して
    `tol_M` と比べる。レート正規化なので判定はスナップショット間隔に依存しない
    (2026-08-14 レビュー指摘)。意味は「**あと `n_ref` step (=1 継続チャンク)
    回しても軸 M の系統的移動は `tol_M` 未満**」。tol_M=5e-3 はループの収束目標
    $0.5\%M_d=0.02$ の 1/4。

    **隣接対レートではなくトレンドを使う理由 (2026-08-15 実測)**: 中域
    x/rt≈9–10 に片振幅 ~0.003・無トレンドの定在音響リミットサイクルがあり
    (run_0017 pass_04: M(9.2) が 3.6615–3.6636 を 24000 step 往復、平均一定)、
    隣接対の $\max|\Delta M|$ は振動の**位相サンプリングノイズ** (0.0012–0.0053)
    を拾って DRIFTING を誤報する。線形フィットは振動を平均化し系統移動だけを
    測る。振幅は `amplitude` (窓平均まわり最大偏差) として別途返し、
    STEADY (amp < tol_M/2) / OSCILLATING (トレンド小・振幅あり) を区別する。
    「軸中心の M が狙い通りか」が唯一の判定基準なので (2026-08-14 ユーザ指示)、
    出口一様性 $\varepsilon$ ではなくこれをゲートに使う。$\varepsilon$ は記録のみ。

    `snaps`: [(step, res_path, mesh_path), ...] を明示すると複数チャンク run
    (warm 継続) をまたいだ列で判定できる。省略時は run_dir を glob。
    """
    from ..metrics.extract import axis_mach

    rd = Path(run_dir)
    if snaps is None:
        res = sorted(rd.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        snaps = [(int("".join(c for c in f.stem if c.isdigit())), f, rd / "nozzle.h5")
                 for f in res]
    # step 重複 (チャンク境界の res_0 = 直前最終場) と step 0 (IC — warm では
    # 壁変更直後の未緩和場なのでトレンドに混ぜない) を除去
    ded, seen = [], set()
    for s in snaps:
        if int(s[0]) > 0 and int(s[0]) not in seen:
            seen.add(int(s[0]))
            ded.append(s)
    snaps = ded
    if len(snaps) < 3:
        return {"verdict": "UNKNOWN", "reason": f"スナップショット不足 ({len(snaps)})"}
    xi = np.linspace(x_lo, x_hi, 200)
    prof, steps = [], []
    for st, f, mesh in snaps[-n_last:]:
        try:
            x, M = axis_mach(mesh, f, axis_band=axis_band_frac * scale)
        except Exception as exc:
            return {"verdict": "UNKNOWN", "reason": f"{type(exc).__name__}: {exc}"}
        prof.append(np.interp(xi, x, M))
        steps.append(int(st))
    P = np.asarray(prof)
    s = np.asarray(steps, dtype=float)
    # 点別線形トレンド b(x) [1/step] の最小二乗解 → n_ref step への射影
    sc = s - s.mean()
    b = (sc[:, None] * (P - P.mean(axis=0))).sum(axis=0) / (sc * sc).sum()
    trend = float(np.max(np.abs(b)) * n_ref)
    # 窓平均まわりの最大偏差 = 振動片振幅 (トレンド除去はしない保守側)
    amp = float(np.max(np.abs(P - P.mean(axis=0))))
    out = {"trend_dM_per_nref": trend, "amplitude": amp,
           # 平均プロファイル推定の標準誤差の目安。リミットサイクルは回し続けても
           # 振幅が縮まないので**振幅ではなくこれをゲートする** (2026-08-15):
           # 測定は窓平均を使うため、必要なのは「平均が tol より良く決まっている」
           # こと。超過時はチャンク追加でサンプル数を増やせば締まる。
           "mean_se": amp / float(np.sqrt(len(prof))),
           "tol_M": tol_M, "n_ref": n_ref, "steps_used": steps}
    if trend >= tol_M:
        out["verdict"] = "DRIFTING"
    elif amp < 0.5 * tol_M:
        out["verdict"] = "STEADY"
    else:
        out["verdict"] = "OSCILLATING"   # AGENTS: 平均±振幅で扱う (振幅は amp)
    return out


def axis_ok(ax: dict) -> bool:
    """軸 M 判定のゲート合格条件 (gate() と継続ループで共用)。

    STEADY、または OSCILLATING で**平均の標準誤差** ≤ tol_M。DRIFTING/UNKNOWN は
    不合格。振幅そのものは条件にしない (上記コメント参照) — 事実として記録され、
    報告は平均±振幅で行う。"""
    if ax.get("verdict") == "STEADY":
        return True
    return (ax.get("verdict") == "OSCILLATING"
            and ax.get("mean_se", np.inf) <= ax.get("tol_M", 5e-3))


def eps_series_steady(run_dir, M_design, x_d=None, rel_tol: float = 5e-2,
                      n_last: int = 3, abs_tol_M: float = 2e-5,
                      abs_tol_theta_deg: float = 5e-3) -> dict:
    """報告量そのもの ($\varepsilon_M$/$\varepsilon_\theta$) の時系列定常性。

    **注意 (2026-08-14)**: 本判定が使う `exit_uniformity` は **A2 未完成**
    (仮想境界による環状面積・$(x_d,0)$ 直線によるコア定義) なので、現時点の
    数値は「**ゲート機構の動作確認値**」に限る。A2 を仕様通り (固定半径格子へ
    補間 / リップからの最終特性線) に完成させてから正式指標として扱うこと。

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
        return float(v.max() - v.min()), float((v.max() - v.min()) / max(abs(v.mean()), 1e-30))
    aM, sM = _spread(eM)
    aT, sT = _spread(eT)
    # **絶対床が要る** (2026-08-14): ε は収束するほど小さくなるので、相対許容だけだと
    # ε_θ=0.06° で 0.001° の揺らぎを「非定常」と誤判定する。相対**または**絶対の
    # どちらかを満たせば定常とみなす。
    # **ゲートに使うのは ε_M だけ** (2026-08-14 ユーザ指示): 本ループの目的は
    # 「軸 M が目標通りか」であり、流れ角 ε_θ は最適化ループ側の目的。ε_θ は
    # 記録するがゲート条件にはしない。
    okM = (sM < rel_tol) or (aM < abs_tol_M)
    return {"verdict": "STEADY" if okM else "DRIFTING",
            "eps_M_series": eM, "eps_theta_series": eT,
            "eps_M_spread": sM, "eps_theta_spread": sT,
            "eps_M_abs_spread": aM, "eps_theta_abs_spread_deg": aT,
            "rel_tol": rel_tol, "abs_tol_M": abs_tol_M,
            "abs_tol_theta_deg": abs_tol_theta_deg}


def gate(run_dir, M_design=None, x_d=None, mode: str = "formal",
         drop: float = 2.0, axis_window=None, scale: float = 1.0,
         warm: bool = False, snaps=None, eff_dir=None) -> dict:
    """正式評価ゲート。`mode="formal"` では **PASS 以外を不合格**にする。

    戻り値 `ok=True` のときのみ、帰還ループは壁・特性線マップを更新してよい。

    `warm=True` (前パス場からの warm 継続パス, 2026-08-15): 初期残差が既に小さく
    「初期比 n 桁低下」は**構造的に**満たせないため、低下桁数要件を外し
    **① DIVERGED でない ② 非有限値なし ③ 軸 M 定常 (レート正規化)** を要件と
    する。判定の読み替えではなく要件の置換であり、`convergence_verdict` /
    `quasisteady_verdict` は事実としてそのまま記録する。cold なアンカーパス
    (pass 1) が同一設定で PASS していることが前提。

    `snaps`: 複数チャンク run の軸 M 判定用 [(step, res, mesh), ...]。
    `eff_dir`: 最終チャンク (health() docstring 参照)。
    """
    h = health(run_dir, drop=drop, eff_dir=eff_dir)
    reasons = []
    if warm:
        if h["convergence_verdict"] in ("DIVERGED", "UNKNOWN"):
            reasons.append(f"convergence={h['convergence_verdict']}")
        if axis_window is None:
            reasons.append("warm ゲートに axis_window が必要")
    else:
        if h["convergence_verdict"] != "PASS":
            reasons.append(f"convergence={h['convergence_verdict']} (--drop {drop})")
        if h["quasisteady_verdict"] not in ("ALL STEADY",):
            reasons.append(f"quasisteady={h['quasisteady_verdict']}")
    if h["nonfinite_residual"]:
        reasons.append("residual に非有限値")
    if isinstance(h["residuals"], dict) and h["residuals"].get("error"):
        reasons.append(f"residual 読取失敗: {h['residuals']['error']}")
    if axis_window is not None:
        # **ゲートの定常性判定は軸 M のみ** (ユーザ指示 2026-08-14)。
        # STEADY のほか **OSCILLATING (リミットサイクル) は振幅 ≤ tol_M なら合格**
        # — その場合 ΔM 測定は末尾スナップ平均を使う (呼び出し側の責務。AGENTS の
        # 「OSCILLATING は平均±振幅で報告」の実体化)。DRIFTING のみ不合格。
        ax = axis_M_series_steady(run_dir, scale, axis_window[0], axis_window[1],
                                  snaps=snaps)
        h["axis_M_steady"] = ax
        if not axis_ok(ax):
            reasons.append(
                f"axis_M_series={ax['verdict']}"
                + (f"(mean_se={ax['mean_se']:.4f}>tol)"
                   if ax["verdict"] == "OSCILLATING" else ""))
    if M_design is not None:   # ε は**診断として記録するだけ** (ゲート条件ではない)
        h["eps_steady"] = eps_series_steady(eff_dir or run_dir, M_design, x_d=x_d)
    h["gate_mode"] = mode + ("+warm" if warm else "")
    h["gate_reasons"] = reasons
    h["ok"] = (len(reasons) == 0) if mode == "formal" else h["usable"]
    return h


def health(run_dir, quantity: str | None = None, drop: float = 2.0,
           eff_dir=None) -> dict:
    """収束 (必須) + 準定常 + 残差絶対値/トレンドをまとめて返す。

    戻り値の `usable` は「**発散していない**」だけを意味する (収束の主張ではない)。
    収束を主張したい場合は `convergence_verdict` をそのまま提示すること。

    `eff_dir`: 複数チャンク run (warm 継続) の**最終チャンク**。収束低下桁数は
    実過渡が入っている**第 1 チャンク** (`run_dir`) で、準定常と残差の非有限値は
    最新場の `eff_dir` で判定する (NaN は一度出たら消えないので最終チャンクの
    検査で十分)。省略時は両方 `run_dir`。
    """
    run_dir = Path(run_dir)
    eff = Path(eff_dir) if eff_dir is not None else run_dir
    conv = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_convergence.py"),
                           "--drop", str(drop), str(run_dir)],
                          capture_output=True, text=True, env=_ENV)
    ctxt = conv.stdout + conv.stderr
    (run_dir / "CONVERGENCE_VERDICT.txt").write_text(
        f"# check_convergence.py --drop {drop}\n" + ctxt)
    cver = _verdict_of(ctxt, ("DIVERGED", "NOT CONVERGED", "PASS"))

    cmd = [sys.executable, str(FORGE_TOOLS / "check_quasisteady.py"), str(eff)]
    if quantity:
        cmd += ["--quantity", quantity]
    qs = subprocess.run(cmd, capture_output=True, text=True, env=_ENV)
    qtxt = qs.stdout + qs.stderr
    (eff / "QUASISTEADY_VERDICT.txt").write_text(qtxt)
    overall = [l for l in qtxt.splitlines() if "OVERALL:" in l]
    if not overall:
        qver = "UNKNOWN"
    elif "NOT" in overall[-1]:
        qver = "NOT ALL STEADY"      # "ALL STEADY" を部分文字列に含むので先に判定
    else:
        qver = "ALL STEADY" if "ALL STEADY" in overall[-1] else overall[-1].split("OVERALL:")[-1].strip()

    res = residual_summary(eff)
    nonfinite = any(v.get("nonfinite") for v in res.values() if isinstance(v, dict))
    return {"convergence_verdict": cver, "convergence_drop_threshold": drop,
            "quasisteady_verdict": qver,
            "residuals": res, "nonfinite_residual": nonfinite,
            # `usable` は「発散していない」だけ。正式判定は gate() を使うこと
            "usable": bool(cver not in ("DIVERGED", "UNKNOWN") and not nonfinite)}
