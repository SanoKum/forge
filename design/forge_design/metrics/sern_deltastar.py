"""⑤ SERN の δ* 抽出 (S5 一発補正の入力)。

node 構造格子 run から、上バンドの各 x station の列 (j 方向) に沿って ρu 分布を取り、
ランプ (j = nj_top−1 側) とカウル内面 (j = 0 側、x ≤ L_cowl) の排除厚さ
    δ* = ∫ (1 − ρu/(ρu)_e) dn,  縁 = 壁から見て最初に edge_frac·max(ρu) に達する点 (探索窓は列の一部)
を計算する。列は鉛直なので dn = dy·cosθ_w で法線長に直す。cell run は非対応 (VALUE 長 ≠ 節点数)。
戻り値は物理単位 [m] の (x, δ*) 表 (runner_sern.prepare の wall_offset にそのまま渡せる)。
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np


def _structured_upper(run_dir: Path):
    info = json.loads((run_dir / "prepare_info.json").read_text())
    m = info["mesh"]; ni, njt, njb, ite = m["ni"], m["nj_top"], m["nj_bot"], m["i_te"]
    with h5py.File(run_dir / "sern.h5") as f:
        cc = f["CELLS/centCoords"][:].reshape(-1, 3)
    N_low = ni * njb; N_up = ni * (njt - 1)

    def up(i, j):
        if j == 0:
            return N_low + N_up + i if i < ite else i * njb + (njb - 1)
        return N_low + i * (njt - 1) + (j - 1)
    idx = np.array([[up(i, j) for j in range(njt)] for i in range(ni)])
    return info, idx, cc


def deltastar_sern(run_dir, res_h5=None, edge_frac: float = 0.99, search_frac: float = 0.5, min_pts: int = 4) -> dict:
    run_dir = Path(run_dir)
    if res_h5 is None:
        res_h5 = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    info, idx, cc = _structured_upper(run_dir)
    with h5py.File(res_h5) as f:
        ro = f["VALUE/ro"][:]; rux = f["VALUE/roUx"][:]
    if len(ro) != len(cc):
        raise ValueError("node run でない (VALUE 長 ≠ 節点数): δ* 抽出は node 専用")
    H = info["H_m"]; ni, njt = idx.shape; ite = info["mesh"]["i_te"]
    X = cc[idx, 0]; Y = cc[idx, 1]; Q = rux[idx]
    out = {"ramp": [], "cowl": []}
    for i in range(1, ni):
        x = float(X[i, 0])
        # --- ramp: 壁 j=njt-1 から下へ
        yw = Y[i, -1]; th = np.arctan2(Y[i, -1] - Y[i - 1, -1], X[i, -1] - X[i - 1, -1])
        n = (yw - Y[i, ::-1]) * np.cos(th); q = Q[i, ::-1]
        out["ramp"].append((x, _dstar(n, q, edge_frac, search_frac, min_pts)))
        # --- cowl 内面: 壁 j=0 から上へ (x ≤ L_cowl のみ)
        if i <= ite and x >= 0.0:
            yw = Y[i, 0]; th = np.arctan2(Y[i, 0] - Y[i - 1, 0], X[i, 0] - X[i - 1, 0])
            n = (Y[i, :] - yw) * np.cos(th); q = Q[i, :]
            out["cowl"].append((x, _dstar(n, q, edge_frac, search_frac, min_pts)))
    res = {k: np.asarray(v) for k, v in out.items()}
    for k in res:
        good = np.isfinite(res[k][:, 1])
        res[k] = res[k][good]
    res["H_m"] = H; res["res_file"] = str(res_h5)
    return res


def _dstar(n, q, edge_frac, search_frac, min_pts) -> float:
    n = np.asarray(n, dtype=float); q = np.asarray(q, dtype=float)
    k = max(int(len(n) * search_frac), min_pts + 1)
    qs = q[:k]
    qe = float(np.max(qs))
    if qe <= 0:
        return np.nan
    ie = int(np.argmax(qs >= edge_frac * qe))
    if ie < min_pts:
        return np.nan
    m = slice(0, ie + 1)
    return float(np.trapezoid(1.0 - q[m] / qe, n[m]))


def smooth_table(tbl, win: int = 9):
    """(x, δ*) 表の移動平均 (端は縮小窓)。"""
    tbl = np.asarray(tbl, dtype=float)
    if len(tbl) < 3:
        return tbl
    d = tbl[:, 1].copy(); out = d.copy()
    h = win // 2
    for i in range(len(d)):
        a, b = max(0, i - h), min(len(d), i + h + 1)
        out[i] = np.mean(d[a:b])
    return np.column_stack([tbl[:, 0], out])


def write_offset_json(res: dict, path, win: int = 9) -> dict:
    wo = {"ramp": smooth_table(res["ramp"], win).tolist(), "cowl": smooth_table(res["cowl"], win).tolist()}
    Path(path).write_text(json.dumps(wo, indent=1))
    return wo
