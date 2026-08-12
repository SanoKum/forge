"""pymoo NSGA-II (サロゲート上の安価な探索) + EHVI infill 提案 (要 .venv-opt)。

infill 候補プール = NSGA-II 最終個体群 (活用) ∪ LHS 探索点 (探査) とし、
全候補を EHVI でランク付けして上位 n_infill 点を返す。決定性: pymoo の
seed と LHS の seed を固定 (同一入力 → 同一提案点列を単体テストで保証)。
"""
from __future__ import annotations

import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize

from .doe import lhs
from .ehvi import ehvi2d, nondominated_mask


class _SurrogateProblem(Problem):
    """サロゲート平均を目的値とする pymoo 問題。"""

    def __init__(self, surro, bounds) -> None:
        bounds = np.asarray(bounds, dtype=float)
        super().__init__(n_var=bounds.shape[0], n_obj=2,
                         xl=bounds[:, 0], xu=bounds[:, 1])
        self._surro = surro

    def _evaluate(self, X, out, *args, **kwargs):
        mu, _ = self._surro.predict(X)
        out["F"] = mu


def propose_infill(surro, X_eval, F_eval, ref, bounds, n_infill: int = 1,
                   seed: int = 0, pop_size: int = 64, n_gen: int = 60,
                   n_explore: int = 1024, min_dist: float = 5e-3,
                   feasible=None) -> np.ndarray:
    """EHVI 最大の未評価点を n_infill 個提案する。

    min_dist: 正規化座標での重複排除半径 (提案同士と評価済み点の両方に適用)。
    feasible: 実スケール x を受けて bool を返す事前フィルタ (幾何整合など)。
    """
    bounds = np.asarray(bounds, dtype=float)
    lo, span = bounds[:, 0], bounds[:, 1] - bounds[:, 0]
    res = minimize(_SurrogateProblem(surro, bounds), NSGA2(pop_size=pop_size),
                   ("n_gen", n_gen), seed=seed, verbose=False)
    cand = np.atleast_2d(res.X)
    cand = np.vstack([cand, lhs(bounds, n_explore, seed=seed + 7)])

    mu, sg = surro.predict(cand)
    F_eval = np.asarray(F_eval, dtype=float)
    front = F_eval[nondominated_mask(F_eval)]
    scores = np.array([ehvi2d(mu[i], sg[i], front, ref) for i in range(len(cand))])

    Xn_eval = (np.atleast_2d(np.asarray(X_eval, dtype=float)) - lo) / span
    picked: list = []
    for i in np.argsort(-scores):
        if feasible is not None and not feasible(cand[i]):
            continue
        xn = (cand[i] - lo) / span
        if np.min(np.linalg.norm(Xn_eval - xn, axis=1)) < min_dist:
            continue
        if picked and min(float(np.linalg.norm(xn - p)) for p in picked) < min_dist:
            continue
        picked.append(xn)
        if len(picked) >= n_infill:
            break
    return lo + span * np.asarray(picked, dtype=float).reshape(-1, bounds.shape[0])
