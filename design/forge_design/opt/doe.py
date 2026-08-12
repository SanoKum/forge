"""DOE: Latin Hypercube (決定的 seed — 同一 seed → 同一点列を単体テストで保証)。"""
from __future__ import annotations

import numpy as np
from scipy.stats import qmc


def lhs(bounds, n: int, seed: int) -> np.ndarray:
    """bounds (d,2) の LHS n 点。戻り (n,d)、実スケール。"""
    bounds = np.asarray(bounds, dtype=float)
    sampler = qmc.LatinHypercube(d=bounds.shape[0], seed=int(seed))
    u = sampler.random(int(n))
    return qmc.scale(u, bounds[:, 0], bounds[:, 1])
