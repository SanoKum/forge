"""SMT KRG ラッパ: 目的ごとに独立 Kriging (要 .venv-opt)。

- 入力は bounds で [0,1]^d に正規化してから学習 (異方性 θ の桁を揃える)。
- SMT KRG は既定で決定的 (単体テストで確認済み)・学習点を補間する。
- predict は (mean, std) を返す (EHVI の入力規約)。
"""
from __future__ import annotations

import numpy as np
from smt.surrogate_models import KRG


class KrigingSet:
    """m 目的ぶんの独立 KRG の束。"""

    def __init__(self, bounds) -> None:
        self.bounds = np.asarray(bounds, dtype=float)
        self.models: list = []

    def _norm(self, X) -> np.ndarray:
        lo, hi = self.bounds[:, 0], self.bounds[:, 1]
        return (np.atleast_2d(np.asarray(X, dtype=float)) - lo) / (hi - lo)

    def fit(self, X, Y) -> "KrigingSet":
        Xn = self._norm(X)
        Y = np.atleast_2d(np.asarray(Y, dtype=float))
        self.models = []
        for j in range(Y.shape[1]):
            m = KRG(print_global=False, theta0=[1e-1] * Xn.shape[1])
            m.set_training_values(Xn, Y[:, j])
            m.train()
            self.models.append(m)
        return self

    def predict(self, X) -> tuple:
        """(mean (n,m), std (n,m))。"""
        Xn = self._norm(X)
        mu = np.column_stack([m.predict_values(Xn).ravel() for m in self.models])
        var = np.column_stack([m.predict_variances(Xn).ravel() for m in self.models])
        return mu, np.sqrt(np.maximum(var, 0.0))
