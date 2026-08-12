"""サロゲート支援 MOO ループ (Phase 2 — 親計画 §4.2 案 A: pymoo + SMT)。

- `ehvi`: 2 目的 EHVI 解析式・hypervolume・非劣解 (numpy/scipy のみ)
- `doe`: LHS DOE (決定的 seed)
- `surrogate`: SMT KRG ラッパ (要 .venv-opt)
- `moo`: pymoo NSGA-II サロゲート探索 + EHVI infill 提案 (要 .venv-opt)
"""
from .doe import lhs
from .ehvi import ehvi2d, hypervolume2d, nondominated_mask

__all__ = ["lhs", "ehvi2d", "hypervolume2d", "nondominated_mask"]
