#!/usr/bin/env python3
"""`moc_diagnostics` の折れ線交差高速化テスト。

通常実行は手作り・固定 seed ランダム折れ線について、AABB sweep と旧全線分対行列の判定を
照合する。`--benchmark` では case/42 相当の A knot 実 MOC 網を 600/1200/2400 点で生成し、
サンプル特性線交差数の一致と時間を測る。
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.geometry.moc_diagnostics import (                          # noqa: E402
    _polylines_cross_aabb, _polylines_cross_vec, _sample_characteristic_crossings)
from forge_design.geometry.moc_inverse import cplus_lines                    # noqa: E402


FAIL = 0


def check(name: str, cond: bool) -> None:
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def legacy_sample_crossings(levels: np.ndarray, n_sample: int) -> int:
    """高速化前と同じ全線分対行列による参照計数。"""
    n = levels.shape[1]
    ms = np.unique(np.linspace(0, n - 1, min(n_sample, n)).astype(int))
    lines = []
    for m in ms:
        line = cplus_lines(levels, int(m))
        if len(line) >= 2:
            lines.append(line[:, :2])
    boxes = [(line[:, 0].min(), line[:, 0].max(), line[:, 1].min(), line[:, 1].max())
             for line in lines]
    count = 0
    for a in range(len(lines)):
        xa0, xa1, ra0, ra1 = boxes[a]
        for b in range(a + 2, len(lines)):
            xb0, xb1, rb0, rb1 = boxes[b]
            if xa1 < xb0 or xb1 < xa0 or ra1 < rb0 or rb1 < ra0:
                continue
            if _polylines_cross_vec(lines[a], lines[b]):
                count += 1
    return count


# 手作りケース: 交差、非交差、端点接触、共線、x 逆向き、x 非単調。
manual = [
    (np.array([[0., 0.], [1., 1.]]), np.array([[0., 1.], [1., 0.]])),
    (np.array([[0., 0.], [1., 0.]]), np.array([[0., 1.], [1., 1.]])),
    (np.array([[0., 0.], [1., 1.]]), np.array([[1., 1.], [2., 0.]])),
    (np.array([[0., 0.], [2., 0.]]), np.array([[1., 0.], [3., 0.]])),
    (np.array([[2., 0.], [1., 1.], [0., 0.]]), np.array([[0., .8], [2., .8]])),
    (np.array([[0., 0.], [2., 1.], [1., 2.], [3., 3.]]),
     np.array([[0., 2.], [2., 1.5], [1., .5], [3., 0.]])),
]
for i, (p, q) in enumerate(manual):
    check(f"手作り折れ線 {i}: 旧実装と一致",
          _polylines_cross_aabb(p, q, chunk_size=2) == _polylines_cross_vec(p, q))

# ランダム walk (非単調) と x 単調折れ線を固定 seed で照合。
rng = np.random.default_rng(20260816)
n_mismatch = 0
for _ in range(500):
    p = np.cumsum(rng.normal(size=(int(rng.integers(2, 80)), 2)), axis=0)
    q = np.cumsum(rng.normal(size=(int(rng.integers(2, 80)), 2)), axis=0)
    n_mismatch += _polylines_cross_aabb(p, q, chunk_size=7) != _polylines_cross_vec(p, q)
for _ in range(500):
    np_, nq = rng.integers(2, 100, size=2)
    p = np.column_stack((np.sort(rng.uniform(-2, 2, np_)), rng.normal(size=np_)))
    q = np.column_stack((np.sort(rng.uniform(-2, 2, nq)), rng.normal(size=nq)))
    if rng.random() < 0.5:
        p = p[::-1]
    n_mismatch += _polylines_cross_aabb(p, q, chunk_size=11) != _polylines_cross_vec(p, q)
check(f"固定 seed ランダム 1000 組: 不一致 {n_mismatch}", n_mismatch == 0)


if "--benchmark" in sys.argv:
    from forge_design.gas import GasSemiPerfect
    from forge_design.geometry.axis_law import KnotQuinticAxisLaw
    from forge_design.geometry.moc_inverse import inverse_design
    from forge_design.geometry.transonic import HallThroat
    from forge_design.geometry.wall_axismach import area_ratio_isentropic

    composition = {"CO2": 0.1677, "H2O": 0.0858, "O2": 0.022, "N2": 0.7245}
    gas = GasSemiPerfect(composition, Tt=1550.0)
    ht = HallThroat(R=3.0, gamma=gas.gamma_star)
    xs_c, *_ = ht.throat_characteristic(n=41)
    x_A = float(xs_c[0])
    M_A, Mp_A, Mpp_A = ht.axis_anchor(x_A)
    law = KnotQuinticAxisLaw(x_A, 45.0, M_A, Mp_A, Mpp_A, 6.0, 2.5)
    r_F = float(np.sqrt(area_ratio_isentropic(6.0, gas)))
    x_end = law.x_E + 2.3 * r_F * np.sqrt(35.0)
    for n_axis in (600, 1200, 2400):
        result = inverse_design(
            ht, lambda x: float(law(max(x, x_A))), x_axis_end=x_end,
            n_axis=n_axis, n_start=41, gamma=gas, M_start=1.05,
            exit_mode="characteristic", x_E=law.x_E, M_d=6.0,
            start_line="throat_char", wall_mode="cplus", axis_dx0=0.03)
        levels = result["levels"]
        t0 = time.perf_counter()
        old = legacy_sample_crossings(levels, 40)
        t_old = time.perf_counter() - t0
        t0 = time.perf_counter()
        new = _sample_characteristic_crossings(levels, n_axis - 1, 40)
        t_new = time.perf_counter() - t0
        check(f"実 MOC n={n_axis}: crossing {old} == {new}", old == new)
        print(f"    old={t_old:.4f}s new={t_new:.4f}s speedup={t_old/max(t_new, 1e-12):.1f}x")


print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
