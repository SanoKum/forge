#!/usr/bin/env python3
"""逆 MOC マーチ (moc_inverse) の単体検証 (Phase 3 §6)。

1. 一様流: M 一定の Cauchy データ → 全点 M 一定・θ=0 (機械精度) + 壁流線が水平
2. 放射源流: 解析場の軸データ + 初期線 → 壁流線が源流レイ (直線) を回復
3. ノズル: Sauer + Bézier 目標 (Md=4, C2 出口) → 出口半径 ≈ 1D 面積比・出口壁角 ≈ 0
4. 決定性
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.geometry.moc_inverse import InverseMOC, inverse_design  # noqa: E402
from forge_design.geometry.moc_kernel import _Pt, pm_nu  # noqa: E402
from forge_design.geometry.bezier import MachBezier  # noqa: E402
from forge_design.geometry.transonic import SauerThroat  # noqa: E402
from forge_design.evaluate.ic import invert_area_ratio  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


G = 1.4

# --- 1. 一様流 -----------------------------------------------------------------
inv = InverseMOC(gamma=G, delta=1.0)
nuU = float(pm_nu(2.0, G))
init = ([_Pt(float(x), 0.0, 0.0, nuU, G) for x in np.linspace(4.0, 0.0, 41)[:-1]]
        + [_Pt(0.0, float(r), 0.0, nuU, G) for r in np.linspace(0.0, 1.0, 21)])
pts = inv.fill(init)
nus = np.array([p.nu for p in pts])
ths = np.array([p.th for p in pts])
check(f"一様流: ν 一定 (max err {np.max(np.abs(nus-nuU)):.1e})",
      float(np.max(np.abs(nus - nuU))) < 1e-12)
check(f"一様流: θ=0 (max {np.max(np.abs(ths)):.1e})", float(np.max(np.abs(ths))) < 1e-12)
wall = inv.wall_streamline(pts, 0.0, 1.0, dx=0.05)
check(f"一様流: 壁流線が水平 (max |r-1| {np.max(np.abs(wall[:,1]-1)):.1e})",
      len(wall) > 20 and float(np.max(np.abs(wall[:, 1] - 1.0))) < 1e-9)

# --- 2. 放射源流 (球状ソース A/A*=(R/R*)^2, 半角 10°) ---------------------------
thc = np.deg2rad(10.0)


def src_M(x, r):
    R = max(np.hypot(x, r), 1.0001)
    return float(invert_area_ratio(np.array([R * R]), np.array([True]), G)[0])


x1 = 1.2
r_top = x1 * np.tan(thc)
# 軸データは回復したい壁端 (x~2.5) の C⁻ 足 (下流) まで延長が必要 — §4.7(c) の実証
init = ([_Pt(float(x), 0.0, 0.0, float(pm_nu(src_M(x, 0.0), G)), G)
         for x in np.linspace(5.5, x1, 141)[:-1]]
        + [_Pt(x1, float(r), float(np.arctan2(r, x1)),
               float(pm_nu(src_M(x1, r), G)), G) for r in np.linspace(0.0, r_top, 25)])
pts = inv.fill(init)
wall = inv.wall_streamline(pts, x1, r_top, dx=0.02)
ray_err = np.abs(wall[:, 1] / wall[:, 0] - np.tan(thc)) / np.tan(thc)
check(f"源流: 壁流線がレイを保持 (max rel err {ray_err.max():.2%}, 到達 x={wall[-1,0]:.2f})",
      wall[-1, 0] > 2.4 and float(ray_err.max()) < 0.01)

# --- 3. ノズル (Sauer + Bézier 目標 Md=4) ---------------------------------------
Md = 4.0
st = SauerThroat(R=2.0, gamma=G)
x0, rr0, MM0, tt0 = st.starting_line(M_start=1.05, n=41)
h = 1e-4
M0 = float(st.mach(x0, 0.0))
dM0 = (st.mach(x0 + h, 0.0) - st.mach(x0 - h, 0.0)) / (2 * h)
xd = 14.0
bz = MachBezier.from_constraints(x0, xd, start=(M0, float(dM0)),
                                 free_cp=[2.0, 2.9, 3.6], end=(Md, 0.0, 0.0))


def target(x):
    if x >= xd:
        return Md
    return float(np.atleast_1d(bz(float(x)))[0])


mu_d = np.arcsin(1.0 / Md)
# 1D 面積比: A/A* = (1/M)((2+(g-1)M^2)/(g+1))^((g+1)/(2(g-1)))
eps_1d = (1.0 / Md) * ((2.0 + (G - 1.0) * Md * Md) / (G + 1.0)) ** ((G + 1.0) / (2.0 * (G - 1.0)))
# 軸目標は壁リップの C⁻ 足まで下流延長が必要 (係数 ~2×: 「リップまで」+「リップの依存域」)
x_end = xd + 2.3 * np.sqrt(eps_1d) / np.tan(mu_d)
res = inverse_design(st, target, x_axis_end=float(x_end), n_axis=420, n_start=41)
wall = res["wall"]
r_exit, th_exit = float(wall[-1, 1]), float(np.rad2deg(wall[-1, 2]))
ratio = r_exit / np.sqrt(eps_1d)  # 期待 √Cd ≈ 0.996-0.998 (R=2 のソニックライン湾曲)
check(f"ノズル: リップまで完成 (x_lip={wall[-1,0]:.2f}, θ_end={th_exit:.3f}°)",
      abs(th_exit) < 0.5 and wall[-1, 0] > xd)
check(f"ノズル: 出口半径/1D√ε = {ratio:.4f} (√Cd 補正込みで ≈1)", 0.985 < ratio < 1.005)
check(f"ノズル: リップ壁 M = {wall[-1,3]:.3f} ≈ Md", abs(wall[-1, 3] - Md) < 0.02)
check("ノズル: 壁が単調拡大", bool(np.all(np.diff(wall[:, 1]) > -1e-9)))
mf = res["mdot_exit"] / res["mdot_start"]
check(f"ノズル: 質量流量整合 ({mf:.4f})", abs(mf - 1.0) < 0.01)

# --- 4. 決定性 ------------------------------------------------------------------
res2 = inverse_design(st, target, x_axis_end=float(x_end), n_axis=420, n_start=41)
check("決定性 (壁 bit 同一)", np.array_equal(res["wall"], res2["wall"]))

# --- 出口一様性メトリクス (A2: 環状面積重み + 有効菱形コア) --------------------
from forge_design.metrics.extract import _annular_areas, test_core_radius  # noqa: E402

r_s = np.linspace(0.1, 0.9, 5)
A_s = _annular_areas(r_s)
check(f"環状面積の総和 = π·r_out² ({A_s.sum():.4f} vs {np.pi:.4f})",
      abs(A_s.sum() - np.pi * 1.0 ** 2) < 1e-12)
# 非一様配置: 隙間・重複がない (テロスコープ) ことを境界から独立に検算
r_ne = np.array([0.05, 0.1, 0.2, 0.4, 0.8])
A_ne = _annular_areas(r_ne)
mid_ne = 0.5 * (r_ne[:-1] + r_ne[1:])
lo0, hiN = max(2 * r_ne[0] - mid_ne[0], 0.0), 2 * r_ne[-1] - mid_ne[-1]
check("環状面積: 非一様配置でもテロスコープ (隙間・重複なし)",
      abs(A_ne.sum() - np.pi * (hiN ** 2 - lo0 ** 2)) < 1e-12)
check("環状面積: 全て正", bool(np.all(A_ne > 0)))
check("コア半径 = (x−x_d)tanμ_d (Md=4, x=24.42, x_d=14 → 2.690)",
      abs(test_core_radius(24.42, 14.0, 4.0, 3.264) - 2.6904) < 1e-3)
check("コア半径は壁半径でクリップ", test_core_radius(60.0, 14.0, 4.0, 3.264) == 3.264)
check("コア半径は x<x_d で 0", test_core_radius(10.0, 14.0, 4.0, 3.264) == 0.0)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
