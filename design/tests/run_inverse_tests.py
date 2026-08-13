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

# --- 出口一様性メトリクス (A2: 実出口面の厳密環状面積 + 有効菱形コア) ---------
from pathlib import Path as _P  # noqa: E402
from forge_design.metrics.extract import (  # noqa: E402
    exit_uniformity, core_radius_straight, core_radius_traced)

_rd = _P("case/41.wind_tunnel_design/run_0001_v1_euler")
if (_rd / "nozzle.h5").exists():
    _res = sorted(_rd.glob("res_[0-9]*.h5"),
                  key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    import h5py as _h5
    with _h5.File(_rd / "nozzle.h5") as _nz:
        _ip = _nz["/BCONDS/2/iPlanes"][:]
        _pc = _nz["/PLANES/centCoords"][:].reshape(-1, 3)
        _rf, _dr = _pc[_ip, 1], _nz["/PLANES/surfArea"][:][_ip]
    _A = 2 * np.pi * _rf * _dr
    _Aex = np.pi * ((_rf + _dr / 2) ** 2 - (_rf - _dr / 2) ** 2)
    check(f"出口面: 2πr·dr = π(r_o²−r_i²) 厳密 (max rel {np.max(np.abs(_A-_Aex)/_Aex):.1e})",
          float(np.max(np.abs(_A - _Aex) / _Aex)) < 1e-4)
    _rw = float((_rf + _dr / 2).max())
    check(f"出口面積の総和 = π·r_wall² ({_A.sum():.6e} vs {np.pi*_rw**2:.6e})",
          abs(_A.sum() - np.pi * _rw ** 2) / (np.pi * _rw ** 2) < 1e-5)
    _e = exit_uniformity(_rd / "nozzle.h5", _res, 4.0, x_d=0.14,
                         core_mode="traced", n_grid=200)
    check(f"真の壁半径 > 最外セル中心 ({_e['r_wall']:.5f})", _e["r_wall"] > _rf.max())
    check(f"εM: 面積積分と固定半径格子が一致 ({_e['eps_M_rms']*100:.3f}% vs "
          f"{_e['eps_M_rms_grid']*100:.3f}%)",
          abs(_e["eps_M_rms"] - _e["eps_M_rms_grid"]) / _e["eps_M_rms"] < 0.02)
    _rt = core_radius_traced(_rd / "nozzle.h5", _res, 0.14, _e["x_plane"])
    _rs = core_radius_straight(_e["x_plane"], 0.14, 4.0, _e["r_wall"])
    check(f"有効菱形: トレース {_rt/_e['r_wall']:.3f} と直線近似 {_rs/_e['r_wall']:.3f} が 2% 内",
          abs(_rt - _rs) / _rt < 0.02)
    check("有効菱形コアは壁半径より小さい (リップまで届かない)", _rt < _e["r_wall"])
else:
    print("skip 出口一様性 (run_0001 が無い)")

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
