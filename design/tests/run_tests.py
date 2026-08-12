#!/usr/bin/env python3
"""forge_design Phase 0 単体テスト (pytest 非依存・素の assert)。"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.geometry.bezier import MachBezier  # noqa: E402
from forge_design.geometry.transonic import SauerThroat  # noqa: E402
from forge_design.geometry.wall import NozzleWall, WallParams  # noqa: E402
from forge_design.meshing.mesh2d import (  # noqa: E402
    Mesh2DParams, generate_axisym_mesh, write_msh41_2d)

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


# --- Bézier: 端点条件・自由度勘定 ---------------------------------------------
bz = MachBezier.from_constraints(
    0.5, 8.0, start=(1.5, 0.8, -0.3), free_cp=[2.2, 2.8, 3.1], end=(4.0, 0.0, 0.0))
check("bezier 次数 = 拘束6+自由3-1 = n=8 (CP 9)", len(bz.cp) == 9)
check("bezier M(x0)", abs(bz(0.5) - 1.5) < 1e-12)
check("bezier M'(x0)", abs(bz.deriv(0.5, 1) - 0.8) < 1e-9)
check("bezier M''(x0)", abs(bz.deriv(0.5, 2) - (-0.3)) < 1e-8)
check("bezier M(x1)", abs(bz(8.0) - 4.0) < 1e-12)
check("bezier M'(x1)=0", abs(bz.deriv(8.0, 1)) < 1e-9)
check("bezier M''(x1)=0", abs(bz.deriv(8.0, 2)) < 1e-8)
bz2 = MachBezier.from_constraints(0.0, 5.0, start=(1.2, 0.5), free_cp=[2.0, 2.5, 3.0])
check("bezier C1 開始 + 端拘束なし", len(bz2.cp) == 5 and abs(bz2.deriv(0.0, 1) - 0.5) < 1e-9)
check("bezier 単調 CP 判定", bz2.cp_monotone())

# --- Sauer 遷音速解 -----------------------------------------------------------
st = SauerThroat(R=2.0, gamma=1.4)
check("sauer 軸ソニック点はスロート下流", st.x_sonic_axis > 0.0)
check("sauer 軸上 M=1 @ソニック点", abs(st.mach(st.x_sonic_axis, 0.0) - 1.0) < 1e-9)
# 壁接線条件 (線形化形 vbar = dr_w/dx — 小擾乱理論の次数で照合)
for x in (0.0, 0.05, 0.15):
    rw = 1.0 + x * x / (2.0 * st.R)
    check(f"sauer 壁接線条件 @x={x}", abs(st.vbar(x, rw) - x / st.R) < 5e-3)
x0, r, M, th = st.starting_line(M_start=1.05, n=21)
check("sauer starting line 全超音速", np.all(M >= 1.05 - 1e-9))
check("sauer starting line 軸で theta=0", abs(th[0]) < 1e-12)
check("sauer Cd 見込み (R=2, 0.9<Cd<1)", 0.9 < st.cd_estimate() < 1.0)

# --- 区分構成壁 ---------------------------------------------------------------
wp = WallParams()
w = NozzleWall(wp)
check("wall validate OK", w.validate() == [])
check("wall throat r=1", abs(w.r(np.array([0.0]))[0] - 1.0) < 1e-12)
check("wall exit r=r_exit", abs(w.r(np.array([wp.L_div]))[0] - wp.r_exit) < 1e-9)
check("wall inlet r=r_inlet", abs(w.r(np.array([w.x_in]))[0] - wp.r_inlet) < 1e-12)
# 接続点の C1 連続 (数値微分)
for xc, name in [(w.x_cs, "pipe/contraction"), (w.x_au, "contraction/Ru"),
                 (w.x_a, "Rd/bell")]:
    h = 1e-6
    dl = (w.r(np.array([xc - h]))[0] - w.r(np.array([xc - 2 * h]))[0]) / h
    dr_ = (w.r(np.array([xc + 2 * h]))[0] - w.r(np.array([xc + h]))[0]) / h
    check(f"wall C1 連続 @{name}", abs(dl - dr_) < 5e-3)
# 収縮の C2 接続 (円弧曲率と一致)
h = 1e-4
xa = w.x_au
d2l = (w.r(np.array([xa - h]))[0] - 2 * w.r(np.array([xa - 2 * h]))[0] + w.r(np.array([xa - 3 * h]))[0]) / h**2
curv = 1.0 / (wp.Ru * np.cos(wp.phi_u) ** 3)
check("wall C2 接続 @contraction/Ru", abs(d2l - curv) / curv < 0.05)

# --- TOP 放物線ベル (Phase 2) --------------------------------------------------
wpt = WallParams(bell_type="top")
wt = NozzleWall(wpt)
check("TOP validate OK", wt.validate() == [])
check("TOP r(P_a)=r_a", abs(wt.r(np.array([wt.x_a]))[0] - wt.r_a) < 1e-12)
check("TOP r(x_e)=r_exit", abs(wt.r(np.array([wpt.L_div]))[0] - wpt.r_exit) < 1e-9)
check("TOP 接線 θn @P_a",
      abs(wt.drdx(np.array([wt.x_a + 1e-12]))[0] - np.tan(wpt.theta_a)) < 1e-6)
check("TOP 接線 θe @出口",
      abs(wt.drdx(np.array([wpt.L_div]))[0] - np.tan(wpt.theta_e)) < 1e-9)
h = 1e-6
dl = (wt.r(np.array([wt.x_a - h]))[0] - wt.r(np.array([wt.x_a - 2 * h]))[0]) / h
dr_ = (wt.r(np.array([wt.x_a + 2 * h]))[0] - wt.r(np.array([wt.x_a + h]))[0]) / h
check("TOP C1 連続 @Rd/bell", abs(dl - dr_) < 5e-3)
xs_t = np.linspace(wt.x_a, wpt.L_div, 500)
check("TOP drdx 単調減少 (放物線の凸性)", np.all(np.diff(wt.drdx(xs_t)) < 1e-12))
# r(x) と数値微分の整合 (t(x) 逆変換の検証)
dnum = (wt.r(xs_t[1:-1] + 1e-6) - wt.r(xs_t[1:-1] - 1e-6)) / 2e-6
check("TOP drdx = 数値微分", np.max(np.abs(wt.drdx(xs_t[1:-1]) - dnum)) < 1e-5)
# 整合チェック: θn < 弦勾配 / θe > 弦勾配 は明示エラー
for bad in (dict(theta_a=np.deg2rad(10.0)), dict(theta_e=np.deg2rad(20.0))):
    try:
        NozzleWall(WallParams(bell_type="top", **bad))
        check(f"TOP 整合 NG 検出 {bad}", False)
    except ValueError:
        check(f"TOP 整合 NG 検出 {list(bad)[0]}", True)

# --- メッシュ: 決定性・整合 ---------------------------------------------------
mp = Mesh2DParams(ni=61, nj=21, scale=0.01)
c1, q1, b1 = generate_axisym_mesh(w, mp)
c2, q2, b2 = generate_axisym_mesh(w, mp)
check("mesh 決定性 (bit 同一)", np.array_equal(c1, c2) and np.array_equal(q1, q2))
check("mesh 軸ノード r=0", np.all(np.abs(c1.reshape(mp.ni, mp.nj, 3)[:, 0, 1]) < 1e-15))
areas = []
for q in q1:
    p = c1[q][:, :2]
    areas.append(0.5 * np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1]))
check("mesh 全 quad 正方向 (CCW)", np.all(np.array(areas) > 0.0))
import hashlib
import tempfile
with tempfile.TemporaryDirectory() as td:
    f1, f2 = Path(td) / "a.msh", Path(td) / "b.msh"
    write_msh41_2d(f1, c1, q1, b1)
    write_msh41_2d(f2, c2, q2, b2)
    check("msh 書き出し決定性 (sha256 同一)",
          hashlib.sha256(f1.read_bytes()).hexdigest() == hashlib.sha256(f2.read_bytes()).hexdigest())

# --- kernel MOC --------------------------------------------------------------
import warnings

from forge_design.geometry.moc_kernel import (  # noqa: E402
    KernelMOC, pm_mach, pm_nu, trace_cminus)
from forge_design.evaluate.ic import invert_area_ratio  # noqa: E402

for M in (1.05, 1.5, 2.5, 4.0):
    check(f"PM 往復 M={M}", abs(pm_mach(pm_nu(M)) - M) < 1e-9)

# 放射源流 (球状ソース A/A*=(R/R*)^2) 照合 — 軸対称源項と符号の検証
g = 1.4
def _fM(x, r):
    R = max(np.hypot(x, r), 1.0001)
    return float(invert_area_ratio(np.array([R ** 2]), np.array([True]), g)[0])
_fth = lambda x, r: float(np.arctan2(r, x))  # noqa: E731
thc = np.deg2rad(10.0)
fr = trace_cminus(_fM, _fth, 1.10, 1.10 * np.tan(thc), 60, g)
ax = KernelMOC(gamma=g).march_wall(
    fr, lambda s: (s, s * np.tan(thc)), lambda s: thc, np.linspace(1.10, 1.9, 41)[1:])
err = max(abs(M - _fM(x, 0)) / M for x, M in ax[1:])
check(f"MOC 放射源流照合 (err={err:.1e})", err < 3e-3)

# 2D 平面自己検証: 円弧 kernel の軸上 ν = ν0 + 2θa (単純波則)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    stp = SauerThroat(R=2.0, gamma=g)
    outp = KernelMOC(gamma=g, delta=0.0).march_arc(stp, np.deg2rad(15.0), n_start=41, n_wall=60)
Mex = pm_mach(pm_nu(outp["axis_M"][0], g) + 2 * np.deg2rad(15.0), g)
check(f"MOC 平面 kernel 単純波則 ({outp['anchor'][0]:.3f} vs {Mex:.3f})",
      abs(outp["anchor"][0] - Mex) / Mex < 0.02)

# 軸対称: 解像度で anchor が動かない (格子収束)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    o1 = KernelMOC(gamma=g).march_arc(SauerThroat(R=2.0), np.deg2rad(15.0), 41, 80)
    o2 = KernelMOC(gamma=g).march_arc(SauerThroat(R=2.0), np.deg2rad(15.0), 81, 160)
check("MOC 軸対称 anchor 格子収束 (<0.5%)",
      abs(o1["anchor"][0] - o2["anchor"][0]) / o2["anchor"][0] < 5e-3)
check("MOC x_k は円弧終端より下流", o1["x_k"] > 2.0 * np.sin(np.deg2rad(15.0)))

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}")
sys.exit(1 if FAIL else 0)
