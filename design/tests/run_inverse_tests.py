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

# --- 5. 充填のベクトル化 ≡ スカラー版 (A8) --------------------------------------
def _fill_scalar(inv_, init_):
    """`InverseMOC.fill` のベクトル化前の実装 (等価性の参照)。"""
    pts_, front = list(init_), list(init_)
    while len(front) >= 2:
        new = []
        for i in range(len(front) - 1):
            B, A = front[i], front[i + 1]
            try:
                P = inv_._interior(A, B)
            except RuntimeError:
                continue
            if not (np.isfinite(P.x) and np.isfinite(P.r)) or P.r < -1e-12:
                continue
            new.append(P)
        front = new
        pts_.extend(new)
    return pts_


x0s, rr0s, MM0s, tt0s = st.starting_line(M_start=1.05, n=31)
_init = ([_Pt(float(x), 0.0, 0.0, float(pm_nu(target(float(x)), G)), G)
          for x in np.linspace(8.0, x0s, 150)[:-1]]
         + [_Pt(float(x0s), float(rr0s[i]), float(tt0s[i]),
                float(pm_nu(float(MM0s[i]), G)), G) for i in range(31)])
_sc = _fill_scalar(inv, _init)
_ve = inv.fill_arrays(_init)
check(f"充填ベクトル化: 点数一致 ({len(_sc)})", len(_sc) == len(_ve))
if len(_sc) == len(_ve):
    _d = np.abs(np.array([[p.x, p.r, p.th, p.nu, p.M] for p in _sc]) - _ve).max(axis=0)
    check(f"充填ベクトル化 ≡ スカラー (max|Δ| x,r,θ,ν,M = {np.max(_d):.1e})",
          float(np.max(_d)) < 1e-12)

# --- 6. スロート特性線を初期値線にした逆設計 (A8) --------------------------------
from forge_design.geometry.transonic import HallThroat  # noqa: E402

_ht = HallThroat(R=2.0, gamma=G)
_xc, _rc, _Mc, _tc = _ht.throat_characteristic(n=41)
_xA = float(_xc[0])
_bz2 = MachBezier.from_constraints(_xA, xd, start=_ht.axis_anchor(_xA)[:2],
                                   free_cp=[2.4, 3.1, 3.7], end=(Md, 0.0, 0.0))


def target_tc(x):
    return Md if x >= xd else float(np.atleast_1d(_bz2(float(x)))[0])


_res_tc = inverse_design(_ht, target_tc, x_axis_end=float(x_end), n_axis=420,
                         n_start=41, start_line="throat_char", exit_mode="lip")
_w = _res_tc["wall"]
check(f"特性線始点: 壁流線がスロートから始まる ((x,r)=({_w[0,0]:.2e},{_w[0,1]:.5f}))",
      abs(_w[0, 0]) < 1e-9 and abs(_w[0, 1] - 1.0) < 1e-9)
check("特性線始点: 壁が単調拡大", bool(np.all(np.diff(_w[:, 1]) > -1e-9)))
_mf_tc = _res_tc["mdot_exit"] / _res_tc["mdot_start"]
# 縦線 (上の §3) との直接比較はしない — 遷音速解 (Sauer/Hall)・目標曲線・x_A が
# 揃っておらず対照になっていないため。同一条件の A/B は plan §9 に記録する。
check(f"特性線始点: 質量流量整合 ({_mf_tc:.4f})", abs(_mf_tc - 1.0) < 0.01)

# --- 7. 流束閉包による壁決定 (A9 — CFD で棄却済み、性質を固定する特性化テスト) ----
# 既定は wall_mode='streamline'。本節は「実装は正しいが壁としては流線に劣る」ことを
# 記録に固定する (plans/archived/tooling-nozzle-moc-flux-closure-wall.md)。
from forge_design.geometry.moc_inverse import (field_interpolator,  # noqa: E402
                                               flux_closure_wall)

# 一様流 (M=2, θ=0) の解析解: ṁ(r) = ρV·πr² なので ṁ* = ρV·π → r_w ≡ 1
inv_u = InverseMOC(gamma=G, delta=1.0)
init_u = ([_Pt(float(x), 0.0, 0.0, nuU, G) for x in np.linspace(4.0, 0.0, 81)[:-1]]
          + [_Pt(0.0, float(r), 0.0, nuU, G) for r in np.linspace(0.0, 1.0, 41)])
pts_u = inv_u.fill(init_u)
itp_u = field_interpolator(pts_u)
mdot_u = inv_u.massflux_across(pts_u, 0.05, 1.0, n=400, itp=itp_u)
xs_u = np.linspace(0.3, 1.5, 25)
rw_u = flux_closure_wall(itp_u, mdot_u, xs_u, np.full_like(xs_u, 1.2), gamma=G)
check(f"流束閉包: 一様流で r_w ≡ 1 (max|r−1| = {np.max(np.abs(rw_u-1)):.1e})",
      bool(np.all(np.isfinite(rw_u))) and float(np.max(np.abs(rw_u - 1.0))) < 2e-3)

# ノズル: 解像度を 2 倍振っても出口半径が動かない (流線壁は動く)
_rs = {}
for _n, _ns in ((250, 31), (500, 41)):
    for _mode in ("streamline", "flux"):
        _r = inverse_design(_ht, target_tc, x_axis_end=float(x_end), n_axis=_n,
                            n_start=_ns, dx_wall=0.02, start_line="throat_char",
                            exit_mode="lip", wall_mode=_mode)
        _rs[(_n, _mode)] = _r
_d_fx = abs(_rs[(500, "flux")]["wall"][-1, 1] - _rs[(250, "flux")]["wall"][-1, 1])
_d_sl = abs(_rs[(500, "streamline")]["wall"][-1, 1] - _rs[(250, "streamline")]["wall"][-1, 1])
check(f"流束閉包: 出口半径が解像度にほぼ依存しない (Δr {_d_fx:.1e} < 流線 {_d_sl:.1e}/3)",
      _d_fx < _d_sl / 3.0)
check(f"流束閉包: 壁始点はスロートのまま ({_rs[(500,'flux')]['wall'][0,1]:.6f})",
      abs(_rs[(500, "flux")]["wall"][0, 1] - 1.0) < 1e-9)
check("流束閉包: 壁が単調拡大",
      bool(np.all(np.diff(_rs[(500, "flux")]["wall"][:, 1]) > -1e-9)))
_wm = _rs[(500, "flux")]["wall_mode"]
check(f"流束閉包: 合成帯で流線壁との差が小さい ({_wm['d_blend_max']:.1e} < 5e-3)",
      _wm["d_blend_max"] < 5e-3)

# **棄却の根拠を固定** (厳密解 = 放射源流のレイ): 流束閉包は解像度に依らない
# **バイアス床**を持ち、流線は細かくすればその床を下回る。粗い場では流束閉包の方が
# 良く見えるので「解像度非依存 = 良い壁」と誤読しやすい — 両解像度で測って固定する。
_e = {}
for _n_ax, _n_st, _dxw in ((141, 25, 0.02), (281, 49, 0.01)):
    _in = ([_Pt(float(x), 0.0, 0.0, float(pm_nu(src_M(x, 0.0), G)), G)
            for x in np.linspace(5.5, x1, _n_ax)[:-1]]
           + [_Pt(x1, float(r), float(np.arctan2(r, x1)),
                  float(pm_nu(src_M(x1, r), G)), G)
              for r in np.linspace(0.0, r_top, _n_st)])
    _p = inv.fill(_in)
    _it = field_interpolator(_p)
    _md = inv.massflux_across(_p, x1 + 0.01, r_top * 1.01, n=800, itp=_it)
    _xq = np.linspace(1.4, 2.4, 21)
    _ex = _xq * np.tan(thc)
    _rq = flux_closure_wall(_it, _md, _xq, _ex * 1.15, gamma=G)
    _ws = inv.wall_streamline(_p, x1, r_top, dx=_dxw, itp=_it)
    _e[_n_ax] = (float(np.max(np.abs(np.interp(_xq, _ws[:, 0], _ws[:, 1]) / _ex - 1.0))),
                 float(np.max(np.abs(_rq / _ex - 1.0))))
    print(f"info 源流 n={_n_ax}: 流線 {_e[_n_ax][0]:.2e} / 流束閉包 {_e[_n_ax][1]:.2e}")
check("源流: 流束閉包はバイアス床を持つ (解像度 2 倍で改善しない)",
      abs(_e[281][1] - _e[141][1]) < 0.5 * _e[141][1])
check(f"源流: 細かい場では流線が流束閉包を下回る "
      f"({_e[281][0]:.2e} < {_e[281][1]:.2e}) — A9 棄却の根拠",
      _e[281][0] < _e[281][1])

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
