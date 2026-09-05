#!/usr/bin/env python3
"""CMC 反応 run の解析: 着火前線 (平均 Y_OH>2e-4 / 条件付き T_Q>1300 K) の時系列、平均 vs 条件付き温度、cmc_dY、
軸上 T・Y_H2O と実験の比較。  python3 analyze_cmc.py run_dir [ref_run ref_step]"""
import sys, glob, re, h5py, numpy as np, pathlib
run = sys.argv[1]; HERE = pathlib.Path(__file__).parent; D = 0.00457
with h5py.File(f"{run}/cabra.h5") as m: xyz = m["MESH/COORD"][:].reshape(-1, 3)
fs = sorted([f for f in glob.glob(f"{run}/res_*.h5") if re.fullmatch(r".*/res_\d+\.h5", f)], key=lambda p: int(re.findall(r"res_(\d+)\.h5", p)[0]))
exc = np.loadtxt(HERE / "exp_centerline.csv", delimiter=",", comments="#")
print("step | Tmax(mean) | ignx/d(OH>2e-4) | TQmax | x/d(TQ>1300) | cmc_dY max | Qdot max | axis T@z/d 14/20/26 (exp 620/1120/1410)")
for f in fs:
    s = int(re.findall(r"res_(\d+)\.h5", f)[0])
    with h5py.File(f) as h:
        V = h["VALUE"]; T = V["T"][:]; oh = V["Y4"][:]; qd = V["chemQdot"][:]
        TQ = V["cmc_TQmax"][:] if "cmc_TQmax" in V else np.zeros_like(T); dY = V["cmc_dY"][:] if "cmc_dY" in V else np.zeros_like(T)
    N = len(T); x = xyz[:N, 0]; y = xyz[:N, 1]
    sel = oh > 2e-4; ign = x[sel].min() / D if sel.any() else float("nan")
    selT = TQ > 1300; ignQ = x[selT].min() / D if selT.any() else float("nan")
    ax = np.where(np.abs(y) < 1e-9)[0]; o = np.argsort(x[ax]); Tax = lambda zd: np.interp(zd * D, x[ax][o], T[ax][o])
    print(f"{s:5d} | {T.max():7.0f} K | {ign:6.1f} | {TQ.max():5.0f} K | {ignQ:6.1f} | {dY.max():.3f} | {qd.max():.2e} | {Tax(14):.0f}/{Tax(20):.0f}/{Tax(26):.0f}")
# 最終 snapshot: 半径 T @ z/d 9, 26 の実験差
exr = np.loadtxt(HERE / "exp_radial_T.csv", delimiter=",", comments="#")
with h5py.File(fs[-1]) as h: T = h["VALUE/T"][:]
N = len(T); x = xyz[:N, 0]; y = xyz[:N, 1]
for zd, col in ((9, 2), (11, 3), (14, 4), (26, 5)):
    m = np.where(np.abs(x - zd * D) < 0.5 * D)[0]; r = y[m] * 1e3; o = np.argsort(r)
    Ti = np.interp(exr[:, 0], r[o], T[m][o]); msk = exr[:, 0] <= 30; d = Ti[msk] - exr[msk, col]
    print(f"final radial T z/d={zd:2d}: mean diff {d.mean():+5.0f} K, max |diff| {np.abs(d).max():4.0f} K (exp r=0: {exr[0,col]:.0f}, forge {Ti[0]:.0f})")
