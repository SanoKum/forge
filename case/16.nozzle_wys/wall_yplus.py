#!/usr/bin/env python3
"""SST run の壁第一層 y+ を収束場から算出・可視化する。

forge は内部で y+ を持つが res に書き出さないため、メッシュ幾何 (壁面 physID とその隣接セル) と
収束場から  y+ = u_tau * y1 / nu_molecular,  u_tau = sqrt(tau_w/rho),  tau_w = mu * u_t / y1
で算出する (線形サブレイヤ近似)。低Re SST では壁第一セルが粘性サブレイヤ内 (mu_turb≈0) なので
分子粘性 mu_lam で tau_w を評価してよい (有効粘性 mu_lam+mu_turb でも y+ は不変)。
y+ の定義中の nu は常に分子粘性。

使い方:
  python3 wall_yplus.py <run_dir> [--wall-physid 3] [--res res_XXXX.h5]
出力: 端末に統計、wall_yplus.png に壁沿い分布。
"""
import argparse, glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("--wall-physid", type=int, default=3, help="no-slip 壁の physID (既定 3)")
ap.add_argument("--res", default=None, help="使用する res ファイル名 (既定: 最新)")
ap.add_argument("-o", "--out", default="wall_yplus.png")
a = ap.parse_args()

def mesh_path(run):
    return [p for p in glob.glob(f"{run}/nozzle_*.h5") if "res_" not in os.path.basename(p)][0]
def latest_res(run):
    fs = [f for f in glob.glob(f"{run}/res_*.h5") if "nan" not in f]
    return max(fs, key=lambda f: int(f.split("res_")[1].split(".")[0]))

m = h5py.File(mesh_path(a.run), "r")
res = f"{a.run}/{a.res}" if a.res else latest_res(a.run)
r = h5py.File(res, "r")["/VALUE"]
b = f"BCONDS/{a.wall_physid}"

iC = np.array(m[f"{b}/iCells"]); iP = np.array(m[f"{b}/iPlanes"])
cc = np.array(m["CELLS/centCoords"]).reshape(-1, 3)[iC]
fc = np.array(m["PLANES/centCoords"]).reshape(-1, 3)[iP]
n = np.array(m["PLANES/surfVect"]).reshape(-1, 3)[iP]
n = n / np.linalg.norm(n, axis=1, keepdims=True)
y1 = np.abs(np.sum((cc - fc) * n, axis=1))                       # 第一セル壁直交距離

U = np.stack([np.array(r["Ux"])[iC], np.array(r["Uy"])[iC], np.array(r["Uz"])[iC]], 1)
ut = np.linalg.norm(U - np.sum(U * n, axis=1, keepdims=True) * n, axis=1)  # 壁接線速度
assert (ut <= np.linalg.norm(U, axis=1) + 1e-6).all(), "u_t > |U|: 法線/インデックスのバグ"

ro = np.array(r["ro"])[iC]; mul = np.array(r["vis_lam"])[iC]; nu = mul / ro
yp = np.sqrt(ut * y1 / nu)                                       # = u_tau*y1/nu (線形サブレイヤ)

print(f"run={a.run}  res={os.path.basename(res)}  wall physID={a.wall_physid}  ({len(yp)} 壁面)")
print(f"  first-cell y1: {y1.min():.2e} .. {y1.max():.2e} m")
print(f"  y+ : min={yp.min():.2f} max={yp.max():.2f} mean={yp.mean():.2f} median={np.median(yp):.2f}")
print(f"  y+<=1: {100*(yp<=1).mean():.1f}%   y+<=2: {100*(yp<=2).mean():.1f}%   y+<=5: {100*(yp<=5).mean():.1f}%")

x = fc[:, 0] * 100.0
o = np.argsort(x)
plt.figure(figsize=(9, 4.5))
plt.plot(x[o], yp[o], "-", lw=1, color="navy")
plt.axhline(1, color="gray", ls="--", lw=1); plt.text(x.max()*0.6, 1.04, "y+=1", color="gray")
plt.xlabel("x [cm]"); plt.ylabel("wall first-cell y+")
plt.ylim(0, max(2.2, yp.max()*1.1)); plt.grid(alpha=0.3)
plt.title(f"Wall y+ along contour  ({os.path.basename(a.run)})")
plt.tight_layout(); plt.savefig(f"{a.run}/{a.out}", dpi=120)
print(f"  saved {a.run}/{a.out}")
