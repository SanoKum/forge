#!/usr/bin/env python3
"""超音速 establish 場に擬似衝撃波を seed する (出口を亜音速化して outlet_statPress が効くように)。

clean な超音速場から outlet_statPress は機能しない (出口超音速 → BC 無視 = lock-in)。
そこで establish 場 (res) を読み、x>X_SEED の下流セルを亜音速・高静圧に置換して
入力 h5 /VALUE/ に書く。これで出口亜音速 → BC が掴み、衝撃波が X_SEED 付近に立ち
背圧 Ps に応じて定位置へ収束する。

使い方:
    python3 seed_shock.py SRC_res.h5 DST_input.h5 --xseed 0.20 --pdown 2.2e6 [--k .. --omega ..]
"""
import argparse, sys, os, h5py, numpy as np
sys.path.insert(0, "/home/sano/work/forge/solver_density_cuda/tools")
from res_h5_to_vtu import parse_conne

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("dst")
ap.add_argument("--xseed", type=float, default=0.20)
ap.add_argument("--pdown", type=float, default=2.2e6)
ap.add_argument("--tdown", type=float, default=288.0)
ap.add_argument("--udown", type=float, default=50.0)
ap.add_argument("--gamma", type=float, default=1.4)
ap.add_argument("--R", type=float, default=287.0)
ap.add_argument("--k", type=float, default=31.547)
ap.add_argument("--omega", type=float, default=3.0018e7)
ap.add_argument("--ncells", type=int, required=True)
a = ap.parse_args()
g = a.gamma

with h5py.File(a.src, "r") as s:
    coord = np.array(s["MESH/COORD"]).reshape(-1, 3)[:, :2]
    conn, offs, _ = parse_conne(np.array(s["MESH/CONNE"]), a.ncells)
    ro = np.array(s["VALUE/ro"]); P = np.array(s["VALUE/P"])
    Ux = np.array(s["VALUE/Ux"]); Uy = np.array(s["VALUE/Uy"]); Uz = np.array(s["VALUE/Uz"])
# cell centroids
cx = np.empty(a.ncells); st = 0
for i, o in enumerate(offs):
    cx[i] = coord[conn[st:o], 0].mean(); st = o

# downstream override -> subsonic high pressure
down = cx > a.xseed
ro_d = a.pdown / (a.R * a.tdown)
ro2 = ro.copy(); Ux2 = Ux.copy(); Uy2 = Uy.copy(); Uz2 = Uz.copy(); P2 = P.copy()
ro2[down] = ro_d; Ux2[down] = a.udown; Uy2[down] = 0.0; Uz2[down] = 0.0; P2[down] = a.pdown

roUx = ro2 * Ux2; roUy = ro2 * Uy2; roUz = ro2 * Uz2
roe = P2 / (g - 1.0) + 0.5 * ro2 * (Ux2**2 + Uy2**2 + Uz2**2)
roK = ro2 * a.k; roOmega = ro2 * a.omega

with h5py.File(a.dst, "r+") as d:
    for name, arr in [("ro", ro2), ("roUx", roUx), ("roUy", roUy), ("roUz", roUz),
                      ("roe", roe), ("roK", roK), ("roOmega", roOmega)]:
        d["VALUE/" + name][...] = arr.astype(d["VALUE/" + name].dtype)
print(f"seeded shock at x={a.xseed*1000:.0f}mm, downstream P={a.pdown/1e6:.1f}MPa "
      f"({down.sum()} of {a.ncells} cells subsonic)")
