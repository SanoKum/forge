#!/usr/bin/env python3
"""res_*.h5 (primitives) を forge 入力 h5 の /VALUE/ (保存量) に焼き直す restart ユーティリティ。

使い方:
    python3 restart_field.py SRC_res.h5 DST_input.h5 [--k K --omega OMEGA]

- SRC の primitives (ro,P,T,Ux,Uy,Uz[,k,omega]) を読み、保存量 (ro,roUx,roUy,roUz,roe,roK,roOmega)
  を計算して DST_input.h5 の /VALUE/ を上書きする。セル順は同一メッシュ前提で一致。
- --k/--omega 指定時は自由流値で roK/roOmega を初期化 (laminar 場 -> SST 起動など)。
  未指定なら SRC の k/omega フィールドをそのまま引き継ぐ (SST -> SST restart)。
"""
import argparse, h5py, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("dst")
ap.add_argument("--gamma", type=float, default=1.4)
ap.add_argument("--k", type=float, default=None)
ap.add_argument("--omega", type=float, default=None)
a = ap.parse_args()
g = a.gamma

with h5py.File(a.src, "r") as s:
    ro = np.array(s["VALUE/ro"]); P = np.array(s["VALUE/P"])
    Ux = np.array(s["VALUE/Ux"]); Uy = np.array(s["VALUE/Uy"]); Uz = np.array(s["VALUE/Uz"])
    has_turb = "VALUE/k" in s and "VALUE/omega" in s
    if a.k is None and has_turb:
        kk = np.array(s["VALUE/k"]); om = np.array(s["VALUE/omega"])
    else:
        kk = np.full_like(ro, a.k if a.k is not None else 1.0)
        om = np.full_like(ro, a.omega if a.omega is not None else 1.0e4)

roUx = ro*Ux; roUy = ro*Uy; roUz = ro*Uz
roe = P/(g-1.0) + 0.5*ro*(Ux**2+Uy**2+Uz**2)
roK = ro*kk; roOmega = ro*om

with h5py.File(a.dst, "r+") as d:
    n = d["VALUE/ro"].shape[0]
    assert n == ro.shape[0], f"cell count mismatch {n} vs {ro.shape[0]}"
    for name, arr in [("ro",ro),("roUx",roUx),("roUy",roUy),("roUz",roUz),
                      ("roe",roe),("roK",roK),("roOmega",roOmega)]:
        d["VALUE/"+name][...] = arr.astype(d["VALUE/"+name].dtype)
print(f"restart: {a.src} -> {a.dst} ({n} cells, "
      f"{'freestream k/omega' if a.k is not None else 'inherited k/omega'})")
