#!/usr/bin/env python3
"""同一メッシュ restart (TP 多成分対応): res_*.h5 の保存量 (ro, roUx, roUy, roUz, roe, roK, roOmega, roY*) を
入力 h5 の /VALUE/ に上書きする (CPG 前提の restart_field.py は TP で使えない)。
  python3 restart_conserved.py SRC_res.h5 DST_input.h5"""
import sys, h5py, numpy as np
src, dst = sys.argv[1], sys.argv[2]
with h5py.File(src) as s, h5py.File(dst, "r+") as d:
    n = d["VALUE/ro"].shape[0]; keys = [k for k in s["VALUE"].keys() if k in ("ro","roUx","roUy","roUz","roe","roK","roOmega") or (k.startswith("roY") and k[3:].isdigit())]
    for k in keys:
        a = np.array(s["VALUE"][k]); assert a.shape[0] == n, (k, a.shape, n)
        if k in d["VALUE"]: d["VALUE"][k][:] = a
        else: d["VALUE"].create_dataset(k, data=a)
    print("restart:", src, "->", dst, keys)
