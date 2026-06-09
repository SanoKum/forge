#!/usr/bin/env python3
# 既存 run の res_*.h5 をリスタート初期場に変換する。
#   - ro,roUx..roe をコピー、roY{s}=Y{s}*ro を再構成 (出力は Y のみ保持のため)
#   - roK,roOmega を freestream 乱流量でシード (SST は k=0 が不動点なので非ゼロ必須)
# 使い方: python3 gen_restart.py <src_res.h5> <mesh.h5> <restart.h5> [nSpecies] [k_ini] [omega_ini]
import sys, numpy as np, h5py, shutil

src_f = sys.argv[1]
mesh_f = sys.argv[2]
dst_f = sys.argv[3]
nsp   = int(sys.argv[4]) if len(sys.argv) > 4 else 2
k_ini = float(sys.argv[5]) if len(sys.argv) > 5 else 10.0
om_ini= float(sys.argv[6]) if len(sys.argv) > 6 else 1.0e4

shutil.copy(mesh_f, dst_f)               # geometry + structure
src = h5py.File(src_f, 'r'); dst = h5py.File(dst_f, 'r+')
ro = src['VALUE/ro'][:]
for k in ['ro','roUx','roUy','roUz','roe']:
    dst['VALUE/'+k][...] = src['VALUE/'+k][:]
# 化学種 roY{s} = Y{s}*ro
for s in range(nsp):
    Y = src['VALUE/Y%d'%s][:]
    nm = 'VALUE/roY%d'%s
    if nm in dst: del dst[nm]
    dst.create_dataset(nm, data=(Y*ro).astype(np.float32))
# 乱流量シード
for nm, arr in [('VALUE/roK', (ro*k_ini).astype(np.float32)),
                ('VALUE/roOmega', (ro*om_ini).astype(np.float32))]:
    if nm in dst: del dst[nm]
    dst.create_dataset(nm, data=arr)
print('restart: ro[%.3f,%.3f] nSpecies=%d  roK=ro*%.1f roOmega=ro*%.0f -> %s'
      % (ro.min(), ro.max(), nsp, k_ini, om_ini, dst_f))
