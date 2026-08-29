#!/usr/bin/env python3
# forge 入力 h5 の VALUE/* を物理単位 N2 衝撃管 IC (cpg) で上書きする。
# node 変換 (VALUE 長 = 節点数) と cell 変換 (VALUE 長 = セル数, CONNE 9列 hex) の両対応。
# 使い方: python3 write_sod_ic.py <input.h5>
import sys, numpy as np, h5py

path = sys.argv[1]
xdia = 0.5
PL, TL = 1.0e6, 2000.0
PR, TR = 1.0e5,  400.0
R, GAM = 8.314462618/0.0280134, 1.4

with h5py.File(path, 'r+') as f:
    co = f['MESH/COORD'][:].reshape(-1, 3)
    n = f['VALUE/ro'].shape[0]
    if n == co.shape[0]:                      # node 変換: DOF = 節点
        x = co[:, 0]
        kind = 'node'
    else:                                     # cell 変換: CONNE (nCells,9) hex
        c = f['MESH/CONNE'][:].reshape(-1, 9)
        assert c.shape[0] == n, f'DOF mismatch: VALUE={n} CONNE rows={c.shape[0]}'
        x = co[c[:, 1:9], 0].mean(1)
        kind = 'cell'
    P = np.where(x < xdia, PL, PR)
    T = np.where(x < xdia, TL, TR)
    ro = P/(R*T)
    roe = P/(GAM-1.0)
    z = np.zeros(n)
    for k, v in [('ro', ro), ('roUx', z), ('roUy', z), ('roUz', z), ('roe', roe),
                 ('P', P), ('T', T), ('Ux', z), ('Uy', z), ('Uz', z)]:
        if 'VALUE/'+k in f:
            f['VALUE/'+k][...] = v.astype(f['VALUE/'+k].dtype)
    print(f'[{kind}] n={n} L(P={PL:.0e},T={TL},ro={ro[x<xdia][0]:.4f}) '
          f'R(P={PR:.0e},T={TR},ro={ro[x>=xdia][0]:.4f}) -> {path}')
