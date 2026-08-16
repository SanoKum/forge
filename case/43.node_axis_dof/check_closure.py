"""r 重み幾何の離散閉性テスト (幾何のみ・run 不要)。

軸対称 (per-radian) の CV では、厳密に
    ∮ r n_y ds = ∫ dA = A_planar,     ∮ r n_x ds = 0,     ∫ r dA = r̄_area·A_planar = V
が成り立つ。forge の r 重みは面ベクトルに面重心半径 r_f を掛けるので、多角形の各辺で
∫ r ds = r_f|S_f| が厳密 → 上の 2 式は**離散でも厳密に成り立つはず**。体積側は r̄ の取り方
(双対重心 = 面積加重重心 = 厳密 / ノード座標 = 誤り) で決まる。本スクリプトはこの 3 つを実測する。
"""
import sys, numpy as np, h5py
from pathlib import Path

def load(h5):
    f = h5py.File(h5)
    st = f["/PLANES/STRUCT"][:]
    sv = f["/PLANES/surfVect"][:].reshape(-1, 3)
    pc = f["/PLANES/centCoords"][:].reshape(-1, 3)
    vol = f["/CELLS/volume"][:]            # r 重み前 (planar 面積 × 単位厚み)
    cc = f["/CELLS/centCoords"][:].reshape(-1, 3)
    nc = f["/MESH/COORD"][:].reshape(-1, 3)
    # PLANES/STRUCT = [nNodes, nodes..., nCells, cells...] の連結
    cells, i, ip = [], 0, 0
    while i < len(st):
        nn = st[i]; i += 1 + nn
        ncl = st[i]; i += 1
        cells.append(st[i:i+ncl].tolist()); i += ncl
        ip += 1
    return sv, pc, vol, cc, nc, cells

def closure(h5, label):
    sv, pc, vol, cc, nc, cells = load(h5)
    n = len(vol)
    rS = sv * np.maximum(pc[:, 1], 0.0)[:, None]     # r 重み面ベクトル (r_f = 面重心半径)
    Sx = np.zeros(n); Sy = np.zeros(n)
    for ip, cl in enumerate(cells):
        if len(cl) > 0 and cl[0] < n: Sx[cl[0]] += rS[ip, 0]; Sy[cl[0]] += rS[ip, 1]
        if len(cl) > 1 and cl[1] < n: Sx[cl[1]] -= rS[ip, 0]; Sy[cl[1]] -= rS[ip, 1]
    A = vol                                            # planar 面積
    ex = np.abs(Sx) / np.maximum(A, 1e-300)            # 0 であるべき (A 正規化)
    ey = np.abs(Sy - A) / np.maximum(A, 1e-300)        # 1 であるべき
    axis = nc[:, 1] < 1e-12
    print(f"== {label}")
    print(f"   |Σ r_f S_f,x| / A_planar : max {ex.max():.3e}  (軸行 {ex[axis].max():.3e})")
    print(f"   |Σ r_f S_f,y − A| / A    : max {ey.max():.3e}  (軸行 {ey[axis].max():.3e})")
    print(f"   V/(r̄_centroid·A) − 1     : max {np.abs(cc[:,1]*A/np.maximum(cc[:,1]*A,1e-300)-1).max():.1e} (定義上 0)")
    print(f"   ノード r で重み付けした場合の体積誤差: 軸行 {np.abs((nc[axis,1]-cc[axis,1])/np.maximum(cc[axis,1],1e-300)).max()*100:.1f}% (=軸で 100% 欠損)")
    return ex, ey

if __name__ == "__main__":
    for h5 in sys.argv[1:]:
        closure(h5, h5)
