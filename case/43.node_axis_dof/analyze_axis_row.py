"""軸行診断: 軸ノード値 vs 第一内点行, 軸 M(x) の run 間比較, 偶関数 fit 外挿との差."""
import sys, json, numpy as np, h5py
from pathlib import Path

def load(rd, step=None):
    rd = Path(rd)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split('_')[1]))
    rf = res[-1] if step is None else rd / f"res_{step}.h5"
    with h5py.File(rd/"nozzle.h5") as nz:
        c = nz["/MESH/COORD"][:].reshape(-1,3)
    with h5py.File(rf) as f:
        d = {k: f["/VALUE/"+k][:] for k in ["ro","P","T","Ux","Uy","sonic"]}
    d["M"] = np.hypot(d["Ux"], d["Uy"])/d["sonic"]
    return c, d, rf.name

def columns(c):
    """構造格子前提: x 座標が同じノードを列にまとめ r 昇順に。"""
    xs = np.round(c[:,0], 9)
    ux = np.unique(xs)
    cols = []
    for x in ux:
        idx = np.where(xs == x)[0]
        idx = idx[np.argsort(c[idx,1])]
        cols.append(idx)
    return ux, cols

def even_fit(r, q, npts=4, deg=2):
    """r>0 の最初の npts 点で q = a0 + a2 r^2 (+ a4 r^4) を最小二乗、a0 を返す."""
    rr, qq = r[1:1+npts], q[1:1+npts]
    A = np.vstack([rr**(2*k) for k in range(deg//2+1)]).T
    a, *_ = np.linalg.lstsq(A, qq, rcond=None)
    return a[0]

def report(rd):
    c, d, name = load(rd)
    ux, cols = columns(c)
    rows = []
    for x, idx in zip(ux, cols):
        if len(idx) < 6: continue
        a, i1 = idx[0], idx[1]
        rows.append((x, c[i1,1], d["ro"][a]/d["ro"][i1], d["P"][a]/d["P"][i1], d["M"][a], d["M"][i1],
                     even_fit(c[idx,1], d["M"][idx], 4, 2), even_fit(c[idx,1], d["M"][idx], 6, 4)))
    r = np.array(rows)
    print(f"== {rd} ({name}) axis nodes {len(r)}")
    print(f"  ro_axis/ro_row1: min {r[:,2].min():.4f} max {r[:,2].max():.4f} | P_axis/P_row1: min {r[:,3].min():.4f} max {r[:,3].max():.4f}")
    print(f"  M_axis - M_row1: max|.| {np.abs(r[:,4]-r[:,5]).max():.5f}")
    print(f"  M_axis - evenfit(4pt,r^2): max|.| {np.abs(r[:,4]-r[:,6]).max():.5f}  (rel {np.abs((r[:,4]-r[:,6])/r[:,4]).max()*100:.3f}%)")
    print(f"  M_axis - evenfit(6pt,r^4): max|.| {np.abs(r[:,4]-r[:,7]).max():.5f}")
    print(f"  P_axis min {d['P'][[i[0] for i in cols]].min():.1f} Pa")
    return r

if __name__ == "__main__":
    R = [report(a) for a in sys.argv[1:]]
    if len(R) == 2:
        # 同一メッシュ前提
        dM = R[1][:,4]-R[0][:,4]
        print(f"\n== run2 - run1: axis M max|Δ| {np.abs(dM).max():.5f} at x={R[0][np.abs(dM).argmax(),0]:.4f}, rms {np.sqrt(np.mean(dM**2)):.5f}")
        dM1 = R[1][:,5]-R[0][:,5]
        print(f"   row1 M max|Δ| {np.abs(dM1).max():.5f}")
        dMf = R[1][:,6]-R[0][:,6]
        print(f"   evenfit M max|Δ| {np.abs(dMf).max():.5f}")
        np.savetxt("axis_compare.csv", np.c_[R[0][:,0], R[0][:,4], R[1][:,4], R[0][:,6], R[1][:,6]], delimiter=",",
                   header="x,M_axis_run1,M_axis_run2,Mfit_run1,Mfit_run2", comments="")
