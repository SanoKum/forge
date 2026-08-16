"""node run の軸量抽出: 軸ノード直読でなく正半径ノードからの偶関数外挿。

背景 (case/42.node_axis_dof): 生産設定 `nodeAxisDirichlet: 1` では軸ノードは解かれず
第一内点のコピー (∂q/∂r=0 の 1 次離散化) なので、軸直読 = 第一内点値 = 真の軸値 + ½ q_rr r₁²
のバイアスを持つ (case/41 で r₁/r_t≈0.036 → ΔM≈0.003 = 0.08% M_d)。軸対称スカラー量は
r の偶関数なので、軸ノードを除外し r>0 の数点で q = a₀ + a₂ r² (+ a₄ r⁴) を最小二乗し a₀ を
軸値とする。`nodeAxisDirichlet: 0` (軸 DOF を解く) でも軸半 CV の離散不整合が残る間は同じ理由で
軸ノードを除外するのが安全。
"""
from pathlib import Path

import numpy as np


def _columns(coords, tol=1e-9):
    """構造格子 (同一 x のノード列) を列にまとめる。戻り値: (x_unique, [idx 配列 (r 昇順)])。"""
    xs = np.round(coords[:, 0], 9)
    ux = np.unique(xs)
    cols = []
    for x in ux:
        idx = np.where(np.abs(xs - x) < tol)[0]
        cols.append(idx[np.argsort(coords[idx, 1])])
    return ux, cols


def even_fit_axis(r, q, npts=4, deg=2, skip_axis=True):
    """1 列の (r, q) から偶関数 fit の軸値 a₀ を返す。r 昇順、r[0] が軸ノード。"""
    i0 = 1 if skip_axis else 0
    rr, qq = np.asarray(r[i0:i0 + npts], float), np.asarray(q[i0:i0 + npts], float)
    if len(rr) < deg // 2 + 1:
        return float(q[0])
    A = np.vstack([rr ** (2 * k) for k in range(deg // 2 + 1)]).T
    a, *_ = np.linalg.lstsq(A, qq, rcond=None)
    return float(a[0])


def axis_profile(run_dir, scale, quantity="M", npts=4, deg=2, nsnap=3, skip_axis=True):
    """node run から軸量の x 分布を偶関数外挿で作る。

    戻り値 dict: x (r_t 単位, 昇順), q_fit (偶関数外挿の軸値), q_axis (軸ノード直読),
    q_row1 (第一内点), r1 (第一内点半径 [m])。quantity は "M" (|U|/c) か VALUE 名。"""
    import h5py
    rd = Path(run_dir)
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    res = [f for f in res if int("".join(c for c in f.stem if c.isdigit())) > 0]
    if not res:
        raise FileNotFoundError(f"{rd} に res_*.h5 が無い")
    with h5py.File(rd / "nozzle.h5") as nz:
        nc = nz["/MESH/COORD"][:].reshape(-1, 3)
    qs = []
    for rf in res[-nsnap:]:
        with h5py.File(rf) as f:
            if quantity == "M":
                q = np.hypot(f["/VALUE/Ux"][:], f["/VALUE/Uy"][:]) / np.maximum(f["/VALUE/sonic"][:], 1e-9)
            else:
                q = f["/VALUE/" + quantity][:]
        if len(q) != len(nc):
            raise ValueError("node run でない (VALUE 長 != 節点数)")
        qs.append(q.astype(np.float64))
    q = np.mean(qs, axis=0)
    ux, cols = _columns(nc)
    x, qf, qa, q1, r1 = [], [], [], [], []
    for xv, idx in zip(ux, cols):
        if len(idx) < npts + 1 or nc[idx[0], 1] > 1e-12:
            continue                      # 軸に届かない列 (壁側の途中列など) は除外
        r = nc[idx, 1]
        x.append(xv / scale)
        qf.append(even_fit_axis(r, q[idx], npts, deg, skip_axis))
        qa.append(float(q[idx[0]])); q1.append(float(q[idx[1]])); r1.append(float(r[1]))
    o = np.argsort(x)
    x = np.asarray(x)[o]
    return {"x": x, "q_fit": np.asarray(qf)[o], "q_axis": np.asarray(qa)[o],
            "q_row1": np.asarray(q1)[o], "r1": np.asarray(r1)[o]}


def axis_curve_evenfit(run_dir, scale, lam=1e-5, npts=4, deg=2):
    """`runner_axismach.axis_curve_node` の偶関数外挿版: 平滑化スプライン (M, M', M'' 同一スプライン)。"""
    from scipy.interpolate import make_smoothing_spline
    p = axis_profile(run_dir, scale, "M", npts, deg)
    spl = make_smoothing_spline(p["x"], p["q_fit"], lam=lam)
    return spl, float(p["x"].min()), float(p["x"].max())
