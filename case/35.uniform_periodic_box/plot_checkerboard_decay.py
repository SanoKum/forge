#!/usr/bin/env python3
"""L1b: res_*.h5 時系列から市松振幅 A_cb(t) = |<P·parity>|/P0 を計算しプロット。

A_cb は圧力場の 2Δ odd-even モードへの射影 (市松そのものの振幅)。初期値 ε=1e-3。
pure KEEP では保持/成長、keepDissType=1 では減衰、が合否 (L1b)。
usage: plot_checkerboard_decay.py RUN_DIR [RUN_DIR ...] [--out PNG]
"""
import argparse, glob, os, re, sys
import numpy as np, h5py
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "solver_density_cuda", "tools"))
from res_h5_to_vtu import parse_conne

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+")
ap.add_argument("--out", default="checkerboard_decay.png")
a = ap.parse_args()

def parity_of(h5path):
    with h5py.File(h5path, "r") as f:
        nc = f["VALUE/ro"].shape[0]
        coord0 = np.array(f["MESH/COORD"]).reshape(-1, 3)
        if coord0.shape[0] == nc:
            # node モード: COORD がノード座標そのもの (CV と 1:1)。centCoords は境界で
            # 双対体積重心にシフトするため使わない (architecture-node-centroid-value-position 未完)。
            cc = coord0
        elif "CELLS/centCoords" in f and f["CELLS/centCoords"].shape[0] == 3*nc:
            cc = np.array(f["CELLS/centCoords"]).reshape(-1, 3)  # cell/node 両対応 (input h5)
        else:
            coord = np.array(f["MESH/COORD"]).reshape(-1, 3)
            conn, offs, _ = parse_conne(np.array(f["MESH/CONNE"]), nc)
            cc = np.zeros((nc, 3)); s = 0
            for i, o in enumerate(offs):
                cc[i] = coord[conn[s:o]].mean(axis=0); s = o
    idx = np.zeros((nc, 3), dtype=np.int64)
    for d in range(3):
        lo, hi = cc[:, d].min(), cc[:, d].max()
        # 一様格子前提: 格子間隔 h を「最小の正の座標差」から推定し rint で整数化 (丸め耐性)
        u = np.unique(np.round(cc[:, d], 9))
        h = np.median(np.diff(u)) if len(u) > 1 else 1.0  # median: ノイズ由来の微小 diff に頑健
        idx[:, d] = np.rint((cc[:, d] - lo) / h).astype(np.int64)
    return np.where((idx.sum(axis=1) % 2) == 0, 1.0, -1.0)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7, 4.5))

P0 = 1.0 / 1.4
for run in a.runs:
    files = sorted(glob.glob(os.path.join(run, "res_*.h5")),
                   key=lambda p: int(re.search(r"res_(\d+)\.h5", p).group(1)))
    files = [p for p in files if "nan" not in p]
    if not files:
        print(f"[skip] {run}: res なし"); continue
    # パリティは run dir の input h5 (CELLS/centCoords 保持) から取る。無ければ最初の res から
    import glob as _g
    inp = [q for q in _g.glob(os.path.join(run, "*.h5"))
           if "res_" not in os.path.basename(q)]
    par = parity_of(inp[0] if inp else files[0])
    steps, amp = [], []
    for p in files:
        with h5py.File(p, "r") as f:
            P = np.array(f["VALUE/P"], dtype=np.float64)
        steps.append(int(re.search(r"res_(\d+)\.h5", p).group(1)))
        # 共分散射影: サブセットで sum(par)!=0 のとき mean(P)*mean(par) の DC バイアスが
        # 乗る (境界 6146 ノードで 3.3e-4 を「残存市松」と誤認した事故の再発防止)。
        amp.append(abs(np.mean(P * par) - np.mean(P) * np.mean(par)) / P0)
    label = os.path.basename(run.rstrip("/"))
    ax.semilogy(steps, np.maximum(amp, 1e-12), "o-", label=label)
    print(f"{label}: A_cb {amp[0]:.3e} -> {amp[-1]:.3e} (step {steps[0]}..{steps[-1]})")

ax.axhline(1e-3, color="gray", ls=":", lw=0.8, label="initial eps")
ax.set_xlabel("step"); ax.set_ylabel(r"$A_{cb}=|\langle P\,\pi\rangle|/P_0$")
ax.set_title("checkerboard amplitude decay (L1b)")
ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(a.out, dpi=130)
print(f"saved {a.out}")
