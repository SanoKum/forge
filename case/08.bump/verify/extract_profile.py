#!/usr/bin/env python3
"""bump 検証用: メッシュ centroid が y≈YLINE のセルを抽出し、
x/L (L=3) でソートして無次元圧力 P/Pt・マッハ数・生値 (ro,u,P) を CSV 出力する。

usage: extract_profile.py --mesh M.h5 --result res.h5 --pt PT --output out.csv [--y 0.3] [--band 0.012]
"""
import argparse, csv, h5py, numpy as np

GAMMA = 1.4
L = 3.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--result", required=True)
    ap.add_argument("--pt", type=float, required=True, help="無次元化に使う inlet 全圧 [Pa]")
    ap.add_argument("--output", required=True)
    ap.add_argument("--y", type=float, default=0.30)
    ap.add_argument("--band", type=float, default=0.012)
    a = ap.parse_args()

    cc = h5py.File(a.mesh, "r")["CELLS"]["centCoords"][:].reshape(-1, 3)
    V = h5py.File(a.result, "r")["VALUE"]
    ro = V["ro"][:]; rx = V["roUx"][:]; ry = V["roUy"][:]; roe = V["roe"][:]
    u2 = (rx*rx + ry*ry) / (ro*ro)
    P = (GAMMA - 1.0) * (roe - 0.5*ro*u2)
    a_snd = np.sqrt(GAMMA * P / ro)
    M = np.sqrt(u2) / a_snd
    u = np.sqrt(u2)

    sel = np.abs(cc[:, 1] - a.y) < a.band
    x = cc[sel, 0] / L
    order = np.argsort(x)
    rows = list(zip(x[order], (P[sel]/a.pt)[order], M[sel][order],
                    ro[sel][order], u[sel][order], P[sel][order]))
    with open(a.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_over_L", "P_over_Pt", "Mach", "ro", "u", "P"])
        for r in rows:
            w.writerow([f"{v:.8e}" for v in r])
    print(f"wrote {a.output}  ({len(rows)} points at y~{a.y})")

if __name__ == "__main__":
    main()
