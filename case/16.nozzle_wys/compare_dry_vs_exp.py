#!/usr/bin/env python3
"""dry run の中心線 p/p0 を Wyslouzil 実験 (wyslouzil_fig3_pp0.csv) と比較する。

dry (凝縮なし) なので実験の isentrope 列 (= 乾燥/等エントロピー圧力) が比較対象。
参考に凝縮あり (1.00kPa 列) も重ねて描く (dry はこの下=isentrope 上に乗るはず)。
3D メッシュは中心軸近傍 (|y|, |z-zc| 小) を抽出する。

usage: compare_dry_vs_exp.py <run_dir> <label> [<run_dir2> <label2> ...]
"""
import sys, glob, os
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

P0 = 59070.0   # 入口全圧 Pt (Wyslouzil 条件)

def latest(run):
    fs = [f for f in glob.glob(f"{run}/res_*.h5") if "nan" not in f]
    return max(fs, key=lambda f: int(f.split("res_")[1].split(".")[0])) if fs else None

def mesh(run):
    return [p for p in glob.glob(f"{run}/nozzle_*.h5") if "res_" not in os.path.basename(p)][0]

def centerline(run):
    res = latest(run)
    if res is None:
        return None
    cc = h5py.File(mesh(run), "r")["/CELLS/centCoords"][:].reshape(-1, 3)
    x, y, z = cc[:, 0], cc[:, 1], cc[:, 2]
    v = h5py.File(res, "r")["/VALUE"]
    P = np.array(v["P"])
    # throat = 最小チャンネル高さの x
    nb = 200; xe = np.linspace(x.min(), x.max(), nb+1); xcb = 0.5*(xe[:-1]+xe[1:])
    h = np.full(nb, np.nan)
    for i in range(nb):
        s = (x >= xe[i]) & (x < xe[i+1])
        if s.sum() > 3: h[i] = y[s].max() - y[s].min()
    ok = np.isfinite(h); xth = xcb[ok][np.argmin(h[ok])]
    # 中心軸近傍 (3D は z 中心も絞る)
    zc = 0.5*(z.min()+z.max()); dz = z.max()-z.min()
    epsy = (y.max()-y.min())*0.02
    epsz = max(dz*0.12, 1e-6)
    sel = (np.abs(y) < epsy) & (np.abs(z - zc) < epsz)
    xs = x[sel]; o = np.argsort(xs)
    return (xs[o]-xth)*100.0, P[sel][o]/P0, os.path.basename(res)

def binavg(x, y, xe):
    xc = 0.5*(xe[:-1]+xe[1:]); yb = np.full(len(xc), np.nan)
    for i in range(len(xc)):
        s = (x >= xe[i]) & (x < xe[i+1])
        if s.sum() > 0: yb[i] = np.median(y[s])
    return xc, yb

def load_exp():
    rows = []
    for ln in open("wyslouzil_fig3_pp0.csv"):
        p = ln.strip().split(",")
        if not p or not p[0] or not p[0][0].isdigit():
            continue
        rows.append([float(p[0]), float(p[1]), float(p[2])])
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2]   # x_cm, isentrope, cond_1kPa

def main(pairs):
    ex, eiso, econd = load_exp()
    xe = np.linspace(-1, 9, 55)
    plt.figure(figsize=(8.5, 6))
    cols = ['navy', 'tab:green', 'tab:orange']
    results = {}
    for i, (run, lab) in enumerate(pairs):
        c = centerline(run)
        if c is None:
            print("MISSING", run); continue
        xb, pb = binavg(c[0], c[1], xe)
        results[lab] = (xb, pb, c[2])
        plt.plot(xb, pb, '-', color=cols[i % len(cols)], lw=2, label=f'forge {lab} (dry)')
    plt.plot(ex, eiso, 's', color='gray', ms=7, label='Wyslouzil isentrope (exp, dry)')
    plt.plot(ex, econd, 'o', color='black', ms=6, mfc='none', label='Wyslouzil cond 1kPa (exp)')
    plt.xlim(-1, 8.5); plt.ylim(0.15, 0.55)
    plt.xlabel('distance from throat [cm]'); plt.ylabel('p / p0')
    plt.grid(alpha=0.3); plt.legend(fontsize=9)
    plt.title('Wyslouzil nozzle: dry centerline p/p0 vs experiment')
    plt.tight_layout()
    out = "dry_vs_exp.png"
    plt.savefig(out, dpi=130)
    print("saved", out)

    print("\n=== p/p0 at experiment stations (dry vs exp isentrope) ===")
    print(" x_cm |  exp_iso | " + " | ".join(f"{lab}" for _, lab in pairs) + " |  diff% (forge-exp)/exp")
    for xq in [0.1, 1.1, 2.1, 3.1, 4.2, 5.2, 6.0, 7.0, 8.0]:
        ei = np.interp(xq, ex, eiso)
        cells = []
        diffs = []
        for run, lab in pairs:
            if lab in results:
                xb, pb, _ = results[lab]; m = np.isfinite(pb)
                fv = np.interp(xq, xb[m], pb[m])
                cells.append(f"{fv:.3f}"); diffs.append(f"{100*(fv-ei)/ei:+.1f}")
            else:
                cells.append("  -  "); diffs.append("  -  ")
        print(f"{xq:5.1f} |  {ei:.3f}  | " + " | ".join(cells) + " |  " + " ".join(diffs))

if __name__ == "__main__":
    args = sys.argv[1:]
    pairs = [(args[i], args[i+1]) for i in range(0, len(args), 2)]
    main(pairs)
