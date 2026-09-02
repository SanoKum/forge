"""forge (node SST) vs SU2 (axisym SST) の境界層厚さ比較 — case/45 M6 CPG。

同一メッシュ (run_0011 の nozzle.msh → nozzle.su2)・同一 BC・同一ガス (CPG γ*=1.27354)。
ステーション x/r_t = {40, 60, 80, 94} で壁法線 (鉛直近似) の ρu プロファイルを取り、
δ99 (0.995 (ρu)_max)・δ* (軸対称質量収支)・θ (運動量厚, 同重み) を比較する。

usage: design/.venv-opt/bin/python compare_bl_su2.py
出力: compare_bl_su2.json / compare_bl_su2.png / 標準出力の表
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/home/sano/work/forge/design")
import h5py  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402
try:
    fm.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
    matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
except Exception:
    pass
from scipy.interpolate import LinearNDInterpolator  # noqa: E402

CASE = Path(__file__).resolve().parent
FORGE_RD = CASE / "run_0011_cpg_ns"
FORGE_PLAIN_RD = CASE / "run_0013_cpg_ns_plainsst"   # dilatation/KL off (SU2 相当)
SU2_RD = CASE / "run_0012_su2_sst"
STATIONS = (40.0, 60.0, 80.0, 94.0)
EDGE_FRAC = 0.995


def load_forge(rd=None):
    rd = rd or FORGE_RD
    S = json.loads((rd / "prepare_info.json").read_text())["scale_m"]
    with h5py.File(rd / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    import re, os
    res = sorted([f for f in os.listdir(rd)
                  if re.fullmatch(r"res_\d+\.h5", f) and f != "res_0.h5"],
                 key=lambda f: int(re.findall(r"\d+", f)[0]))[-1]
    with h5py.File(rd / res) as f:
        q = f["/VALUE/ro"][:] * f["/VALUE/Ux"][:]
    return nc[:, 0] / S, nc[:, 1] / S, q, S, res


def load_su2(S):
    """RESTART_ASCII (restart_flow.csv): 列名から x/y/Momentum_x (=ρu) を拾う。"""
    path = SU2_RD / "restart_flow.csv"
    with open(path) as f:
        rdr = csv.reader(f)
        hdr = [h.strip().strip('"') for h in next(rdr)]
        data = np.loadtxt(f, delimiter=",")
    col = {h: i for i, h in enumerate(hdr)}
    x = data[:, col["x"]] / S
    y = data[:, col["y"]] / S
    if "Momentum_x" in col:
        q = data[:, col["Momentum_x"]]
    else:  # 保険: Density × Velocity_x
        q = data[:, col["Density"]] * data[:, col["Velocity_x"]]
    return x, y, q


def bl_metrics(itp, xs, rw, n=600, window=3.0):
    y = np.linspace(0.0, min(window, 0.95 * rw), n)
    q = itp(np.full(n, xs), rw - y)
    ok = np.isfinite(q)
    q = np.where(ok, q, 0.0)
    qe = float(np.max(q))
    ie = int(np.argmax(q >= EDGE_FRAC * qe))
    ye = float(y[ie])
    m = y <= ye
    d = 1.0 - q[m] / max(qe, 1e-30)
    w = 1.0 - y[m] / rw
    I = float(np.trapezoid(d * w, y[m]))
    dstar = float(rw * (1.0 - np.sqrt(max(1.0 - 2.0 * I / rw, 0.0))))
    # 運動量厚 θ (同じ円環重み): ∫ (ρu/ρeue)(1-u/ue) ~ ρu 近似で (q/qe)(1-q/qe)
    It = float(np.trapezoid((q[m] / qe) * (1.0 - q[m] / qe) * w, y[m]))
    theta = float(rw * (1.0 - np.sqrt(max(1.0 - 2.0 * It / rw, 0.0))))
    return {"y_edge": ye, "dstar": dstar, "theta": theta,
            "profile_y": y, "profile_q": q / max(qe, 1e-30)}


def main():
    fx, fr, fq, S, res = load_forge()
    print(f"forge: {res}  scale={S}")
    px, pr, pq, _, pres = load_forge(FORGE_PLAIN_RD)
    print(f"forge plain: {pres}")
    sx, sy, sq = load_su2(S)
    print(f"su2: {len(sx)} pts")
    itp_f = LinearNDInterpolator(np.c_[fx, fr], fq)
    itp_p = LinearNDInterpolator(np.c_[px, pr], pq)
    itp_s = LinearNDInterpolator(np.c_[sx, sy], sq)
    wp = np.loadtxt(FORGE_RD / "wall_physical.csv", delimiter=",", skiprows=1) / S

    rows, out = [], {"stations": []}
    fig, axes = plt.subplots(1, len(STATIONS), figsize=(4 * len(STATIONS), 4.4),
                             sharey=False)
    for ax, xs in zip(np.atleast_1d(axes), STATIONS):
        rw = float(np.interp(xs, wp[:, 0], wp[:, 1]))
        bf = bl_metrics(itp_f, xs, rw)
        bp = bl_metrics(itp_p, xs, rw)
        bs = bl_metrics(itp_s, xs, rw)
        rows.append((xs, bf, bs, bp))
        out["stations"].append({
            "x_rt": xs, "rw": rw,
            "forge": {k: bf[k] for k in ("y_edge", "dstar", "theta")},
            "forge_plain": {k: bp[k] for k in ("y_edge", "dstar", "theta")},
            "su2": {k: bs[k] for k in ("y_edge", "dstar", "theta")},
            "ratio_d99_plain": bs["y_edge"] / bp["y_edge"] if bp["y_edge"] else None,
            "ratio_theta_plain": bs["theta"] / bp["theta"] if bp["theta"] else None,
            "ratio_dstar": bs["dstar"] / bf["dstar"] if bf["dstar"] else None,
            "ratio_d99": bs["y_edge"] / bf["y_edge"] if bf["y_edge"] else None})
        ax.plot(bf["profile_q"], bf["profile_y"], color="#1f77b4", lw=1.5,
                label="forge SST (生産設定: 圧縮性生産項 ON)")
        ax.plot(bp["profile_q"], bp["profile_y"], color="#2ca02c", lw=1.5,
                label="forge SST (素: 補正 OFF)")
        ax.plot(bs["profile_q"], bs["profile_y"], color="#d62728", lw=1.5,
                ls="--", label="SU2 SST V2003m")
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, max(bf["y_edge"], bs["y_edge"]) * 1.8)
        ax.set_title(f"x = {xs:g} r_t")
        ax.set_xlabel("ρu / (ρu)_e")
        ax.grid(alpha=.3)
    np.atleast_1d(axes)[0].set_ylabel("壁からの距離 y [r_t]")
    np.atleast_1d(axes)[0].legend(fontsize=8)
    fig.suptitle("M6 d1.55 CPG ノズル境界層: forge vs SU2 (同一メッシュ・同一BC)")
    fig.tight_layout()
    fig.savefig(CASE / "compare_bl_su2.png", dpi=120)

    print(f"\n{'x':>5} {'δ99 f':>8} {'δ99 s':>8} {'比':>6} "
          f"{'δ* f':>8} {'δ* s':>8} {'比':>6} {'θ f':>8} {'θ s':>8} {'比':>6}")
    for xs, bf, bs, bp in rows:
        print(f"{xs:5.0f} {bf['y_edge']:8.4f} {bs['y_edge']:8.4f} "
              f"{bs['y_edge']/bf['y_edge']:6.3f} "
              f"{bf['dstar']:8.4f} {bs['dstar']:8.4f} {bs['dstar']/bf['dstar']:6.3f} "
              f"{bf['theta']:8.4f} {bs['theta']:8.4f} {bs['theta']/bf['theta']:6.3f}")
    print("\n(plain vs SU2)")
    for xs, bf, bs, bp in rows:
        print(f"{xs:5.0f} d99 {bs['y_edge']/bp['y_edge']:6.3f}  "
              f"dstar {bs['dstar']/bp['dstar']:6.3f}  theta {bs['theta']/bp['theta']:6.3f}")
    (CASE / "compare_bl_su2.json").write_text(json.dumps(out, indent=1))
    print("\nsaved: compare_bl_su2.{json,png}")


if __name__ == "__main__":
    main()
