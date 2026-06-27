#!/usr/bin/env python3
"""
準定常 (quasi-steady) 判定ツール (AGENTS.md「準定常確認 (必須)」の実体化)。

`check_convergence.py` が **残差** の収束を見るのに対し、本ツールは **報告する派生量そのもの**
(衝撃位置 shock / 上下非対称 asym / 最大マッハ machmax / 最大静圧 pmax 等) を、run ディレクトリ内の
全 `res_*.h5` スナップショット時系列で評価し、その量が **頭打ち (定常化) したか** を判定する。

目的: 「残差プラトー = 報告していい」と誤認して、**過渡ピーク / ドリフト中の量を定常値として報告する**ことを
防ぐ (本ツールが無かったため、過渡 0.25 の非対称を定常偏りと誤報告した事例があった)。衝撃位置・非対称・
CL/CD・massflux・推力・peak μt 等を「○○だ」と報告する応答は、必ず本ツールの VERDICT を引用すること。

使い方:
  python3 tools/check_quasisteady.py <run_dir> [--quantity shock,asym] [--mesh mesh.h5]
  python3 tools/check_quasisteady.py <run_dir> --tail 0.4 --drift 0.05 --osc 0.10
終了コード: 全量が STEADY なら 0、1つでも DRIFTING/UNSETTLED があれば 1 (CI/スクリプトで使える)。
OSCILLATING (リミットサイクル) は 0 扱いだが「平均±振幅」で報告すること (瞬時値で報告しない)。

判定 (各量の末尾 tail について):
  DRIFTING            : 末尾が単調トレンドで |傾き×tail幅|/平均 > --drift (まだ動いている。run を伸ばす)
  OSCILLATING         : 末尾の振れ (max-min)/平均 > --osc だがトレンド小 (リミットサイクル → 平均±振幅で報告)
  TRANSIENT-UNSETTLED : スナップショットが少なすぎる / 全系列の極値が末尾にある (過渡が減衰しきっていない)
  STEADY              : 末尾が許容内で平坦
"""
import sys, os, glob, argparse, math
import numpy as np
import h5py

GAMMA = 1.4


def find_mesh(run_dir, explicit):
    if explicit:
        return explicit if os.path.isabs(explicit) else os.path.join(run_dir, explicit)
    # 入力メッシュ = res_ でない .h5 で /CELLS/centCoords を持つもの
    for f in sorted(glob.glob(os.path.join(run_dir, '*.h5'))):
        if os.path.basename(f).startswith('res_'):
            continue
        try:
            with h5py.File(f, 'r') as h:
                if '/CELLS/centCoords' in h:
                    return f
        except Exception:
            pass
    return None


def res_files(run_dir):
    fs = glob.glob(os.path.join(run_dir, 'res_*.h5'))
    fs = [f for f in fs if 'wall' not in os.path.basename(f)]
    def step(f):
        d = ''.join(ch for ch in os.path.basename(f) if ch.isdigit())
        return int(d) if d else -1
    return sorted(fs, key=step), [step(f) for f in sorted(fs, key=step)]


def centroids(mesh):
    return h5py.File(mesh, 'r')['/CELLS/centCoords'][:].reshape(-1, 3)


def mirror_index(cc):
    try:
        from scipy.spatial import cKDTree
    except Exception:
        return None, None
    x, y = cc[:, 0], cc[:, 1]
    t = cKDTree(np.column_stack([x, y]))
    d, idx = t.query(np.column_stack([x, -y]), k=1)
    return idx, d.max()


# ---- 量の抽出子: (mesh cc, res VALUE) -> scalar ----
def q_shock(cc, V, yband=1e-3):
    x, y = cc[:, 0], cc[:, 1]
    P = V['P'][:]
    m = np.abs(y) < yband
    o = np.argsort(x[m]); xx = x[m][o]; PP = P[m][o]
    if len(xx) < 10:
        return float('nan')
    p0 = np.median(PP[xx < 0.05]) if np.any(xx < 0.05) else PP[0]
    thr = p0 + 0.15 * (PP.max() - p0)
    idx = np.where(PP > thr)[0]
    return xx[idx[0]] * 1e3 if len(idx) else float('nan')  # mm


def q_machmax(cc, V):
    return float((np.sqrt(V['Ux'][:]**2 + V['Uy'][:]**2 + V['Uz'][:]**2) / V['sonic'][:]).max())


def q_pmax(cc, V):
    return float(V['P'][:].max())


def make_q_asym(cc):
    idx, dmax = mirror_index(cc)
    if idx is None or dmax > 1e-5:   # scipy 無し or 非対称メッシュ
        return None
    def f(cc_, V):
        M = np.sqrt(V['Ux'][:]**2 + V['Uy'][:]**2) / V['sonic'][:]
        return float(np.linalg.norm(M - M[idx]) / max(np.linalg.norm(M), 1e-30))
    return f


def classify(steps, vals, tail_frac, drift_tol, osc_tol, min_snaps):
    s = np.array(steps, float); v = np.array(vals, float)
    good = np.isfinite(v)
    s, v = s[good], v[good]
    n = len(v)
    if n < min_snaps:
        return 'TRANSIENT-UNSETTLED', f"only {n} snapshot(s) (<{min_snaps})", None
    k = max(3, int(math.ceil(tail_frac * n)))
    st, vt = s[-k:], v[-k:]
    mean = float(np.mean(vt)); scale = max(abs(mean), 1e-30)
    span = float(vt.max() - vt.min())
    fluct = span / scale
    # 末尾の線形トレンド (傾き×幅 / 平均)
    if st.max() > st.min():
        slope = np.polyfit(st, vt, 1)[0]
        drift = abs(slope * (st.max() - st.min())) / scale
    else:
        drift = 0.0
    # 全系列の極値が末尾末端にある = まだ成長/減衰中
    extremum_at_end = (np.argmax(v) >= n - 2) or (np.argmin(v) >= n - 2)
    detail = f"tail mean={mean:.4g}  drift={drift*100:.1f}%/tail  fluct={fluct*100:.1f}%"
    amp = span / 2.0
    if drift > drift_tol:
        return 'DRIFTING', detail + ("  (extremum at tail-end)" if extremum_at_end else ""), (mean, amp)
    if fluct > osc_tol:
        return 'OSCILLATING', detail + f"  -> report {mean:.4g} +/- {amp:.2g}", (mean, amp)
    if extremum_at_end and drift > drift_tol * 0.5:
        return 'TRANSIENT-UNSETTLED', detail + "  (still trending at tail-end)", (mean, amp)
    return 'STEADY', detail, (mean, amp)


SEV = {'STEADY': 0, 'OSCILLATING': 1, 'TRANSIENT-UNSETTLED': 2, 'DRIFTING': 3}
BUILTIN = ['shock', 'asym', 'machmax', 'pmax']


def analyze(run_dir, want, tail_frac, drift_tol, osc_tol, min_snaps, mesh_arg):
    mesh = find_mesh(run_dir, mesh_arg)
    if mesh is None:
        print(f"\n=== {run_dir}  -> NO mesh (.h5 with /CELLS/centCoords) ==="); return 3
    fs, steps = res_files(run_dir)
    if len(fs) < 2:
        print(f"\n=== {run_dir}  -> TRANSIENT-UNSETTLED (only {len(fs)} res_*.h5) ==="); return 2
    cc = centroids(mesh)
    extractors = {'shock': q_shock, 'machmax': q_machmax, 'pmax': q_pmax}
    qa = make_q_asym(cc)
    if qa:
        extractors['asym'] = qa
    quantities = [q for q in want if q in extractors]
    series = {q: [] for q in quantities}
    for f in fs:
        V = h5py.File(f, 'r')['VALUE']
        for q in quantities:
            try:
                series[q].append(extractors[q](cc, V))
            except Exception:
                series[q].append(float('nan'))
    worst = 0
    lines = []
    for q in quantities:
        verdict, detail, _ = classify(steps, series[q], tail_frac, drift_tol, osc_tol, min_snaps)
        worst = max(worst, SEV[verdict])
        lines.append(f"  {q:9s}: {detail:55s} {verdict}")
    overall = [k for k, vv in SEV.items() if vv == worst][0]
    print(f"\n=== {run_dir}  [{len(fs)} snapshots, steps {steps[0]}..{steps[-1]}]  -> {overall} ===")
    for l in lines:
        print(l)
    if 'asym' not in quantities and 'asym' in want:
        print("  (asym skipped: mesh not up-down symmetric or scipy missing)")
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dirs', nargs='+')
    ap.add_argument('--quantity', default=','.join(BUILTIN),
                    help='comma list of: ' + ','.join(BUILTIN) + ' (default all applicable)')
    ap.add_argument('--mesh', default=None, help='input mesh h5 (auto-detected if omitted)')
    ap.add_argument('--tail', type=float, default=0.4, help='tail fraction of snapshots for the steadiness check')
    ap.add_argument('--drift', type=float, default=0.05, help='max fractional trend across tail for STEADY')
    ap.add_argument('--osc', type=float, default=0.10, help='max fractional fluctuation across tail for STEADY')
    ap.add_argument('--min-snaps', type=int, default=4, help='min snapshots required to judge')
    args = ap.parse_args()
    want = [q.strip() for q in args.quantity.split(',') if q.strip()]
    worst = 0
    for rd in args.run_dirs:
        worst = max(worst, analyze(rd, want, args.tail, args.drift, args.osc, args.min_snaps, args.mesh))
    print(f"\nOVERALL: {'ALL STEADY' if worst == 0 else 'NOT ALL STEADY (see above)'}")
    sys.exit(0 if worst == 0 else 1)


if __name__ == '__main__':
    main()
