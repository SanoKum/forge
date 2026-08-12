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
    # 主スナップショット res_<step>.h5 のみ対象 (res_wall_*/res_outlet_*/res_nan_* 等の境界・
    # 診断ファイルは除外)。旧実装は全数字連結で step を作っており、拡張子の「5」や bcond id まで
    # step に混入していた (res_0.h5→5, res_outlet_2_12000.h5→2120005; 2026-08-11 レビュー指摘)。
    import re
    pat = re.compile(r'^res_(\d+)\.h5$')
    pairs = []
    for f in glob.glob(os.path.join(run_dir, 'res_*.h5')):
        m = pat.match(os.path.basename(f))
        if m:
            pairs.append((int(m.group(1)), f))
    pairs.sort()
    return [f for _, f in pairs], [s for s, _ in pairs]


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


# ---- 平板専用: 運動量厚さ theta と Cf/Cf_KS (Karman-Schoenherr 外部相関) ----
# 前提: 壁は y=0 で x>0 が平板 (case/26 系)。x<0 は slip の助走区間なので除外する。
# Cf は「壁関数/壁解像を問わず同じ定義」にするため Reichardt 則の逆解きで u_tau を求める
# (壁解像 y+<1 でも Reichardt は u+=y+ に縮退するので同一定義で扱える)。
# theta は積分核が十分ゼロになる外部流まで積分する (delta99 で打ち切らない)。
KAPPA_WL = 0.41


def _reichardt_uplus(yp):
    return np.log(1.0 + KAPPA_WL * yp) / KAPPA_WL + 7.8 * (
        1.0 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-yp / 3.0))


def _solve_utau(Ut, y, nu):
    utau = np.sqrt(max(nu * Ut / max(y, 1e-30), 1e-30))
    for _ in range(80):
        utau = max(utau, 1e-12)
        f = Ut / utau - _reichardt_uplus(utau * y / nu)
        d = 1e-6 * utau
        df = ((Ut / (utau + d) - _reichardt_uplus((utau + d) * y / nu)) - f) / d
        if abs(df) < 1e-20:
            break
        utau -= f / df
    return max(utau, 0.0)


def cf_karman_schoenherr(re_theta):
    """NASA TMR flat-plate validation と同じ Karman-Schoenherr 相関。"""
    if not np.isfinite(re_theta) or re_theta <= 1.0:
        return float('nan')
    l = np.log10(re_theta)
    return 1.0 / (17.08 * l * l + 25.11 * l + 6.012)


def _plate_column(cc, V, xs):
    """平板上 x=xs 最近傍の壁法線カラムを壁から昇順で返す。"""
    x, y = cc[:, 0], cc[:, 1]
    on = x > 1e-6
    if not np.any(on):
        return None
    xa = np.unique(np.round(x[on], 6))
    xc = xa[np.argmin(np.abs(xa - xs))]
    idx = np.sort(np.where(np.abs(x - xc) < 1e-4)[0])
    o = np.argsort(y[idx])
    return idx[o]


def _theta_and_utau(cc, V, xs, ytop=None):
    col = _plate_column(cc, V, xs)
    if col is None or len(col) < 5:
        return float('nan'), float('nan'), float('nan')
    yy = cc[col, 1]
    ro = V['ro'][:][col]
    ux = V['roUx'][:][col] / ro
    mu = V['vis_lam'][:][col]
    # 壁点 (y=0, u=0) を必ず含める: node は既にある / cell は補う
    if yy[0] > 1e-12:
        yy = np.concatenate(([0.0], yy))
        ux = np.concatenate(([0.0], ux))
        ro = np.concatenate((ro[:1], ro))
        mu = np.concatenate((mu[:1], mu))
        first = 1
    else:
        first = 1                      # 壁点の次が第一 DOF
    m = np.ones(len(yy), dtype=bool) if ytop is None else (yy <= ytop)
    ue = ux[m].max()
    roe = ro[m][int(np.argmax(ux[m]))]
    core = (ro[m] / roe) * (ux[m] / ue) * (1.0 - ux[m] / ue)
    theta = float(np.trapz(core, yy[m]))
    nu1 = mu[first] / ro[first]
    utau = _solve_utau(ux[first], yy[first], nu1)
    return theta, utau, ue


def make_q_flatplate(xs, ytop):
    def q_theta(cc, V):
        return _theta_and_utau(cc, V, xs, ytop)[0]

    def q_cf_retheta(cc, V):
        theta, utau, ue = _theta_and_utau(cc, V, xs, ytop)
        col = _plate_column(cc, V, xs)
        if col is None:
            return float('nan')
        ro = V['ro'][:][col]
        mu = V['vis_lam'][:][col]
        ux = V['roUx'][:][col] / ro
        j = int(np.argmax(ux))
        roe, mue = ro[j], mu[j]
        re_theta = roe * ue * theta / mue
        cf = 2.0 * (utau / ue) ** 2
        ks = cf_karman_schoenherr(re_theta)
        return float(cf / ks) if np.isfinite(ks) and ks > 0 else float('nan')

    return q_theta, q_cf_retheta


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
# 平板専用 (明示指定のときだけ有効): --quantity theta,cf_retheta
FLATPLATE = ['theta', 'cf_retheta']
# 壁モデル応力 (res_wall_<physID>_<step>.h5 の utau/ro から tau = rho*utau^2 の面平均)。
# **壁関数 run 専用** — wallTreatmentSST=0 (壁解像) では utau が出力されず全ゼロになるので
# NOT APPLICABLE を返す。壁解像の壁応力が要るなら接線 traction を別途算出する必要がある。
# 範囲は --wall-xmin/--wall-xmax、対象壁は --wall-phys-id で指定する (ケース固有値をハードコードしない)。
WALL = ['wall_model_tau']


def analyze(run_dir, want, tail_frac, drift_tol, osc_tol, min_snaps, mesh_arg, cf_x=0.6, cf_ytop=None,
            wall_phys_id=None, wall_xmin=None, wall_xmax=None):
    mesh = find_mesh(run_dir, mesh_arg)
    if mesh is None:
        print(f"\n=== {run_dir}  -> NO mesh (.h5 with /CELLS/centCoords) ==="); return 3
    fs, steps = res_files(run_dir)
    if len(fs) < 2:
        print(f"\n=== {run_dir}  -> TRANSIENT-UNSETTLED (only {len(fs)} res_*.h5) ==="); return 2
    cc = centroids(mesh)
    extractors = {'shock': q_shock, 'machmax': q_machmax, 'pmax': q_pmax}
    if any(q in want for q in FLATPLATE):
        qt, qc = make_q_flatplate(cf_x, cf_ytop)
        extractors['theta'] = qt
        extractors['cf_retheta'] = qc
    qa = make_q_asym(cc)
    if qa:
        extractors['asym'] = qa
    # 未知の量を黙って無視すると「評価していないのに ALL STEADY」になるので必ずエラーにする。
    # ただし判定は **固定集合** に対して行う (run 依存で適用不能な asym 等を unknown 扱いしない)。
    known = set(BUILTIN) | set(FLATPLATE) | set(WALL)
    unknown = [q for q in want if q not in known]
    if unknown:
        print(f"\n=== {run_dir}  -> ERROR: unknown quantity {unknown} "
              f"(available: {sorted(known)}) ===")
        return 3
    # 適用不能な量 (非対称メッシュの asym 等) は従来どおり skip する
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
    if 'wall_model_tau' in want:
        import re as _re
        pat_w = _re.compile(r'^res_wall_(\d+)_(\d+)\.h5$')
        by_step = {}
        for f in glob.glob(os.path.join(run_dir, 'res_wall_*.h5')):
            m = pat_w.match(os.path.basename(f))
            if not m:
                continue
            pid, st = int(m.group(1)), int(m.group(2))
            if wall_phys_id is not None and pid != wall_phys_id:
                continue
            by_step.setdefault(st, []).append(f)
        wsteps, wvals = [], []
        allzero = True
        for st in sorted(by_step):
            num = den = 0.0                    # 同一 step の複数壁面は面積(点数)重みで集約
            for f in by_step[st]:
                with h5py.File(f, 'r') as hw:
                    if 'VALUE/utau' not in hw or 'VALUE/ro' not in hw:
                        continue
                    xw = hw['MESH/COORD'][:].reshape(-1, 3)[:, 0]
                    tau = hw['VALUE/ro'][:] * hw['VALUE/utau'][:] ** 2
                    mm = np.ones_like(xw, bool)
                    if wall_xmin is not None:
                        mm &= (xw >= wall_xmin)
                    if wall_xmax is not None:
                        mm &= (xw <= wall_xmax)
                    if mm.sum() == 0:
                        continue
                    if np.any(hw['VALUE/utau'][:][mm] > 0.0):
                        allzero = False
                    num += float(tau[mm].sum()); den += float(mm.sum())
            if den > 0:
                wsteps.append(st); wvals.append(num / den)
        if allzero and wsteps:
            lines.append(f"  {'wall_model_tau':14s}: utau is all zero -> NOT APPLICABLE "
                         f"(wall-function run only; wallTreatmentSST=0 does not output utau)")
            worst = max(worst, 3)
        elif len(wsteps) < 2:
            lines.append(f"  {'wall_model_tau':14s}: only {len(wsteps)} usable res_wall_*.h5 "
                         f"-> TRANSIENT-UNSETTLED")
            worst = max(worst, 2)
        else:
            v, d, _ = classify(wsteps, wvals, tail_frac, drift_tol, osc_tol, min_snaps)
            worst = max(worst, SEV[v])
            lines.append(f"  {'wall_model_tau':14s}: {d:50s} {v}")
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
                    help='comma list of: ' + ','.join(BUILTIN + FLATPLATE + WALL) +
                         ' (default all applicable; theta/cf_retheta are flat-plate specific '
                         'and assume the wall is y=0 with the plate at x>0)')
    ap.add_argument('--cf-x', type=float, default=0.6,
                    help='streamwise station [m] for theta / cf_retheta (default 0.6)')
    ap.add_argument('--wall-phys-id', type=int, default=None,
                    help='wall_model_tau: 対象壁の physID (res_wall_<physID>_*.h5)。既定は全壁を集約')
    ap.add_argument('--wall-xmin', type=float, default=None,
                    help='wall_model_tau: 集約する x の下限 [m] (前縁・淀み域を除くのに使う)')
    ap.add_argument('--wall-xmax', type=float, default=None, help='wall_model_tau: x の上限 [m]')
    ap.add_argument('--cf-ytop', type=float, default=None,
                    help='upper limit [m] of the theta integral (default: full column). '
                         'Use to check sensitivity of theta to the truncation height.')
    ap.add_argument('--mesh', default=None, help='input mesh h5 (auto-detected if omitted)')
    ap.add_argument('--tail', type=float, default=0.4, help='tail fraction of snapshots for the steadiness check')
    ap.add_argument('--drift', type=float, default=0.05, help='max fractional trend across tail for STEADY')
    ap.add_argument('--osc', type=float, default=0.10, help='max fractional fluctuation across tail for STEADY')
    ap.add_argument('--min-snaps', type=int, default=4, help='min snapshots required to judge')
    args = ap.parse_args()
    want = [q.strip() for q in args.quantity.split(',') if q.strip()]
    worst = 0
    for rd in args.run_dirs:
        worst = max(worst, analyze(rd, want, args.tail, args.drift, args.osc, args.min_snaps,
                                   args.mesh, args.cf_x, args.cf_ytop,
                                   args.wall_phys_id, args.wall_xmin, args.wall_xmax))
    print(f"\nOVERALL: {'ALL STEADY' if worst == 0 else 'NOT ALL STEADY (see above)'}")
    sys.exit(0 if worst == 0 else 1)


if __name__ == '__main__':
    main()
