#!/usr/bin/env python3
"""限界サイクル振幅の比較解析 (rho-Y common limiter 診断用)。

限界サイクル振幅 = settled スナップショット間の時間差の混合層内 std。
  A_phi = mean_i std_{mixing}( phi(t_i) - phi(t_{i-1}) )
混合層 = 0.05 < Y_He < 0.95 のセル。定常背景を差し引き純粋な振動成分を抽出する。
rms_ro floor = residual_history.csv の rms_ro を settled 区間で平均。
"""
import sys, glob, os, re
import numpy as np
import h5py

def settled_snaps(run, n=6):
    fs = glob.glob(os.path.join(run, "res_*.h5"))
    def step(f):
        m = re.search(r"res_(\d+)\.h5", f); return int(m.group(1)) if m else -1
    fs = sorted([f for f in fs if step(f) > 0], key=step)
    return fs[-n:]

def rms_ro_floor(run, frac=0.5):
    csv = os.path.join(run, "residual_history.csv")
    if not os.path.exists(csv): return float('nan')
    import csv as _c
    rows=[]
    with open(csv) as fh:
        r=_c.DictReader(fh)
        for row in r:
            if row.get('phase')=='outer_end' or row.get('inner')=='-1':
                try: rows.append(float(row['rms_ro']))
                except: pass
    if not rows: return float('nan')
    rows=np.array(rows); tail=rows[int(len(rows)*frac):]
    tail=tail[np.isfinite(tail)]
    return float(np.mean(tail)) if len(tail) else float('nan')

def analyze(run):
    snaps = settled_snaps(run)
    if len(snaps) < 2:
        return None
    P=[]; Ur=[]; Y0=[]; T=[]
    for f in snaps:
        V=h5py.File(f,'r')['VALUE']
        P.append(V['P'][()]); Ur.append(V['Uy'][()]); Y0.append(V['Y0'][()]); T.append(V['T'][()])
    P=np.array(P); Ur=np.array(Ur); Y0=np.array(Y0); T=np.array(T)
    if not np.isfinite(P).all():
        return {'diverged':True}
    # 混合層マスク (各スナップで Y_He 中間域)。全スナップ平均 Y で固定。
    Ybar=Y0.mean(axis=0)
    mask=(Ybar>0.05)&(Ybar<0.95)
    def amp(arr):
        d=np.diff(arr,axis=0)          # 時間差
        return float(np.mean([np.std(dd[mask]) for dd in d]))
    return {
        'diverged':False,
        'nmix':int(mask.sum()),
        'A_P':amp(P), 'A_ur':amp(Ur), 'A_YHe':amp(Y0), 'A_T':amp(T),
        'rms_ro':rms_ro_floor(run),
        # 過剰拡散チェック: 混合層厚さ (Y_He 0.1-0.9 を含むセル数) と |∇Y_He| の代理 (Y std)
        'mix_thick_cells':int(((Ybar>0.1)&(Ybar<0.9)).sum()),
        'YHe_max':float(Ybar.max()),
    }

if __name__=='__main__':
    print(f"{'run':<40} {'rms_ro':>10} {'A_P':>9} {'A_ur':>8} {'A_YHe':>9} {'A_T':>8} {'mixcells':>8} {'thick':>6}")
    for run in sys.argv[1:]:
        r=analyze(run)
        name=os.path.basename(run.rstrip('/'))
        if r is None:
            print(f"{name:<40} (insufficient snapshots)")
        elif r.get('diverged'):
            print(f"{name:<40} DIVERGED (NaN in snapshot)")
        else:
            print(f"{name:<40} {r['rms_ro']:>10.3e} {r['A_P']:>9.2f} {r['A_ur']:>8.3f} {r['A_YHe']:>9.3e} {r['A_T']:>8.3f} {r['nmix']:>8} {r['mix_thick_cells']:>6}")
