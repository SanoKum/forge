#!/usr/bin/env python3
"""
収束判定ツール (AGENTS.md「収束確認 (必須)」の実体化)。

forge の run ディレクトリの residual_history.csv を読み、**全保存量の残差列**
(rms_ro, rms_roUx, rms_roUy, rms_roUz, rms_roe, RANS時 rms_roK/rms_roOmega) について
初期値・最終値・低下桁数・末尾トレンド (falling/flat/rising) を出し、明確な VERDICT を返す。

目的: 「rms_ro と NaN だけ見て収束と判断する」ことを防ぐ (AGENTS.md 違反の常習を防止)。
結果を「収束した」「一致した」と報告する前に必ず本ツールを通すこと。

使い方:
  python3 tools/check_convergence.py <run_dir> [<run_dir2> ...]
  python3 tools/check_convergence.py --drop 4 --tail 0.2 run_a run_b   # 判定閾値を調整
終了コード: 全 run が PASS なら 0、1つでも未収束/NaN があれば 1 (CI/スクリプトで使える)。

判定基準 (既定):
  - NaN/Inf が無い。
  - 各非ゼロ保存量残差が初期から >= --drop 桁 (既定 3) 低下している。
  - 末尾 (--tail, 既定 20%) が rising でない (flat/falling)。
"all-zero" 列 (例: 2D の rms_roUz) は判定から除外する。
"""
import csv, math, sys, argparse, os

CONSERVED = ['rms_ro', 'rms_roUx', 'rms_roUy', 'rms_roUz', 'rms_roe',
             'rms_roK', 'rms_roOmega']


def load_series(path):
    rows = [r for r in csv.DictReader(open(path)) if r.get('phase') == 'outer_end']
    if not rows:
        # outer_end が無い構成 (純 explicit) では最終 inner / 全行を使う
        rows = list(csv.DictReader(open(path)))
    cols = {}
    for c in CONSERVED:
        if rows and c in rows[0]:
            cols[c] = [float(r[c]) for r in rows]
    return rows, cols


def analyze(path, min_drop, tail_frac):
    rows, cols = load_series(path)
    if not rows:
        return None
    laststep = rows[-1].get('step', '?')
    report = {}
    ok = True
    any_nan = False
    any_stalled = False
    any_converging = False
    for c, ser in cols.items():
        nz = [v for v in ser if v != 0.0]
        if not nz:
            report[c] = ('all-zero (inactive, skip)', True)
            continue
        if any(math.isnan(v) or math.isinf(v) for v in ser):
            report[c] = ('NaN/Inf present  <-- DIVERGED', False)
            ok = False; any_nan = True
            continue
        init, fin = ser[0], ser[-1]
        n = len(ser)
        a = ser[int(n * (1 - tail_frac)):]
        b = ser[int(n * (1 - 2 * tail_frac)):int(n * (1 - tail_frac))] or a
        ma = sum(abs(x) for x in a) / len(a)
        mb = sum(abs(x) for x in b) / len(b)
        # rising = 窓平均が 5% 増 **かつ** プラトーを実際に離脱している (系列最小の
        # 2 倍超)。後者のガードが無いと、深く収束したプラトー (例: 3.4 桁低下後の
        # 1e-7 台) のリミットサイクル呼吸 ±5% を「発散傾向」と誤判定する
        # (run 毎の非決定性で同一設定が PASS/NOT CONVERGED に割れる実害が出た
        # 2026-08-15)。本物のリバウンド発散は最小値の 2 倍を速やかに超えるので
        # 検出力は保たれる。
        smin = min(abs(x) for x in ser if x != 0.0) if any(x != 0.0 for x in ser) else 0.0
        trend = ('rising' if (ma > mb * 1.05 and ma > 2.0 * smin)
                 else ('flat' if ma > mb * 0.9 else 'falling'))

        # init==0 (例: アライン格子で Uy が初期厳密 0): 相対低下は無意味。最大初期残差から
        # 立ち上がっていなければ inactive 扱い (rising でなければ収束済みとみなす)。
        if init == 0.0:
            col_ok = trend != 'rising'
            report[c] = (f"init=0 fin={fin:.2e} (near-zero from start) {trend:7s}"
                         f"{'' if col_ok else '  <-- rising!'}", col_ok)
            ok = ok and col_ok
            continue

        drop = math.log10(abs(init) / abs(fin)) if fin != 0 else float('inf')
        col_ok = drop >= min_drop and trend != 'rising'
        ok = ok and col_ok
        # 未達の理由を区別: falling=収束途中(あと steps)、flat=停滞、rising=発散傾向
        status = ''
        if not col_ok:
            if trend == 'rising':   status = '  <-- RISING (divergent)'
            elif trend == 'flat':   status = '  <-- STALLED (plateau)'; any_stalled = True
            else:                   status = '  <-- still converging'; any_converging = True
        report[c] = (f"init={init:.2e} fin={fin:.2e} drop={drop:4.1f}dec {trend:7s}{status}", col_ok)
    return laststep, report, ok, any_nan, any_stalled, any_converging


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dirs', nargs='+')
    ap.add_argument('--drop', type=float, default=3.0, help='required orders-of-magnitude drop')
    ap.add_argument('--tail', type=float, default=0.2, help='tail fraction for trend check')
    args = ap.parse_args()

    all_pass = True
    for rd in args.run_dirs:
        path = rd if rd.endswith('.csv') else os.path.join(rd, 'residual_history.csv')
        if not os.path.exists(path):
            print(f"[{rd}] NO residual_history.csv"); all_pass = False; continue
        res = analyze(path, args.drop, args.tail)
        if res is None:
            print(f"[{rd}] empty residual file"); all_pass = False; continue
        laststep, report, ok, any_nan, any_stalled, any_converging = res
        if ok:
            verdict = 'PASS (converged)'
        elif any_nan:
            verdict = 'DIVERGED (NaN/Inf)'
        elif any_stalled:
            verdict = 'NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)'
        elif any_converging:
            verdict = 'NOT CONVERGED (still converging — run more steps)'
        else:
            verdict = 'NOT CONVERGED'
        print(f"\n=== {rd}  [last step {laststep}]  -> {verdict} ===")
        for c, (msg, _) in report.items():
            print(f"  {c:12s}: {msg}")
        all_pass = all_pass and ok
    print(f"\nOVERALL: {'ALL PASS' if all_pass else 'CHECK FAILURES ABOVE'}")
    sys.exit(0 if all_pass else 1)


if __name__ == '__main__':
    main()
