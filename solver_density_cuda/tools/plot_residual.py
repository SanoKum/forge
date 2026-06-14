#!/usr/bin/env python3
"""forge の residual_history.csv を semilog でプロットする(堅牢版)。

存在する rms_* 列を自動検出して描く。nan/inf は欠損として扱う。
plot_implicit_residuals.py が特定の列(abs_*/dq_*)を要求して落ちる CSV でも動く。

使い方:
    plot_residual.py [residual_history.csv] [-o out.png] [--phase outer_end] [--all]
"""
import argparse
import csv
import math
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", default="residual_history.csv")
    ap.add_argument("-o", "--output", default=None, help="出力 PNG (既定: 入力と同名 .png)")
    ap.add_argument("--phase", default="outer_end",
                    help="この phase 行のみ描画 (既定 outer_end, --all で全行)")
    ap.add_argument("--all", action="store_true", help="phase で絞らず全行")
    ap.add_argument("--cols", nargs="*", default=None,
                    help="描画する列を明示 (既定: 存在する rms_* 全部)")
    ap.add_argument("--abs", action="store_true",
                    help="絶対残差を描く (既定は初手残差を基準とした相対残差)")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = Path(args.input)
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f"{path}: 行がありません")

    if not args.all and "phase" in rows[0]:
        rows = [r for r in rows if r.get("phase") == args.phase] or rows

    cols = args.cols or [c for c in rows[0]
                         if c.startswith("rms_") and not c.startswith("rms_dq")]
    if not cols:
        cols = [c for c in rows[0] if c.startswith("rms_")]

    def v(s):
        try:
            x = float(s)
            return abs(x) if math.isfinite(x) else float("nan")
        except (TypeError, ValueError):
            return float("nan")

    xs = list(range(len(rows)))
    if "step" in rows[0]:
        xs = [int(float(r["step"])) for r in rows]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for c in cols:
        ys = [v(r.get(c)) for r in rows]
        if not args.abs:
            # 相対残差: 初手 (最初の有限・非ゼロ値) を基準に正規化
            ref = next((y for y in ys if math.isfinite(y) and y > 0.0), None)
            if ref:
                ys = [y / ref for y in ys]
        ax.semilogy(xs, [y + 1e-300 for y in ys], marker=".", ms=3, label=c)
    ax.set_xlabel("step")
    ax.set_ylabel("rms residual" if args.abs else "rms residual / initial (相対残差)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    ax.set_title(str(path.resolve().parent.name))
    fig.tight_layout()

    out = args.output or str(path.with_suffix(".png"))
    fig.savefig(out, dpi=120)
    # 末尾の値も表示しておくと便利
    last = rows[-1]
    summary = "  ".join(f"{c}={last.get(c)}" for c in cols)
    print(f"-> {out}")
    print(f"   last (step {last.get('step','?')}): {summary}")


if __name__ == "__main__":
    main()
