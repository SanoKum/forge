#!/usr/bin/env python3
"""画像から複数系列を手動 digitize する軽量ツール。

ginput の「Enter で終わらない / 右クリックで点が消える / 系列ごとに較正がずれる」を解消:
  - 軸較正は **起動時に一度だけ**。以降の全系列で同じ変換を共有 (系列間で軸がずれない)。
  - 操作は画面下の **ボタン** (Undo / New series / Finish)。focus に依存せず確実に動く。
  - 右クリックは無視 (点を消さない)。左クリックでのみ点追加。打った点は即マーカー描画され消えない。

使い方:
  python3 digitize.py IMG --xcal X0 X1 --ycal Y0 Y1 [-o out.csv] [--names n1 n2 ...] [--logx] [--logy]
手順:
  1. 起動後タイトルの指示どおり 4 点クリック: x=X0 → x=X1 → y=Y0 → y=Y1 の各位置 (目盛上)。
     値は --xcal/--ycal で渡すのでクリックのみ (端末入力なし=focus が逃げない)。
  2. データ点を左クリックで追加。
  3. [Undo]=直前取消 / [New series]=次系列へ / [Finish]=保存して終了。
     (キーでも可: u=Undo, n=New series, f or Enter=Finish)
出力 (既定 --format wide): 各系列を共通 x グリッドに内挿し
  `<xcol>,<name1>,<name2>,...` 形式で直接出力 → compare_fig3.py がそのまま読める (マージ不要)。
  (--format long で series,x,y の生クリックも可)。確認 PNG も併せて出力。

例 (この Fig.3 p/p0 図, x:0-8cm, y:0.2-0.5。出力は既存 wyslouzil_fig3_1kPa.csv と同じ列):
  python3 digitize.py wyslouzil_fig3_pp0.png --xcal 0 8 --ycal 0.2 0.5 \
          --xcol distance_cm --names pp0_isentrope pp0_cond_1kPa -o wyslouzil_fig3_1kPa_redo.csv
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("TkAgg")          # 失敗時は環境変数 MPLBACKEND=QtAgg 等を試す
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Button

ap = argparse.ArgumentParser()
ap.add_argument("image")
ap.add_argument("--xcal", nargs=2, type=float, required=True, metavar=("X0", "X1"),
                help="x 軸較正 2 点の既知値 (この順にクリック)")
ap.add_argument("--ycal", nargs=2, type=float, required=True, metavar=("Y0", "Y1"),
                help="y 軸較正 2 点の既知値")
ap.add_argument("-o", "--out", default="digitized.csv")
ap.add_argument("--names", nargs="*", default=None, help="系列名 (=ワイド時の y 列名)")
ap.add_argument("--xcol", default="x", help="ワイド出力の x 列名 (例 distance_cm)")
ap.add_argument("--format", choices=["wide", "long"], default="wide",
                help="wide=各系列を共通 x グリッドに内挿し distance,col1,col2,... 形式で直接出力 "
                     "(compare_fig3.py がそのまま読める)。long=series,x,y の生クリック")
ap.add_argument("--logx", action="store_true")
ap.add_argument("--logy", action="store_true")
a = ap.parse_args()

im = mpimg.imread(a.image)
fig, ax = plt.subplots(figsize=(11, 9))
plt.subplots_adjust(bottom=0.12)
ax.imshow(im)

S = {"phase": "calib", "calib": [], "series": [[]], "smarkers": [[]]}

def name_of(i):
    if a.names and i < len(a.names):
        return a.names[i]
    return f"series_{i}"

title = ax.set_title(f"calib 1/4: click the position of x = {a.xcal[0]}")

def build_transform():
    (px0, _), (px1, _), (_, py0), (_, py1) = S["calib"]
    def mk(p0, v0, p1, v1, log):
        if log:
            v0, v1 = np.log10(v0), np.log10(v1)
        m = (v1 - v0) / (p1 - p0)
        b = v0 - m * p0
        return (lambda p: 10 ** (m * p + b)) if log else (lambda p: m * p + b)
    return mk(px0, a.xcal[0], px1, a.xcal[1], a.logx), mk(py0, a.ycal[0], py1, a.ycal[1], a.logy)

def on_click(event):
    if event.inaxes is not ax or event.button != 1:
        return  # 主軸内の左クリックのみ。右/中クリックは無視 → 点を消さない
    x, y = event.xdata, event.ydata
    if S["phase"] == "calib":
        S["calib"].append((x, y))
        ax.plot(x, y, "g+", ms=13, mew=2)
        msgs = [f"calib 2/4: click the position of x = {a.xcal[1]}",
                f"calib 3/4: click the position of y = {a.ycal[0]}",
                f"calib 4/4: click the position of y = {a.ycal[1]}"]
        n = len(S["calib"])
        if n < 4:
            title.set_text(msgs[n - 1])
        else:
            S["phase"] = "data"
            S["fx"], S["fy"] = build_transform()
            title.set_text(f"series '{name_of(0)}': left-click to add points  (Undo / New series / Finish)")
    else:
        S["series"][-1].append((x, y))
        mk, = ax.plot(x, y, "r.", ms=8)
        S["smarkers"][-1].append(mk)
    fig.canvas.draw()

def undo(_=None):
    # 末尾から、点を持つ最後の系列を探して 1 点取り消す (系列境界も跨ぐ)。即時再描画。
    if S["phase"] != "data":
        return
    for k in range(len(S["series"]) - 1, -1, -1):
        if S["series"][k]:
            S["series"][k].pop()
            S["smarkers"][k].pop().remove()
            left = sum(len(s) for s in S["series"])
            print(f"undo: 1 点取消 (残り {left} 点)")
            fig.canvas.draw()
            return
    print("undo: 取り消す点がありません")

def new_series(_=None):
    if S["phase"] != "data" or not S["series"][-1]:
        return
    S["series"].append([])
    S["smarkers"].append([])
    i = len(S["series"]) - 1
    title.set_text(f"series '{name_of(i)}': left-click to add points")
    fig.canvas.draw()

def save_and_close(_=None):
    save()
    plt.close(fig)

def save():
    fx, fy = S.get("fx"), S.get("fy")
    if fx is None:
        print("未較正のため保存しません")
        return
    # 各系列を (data 座標, x 昇順) に変換
    series = []  # (name, xs(sorted), ys)
    for i, s in enumerate(S["series"]):
        if not s:
            continue
        xy = sorted(((fx(px), fy(py)) for px, py in s), key=lambda r: r[0])
        series.append((name_of(i), np.array([p[0] for p in xy]), np.array([p[1] for p in xy])))
    if not series:
        print("点がありません"); return

    if a.format == "long":
        with open(a.out, "w") as f:
            f.write("series,x,y\n")
            for nm, xs, ys in series:
                for x, y in zip(xs, ys):
                    f.write(f"{nm},{x:.4f},{y:.4f}\n")
    else:
        # 共通 x グリッド = 全系列クリック x の和集合を、全系列の **重なり範囲**
        # [max(各min), min(各max)] にトリム。これで各系列を必ず内挿のみ (外挿なし) で評価でき、
        # 空欄のない矩形 CSV になる (どの parser でも読める。外挿の嘘も入らない)。
        lo = max(xs.min() for _, xs, _ in series)
        hi = min(xs.max() for _, xs, _ in series)
        grid = np.unique(np.round(np.concatenate([xs for _, xs, _ in series]), 4))
        grid = grid[(grid >= lo - 1e-9) & (grid <= hi + 1e-9)]
        if grid.size == 0:
            print("系列の x 範囲が重なっていません。--format long で出力してください。"); return
        cols = [(nm, np.interp(grid, xs, ys)) for nm, xs, ys in series]
        with open(a.out, "w") as f:
            f.write(a.xcol + "," + ",".join(nm for nm, _ in cols) + "\n")
            for j, xg in enumerate(grid):
                f.write(f"{xg:.4f}," + ",".join(f"{c[j]:.4f}" for _, c in cols) + "\n")

    npts = sum(len(xs) for _, xs, _ in series)
    print(f"\n{npts} 点 / {len(series)} 系列を {a.out} ({a.format}) に保存:")
    for nm, xs, ys in series:
        print(f"  {nm}: {len(xs)} 点  x[{xs.min():.2f}..{xs.max():.2f}] y[{ys.min():.4f}..{ys.max():.4f}]")
    chk = a.out.replace(".csv", "_check.png")
    fig.savefig(chk, dpi=100)
    print("確認図:", chk)

def on_key(event):
    if event.key in ("enter", "f"):
        save_and_close()
    elif event.key == "u":
        undo()
    elif event.key == "n":
        new_series()

fig.canvas.mpl_connect("button_press_event", on_click)
fig.canvas.mpl_connect("key_press_event", on_key)
bu = Button(fig.add_axes([0.30, 0.02, 0.12, 0.06]), "Undo")
bn = Button(fig.add_axes([0.44, 0.02, 0.16, 0.06]), "New series")
bf = Button(fig.add_axes([0.62, 0.02, 0.12, 0.06]), "Finish")
bu.on_clicked(undo)
bn.on_clicked(new_series)
bf.on_clicked(save_and_close)

print(__doc__)
plt.show()
