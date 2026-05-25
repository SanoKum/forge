# Point Entry Flow

写真起こしした子午面点列は、まず 1 枚の CSV に打ち込む。

## CSV Format

ヘッダは `surface,x,r` とする。

- `surface`: `wall` または `needle`
- `x`: 軸方向座標
- `r`: 半径方向座標

点列は厳密に昇順でなくてもよいが、同じ `surface` の中で同じ `x` を重複させないこと。

## Entry Order

1. `wall` の点を入口側から出口側へ打つ
2. `needle` の点を入口側から出口側へ打つ
3. `r < 0` を入れない
4. 接触していそうな箇所ほど点を密にする

## Minimal Workflow

1. [examples/point_entry_template.csv](/home/sano/work/forge/case/22.choke_valve/examples/point_entry_template.csv) をコピーして新しい CSV を作る
2. `wall` と `needle` の点列を実測値で置き換える
3. 次のように CLI を実行する

```bash
PYTHONPATH=src python3 -m choke_valve.cli your_points.csv \
  --x-min 0.0 \
  --x-max 20.0 \
  --dx 0.1 \
  --p0 1000000 \
  --t0 300 \
  --pb 700000 \
  --plot outputs/your_points.png
```

## Notes

- `--x-min`, `--x-max`, `--dx` は CSV 入力時のサンプリング条件
- まずは図を見て `r_wall`, `r_needle`, `Area` が不自然でないか確認する
- 幾何が不自然なら、流れ場結果ではなく点列を先に直す