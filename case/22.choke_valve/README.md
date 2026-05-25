# Choke Valve Tools

粗い子午面点列から軸対称チョークバルブの断面積分布 $A(x;\delta)$ を作り、その分布に基づいて CH4 の準一次元圧縮性流れを計算するための最小ツールです。

初版の前提は次のとおりです。

- 流体は CH4
- 理想気体
- 定比熱、定数 $\gamma$、定数 $R$
- 定常、断熱、無摩擦、準一次元
- 正常衝撃波は 1 本まで

## できること

- 写真起こしした `wall` と `needle` の点列から断面積分布を作る
- 開度 `offset` ごとの喉部位置、喉部面積、最小 gap を見る
- 入口全圧 `p0`、入口全温 `t0`、背圧 `pb` を与えて 1 ケースの流れ場を解く
- 開度 sweep に対して流量、出口 Mach、shock 位置、regime をまとめる

## セットアップ

Python 3.12 以上を想定しています。

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

インストール後は次の CLI が使えます。

- `choke-valve`
- `choke-valve-sweep`
- `choke-valve-flow-sweep`

開発中にインストールせず使う場合は `PYTHONPATH=src python3 -m choke_valve...` でも実行できます。

## 入力データ

CSV は `surface,x,r` ヘッダを持ちます。

- `surface`: `wall` または `needle`
- `x`: 軸方向座標
- `r`: 半径方向座標

例は [examples/point_entry_template.csv](examples/point_entry_template.csv) と [examples/point_syashin.csv](examples/point_syashin.csv) です。

```csv
surface,x,r
wall,0.0,80.0
wall,120.0,40.0
wall,600.0,60.0
needle,0.0,33.5
needle,75.0,45.5
needle,220.0,0.0
```

## 単位

- geometry の `x`, `r` は mm
- 面積出力は mm^2
- 流量出力は kg/s
- 圧力入力 `p0`, `pb` は Pa
- 温度入力 `t0` は K

注意:

- solver に与える圧力は絶対圧として扱います
- `pb` は必ず `p0` より小さくしてください
- geometry 側は mm 入力ですが、流量計算では内部で m^2 に変換しています

## サンプリング

CSV から geometry を作るときは、`x_min`, `x_max`, `dx` で等間隔サンプリングします。

- `x_min`: サンプリング開始位置
- `x_max`: サンプリング終了位置
- `dx`: 刻み幅

たとえば `x_min=0`, `x_max=600`, `dx=1.0` なら、ほぼ 601 点で `A(x)` を作ります。

## 1 ケース解く

単発の geometry と flow を見る基本コマンドです。

```bash
choke-valve examples/point_syashin.csv \
  --x-min 0 \
  --x-max 600 \
  --dx 1.0 \
  --offset -100 \
  --p0 5900000 \
  --t0 623.15 \
  --pb 5000000 \
  --plot outputs/example_case.png \
  --output outputs/example_case.csv
```

このコマンドで次を得ます。

- 断面積分布 CSV
- geometry / area / Mach / pressure の 4 段プロット
- 標準エラー出力に `throat_x`, `throat_area`, `regime`, `mass_flow_kg_s`, `shock_x` など

## geometry だけ sweep する

開度に対する喉部位置、喉部面積、最小 gap を見たい場合です。

```bash
choke-valve-sweep examples/point_syashin.csv \
  --x-min 0 \
  --x-max 600 \
  --dx 1.0 \
  --offsets 0 -10 -40 -80 -100 -120 \
  --plot outputs/offset_sweep/summary.png \
  --csv outputs/offset_sweep/summary.csv
```

## flow を開度 sweep する

背圧を固定して、開度ごとの流量や regime を比較する場合です。

```bash
choke-valve-flow-sweep examples/point_syashin.csv \
  --x-min 0 \
  --x-max 600 \
  --dx 1.0 \
  --offsets 0 -10 -40 -80 -100 -120 \
  --p0 5900000 \
  --t0 623.15 \
  --pb 5000000 \
  --plot outputs/opening_sweep/summary.png \
  --mass-flow-plot outputs/opening_sweep/mass_flow_vs_opening.png \
  --csv outputs/opening_sweep/summary.csv
```

この summary では次を出します。

- 流量 vs 開度
- 出口 Mach vs 開度
- shock 位置 vs 開度
- regime vs 開度
- 追加で `開度-流量` の単独図

## 背圧ごとにフォルダを切って一括生成する

`outputs/opening_sweep_by_pb` のように、同じ背圧の結果を同じフォルダへまとめたい場合は次のコマンドを使います。

```bash
mkdir -p outputs/opening_sweep_by_pb && \
pbs='2.9 4.0 5.0 5.5 5.8' && \
offsets='0:-0 m10:-10 m40:-40 m80:-80 m100:-100 m120:-120' && \
for pb in $pbs; do
  pb_tag=${pb/./p}
  dir=outputs/opening_sweep_by_pb/pb${pb_tag}
  mkdir -p "$dir"
  pb_pa=$(awk "BEGIN { printf \"%.0f\", $pb*1000000 }")

  for spec in $offsets; do
    label=${spec%%:*}
    offset=${spec##*:}
    PYTHONPATH=src python3 -m choke_valve.cli examples/point_syashin.csv \
      --x-min 0 \
      --x-max 600 \
      --dx 1.0 \
      --offset "$offset" \
      --p0 5900000 \
      --t0 623.15 \
      --pb "$pb_pa" \
      --plot "$dir/point_syashin_${label}_pb${pb_tag}.png" \
      > "$dir/point_syashin_${label}_pb${pb_tag}.csv" \
      2> "$dir/point_syashin_${label}_pb${pb_tag}.txt" || exit 1
  done

  PYTHONPATH=src python3 -m choke_valve.flow_sweep_cli examples/point_syashin.csv \
    --x-min 0 \
    --x-max 600 \
    --dx 1.0 \
    --offsets 0 -10 -40 -80 -100 -120 \
    --p0 5900000 \
    --t0 623.15 \
    --pb "$pb_pa" \
    --plot "$dir/point_syashin_pb${pb_tag}_opening_summary.png" \
    --mass-flow-plot "$dir/point_syashin_pb${pb_tag}_mass_flow_vs_opening.png" \
    --csv "$dir/point_syashin_pb${pb_tag}_opening_summary.csv" || exit 1
done
```

各フォルダには次が入ります。

- 各開度の個別 4 段図 `point_syashin_...png`
- 各開度の面積分布 CSV `point_syashin_...csv`
- 各開度のログ `point_syashin_...txt`
- 開度 sweep の summary 図 `point_syashin_pb..._opening_summary.png`
- 開度-流量図 `point_syashin_pb..._mass_flow_vs_opening.png`
- 開度 sweep の summary 表 `point_syashin_pb..._opening_summary.csv`

## 出力の読み方

`regime` は次のいずれかです。

- `subsonic`: 全域亜音速
- `choked_supersonic`: チョーク後、内部 shock なし
- `choked_internal_shock`: チョーク後、内部 normal shock あり
- `choked_external_shock`: shock は出口外側
- `choked_internal_shock_unresolved`: 内部 shock を想定したが位置が確定しなかった

主な出力値の意味:

- `throat_x`: 喉部位置 [mm]
- `throat_area`: 喉部面積 [mm^2]
- `minimum_gap`: 最小半径 gap [mm]
- `mass_flow_kg_s`: 質量流量 [kg/s]
- `critical_back_pressure`: チョーク境界の背圧 [Pa]
- `shock_x`: 内部 shock 位置 [mm]

## 既存の例

このリポジトリには、既にいくつかの出力例があります。

- [outputs/opening_sweep_by_pb](outputs/opening_sweep_by_pb): 背圧ごとのフォルダに、各開度の個別図と summary を整理したもの
- [outputs/pb_sweep_by_opening](outputs/pb_sweep_by_opening): 開度ごとのフォルダに、各背圧の個別図を整理したもの

## テスト

```bash
PYTHONPATH=src python3 -m unittest discover -s tests
```

## 注意

- geometry が不自然なら、flow 結果より先に点列を修正してください
- 喉部位置と喉部面積は `dx` に依存します
- 初版では粘性、熱伝達、非理想気体、境界層閉塞は扱っていません

## 関連ファイル

- [plan.md](plan.md): 実装方針
- [point_entry_flow.md](point_entry_flow.md): 点列入力の最小手順
- [src/choke_valve/geometry.py](src/choke_valve/geometry.py): geometry 前処理
- [src/choke_valve/flow.py](src/choke_valve/flow.py): 準一次元 flow solver