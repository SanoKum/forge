# 08.bump

## 位置づけ

`08.bump` は、既存の line profile 比較フローを使った回帰確認に向くケースである。数値結果を変えないはずの変更を確認するときに優先度が高い。

代表の実行先は次のディレクトリ:

- `case/08.bump/run_slau_mach1.65_4pctHt`

## このケースを使う場面

- 変更前後で密度、圧力、速度分布が変わっていないかを確認したいとき
- 単に完走するかだけでなく、既存ベースラインとの差分も見たいとき
- `20.naca_ml` の確認だけでは不十分なとき

## 現状の設定目安

`solverConfig.yaml` の現状設定では次が基準になる。

- solver: `SLAU`
- `nStep: 2000`
- `outStepInterval: 20`
- 結果比較に使う主なファイル: `res_2000.h5`

## 基本の確認手順

1. `case/08.bump/run_slau_mach1.65_4pctHt` で計算を最後まで流す。
2. `solver_density_cuda/tools/extract_line_profile.py` を使って、無次元 `y = 0.5` の line profile を CSV に書き出す。
3. 同じ run ディレクトリにある保存済みベースライン CSV と比較する。

実行例:

```bash
cd /home/sano/work/forge
python3 solver_density_cuda/tools/extract_line_profile.py \
  --mesh case/08.bump/run_slau_mach1.65_4pctHt/bump_4pct.h5 \
  --result case/08.bump/run_slau_mach1.65_4pctHt/res_2000.h5 \
  --y 0.5 \
  --output case/08.bump/run_slau_mach1.65_4pctHt/line_y0.5_res_2000_latest.csv \
  --compare case/08.bump/run_slau_mach1.65_4pctHt/line_y0.5_res_2000_baseline.csv
```

## 比較ルール

- 比較対象は主に velocity、density、pressure。
- スクリプト既定値の `abs_tol=1e-2`、`rel_tol=1e-4` をまず使う。
- 数値結果を変えない想定の変更なら、この比較が既定許容差内に収まることを期待する。
- アルゴリズムを変える変更なら、差分を失敗扱いせず、どの量がどれだけ変わったかを報告する。

## 既存ファイル

この run ディレクトリには、比較用のファイルが既に置かれている。

- `line_y0.5_res_2000_baseline.csv`
- `line_y0.5_res_2000_latest.csv`
- `check.pvsm`

必要に応じて ParaView で可視化確認も行う。
