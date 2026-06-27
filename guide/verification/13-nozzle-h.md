# 13.nozzle_H

## 位置づけ

`13.nozzle_H` は、ノズル系の圧縮性流れを追加で確認したいときの標準検証ケースの 1 つである。

この文書では、代表 run として次を使う。

- `case/13.nozzle_H/run_ROE`

`run_ROE` は実際に出力ファイル群が揃っており、確認用の既定 run として扱いやすい。

## 構成

`13.nozzle_H` には複数の run ディレクトリがある。

- `run_AUSM+`
- `run_KEEP_AUSM`
- `run_ROE`
- `run_slau`

明示がない限り、この文書では `run_ROE` を参照する。

## 現状の設定目安

`run_ROE/solverConfig.yaml` の現状設定では次が基準になる。

- solver: `ROE`
- `nStep: 100000`
- `outStepInterval: 2000`
- 出力ファイル: `nozzle_H.h5`

## 確認観点

- `forge` が `run_ROE` から正常に起動し、結果出力を継続できるか。
- `res_0.h5`、`res_2000.h5`、`res_100000.h5` などの出力が揃うか。
- ノズル系ケースでのみ出る設定依存の不具合がないか。

## 実行の考え方

このケースは `20.naca_ml` のような専用 `run_case.sh` を前提にせず、既に用意されたメッシュと設定を使って run ディレクトリから `forge` を実行する確認に向いている。

Docker 例:

```bash
cd /home/sano/work/forge
docker run --rm --gpus all \
  -v /home/sano/work/forge:/workspace \
  -w /workspace/case/13.nozzle_H/run_ROE \
  forge-solver:cuda-dev \
  bash -lc '/workspace/solver_density_cuda/build/forge'
```

実際の確認では、Docker イメージやビルド成果物が既に用意されている前提で、対象 run ディレクトリから `forge` を呼ぶ。

## 使いどころ

- `20.naca_ml` の確認後に、別系統の流れ場でも見たいとき
- `08.bump` とは異なるケースで圧縮性ソルバーの挙動を見たいとき
- ノズル関連の設定や境界条件変更の影響を確認したいとき
