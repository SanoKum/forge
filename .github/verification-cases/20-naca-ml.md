# 20.naca_ml

## 位置づけ

`20.naca_ml` は、コード変更時の既定検証ケースとして扱う。ユーザーから別指定がない限り、まずこのケースで確認する。

既定の実行先は次のディレクトリ:

- `case/20.naca_ml/001.test/run_slau`

`001.test` は実計算用のサブケースで、`000.template` ではなくこちらを優先する。

## 構成

- `case/20.naca_ml/001.test/mesh`: 形状生成とメッシュ生成を行うディレクトリ
- `case/20.naca_ml/001.test/run_slau`: `solverConfig.yaml` と `run_case.sh` がある実行ディレクトリ

## 推奨手順

このケースでは `run_case.sh` を優先する。

```bash
cd /home/sano/work/forge/case/20.naca_ml/001.test/run_slau
./run_case.sh
```

`run_case.sh` は次を順に実行する。

1. `mesh/` 側で NACA 形状入力を生成する。
2. Gmsh で `naca.msh` を作る。
3. `run_slau` に `naca.msh` を配置する。
4. `convertGmshToForge` で `naca.h5` に変換する。
5. `forge` を起動する。

## 目安となる設定

`solverConfig.yaml` の現状設定では次を確認しやすい。

- solver: `SLAU`
- `nStep: 4000`
- `outStepInterval: 1000`
- 出力ファイル: `naca.h5`

## 確認観点

- `run_case.sh` が最後まで通るか。
- `res_0.h5`、`res_1000.h5`、`res_2000.h5` のような途中結果が生成されるか。
- 明らかなクラッシュ、設定読み込み失敗、メッシュ変換失敗が出ていないか。

まずはこのケースで最低限の動作確認を行い、必要があれば `08.bump` や `13.nozzle_H` へ広げる。
