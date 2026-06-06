# forge テスト

`forge` の回帰確認用ハーネス。**GPU 必須**。CI は未整備のため、現状はローカル
(コミット前 / マージ前) で手動実行する位置づけ。

## 構成

```
tests/regression/
├── run_regression.py            # ハーネス本体
├── cases/<name>.json            # ケース設定 (対象ケース・縮小ステップ・許容差)
└── baselines/<name>/residual_history.csv   # committed golden (outer_begin 行)
```

回帰の判定は `residual_history.csv` の `outer_begin` 行 (step ごと 1 行) を
baseline と相対許容差 (`rtol`) で比較して行う。GPU の atomicAdd 由来の浮動小数
非決定性があるため完全一致は求めず、桁違いの回帰・発散・NaN を捕捉する設計。

## 前提

- GPU が利用可能であること。
- forge がビルド済みであること (`tools/build_native_wsl.sh` または Docker の `tools/build.sh`)。
- 対象ケースの雛型 run にメッシュ (`naca.h5`) があること。無ければ先に
  `case/20.naca_ml/001.test/run_slau/run_case.sh` 等でメッシュを生成する。
- Python に `numpy` が入っていること。

実行は AGENTS.md の共通ルールに従い、既存 run を直接使わない。ハーネスは
ケース内に一時 run (`run_regression_*`) を複製して実行し、既定では実行後に削除する
(残すには `--keep`)。`run_regression_*` は `.gitignore` 済み。

## 使い方

```bash
cd /home/sano/work/forge/solver_density_cuda

# 1) スモーク (短ステップ、有限性のみ。baseline 不要・高速)
python3 tests/regression/run_regression.py --smoke

# 2) 回帰比較 (baseline と比較し PASS/FAIL)
python3 tests/regression/run_regression.py

# 3) baseline 更新 (数値を変える意図的変更の後に実行)
python3 tests/regression/run_regression.py --update-baseline

# 4) 実行せず既存 CSV を比較 (GPU 不要・比較ロジック確認用)
python3 tests/regression/run_regression.py --compare-only path/to/residual_history.csv
```

runner は既定 `auto` (native バイナリがあれば native、無ければ Docker)。
明示するには `--runner native` / `--runner docker`。Docker イメージ名は
`--image` で指定 (既定 `forge-solver:cuda-dev`)。

終了コード: `0`=PASS、`1`=回帰検出、`2`=セットアップ/実行エラー。

## ケースの追加

1. `cases/<name>.json` を作る (既存の `naca_slau.json` を参照)。
2. `python3 tests/regression/run_regression.py --case <name> --update-baseline` で
   golden を生成し、`baselines/<name>/residual_history.csv` をコミットする。

## 既知の制約・今後の課題

- GPU 必須のため自動 CI 未対応 (self-hosted GPU runner または build-only CI は backlog)。
- 回帰指標は残差履歴のみ。壁面分布や line profile の回帰比較は未対応 (backlog)。
