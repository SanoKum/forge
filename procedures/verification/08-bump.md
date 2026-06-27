# 08.bump

## 位置づけ

`08.bump` は、既存の line profile 比較フローを使った回帰確認に向くケースである。数値結果を変えないはずの変更を確認するときに優先度が高い。

低マッハ (亜音速 M~0.5, `bump.h5`, 粘性) と高マッハ (超音速 M=1.65, `bump_4pct.h5`, 非粘性) の 2 構成を、SLAU の陽解法 (RK3, `timeIntegration:3`) と陰解法 (block DPLUR, `timeIntegration:11`+`blockDPLUR:1`) の両方で回帰確認できる。**「y≈0.3 の数値解が基準とずれていないか」と「収束度合いが悪化していないか」を自動判定する検証一式を `case/08.bump/verify/` に用意してある (推奨。下記「低M・高M 回帰検証 (y≈0.3 + 収束度)」節)。**

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

## 低M・高M 回帰検証 (y≈0.3 + 収束度)

コード改良時の回帰確認として、低マッハ・高マッハの 2 構成を SLAU 陽解法/陰解法で流し、次の 2 点を `case/08.bump/verify/` の committed baseline と自動比較する。

1. **数値解の一致**: y≈0.3 のラインに沿った無次元圧力 `P/Pt` とマッハ数 `Mach` が baseline からずれていないか (既定許容: 最大相対差 0.5%)。
2. **収束度合いの非劣化**: 残差 `rms_ro` の floor が baseline の 3 倍以内か、`rms_ro<1e-4` への到達ステップが baseline の 3 倍以内か、発散 (NaN) していないか。

### 構成ファイル (`case/08.bump/verify/`)

- `run_verification.sh` — ランナー。新しい `_work/run_<case>_<method>` を作って `forge` を流し、プロファイル抽出→比較まで一括実行する。
- `gen_config.sh` — 低M/高M × 陽/陰 の `solverConfig.yaml` 生成 (mesh・bcond は正準 run から複製)。
- `extract_profile.py` — y≈0.3 の `P/Pt`・`Mach`・生値を CSV 出力。
- `compare.py` — latest と baseline を比較し PASS/FAIL を表示、差分図 `verify_profiles.png` を生成。
- `baseline_{loM,hiM}_{exp,imp}_y03.csv` — y≈0.3 数値解の基準 (2026-06 時点の収束解)。
- `baseline_convergence.csv` — 残差 floor・`<1e-4` 到達ステップの基準。

### 実行手順

```bash
# 1. native ビルドが stale でないことを確認 (AGENTS.md「実行前のバイナリ鮮度チェック」)
# 2. 検証を流す
cd /home/sano/work/forge/case/08.bump
bash verify/run_verification.sh quick   # 陰解法のみ (loM+hiM, ~70s)。既定
bash verify/run_verification.sh full    # 陽+陰の4本 (loM exp が ~4分)
```

`_work/` は実行ごとに再生成される一時ディレクトリ。結果図は `_work/verify_profiles.png`、各 run の残差は `_work/run_*/residual_history.png`。

### 判定の考え方

- 数値結果を変えない想定の変更なら、全 run が `PASS` (ΔP・ΔM < 0.5%、収束悪化なし) になることを期待する。GPU の atomic 集約順序による run 間ばらつきは ~0.05% 程度で、許容差 0.5% はこのノイズより十分大きく実回帰より小さい。
- アルゴリズムを意図的に変える変更で差分が出る場合は、失敗と決めつけず「どの量がどれだけ変わったか」を報告し、必要なら baseline を更新する。
- baseline 更新時は、変更が妥当 (収束解・無発散・物理的) であることを確認したうえで `baseline_*_y03.csv` と `baseline_convergence.csv` を取り直す。

### 基準値 (2026-06, native build, RTX 3060)

| 構成 | 残差 floor `rms_ro` | `<1e-4` 到達 step |
| --- | --- | --- |
| loM 陽 (cfl 1.0, 40000 step) | 4.8e-6 | 1109 |
| loM 陰 (cfl_pseudo 10, 8000 step) | 2.3e-6 | 481 |
| hiM 陽 (cfl 1.0, 8000 step) | 4.8e-6 | 13 |
| hiM 陰 (cfl_pseudo 10, 3000 step) | 2.1e-6 | 1 |

陽解法 RK3 は安定上限が `cfl_pseudo≈1` (1.5 以上で即発散)、陰解法 block DPLUR は `cfl_pseudo=50` まで安定で陽解法と同一解に約 5 倍速く収束する (検証 2026-06)。
