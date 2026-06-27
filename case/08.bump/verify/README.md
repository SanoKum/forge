# 08.bump 低M・高M 回帰検証

コード改良時に、低マッハ (亜音速 M~0.5) と高マッハ (超音速 M=1.65) の bump 流れで
**y≈0.3 の数値解 (P/Pt, Mach) が基準とずれていないか**、
**収束度合い (残差 floor・収束ステップ・発散有無) が悪化していないか** を自動判定する。

正本の手順・判定基準は [`procedures/verification/08-bump.md`](../../../procedures/verification/08-bump.md) を参照。

## 使い方

```bash
cd /home/sano/work/forge/case/08.bump
bash verify/run_verification.sh quick   # 陰解法のみ (~70s), 既定
bash verify/run_verification.sh full    # 陽+陰の4本 (~5分)
```

## 構成

| ファイル | 役割 | commit |
| --- | --- | --- |
| `run_verification.sh` | ランナー (forge 実行→抽出→比較) | ✓ |
| `gen_config.sh` | solverConfig.yaml 生成 | ✓ |
| `extract_profile.py` | y≈0.3 プロファイル抽出 | ✓ |
| `compare.py` | baseline 比較・PASS/FAIL | ✓ |
| `baseline_*_y03.csv` | y≈0.3 数値解の基準 | ✓ |
| `baseline_convergence.csv` | 残差 floor・収束 step の基準 | ✓ |
| `_work/` | 実行ごとに再生成される一時 run・結果図 | ✗ |
