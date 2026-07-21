# limiter=-1 (off) が Venkatakrishnan フル計算に落ちるバグの修正

## メタ

- **area**: `architecture`
- **status**: `done`
- **related_docs**:
  - `methods/architecture/overview.md`
  - `plans/accepted/architecture-perphase-profiling-hotspot.md` (前段のプロファイリング手法)
- **created**: `2026-07-22`
- **owner**: `CFD Dev`

## 1. 目的・発見経緯

case/39 周期丘 DDES (KEEP スタック) の速度計測 (`FORGE_PROFILE=1`) で `limiter` フェーズが
**全体の 21.8%** を占めていることが判明。KEEP の流束カーネル (`KEEP_d`, convectiveFlux_d.cu) は
`LimiterFields` を一切受け取らず、`limiter_*` は完全に不要 (config も `limiter: -1` で「不関与」と明記)。
[`limiter_d.cu`](../../solver_density_cuda/cuda_forge/limiter_d.cu) の早期 return 判定が
`cfg.limiter == 0` のみで `-1` を素通りさせ、`venkata_limiter` 分岐 (`limiter: 2` と同一コスト) を
全 cell・全変数で計算していた。`solverConfig.cpp` の検証は `limiter` を `{0,1,2,-1}` のみ受理し、
エラーメッセージも "0, 1, 2, or -1" と `-1` を正式な値として扱っている。

## 2. 修正

`limiter_d_wrapper` の早期 return を `cfg.limiter == 0` → `cfg.limiter <= 0` に変更 (1 行)。
`limiter_*` 配列は既存どおり `fill_limiter_d` で 1.0 埋め (= 無制限 2 次) されるため、
これを読む SLAU/HLLE/ROE + `limiter: -1` の組み合わせがあれば数値挙動が変わりうるが、
リポジトリ全体を grep した結果 **`limiter: -1` は KEEP と排他的に使われており** (唯一の例外
`case/38.channel_wmles/run_diag_slau` は solver をコピペで書き換えた際コメントが KEEP のまま残った
明らかな取り違えの使い捨て診断 run で、既に実行済み・再実行対象でない)、実質的な回帰リスクはない。

## 3. 検証

case/39 `run_diag_profile` (xcluster 粗格子, KEEP+dual-time+DDES, 500 step, GPU 単独):

- **速度**: limiter フェーズ 21.8%→0.58%、総時間 **82.97s→68.82s (−17.0%)**。
- **正しさ**: 修正前後の `residual_history.csv` は最初の行から ~1e-7 相対で揺れるが、
  **同一 (修正後) バイナリを 2 回実行した run-to-run 差も同じオーダー** (GPU reduction の
  非決定性、[[cell-atomicadd-nondeterminism]] と同種の node 版ノイズ床)。500 step 後に
  最大 10% まで開くのは DDES の乱流カオス的鋭敏性による増幅で、修正前後の系統差ではない
  (old-vs-new と new-vs-new の初期差が同オーダーであることで確認)。
  → **修正は KEEP の物理に対して不変**、既存の GPU 非決定性ノイズ床の範囲内。

## 4. 影響範囲

`limiter: -1` を使う既存 run 全て (case/38.channel_wmles の KEEP-LES キャンペーン全体、
case/39.periodic_hills の DDES 系列) が今後 ~20% 高速化する。挙動は不変 (§3)。

## 変更ログ

- 2026-07-22: 発見・修正・検証・完了。
