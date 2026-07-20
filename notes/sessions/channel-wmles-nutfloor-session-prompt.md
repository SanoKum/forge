# 引き継ぎプロンプト: WMLES チャネル較正 → νt 床の判別実験 (2026-07-20 セッション末)

> 新セッションの冒頭にこのファイルを読ませる。前身は
> [iddes-phase2-session-prompt.md](iddes-phase2-session-prompt.md)。
> 本セッションで「等温壁検証・体積力 config・周期チャネル Reτ550 較正キャンペーン
> (純 LES → IDDES → 2 倍細分判別)」まで完了し、律速が shield/界面と確定した。
> **次タスクは νt 床の判別実験 — ユーザー指示により「実験のみ」。恒久実装は結果を見て
> ユーザーと判断するまで禁止。**

## 完了済み (再調査・再実装しないこと)

主要 commit (feature/median-dual-3d):
- `ef3a8d1` cell 等温壁 ghost T バグ修正 (鏡像外挿, 流束 1/2 バグ)。純伝導厳密解で機械精度検証
- `8f82826` node 等温壁の壁ノード T ピン + **block-DPLUR エネルギー行 decouple (`iso_wall_flag`)**
  — decouple 無しの状態ピンは implicit 数 step 発散 (実測)。implicit 実用 cfl_pseudo 上限 ~5
- `96a93d6` 移動等温壁 (ghost 2U_w−U_c) + Couette 厳密解検証 / SST EWT×等温壁の熱壁関数
  ギャップ実測 (**q_w +45% 過大** — Kader 移植は未着手の別タスク)
- `10da5d7`/`7eda32f` nodeWallDirichlet_d 集約リファクタ + thermo_state_at_T 共通化 (挙動不変)
- `9439fb8` **体積力 config** (`bodyForce: [fx,fy,fz]`, Poiseuille 厳密解検証済)
- `fb955c1`〜 **case/38 周期チャネル Reτ550** (node, KEEP+ES precond σ0.02, WALE→IDDES, 体積力駆動)
- `cbd42f5` **node k-Dirichlet の explicit 経路バグ修正** (状態ピンを applyBconds 位相へ) +
  `nodeOmegaWfDirichlet` (第一内層 ω ピン, **opt-in 既定 OFF** — §3.7 凹コーナー race 懸念のため)
- `51ac11e` wf ピンノードの **SST Bradshaw リミッタ迂回** (νt=ρk/ω)。ジャンプ→S 大→νt 頭打ち
  →ジャンプ維持の正帰還ロックを切断 (log 則張替でも再形成される構造不安定を run_0012 で実証済)

## チャネル較正の確定事実 (case/38 README run 表が一次情報)

| 構成 | log 層誤差 (基準 \|5\|%) |
| --- | --- |
| 純 LES (WALE) + 代数壁応力, σ0.05 | +42% |
| 同, σ0.02 (limiter は KEEP スタックで不関与) | +30% |
| SST-IDDES + k/ω ピン + リミッタ迂回 (粗格子 Δ⁺46-72) | **+12.8%** (`run_0014`, 図 `channel_ret550_analysis.png`) |
| 同・2 倍細分 (Δ⁺23-36) | **+22.2%** (`run_0015`) — **悪化** |

- **解像応力の引継ぎ高さ y_to ≈ 2.4Δy** (粗 110⁺/細 56⁺で格子比例を実証)。実効フィルタ ≈2Δ (スペクトル)
- 細分で悪化した理由: h_max=Δx が半減 → **f_B の RANS 帯 (0.53h_max) が第 1 内点にすら届かず
  shield 消滅**。一様格子 IDDES の自己矛盾 (細分するほど shield が薄る)。IDDES WMLES 分岐の
  設計想定は壁法線解像格子 (y₁⁺~1) であり、我々の一様格子運用は設計封筒外
- 誤差は y_to より下の欠損帯で生成され、上は log 則と平行 (コアの LES は健全: 解像率 ~100%,
  モデル k/解像 TKE ≈18%, コア νt ≈ 理論 SGS 値)
- u_τ は全構成で 3.77-3.81 に静定 (目標 3.85, 大域収支は壁モデル経由で閉じる)

## 次タスク: νt 床の判別実験 (実験のみ!)

**目的**: 「解像できない帯 (d < y_to) に壁法則整合 νt を敷けば log 層が閉じるか」の白黒。
**ユーザー合意事項: 恒久実装はしない。実験用の一時変更 (フラグ/ローカルパッチ) で判定し、
結果を見てから設計を相談する。**

- 対象: 粗格子 (`run_0014` の続き)。床 = 内点第 1-2 ノード (一様格子で 2.4Δy より下は常に 2 ノード。
  「2」はハードコードでなく d < C·h_max, C≈1.55 の判定式で。C はアスペクト比依存の実測値)
- 床の値: mixing-length νt = ρ(κd)²|S| または ρκd·u_τ (wf の u_τ 利用可)。既存 wf ピン機構
  (roK_wf/roOmega_wf を第 2 ノードにも書く) の一時拡張が最小実装
- **合否指紋**: 床適用後の応力分解で y_to が 2.4Δy に**留まれば成功** (log 誤差は界面 1 セル分の
  数 % に縮むはず) / y_to が**上へ逃げれば追いかけっこ** (界面抑制 — その場合は界面強制系の
  研究課題になる。深追いせずユーザーへ報告)
- 対抗馬 (基準線): 壁法線ストレッチ格子 (y₁⁺~1-5) の正統 IDDES — 文献比較可能な基準線として価値あり

## 罠 (今日ハマったもの・恒久ルール)

- **u_τ・派生量はプラトー到達までに評価・診断しない** (ユーザーから 2 回指摘。発展途中の
  プロファイル形状からスキーム診断をしない)。チャネルのスピンアップは外側時定数 ~100FT 級
- **node SST の res 出力は壁ノード roOmega に inf/1e32 ゴミを含む** (ピン前位相の出力)。
  restart/interp の seed にすると step1 で即 NaN。**k/ω は seed せず一様 IC (k=1, ω=5000) で
  再平衡させるのが安定** (速度場だけ移植)。[[forge-sst-restart-nonfidelity]] に追記済み
- **壁 bvar 出力 (res_<wall>_*.h5) の ro/Psb は field と別位相で信用不可** (τ_w 診断を 2 回汚した)。
  τ_w は field の壁ノード状態 + utau から組む。bvar 出力の位相整合は未修正の小バグ (残タスク)
- 新 run dir 作成時: seed 読取後に旧 res_*/log を**必ず rm** [[run-dir-copy-clean-res]]
- interp_field.py は node の VIZMESH 形式で落ちる → 一様直方体格子なら NN 移植を直接書く
- forge の kill は `pkill -x forge` (`-f` は自分のシェルにマッチして自爆する)
- 回帰は `tests/regression/run_regression.py --case all` (4 ケース)。ベースラインは現行
- 大域運動量は必ず釣り合う (準定常なら)。「不均衡」に見えたら測り方 (bvar 位相・蓄積項) を疑う
- KEEP スタックでは limiter 設定は流束に不関与 (ES recon-jump はリミタ無し線形再構成)

## 他の残タスク (優先度は νt 床実験の後)

1. SST EWT の熱壁関数 (Kader 移植, wallLaw_d.cuh 部品流用) — 等温壁 q_w +45% ギャップの解消
2. node slip 接線密度勾配バグ (ghostless plan の slip 再実装に移管済, 回帰ケース指定済)
3. 壁 bvar 出力の位相整合
4. `nodeOmegaWfDirichlet` の凹コーナー再検証 (backstep) — default 化するなら必須
5. iddes plan §4.6 の記録同期 (チャネル結果を IDDES 側 plan にも反映)

## 関連 memory (自動ロードされる)

develop-flow-before-reporting (再発事例追記済) / forge-sst-restart-nonfidelity (roOmega inf 追記済) /
run-dir-copy-clean-res / node-slip-spurious-flow / user-prefers-node-base / verify-node-and-cell-both
