# SST 壁関数のエネルギー流束モデル置換 (Kader q_w 壁法則 — 等温壁対応と T_aw 強閉包)

## メタ

- **area**: `turbulence / boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5(f) (弱閉包と「状態適用は暴走」の記録) / §10 (WMLES Kader)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) §3.7
- **related_plans**:
  - [`../accepted/turbulence-sst-thermal-wall-function.md`](../accepted/turbulence-sst-thermal-wall-function.md) (前段: 断熱壁の弱閉包。本 plan はその「やらない」に切り出した将来課題)
  - [`turbulence-wmles-wall-stress.md`](turbulence-wmles-wall-stress.md) (Kader 温度壁法則の既存実装元)
- **created**: `2026-08-11`
- **owner**: `sano`

## 1. 目的

SST automatic wall treatment (`wallTreatmentSST=1`) のエネルギー側を**壁面値の弱閉包**
(accepted plan) から**壁隣接伝導流束のモデル置換**へ拡張する。得られるもの:

1. **等温壁 × 壁関数メッシュの正しい壁熱流束 q_w** — 冷却ノズル壁の熱負荷予測
   (設計ツールの本命量)。現状の等温壁は解像伝導 (μ∂T/∂n) で q_w を評価するため、
   y+≫1 メッシュでは sublayer の温度勾配が解像されず q_w を大きく誤る
   (τ_w の解像勾配過小と同型の欠陥のエネルギー版)。
2. **断熱壁の T_aw「強閉包」**(optional) — T_aw を壁状態として保持しても暴走しない
   構造。弱閉包で壁温出力は足りているため優先度は 1. より低いが、BL 内温度場まで
   回復整合にしたい場合の上位互換。

## 2. 背景 (accepted plan の教訓 — なぜ流束置換が必須か)

T_aw を状態 (node 温度ピン / cell ghost) として課すと、**壁隣接の解像伝導が「壁と第一内点の
恒常的回復オフセット ΔT=r·U_t²/2cp」を実勾配とみなして毎ステップ BL を加熱**し、
T_aw=T_rep+Δ が追随して正帰還で暴走する (実測: node 壁温 1832 K / cell Tt 飽和,
case/40 run_0038/0039 旧版)。SU2 が状態設定で安定なのは、壁関数使用時に**壁隣接の
粘性/伝導流束を壁モデル値で置換**していて解像スキームがピンのジャンプを見ないため。
つまり「状態を持つなら流束もモデル化する」が構造的要件。

## 3. スコープ

- **やる**: `wallTreatmentSST==1` の壁 (まず `wall_isothermal`、次に opt-in で `wall`) の
  エネルギー流束を Kader 型壁法則 q_w に置換する経路。node (W–I 双対面) と cell (壁面
  flux) の両系。config は `sstEnergyWallFunction` (仮, 既定 0) で opt-in。
- **やらない**: WMLES 経路の変更 (共有化はするが挙動不変)。low-Re (`wallTreatmentSST=0`)。
  化学種拡散の壁モデル化。

## 4. 設計方針

**既存資産を最大限流用する** — 必要な機構は全て WMLES 側に実装済みで、SST への配線が主作業:

- **壁法則 (モデル層)**: `wmlesWallModel_d.cu` の Kader 温度壁法則
  (T⁺(y⁺, Pr, Pr_t) ブレンド) をデバイス関数として共通化し、SST 壁関数
  (`ransWallFunction_d.cu`, u_τ は Reichardt Newton で既得) から呼で q_w を得る:
  - 等温壁: q_w = ρ cp u_τ (T_w − T_aw,rep)/T⁺(y⁺_rep) — T_w は BC 指定値、
    T_aw,rep = T_rep + r·U_t²/2cp (圧縮性回復補正込み。Kader+Crocco の標準形)。
  - 断熱壁 (強閉包 opt): q_w = 0 を課しつつ壁状態 T=T_aw を保持 (等温壁ピン機構を
    q_w=0 モデル置換と併用 — accepted plan で暴走した構成に流束置換を足した形)。
- **流束層 (置換機構)**: 既存の `Qw_Wall` マーカ機構をそのまま使う —
  `viscousFlux_d` の W–I 面 xor 置換 (AddQWall: 解像伝導 → q_w·S) と
  `viscousFlux_wall_d` の cell 壁面版は WMLES 等温壁で実装・検証済み。SST 壁関数から
  同じマーカに q_w を書き込むだけで node/cell とも流束置換が効く。
- **状態層**: node 等温壁は既存の `pin_wall_node_temperature_d` + `res_roe` 0 化 +
  DPLUR roe 行 decouple (`iso_wall_flag`) をそのまま使用 (変更なし)。
- **責務 3 層 (nodeWallDirichlet_d.cuh 冒頭の設計) に沿う**: モデル層 = SST 壁関数が
  q_w を計算し `Qw_Wall` を書く / 流束層 = 既存 xor 置換 / 状態層 = 既存ピン。
  新設カーネルは「SST 版 q_w 計算」のみで済む見込み。

**注意点 (今夜の教訓の反映)**:

- W–I 面の置換は**加算でなく置換** (xor 機構) であること — 置換漏れの面が 1 枚でも
  残ると回復オフセットをその面が汲んで暴走が再発する。
- 淀み点 (u_τ→0) で T⁺ 評価が退化する — WMLES 実装のガードを踏襲し、u_τ 床でなく
  q_w→解像値へ連続に退避させる。
- コーナー (壁∩出口等) の複数マーカ半割面での q_w 二重計上に注意 (マルチマーカ方式)。
- explicit / implicit 両対応 (implicit は roe 行 decouple が等温壁機構で既に整合)。

## 5. 実装ステップ

1. methods/turbulence §6.5 に (g) として q_w 壁法則を追記 (Kader 式・適用条件・
   WMLES §10 との共有関係)。
2. Kader T⁺ デバイス関数の共通化 (wmles → 共有ヘッダ)。
3. `ransWallFunction_d.cu`: `wallTreatmentSST==1 && sstEnergyWallFunction==1` の
   `wall_isothermal` bplane で q_w を計算し `Qw_Wall`/bvar (`qwall`) を書く。
4. case 検証 (§6) → 断熱強閉包 (optional) は等温壁が通ってから。

## 6. 検証

- [ ] **平板等温壁の y+ 掃引** (case/26 系): y+ ≈ 1 / 30 / 100 メッシュで q_w(x) が
  相関式 (Colburn/Kader) と low-Re 基準解に整合し、y+ 依存が消えること。SU2
  (壁関数+等温壁) との同一メッシュ比較。
- [ ] **case/40 等温壁変種** (T_w = 1000 K 等): 壁熱流束積分が y+1 low-Re 基準と
  壁関数メッシュで ~数 % 一致 (熱負荷予測の生産検証)。
- [ ] **断熱強閉包** (optional): run_0038 旧版の状態ピン構成 + 流束置換で暴走が消え、
  壁温 = T_aw・BL 内温度場が y+1 基準の回復分布に近づくこと。
- [ ] OFF 回帰: 既定 0 でビット不変。WMLES 経路もビット不変。

## 7. 変更ログ

- `2026-08-11` — 起票 (accepted turbulence-sst-thermal-wall-function の切り出し。
  設計のみ、実装未着手)。
