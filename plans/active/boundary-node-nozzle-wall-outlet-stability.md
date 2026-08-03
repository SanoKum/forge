# node (median-dual) 2D の傾斜壁・出口コーナー不安定 (ノズル系で初手〜百 step 発散)

## メタ

- **area**: `boundary / discretization`
- **status**: `draft`  <!-- 調査完了・修正未着手 -->
- **related_docs**:
  - `methods/discretization.md` / `methods/boundary.md`
- **related_plans**:
  - [`tooling-nozzle-phase0-foundation.md`](tooling-nozzle-phase0-foundation.md) (発見元: node での再計算が全滅)
  - [`../accepted/turbulence-node-inlet-dirichlet-conserved.md`](../accepted/turbulence-node-inlet-dirichlet-conserved.md) (入口側の類似修正 — 出口版が欠けている疑い)
  - [`../accepted/diffusion-node-boundary-real-distance.md`](../accepted/diffusion-node-boundary-real-distance.md) / [`../accepted/turbulence-node-sst-wallfunction.md`](../accepted/turbulence-node-sst-wallfunction.md)
- **created**: `2026-08-03`
- **owner**: `sano` (調査: Claude 自走)

## 1. 目的 (現象)

forge_design の③ベルノズルメッシュ (構造化 quad 221×65、壁クラスタ AR≦19、品質 PASS、
**cell モードでは同一 msh が完走**) を node で回すと、あらゆる設定で序盤発散する。
2D node の検証実績 (平板 case/26・backstep case/18・case/36) は全て**軸平行壁**であり、
傾斜壁 + 超音速沿い流れ + 超音速流出コーナーを持つノズル系は node 初適用で未踏だった。

## 2. 調査結果 (2026-08-03, 全て 400 step 以下の再現 run at scratchpad)

### 切り分けマトリクス (単独では全て無効)

| 変更 | 結果 |
| --- | --- |
| SST → laminar / convMethod 1→0 / cfl 4→0.5 / axisym→平面 | いずれも step 2–4 で NaN |
| outlet_statPress→outflow / wall→slip / inlet_Pressure→速度入口 | 同上 (残差爆発値がビット同一 = 種は境界種別に非依存) |
| bndFirstOrder 0→1 / implicitRelax 1.0→0.3 / 壁ノード IC 速度ゼロ化 | 同上 |
| 双対幾何検査 (体積・面積・法線外向き) | 異常なし |

### 段階起動 (stage A: convMethod 0 + cfl_pseudo 3 + implicitRelax 0.5) で層が分離

| 構成 | 発散 step | 種の位置 |
| --- | --- | --- |
| 既定設定 (2次, relax 1, cfl 4) | 2–4 | 傾斜壁 (収縮部/ベル) の壁・壁隣接ノード。step 1 で壁 slip ノードの ΔroUx ~ O(10³) |
| stage A + SST (wf 0/1 とも) | 11–15 | **出口∩壁コーナーの outlet ノードで roOmega → 1e19–1e20** (流れ残差は健全 1.7e-4) |
| stage A + laminar | 128 | ベル中腹の壁隣接帯 (x≈47mm)。それまで残差は順調に低下 (1e-4) |
| 参考: 準直線壁 (勾配≦5°) 対照 | 5–6 (SST) | 壁隣接。傾斜が緩いと遅延 → 勾配依存性あり |

### 結論 (二層の独立した node 問題)

1. **壁境界の運動量不安定 (傾斜壁×沿い超音速流)**: 設定を柔らかくすると遅延するが根治しない
   (2→15→128 step)。壁隣接ノードで成長するモードで、slip でも no-slip でも発生。既知の
   未修正バグ「node slip + 接線密度勾配の市松スプリアス流」(node-slip-spurious-flow) と
   同族の可能性。壁半割面流束の傾斜壁での取り扱い (弱形式置換・実距離 over-relax の適用域)
   を疑う。
2. **出口境界の SST スカラー爆発 (コーナー)**: roK/roOmega が出口∩壁コーナーから成長。
   入口側は [turbulence-node-inlet-dirichlet-conserved] で「Dirichlet 保存量整合+残差除外」
   が入っているが、**出口側に同等の処置が無い**。

## 3. 修正方針 (次セッションの作業項目)

1. 最小再現の恒久化: forge_design で 30×15 程度の小メッシュ問題 YAML を作り、数十秒で
   回る再現ケースを `case/40` に置く (再現 run は本 plan 起票時点では scratchpad のみ)。
2. 問題 1: 壁半割面の運動量流束を傾斜壁で点検 (境界一次化の適用範囲、mirror-ghost 退化、
   mesh.cpp/viscousFlux の境界パス)。step 1 の ΔroUx 分布 (壁 slip ノードに集中) を
   指紋として単体切り分け。
3. 問題 2: 出口 (outlet_statPress / outflow) の node スカラー処置を入口修正と同型に
   (保存量整合 or 残差除外 or コーナーノードの優先順位規則)。
4. 検証: 本ノズルメッシュ node SST 12000 step 完走 + cell (run_0002) と η_CF 一致 ~1%。

## 4. 影響・当面の運用

- **設計ツール Phase 0/1 は当面 cell モードで進める** (node 既定の方針 [user-prefers-node-base]
  は維持し、本修正後に node へ切替)。tool 側は `mesh.discretization` で両対応済み。
- 発見元 run: `case/40.nozzle_design_tool/run_0005_bell_smoke_node` (全 step NaN、記録として保持)。

## 5. 変更ログ

- `2026-08-03` — 起票。調査マトリクス完了 (上表)。修正は未着手。
