# forge 実装計画インデックス

`.github/plans/` 配下の実装計画 (`.md`) の一覧と状態。
新規 plan を追加・状態を変えた場合は、本表を必ず同期させる。

計画の作り方は [`_template.md`](_template.md) を起点に、命名規約 `<area>-<short-slug>.md`
(長期テーマは `<area>-<theme>/README.md`) で配置する。
運用フロー全体は [`../../AGENTS.md`](../../AGENTS.md) の
「## 開発フロー」節を参照。

## 一覧

| Plan | area | status | related_docs | 概要 |
| --- | --- | --- | --- | --- |
| [diffusion-viscous-shear-flux.md](diffusion-viscous-shear-flux.md) | diffusion | done | [`docs/diffusion/`](../../docs/diffusion/) | 粘性せん断フラックスの修正 (内部面の法線項 `delta_x→delta` + 転置項 + 壁面 no-slip 抗力、軸平行格子で横方向粘性が落ちる不具合、SU2 検証済) |
| [gpu-implicit-plan.md](gpu-implicit-plan.md) | time_integration | in_progress | [`docs/time_integration/`](../../docs/time_integration/) | GPU pseudo-time implicit (matrix-free + 近似 Jacobian + block DPLUR)。制御フロー整理・古典 DPLUR 化・dual-time 受け入れ構造（§0 改訂メモ）。dual-time 本体は後続 |
| [convection-slau2-lowmach.md](convection-slau2-lowmach.md) | convection | done | [`docs/convection/`](../../docs/convection/) | `solver: SLAU2` を追加 (圧力束第 3 項のみ差し替え、質量流束不変)。実装・検証完了。**検証所見**: 低マッハのエネルギー残差フロアは SLAU2 では解消せず (カップリングは共通の質量流束項。SLAU2 は衝撃波ロバスト性向け)。根治は前処理が必要 |
| [time_integration-lowmach-preconditioning.md](time_integration-lowmach-preconditioning.md) | time_integration | done | [`docs/convection/`](../../docs/convection/), [`docs/time_integration/`](../../docs/time_integration/) | Weiss–Smith 低マッハ前処理 (`lowMachPrecond`, 既定 0)。**フラックス散逸前処理 `c_hat→c'` のみ採用** — `case/23.axi_nozzle` 3D ノズルの低マッハ自励振動 (limit cycle) を 20k step 収束ベースで ε=0.15 −32% 減衰 (超音速域不変)。**LHS 固有値前処理・setDT 前処理は block DPLUR で対角優位性を崩し有害→棄却**。設計分析 (§9 `2026-06-08`) で「振動=RHS 駆動」と確定し **Phase 3 (b) Thornber 再構成補正 (`lowMachThornber`) を実装・20k 検証 → 負の結果** (速度ジャンプを減らす=散逸を減らす逆符号の補正で、under-damping のこの症状には無効〜僅かに悪化。機能は opt-in 残置・根治不採用)。根治レバーは Phase1 の圧力散逸 `c'` 側に残り、ε 安定限界超えは完全 (a) に帰着。**Phase 4 (a) 完全 $\Gamma^{-1}A$ 前処理 (`lowMachPrecond=2`) を実装・検証 → 低マッハ自励振動を根治 (採用)**: $\Gamma_c$ ランク1閉形の時間項+物理フラックス FVS+前処理 $\Delta\tau'$、対角 5×5 を倍精度組み立て・解。Phase1 では発散した ε=0.05 を前処理 LHS が安定化し、chamber 振幅 0.882%→**0.087%・定常収束**(超音速域不変)。収束加速は副次でブレークイーブン。plan **done** |
| [architecture-axisymmetric.md](architecture-axisymmetric.md) | architecture | done | [`docs/axisymmetric/`](../../docs/axisymmetric/) | 軸対称モード Phase 1 (B 流儀 r 重み付け + 圧力ソース + 軸 BC、非粘性) |
| [architecture-axisym-nozzle-geometry.md](architecture-axisym-nozzle-geometry.md) | architecture | done | [`docs/axisymmetric/`](../../docs/axisymmetric/) | 軸対称検証ノズルの喉部直後を 5 次多項式で滑らかに接続する geometry 改善 |
| [architecture-rans-sst.md](architecture-rans-sst.md) | architecture | done | [`docs/turbulence/`](../../docs/turbulence/) | Menter SST (k-ω) を explicit 軸対称ノズル (run_0087〜0090) で 4 段階検証完了。advection・diffusion・source (F1/F2 ブレンド)・渦粘性すべて実装済。軸対称 geometric source は子 plan へ |
| [architecture-axisym-sst.md](architecture-axisym-sst.md) | architecture | done | [`docs/turbulence/`](../../docs/turbulence/) | 親 SST の子 plan。フープひずみ $2(u_r/r)^2$ (run_0091) + 生産項の圧縮性補正 `dilatationCorrection` (0:off/1:deviatoric/2:+等方項) を実装・段階検証。run_0093(A+B) で k −12%・vis_turb −14% の膨張減衰を確認。既定値 2 (全 SST ケースに適用) |
| [time_integration-explicit-pointimplicit-sst.md](time_integration-explicit-pointimplicit-sst.md) | time_integration | done | [`docs/time_integration/`](../../docs/time_integration/) | 陽解法 RK (`timeIntegration==1/3`) のスカラー (k/ω) 源項を消散ヤコビアン (`src_jac_k`/`src_jac_omega`) で point-implicit 減衰。block 陰解法と同じ対角を残差増分に適用し、純陽的だと壁近傍 ω で発散する RANS を陽解法 RK で安定化。LES (src_jac=0) は無影響。`case/18.backstep` 2D で検証 |
| [time_integration-implicit-stable-cfl.md](time_integration-implicit-stable-cfl.md) | time_integration | done | [`docs/time_integration/`](../../docs/time_integration/) | 陰解法 (block DPLUR) の安定 `cfl_pseudo` 引き上げ。診断で律速は **k/ω 輸送項の陽的扱い**と確定 (平均流は静止し k/ω のみ発散・`scalarDiffusion=0` で消失)。k/ω 移流+拡散のスペクトル半径を point-implicit 対角 `transport_diag` に陰化 (defect-correction で定常解不変)。`case/26.flat_plate_sst` で安定 `cfl_pseudo` 5–6→**120** (約20倍)、壁法則・Cf は <0.1% 不変。生産項陰化 (#1) は ceiling 不変で棄却、平均流粘性 (#2) は律速でなく無効 |
| [thermophysics-multicomponent-tpgas.md](thermophysics-multicomponent-tpgas.md) | thermophysics | done | [`docs/thermophysics/`](../../docs/thermophysics/) | 多成分 thermally-perfect gas (非反応)。**M1 NASA-9熱力学+SLAU/ROE TP整合 / M2 多成分質量分率輸送 / M3 kinetic theory輸送係数 / M4 化学種拡散+エネルギー結合 を実装・検証完了** (衝撃管CEA照合・拡散係数0.987一致・27.axi_nozzle_plumeでTP-N2がNASA等エントロピーに0.04%一致)。TP境界条件整合化・陰解法Jacobianのγ[ic]化も完了。後続: 組成依存BC・RANS+TP検証 (plan §10) |

## 状態の意味

- `draft` — 設計検討中、実装未着手。
- `in_progress` — 実装または検証進行中。
- `done` — 関連 docs 更新・実装・検証がすべて完了。
- `superseded` — 別 plan に置き換えられた。後継 plan を `related_plans` に明記すること。
