# 粘性せん断フラックスの修正 (内部面の法線項 + 壁面 no-slip 抗力)

## メタ

- **area**: `diffusion`
- **status**: `done`
- **related_docs**:
  - `docs/diffusion/theory.md`
  - `docs/diffusion/implementation.md`
  - `docs/boundary/implementation.md`
- **created**: `2026-06-06`
- **owner**: `未定`

## 1. 目的

粘性せん断フラックスが軸平行格子で法線方向の物理粘性を過小評価する不具合を修正する。
原因は同一クラスの2箇所:(A) 内部面カーネル `viscousFlux_d` が法線項に**成分** `delta_x`
を使い(スカラー `delta` が正)、(B) 壁面カーネル `viscousFlux_wall_d` が `*sxx` を使う。
いずれも軸平行 $y$ 法線面で流れ方向運動量の横方向拡散 $\mu\,\partial u_x/\partial y$ が落ちる。
完了時には、平板チャネル (case 24) で SU2 laminar 参照解と一致する放物線境界層が形成され、
壁面せん断 `twall_x` が物理的に非ゼロとなる。

## 2. スコープ

- **やる**:
  - 内部面 `viscousFlux_d` の運動量法線項を `delta_x/delta_y/delta_z` → スカラー `delta`
    に統一(熱伝導束と同型)。Laplacian 接線補正の勾配添字も同一成分 `dUxd{x,y,z}f` へ是正。
  - 完全な Newton 応力にするため**転置項** $\mu\,\partial u_j/\partial x_i\,S_j$ を内部面・壁面とも有効化
    (面平均/セル中心勾配にフル `S` を内積)。
  - 壁面 `viscousFlux_wall_d` の運動量せん断評価を内部面と同じ $\tau_{ij}S_j$ 形へ修正
    (`*sxx` → 面積 `sss`)。壁面出力 `twall_*` の整合。
  - case 24 (平板) で SU2 参照との再検証。
- **やらない**:
  - case 25 (軸対称パイプ) の SU2 再検証、および RANS/LES 検証ケースの再ベースライン
    (フォローアップ。本 plan では平板チャネルでの妥当性確認まで)。
  - 温度依存粘性 (Sutherland) の導入 (別 TODO)。
  - 低マッハ前処理・収束改善 (別 plan)。

## 3. 関連 docs と前提

- 正しい離散化: [`docs/diffusion/theory.md`](../../docs/diffusion/theory.md)「Over-relaxed 離散化」「壁面寄与」。
- 現状実装と不具合: [`docs/diffusion/implementation.md`](../../docs/diffusion/implementation.md) `viscousFlux_wall_d` 節。
- 壁 BC のゴースト構成 (速度反転ミラー): [`docs/boundary/implementation.md`](../../docs/boundary/implementation.md)、
  `cuda_forge/boundaryCond_d.cu` の `wall_d` (`Ux[ig]=-Ux[ic]`)。対流側は no-penetration の
  みを保証し、接線 no-slip は粘性束に依存する点が前提。

## 4. 設計方針

正本は熱伝導束 ([`L123-124`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L123))。これは
スカラー `delta` の法線項 + 同成分勾配の接線補正 `k` で正しく組まれている。運動量も同型に揃える。

### (A) 内部面 `viscousFlux_d` ([`L105-121`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L105))

現状 → 修正:

```cpp
// 現状 (tau_x)
tau_x  = mu*((Ux1-Ux0)/dcc)*delta_x;                  // 成分 delta_x: y法線面で 0
tau_x += mu*(dUxdxf*k_x + dUydxf*k_y + dUzdxf*k_z);   // 転置勾配 ∂u_j/∂x (k で誤用)
// 修正後 (完全 Newton 応力 tau_ij S_j)
tau_x  = mu*((Ux1-Ux0)/dcc)*delta;                    // (1) Laplacian 法線 (スカラー delta)
tau_x += mu*(dUxdxf*k_x + dUxdyf*k_y + dUxdzf*k_z);   // (2) Laplacian 接線 (同成分 ∂u_x/∂x_j)
tau_x += mu*(dUxdxf*sxx + dUydxf*syy + dUzdxf*szz);   // (3) 転置 ∂u_j/∂x にフル S
tau_x += -mu*2.0/3.0*divu*sxx;                        // (4) 発散項
```

`tau_y`/`tau_z` は成分を入れ替えて同型。直交格子では `k=0` なので (2) は無影響だが、
法線項 `delta_x→delta` が本質(`delta = delta_x+delta_y+delta_z` は直交のみ成立)。
(3) の転置項は純せん断 $u_x=\gamma y$ の x 法線面 $F_y=\mu\gamma S$ を回復するために必要。

### (B) 壁面 `viscousFlux_wall_d` ([`L256-266`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L256))

```cpp
// 現状: x運動量に sxx → 水平壁 (sxx≈0) で no-slip 抗力ゼロ
tau_x = mu*((Ux[ig]-Ux[ic])/dcc)*sxx;
// 修正後: ミラーゴースト u^g=-u^c で d∥S なので delta=sss、k=0
tau_x  = mu*((Ux[ig]-Ux[ic])/dcc)*sss;               // 法線項 (面積 sss)
tau_x += mu*(dUxdxf*sxx + dUydxf*syy + dUzdxf*szz);  // 転置項 (セル中心勾配にフル S)
tau_x += -mu*2.0/3.0*divu*sxx;                       // 発散項
```

発散項 $-\tfrac{2}{3}\mu(\nabla\!\cdot\!\mathbf{u})S_i$ は成分 `s*x/s*y/s*z` で残す。
`twall_x_b/twall_y_b/twall_z_b` は修正後の応力ベクトルから再計算し、ストリーム方向せん断が
`twall_x` に正しく現れることを確認する。

数式の正本は theory.md「Over-relaxed 離散化」「壁面寄与」を参照し、本 plan には差分のみ記載する。

## 5. 実装ステップ

1. `cuda_forge/viscousFlux_d.cu` の内部面 `viscousFlux_d` (L105-121) の法線項を
   `delta` に、接線項を同成分勾配に是正 (§4-A)。`tau_y/tau_z` も同様。
2. 同ファイルの壁面 `viscousFlux_wall_d` (L260-268) を §4-B の形へ修正。`twall_*_b` も整合。
3. `wall` / `wall_isothermal` 共通カーネルのため、熱伝導束 (`heatflux`, Ts[ig]) には
   手を入れず運動量項のみ変更することを確認。
4. ビルド (`solver_density_cuda/build`、`docker run ... forge` 経路) が通ることを確認。

## 6. 検証

- **単体 / ビルド**: `forge` が Docker (`forge-solver:cuda-dev`, `tools/build.sh`) でエラーなく完走 ✓。
- **検証ケース**: `case/24.laminar_channel_bl/run_0003_slau_laminar_viscfix/` で再計算 ✓。
  参照は SU2 laminar 解 [`run_0002_su2_laminar/`](../../case/24.laminar_channel_bl/run_0002_su2_laminar/)
  (`compare_outlet.py` / `outlet_compare.png`)。
- **判定基準と結果** (修正前 → 修正後 / SU2):
  - 壁隣接セル Ux: 17.98 → **0.24** m/s (≈ no-slip) ✓
  - 中心/平均速度比: 1.18 → **1.53** (SU2 1.51、放物線) ✓
  - 質量流量: 0.312 → **0.141** (SU2 0.155、差 約9%、修正前は約2倍) ✓
  - `twall_x`: ≡0 → **平均 -29.9 Pa** (非ゼロ・ストリーム方向)、`twall_y` 平均 ≈0 ✓
  - `residual_history.png` 生成済み ✓
- **長時間収束 + 摩擦則検証** (`run_0004_slau_laminar_viscfix_long`, 150k step):
  流量は 40k step で 0.14339 に完全収束。完全発達摩擦則 $\bar u=(H^2/12\mu)(-dp/dx)$ に対し
  forge は **0.5% 適合**(ū=12.28 vs 予測 12.22)。ROE スキーム (`run_0005`) でも 0.14349 と
  一致し、対流スキーム非依存。
- **SU2 リファレンスの再収束**: 当初 SU2 (`run_0002`, 前処理なし 20k iter) は 0.155 で摩擦則を
  14% 破る収束不足だった。`LOW_MACH_PREC=YES` + 長時間 (`run_0006`, 60000 iter 完走) で
  **0.14321 に収束**し、Ucenter/Uavg=1.501(放物線理論1.5と0.1%)。
  **forge 0.14339 と収束 SU2 0.14321 は 0.124% 一致** → 修正は正しい。
  (rms[RhoE]=-1.31 は低マッハ密度ベースの限界で運動量収束には無影響。)
- **case 25 (軸対称パイプ) 検証**: `run_0002_slau_laminar_viscfix` (60k step) で Ux_wall=0.028 m/s
  (no-slip 復活)、HP 壁面摩擦則 τ_wall=R/2·(-dP/dx) を **0.18% 適合**、mdot=0.000502 kg/s。
  SU2 参照 (`run_0003_su2_laminar`, LOW_MACH_PREC, 60k) は mdot=0.000509、Uc/Ua=2.007。
  **forge vs SU2 mdot 差: 1.42%**。forge Uc/Ua=2.11 (HP理論2.0 から5%) は run_0001 プラグ流
  起点からの収束過渡である可能性。RANS/LES の再ベースラインは未実施 (別フォローアップ)。

## 7. 影響範囲

- 触るファイル: `solver_density_cuda/cuda_forge/viscousFlux_d.cu` (内部面カーネル + 壁面カーネル)。
- 既存ケースへの影響: **粘性を含む全ケース**で法線方向の物理粘性が増える(内部面修正)。
  特に軸平行格子のケースで挙動が変わる。加えて全粘性壁ケースの壁摩擦・$y^+$・抗力係数が
  動く(壁面修正)。RANS/LES 含む検証ケースは広く再ベースラインが必要。
- docs: `docs/diffusion/theory.md`・`docs/diffusion/implementation.md` 更新済み、
  実装完了後に `docs/index.md` との整合を確認。

## 8. 完了条件

- [x] 関連 `docs/diffusion/theory.md` 更新済み
- [x] 関連 `docs/diffusion/implementation.md` 更新済み
- [x] 実装・検証完了 (本 plan の §6、case 24 で SU2 と一致を確認)
- [x] `.github/plans/README.md` の状態を `done` に更新
- [x] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `2026-06-06` — 初稿。case 24 を SU2 8.5.0 laminar NS で再計算し、forge 出口分布が
  プラグ流 (壁すべり ~18 m/s・流量約2倍・Mach 0.090 vs 0.057) になることを確認。
  原因を `viscousFlux_wall_d` の成分別 `*sxx` 射影に特定 (`twall_x≡0` で裏付け)。
- `2026-06-06` — スコープ拡大。内部面 `viscousFlux_d` の法線項も成分 `delta_x`
  (本来スカラー `delta`) を使っており、軸平行 $y$ 法線面で横方向運動量拡散が落ちる
  同種の不具合と判明。純せん断 $u_x=\gamma y$ のテストと熱伝導束 (正しく `delta` 使用) との
  対比で確認。plan を内部面+壁面の両方を対象に改題・改訂し、ファイル名を
  `diffusion-wall-viscous-flux.md` → `diffusion-viscous-shear-flux.md` にリネーム。
- `2026-06-06` — 実装・検証完了。`viscousFlux_d.cu` の内部面・壁面カーネルを完全な
  Newton 応力 (Laplacian 法線スカラー `delta` + 同成分接線 `k` + 転置項フル `S` + 発散項) に
  修正。転置項 $\mu\,\partial u_j/\partial x_i\,S_j$ はレビュー指摘で有効化 (元コードでコメント
  アウトされていたもの)。`run_0003_slau_laminar_viscfix` で SU2 と一致する放物線・no-slip・
  `twall_x` 非ゼロを確認 (§6)。`status` を `done` に更新。
  残り: case 25 / RANS・LES の再ベースラインはフォローアップ。
- `2026-06-06` — 流量の定量検証。長時間 (150k step) で forge 流量は 0.14339 に完全収束し
  完全発達摩擦則を 0.5% で満たすことを確認 (ROE でも 0.14349、対流スキーム非依存)。当初の
  SU2 参照 (0.155) は低マッハ収束不足だったと判明し、`LOW_MACH_PREC=YES` + 長時間で
  再収束させると 0.144 (摩擦則適合) となり **forge と 0.4% 一致**。修正の妥当性を二重に裏付け。
- `2026-06-06` — SU2 `run_0006` 60000iter 完走確認。最終流量 0.14321、Ucenter/Uavg=1.501
  (放物線理論との差 0.1%)。**forge 0.14339 との差は 0.124%** に更に絞られた。
  rms[RhoE]=-1.31 に留まるが低マッハ密度ベースの既知限界で運動量収束には無影響。
- `2026-06-06` — case 25 (軸対称パイプ) 検証完了。`run_0002_slau_laminar_viscfix` (60k step):
  Ux_wall=0.028 m/s (no-slip)、Uc/Ua=2.11、HP 摩擦則 0.18% 適合、mdot=0.000502 kg/s。
  SU2 参照 (`run_0003_su2_laminar`, 60k): mdot=0.000509、Uc/Ua=2.007。**差 1.42%**。
  commit `abd833d` (branch `fix/viscous-shear-flux`) でソース修正を分離コミット済み。
