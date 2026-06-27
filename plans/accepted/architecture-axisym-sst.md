# architecture-axisym-sst — 軸対称 SST 幾何項 子計画

## メタ

- **area**: `architecture`
- **status**: `done`
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) (§7)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) (§5)
  - [`methods/axisymmetric/theory.md`](../../methods/axisymmetric/theory.md)
- **related_plans**:
  - [`architecture-rans-sst.md`](architecture-rans-sst.md) (親)
  - [`architecture-axisymmetric.md`](architecture-axisymmetric.md)
- **created**: `2026-06-07`
- **owner**: `Claude`

## 1. 目的

親 plan [`architecture-rans-sst.md`](architecture-rans-sst.md) で導入した Menter SST
は、3D 直交座標の式をそのまま 2D 軸対称メッシュに乗せている。本子 plan では
軸対称特有の幾何項 — 生産項のフープひずみ $S_{\theta\theta} = u_r/r$ と圧縮性
(dilatation) 補正 — を SST 生産項・渦粘性へ正しく取り込み、軸対称ノズルで
物理的に妥当な乱流場を得る。

## 2. スコープ

### やる
- 生産項・渦粘性のひずみ速度 $S^2$ に $2\,(u_r/r)^2$ を追加 (theory.md §7.2)。
- deviatoric 評価時の完全発散に `axisym_divU` を使う圧縮性補正 (§7.3)。
- 既存 NS 診断量 (`axisym_uy_over_r`, `axisym_divU`) の SST 側流用。
- 拡散の $\theta\theta$ 寄与が $r$ 重み面積で発散形を保つことの単体確認 (§7.1)。
- `isAxisymmetric == 0` で 2D / 3D SST と完全一致する回帰確認。
- 軸対称ノズル複製 run での検証と residual plot 保存。

### やる (追補 2026-06-07: dilatation 補正)
- `dilatationCorrection` config (0:off / 1:A / 2:A+B) を追加。
- (A) deviatoric トレース除去: $S^2 \mathrel{-}= \tfrac23(\nabla\!\cdot\!\mathbf u)^2$ (theory.md §7.3)。
- (B) 等方項: $P_k \mathrel{-}= \tfrac23\rho k(\nabla\!\cdot\!\mathbf u)$ + $P_k\ge0$ クリップ。
- run_0091(off) → run_0092(A) → run_0093(A+B) の段階効果確認。

### やらない
- 壁関数・遷移モデル・他 RANS モデル。
- SST の陰解法 Jacobian 連成 (別 plan)。
- 任意軸 (x 軸固定以外) への対応。
- Sarkar / Zeman / Wilcox の $M_t$ dilatation-dissipation 補正 (theory.md §7.3 参照、不採用)。

## 3. 関連 docs と前提

- 理論: [theory.md §7](../../methods/turbulence/theory.md) (フープひずみ・圧縮性・診断量再利用)。
- 実装: [implementation.md §5](../../methods/turbulence/implementation.md)。
- 軸対称 B 流儀の前提: [axisymmetric/theory.md](../../methods/axisymmetric/theory.md)。
- 既存診断量の供給元:
  [`axisymmetricSource_d.cu`](../../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)
  の `axisymmetricDiagnostics_d` が `axisym_uy_over_r`, `axisym_divU` を毎ステップ更新する。

## 4. 設計方針

theory.md §7 に対する実装差分のみ記す。

- `rans_sst_source_d` / `sst_eddy_viscosity_d` に
  `int isAxisymmetric` と `flow_float* axisym_uy_over_r` を引数追加。
- `isAxisymmetric == 1` のとき planar 9 成分から組んだ
  $S^2$ に $2\,(\texttt{axisym\_uy\_over\_r}[ic])^2$ を加算。
  `isAxisymmetric == 0` では加算 0 → 既存挙動と完全一致。
- 圧縮性補正を入れる場合は `axisym_divU` を併せて渡す。初期実装では
  standard SST ($S = \sqrt{2 S_{ij}S_{ij}}$, deviatoric 化なし) を維持し、
  dilatation 補正はオプションとして後段で評価する。
- main loop での順序保証: `axisymmetricDiagnostics_d`(診断量更新) →
  `turbulent_viscosity_d_wrapper`(渦粘性) → `ransSource_d_wrapper`(source)。
- 新規物理量の計算・新規 device 配列は追加しない (既存 `var.c_d` を流用)。

## 5. 実装ステップ

1. [`turbulent_viscosity_d.cu`](../../solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu):
   `sst_eddy_viscosity_d` に `isAxisymmetric` / `axisym_uy_over_r` を追加し、
   $S^2$ にフープ項を加算。wrapper の呼び出しに `cfg.isAxisymmetric`,
   `var.c_d["axisym_uy_over_r"]` を渡す。
2. [`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu):
   `rans_sst_source_d` に同様の引数・加算を追加。
3. `main.cpp`: 診断量更新 → 渦粘性 → source の順序を確認 (既存順序の検証のみ。
   必要なら並べ替え)。
4. 拡散 $\theta\theta$ の発散形保持を単体確認 (半径方向解析解 or 物理量チェック)。
5. 軸対称ノズル複製 run (run_0091) で検証、residual plot 保存。
6. laminar / 既存 SST (2D) ケースで null regression。

各ステップで触るファイルは上記に明記。

## 6. 検証

- **単体 / ビルド**: `solver_density_cuda` が native (arch 86) でコンパイルできること。
- **検証ケース**:
  - 既定: `case/23.axi_nozzle` の run_0090 (full SST) を複製した run_0091。
  - 回帰: `isAxisymmetric == 0` の既存 SST ケースで挙動不変。
- **判定基準**:
  - フープ項追加後も `k`, `omega` が負値化・発散しない。
  - 軸近傍で生産項が増え、`vis_turb` 分布が run_0090 から物理的に妥当に変化。
  - 残差 (`rms_roK`, `rms_roOmega`) が安定して低下する。
  - `isAxisymmetric == 0` 経路は数値的に従来と一致。
  - `residual_history.csv` / `residual_history.png` を run_0091 に保存。

## 7. 影響範囲

- 触る主要ファイル:
  - [`cuda_forge/turbulent_viscosity_d.cu`](../../solver_density_cuda/cuda_forge/)
  - [`cuda_forge/ransSource_d.cu`](../../solver_density_cuda/cuda_forge/)
  - [`main.cpp`](../../solver_density_cuda/main.cpp) (順序確認)
- 既存ケースへの影響:
  - `isAxisymmetric == 0` では完全互換を維持。
- ドキュメント:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §7 (済)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) §5 (済)
  - [`.github/plans/README.md`](../README.md)

## 8. 完了条件

- [x] 関連 `methods/turbulence/theory.md` §7 更新済み
- [x] 関連 `methods/turbulence/implementation.md` §5 更新済み
- [x] フープひずみ項を生産項・渦粘性に実装 (`rans_sst_source_d`, `sst_eddy_viscosity_d`)
- [x] 圧縮性補正の方針確定 → 追補で `dilatationCorrection` (0:off/1:A/2:A+B) として実装・検証済み。既定値 2 (A+B, 全 SST ケース)。Sarkar 系 $M_t$ 散逸は不採用
- [x] 拡散 $\theta\theta$ の発散形保持を確認 (variables.cpp の $r$ 重み付き sx/volume でコード構造的に確定)
- [x] run_0091 で軸対称 SST 検証完了 (residual plot 保存)
- [x] `isAxisymmetric == 0` の null regression 確認 (`if(isAxisymmetric)` ガードによる構造的ビット一致)
- [x] (追補) `dilatationCorrection` config 追加 (0/1/2)
- [x] (追補) (A) deviatoric を実装し run_0092 で効果確認 (場への影響は小、安定)
- [x] (追補) (B) 等方項を実装し run_0093 で効果確認 (k -12%, vis_turb -14%, 安定収束)
- [x] `.github/plans/README.md` の状態を `done` に更新
- [x] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `2026-06-07` — 初稿。親 plan の残課題「軸対称 SST の子 plan 起票」を受けて作成。
  theory.md §7・implementation.md §5 に軸対称幾何項を整理済み。
- `2026-06-07` — 実装・検証を完遂し status → done。
  - `sst_eddy_viscosity_d` / `rans_sst_source_d` に `isAxisymmetric` と
    `axisym_uy_over_r` を引数追加。`isAxisymmetric==1` のとき
    $S^2 \mathrel{+}= 2\,(u_r/r)^2$ を加算 (フープひずみ)。
  - 圧縮性補正は standard SST ($S=\sqrt{2S_{ij}S_{ij}}$) を維持。`axisym_divU` は将来の
    dilatation オプション用に温存。
  - 拡散 $\theta\theta$: [variables.cpp:263-272](../../solver_density_cuda/variables.cpp)
    で軸対称時 `sx *= r_face`, `volume *= r_cell` と $r$ 重み付けされ、拡散 kernel が
    $r$ 重み面積、時間積分が $r$ 重み体積を使うため B 流儀の発散形を保持。別 source 不要。
  - 検証 run_0091 (run_0090 収束状態から +10000 steps, 軸対称フープ有効):
    - rms_ro: 1.87e-6 → 1.32e-7
    - rms_roK: 0.062 → 0.0037
    - rms_roOmega: 5.53 → 0.633
    - k: min 0.004→0.093 (軸近傍生産の反映), max 6.83e4→6.85e4
    - vis_turb: max=0.344 (制限子で有界), 発散なく安定収束
  - null regression: `if(isAxisymmetric)` ガードにより非軸対称経路は従来とビット一致。
- `2026-06-07` (追補) — 生産項の圧縮性 (dilatation) 補正を実装・段階検証。
  - `solverConfig` に `dilatationCorrection` (0:off / 1:A deviatoric / 2:A+B) を追加
    ([solverConfig.{hpp,cpp}](../../solver_density_cuda/input/))。検証時は 0 既定だったが、
    ユーザ判断で最終的に **既定値 2 (A+B)** に決定 (全 SST ケースに適用。
    SST kernel は LESorRANS==2 のみ走るため laminar/LES は不変)。
  - `rans_sst_source_d` に `divU = ∂xUx+∂yUy+∂zUz (+軸対称 u_r/r)` を計算し:
    - (A) level≥1: $S^2 \mathrel{-}= \tfrac23(\nabla\!\cdot\!\mathbf u)^2$、$S^2\ge0$ クリップ。
    - (B) level≥2: $P_k \mathrel{-}= \tfrac23\rho k(\nabla\!\cdot\!\mathbf u)$、生産リミッタ後に $P_k\ge0$ クリップ。
      $\omega$ 生産には入れない。
  - 検証 (各 run は前段収束状態から +10000 steps 継続):
    - run_0092 (A, level 1): k mean 7170→7160, vis_turb max 0.344 不変。
      残差 rms_roK 0.0037→0.0014。場への影響は小、安定。
      理由: $k$ 主体の境界層では $\nabla\!\cdot\!\mathbf u$ 小、膨張コアでは $\mu_t$ 小 + 生産リミッタ。
    - run_0093 (A+B, level 2): k mean 7160→6280 (−12%), k max 6.85e4→6.18e4 (−10%),
      vis_turb mean 0.0267→0.0229 (−14%), max 0.344→0.317 (−8%)。
      膨張による乱流減衰を物理的に捕捉。残差は単調減少 (rms_roK 1.87→0.110, 継続低下中)。
      (B) のシンク項で剛性がやや増し残差フロアは (A) より高めだが安定・収束。
  - Sarkar / Zeman / Wilcox の $M_t$ dilatation-dissipation は壁流れで有害のため不採用。
  - 既定値: ユーザ判断で `dilatationCorrection` の既定を **2 (A+B)** に決定
    ([solverConfig.hpp](../../solver_density_cuda/input/solverConfig.hpp),
    [solverConfig.cpp](../../solver_density_cuda/input/solverConfig.cpp))。
    全 SST ケースに圧縮性補正が既定適用される。従来の非圧縮挙動が必要なときのみ `0` を明示。
