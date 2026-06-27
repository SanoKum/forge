# 別セッション用プロンプト: backstep の圧力・密度振動の切り分けと打ち手検証

> このファイルは「`case/18.backstep` で node/cell を問わず観測される **圧力 `P`・密度 `ρ` の振動**の
> 原因を切り分け、打ち手を検証する」作業を**別セッションで進めるための指示プロンプト**である。
> 新しいセッションの最初の入力としてこの内容を貼り付けて使う。
> コード改修に入る場合は AGENTS.md の開発フロー (methods→plan→実装) を必ず踏むこと。
> **検証ルール厳守**: forge を回したら `solver_density_cuda/tools/check_convergence.py` と
> `check_quasisteady.py` の VERDICT を必ず引用する。残差プラトー≠収束、過渡ピーク≠定常。
> 各 run は新しい `run_*` で実行し `case/18.backstep/README.md` の run 一覧表を同期する。

---

## 背景 / 問題

forge は密度ベース圧縮性ソルバ (`solver_density_cuda/`)。`case/18.backstep` は低マッハ
(M≈0.11、Pt=100000/Ps=99118、U≈38 m/s) の後退段 RANS-SST を **block-DPLUR 陰解法**
(`timeIntegration:11`, `blockDPLUR:1`, `nStepInner:5`, `cfl_pseudo:1`) で回している。
直近セッションで **node-centered SST の近壁過剰乱流は修正済み** (plan
`turbulence-node-wall-function-coverage`、生産被覆+k Dirichlet)。場平均 vis_turb・x_R は
cell/SU2 と整合するようになった。

**残る問題 = 圧力 `P`・密度 `ρ` の振動が激しい (node・cell の両方)。** discretization に依らないので、
**低マッハの圧力-速度結合 / 数値散逸スケール / リミッタ / リミットサイクル**といったスキーム側が疑わしい。

### 既知の状況 (直近セッションの観測)

- **全 run が残差プラトー (`NOT CONVERGED`)**。cell/node/SU2 いずれも低マッハ剥離で `rms_ro`/`rms[Rho]` が
  -3 付近で下げ止まる (SU2 v8 でも `rms[Rho]≈-2.95`)。`check_quasisteady` では `machmax`/`pmax` は `STEADY`
  (drift~0) だが、残差は落ちない = リミットサイクル的。**「振動」がこの残差プラトーと同根か、別の空間振動かを
  まず切り分ける**こと。
- 代表 run: `case/18.backstep/run_0067_node_wfcov_long` (node)・`run_0057_cell_autowall_imp` (cell)・
  `run_0058_su2_sst` (SU2 基準)。設定は `solver:SLAU`, `convMethod:2 (MUSCL)`, `limiter:2 (venkata)`,
  `cfl_pseudo:1`, `implicitRelax:0.5`。
- メッシュ `mesh/backstep_planar.h5` (10997 node / 平面2D)。`check_mesh_quality.py` は PASS 域。
- 過去 (run_0039 系) で node 低マッハ淀み域に **odd-even チェッカーボード**が出ていた記録あり
  (memory [[node-mode-periodic-and-backstep-status]])。これが cell でも出ているかは要確認。

## まずやること: 振動の「正体」を切り分ける (打ち手の前に必須)

憶測で対策を当てず、**振動が次のどれかをデータで判定**する:

1. **空間 odd-even (チェッカーボード)**: 隣接 CV で `P`/`ρ` が市松に振動。collocated 低マッハの
   圧力-速度デカップリング。→ 隣接セル偏差指標 (`|P_i - 近傍平均|`) を場で可視化・定量。
2. **時間リミットサイクル**: 残差が下がらず一定振幅で振動。プローブ `P(t)` の時系列を取り FFT/RMS。
   `probe.yaml` に段差後流・せん断層・出口手前の点を置き `outStepInterval` 細かめで時系列収集。
3. **リミッタ・チャタリング**: venkatakrishnan リミッタが step ごとに ON/OFF して `P` がパタつく。
   `limiter_P` / `ducros` の出力を時系列で見る。1次化 (`convMethod:0`) で消えるか。
4. **発達途中の物理的非定常**: 剥離せん断層の渦放出 (本質的 URANS)。SU2 でも出るなら物理寄り。
5. **単精度 (float32) の桁落ち**: forge 既定は float32。低マッハでは `P=(γ-1)(roe-ρe_k)` が
   大きな全エネルギー (roe≈2.5e5) と運動エネルギーの**差**で評価され、`P` の有効桁が落ちて
   round-off が振動に化ける (catastrophic cancellation)。`ρ` も同様に状態方程式経由で汚れる。
   → **`build-double` (倍精度ビルド) で同一 run を回し振動が減るか**で精度寄与を判定 (下記 D)。

判定の道具:
- プローブ時系列: `probe.yaml` に点を設定し `P`,`ro`,`Ux` を収集 → 振幅・周波数・空間相関。
- 空間 odd-even 指標: `res_*.h5` の `P`/`ro` で各 CV と双対/セル隣接平均の差を計算 (node は
  `CELLS/centCoords`+`PLANES`、cell は CONNE 隣接)。
- `check_quasisteady.py <run> --quantity pmax,machmax` (派生量の定常性)、`check_convergence.py` (残差)。
- **node/cell/SU2 三者比較** (`procedures/su2-cross-check.md`): 振動が三者共通なら物理/低マッハ一般、
  forge 二者のみなら forge スキーム固有、node のみなら median-dual 固有。
- **単精度 vs 倍精度比較**: 既定 `build-native` (float32) と `build-double` (float64) で同一 run を回し、
  振動振幅 (P プローブ RMS・空間 odd-even 指標) を比較。倍精度で大きく減れば桁落ち寄与が確定。
  SU2 (倍精度) が振動少なめなら傍証になる。

## 打ち手の候補 (切り分け結果に応じて選ぶ)

> **重要**: 1 回に 1 レバーだけ変えて A/B する。既存の正しい既定 (`wallTreatmentSST:1`,
> `nodeWallDirichlet:1`, `nodeKwfDirichlet:1`) は維持したまま比較する。

### A. 低マッハ圧力-速度結合 / 散逸スケール (本命)
- **SLAU2 圧力流束** (`solver:` を SLAU→SLAU2 相当、plan [`convection-slau2-lowmach`](../../plans/accepted/convection-slau2-lowmach.md))。
  SLAU の低マッハ圧力散逸が不足し市松が減衰しない可能性。SLAU2 の圧力項で改善するか。
- **低マッハ前処理** (`lowMachPrecond`, plan [`time_integration-lowmach-preconditioning`](../../plans/accepted/time_integration-lowmach-preconditioning.md))。
  密度ベースは M→0 で音速項が散逸を支配し圧力-速度が緩む。Weiss–Smith 前処理で散逸を対流スケールに
  そろえる。**注意**: 過去に node + `lowMachPrecond=2` は step0 で NaN になった記録あり (memory)。node で
  使うなら前処理の node 対応を先に確認。
- **free-stream 保存流束** (plan [`convection-freestream-preserving-flux`](../../plans/active/convection-freestream-preserving-flux.md))。
  基準静圧差分で低マッハの圧力桁落ちを抑える。振動が圧力桁落ち由来なら効く。

### B. 勾配・リミッタ (空間振動向け)
- **リミッタ**: venkata→barth、またはリミッタ凍結 (収束後固定)。`limiter:1`/freeze で振動が減るか。
- **1次化テスト** (`convMethod:0`): 振動が消えれば 2 次再構成 (MUSCL) が震源。診断用 (恒久解でない)。
- **LSQ 勾配** (node, plan [`discretization-lsq-gradient`](../../plans/active/discretization-lsq-gradient.md))・
  **Ducros センサ**: 既定では SLAU 散逸が支配で効果薄の記録 (memory [[ducros-limiter-inert-with-slau]])。
  ただし振動が勾配起因なら再評価。

### C. 時間積分 / リミットサイクル向け
- `cfl_pseudo` / `nStepInner` / `implicitRelax` の感度 (リミットサイクルが数値時間積分由来か)。
- **dual-time URANS**: 物理的非定常 (せん断層渦放出) なら定常を求めず dual-time で時間平均評価。
  SU2 cross-check も同様に固定 CFL / URANS で挙動を見る (`procedures/su2-cross-check.md` の
  「CFL_ADAPT 共振」節)。
- cell の **atomicAdd 非決定性**は残差を ~6e-4 で底打ちさせる (memory [[cell-atomicadd-nondeterminism]])。
  cell の振動下限がこれで決まっていないか (node はほぼ決定的なので node 側で純スキーム挙動を見る)。

### D. 数値精度 (single/double, float32 の桁落ち)
- **倍精度ビルドで A/B**: `solver_density_cuda/build-double/forge` で同一 run を回し、振動が単精度 (`build-native`)
  比でどれだけ減るか定量。減れば低マッハ `P=(γ-1)(roe-ρe_k)` の桁落ちが (少なくとも一因) と確定。
  - **注意**: 「精度のせいに見えて実は別物」の前例がある — 軸対称近軸 k スパイクは float32 でなく block-DPLUR の
    粘性対角幾何が真因で、倍精度は対症療法だった (memory [[mixed-precision-axisym-refuted]])。**倍精度で消えても
    『精度が真因』と即断せず**、桁落ちが本質か (= スキーム再定式化で単精度のまま直せるか) を見る。
  - 一方 float32 が低マッハ・非直交で free-stream を保持できず発散する実例もある (memory
    [[forge-freestream-nonorthogonal]])。backstep_planar は構造的で非直交は小さいが、低マッハ桁落ちは別問題。
- **桁落ちを断つ再定式化 (単精度を保ちたい場合)**: 圧力を基準値からのゲージ圧 `P' = P - P_ref` で持つ/
  流束を基準静圧差分で組む (plan [`convection-freestream-preserving-flux`](../../plans/active/convection-freestream-preserving-flux.md))、
  `roe` から KE を引く順序・スケーリングの見直し等。倍精度全面化はメモリ・速度コスト大なので、
  **倍精度は「精度寄与の切り分け診断」に使い、恒久解は再定式化 or 低マッハ前処理 (A) を優先**するのが筋。
- 切り分けの組み合わせ: 「単精度 SLAU」「倍精度 SLAU」「単精度 SLAU2/前処理」を比べ、
  *精度*と*スキーム散逸*のどちらが効くか分離する。

## 判定基準 / 成果物

- 振動の正体 (空間 odd-even / 時間リミットサイクル / リミッタ / 物理) を**データで特定**し notes に記録。
- 有効な打ち手を 1 つ以上、A/B で**定量的に実証** (`P` プローブ RMS・空間 odd-even 指標の before/after、
  `check_convergence`/`check_quasisteady` VERDICT 付き)。三者比較で forge 固有か低マッハ一般かを述べる。
- スキーム改修に至る場合は methods→plan→実装の順で。打ち手が config レベルで足りるなら設定推奨を
  `procedures/solver-settings.md` / `divergence-and-startup.md` に反映。
- `case/18.backstep/README.md` の run 一覧を同期。結論と根拠 run パスを明示。

## 参照
- 関連 plan: `convection-slau2-lowmach`, `time_integration-lowmach-preconditioning`,
  `convection-freestream-preserving-flux`, `discretization-lsq-gradient`, `discretization-median-dual`。
- methods: `methods/convection/`, `methods/time_integration/`, `methods/discretization.md`。
- 手順: `procedures/su2-cross-check.md`, `procedures/solver-settings.md`,
  `procedures/divergence-and-startup.md`。
- memory: node-mode-periodic-and-backstep-status, ducros-limiter-inert-with-slau,
  cell-atomicadd-nondeterminism, forge-freestream-nonorthogonal,
  mixed-precision-axisym-refuted (精度のせいに見えて幾何だった前例), native-build-now-possible。
- ビルド: `build-native` (float32, 既定)・`build-double` (float64, 精度切り分け用)。
  精度評価・ベンチは native を既定 (procedures/development-environment.md)。
