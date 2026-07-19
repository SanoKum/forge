# WMLES 用代数壁応力モデル (Reichardt + Kader)

## メタ

- **area**: `turbulence / boundary`
- **status**: `in_progress`
- **related_docs**:
  - `methods/turbulence/theory.md` / `methods/turbulence/implementation.md` (壁法則・SST 壁関数の現在仕様)
  - `methods/boundary.md` (壁 BC の現在仕様)
- **related_plans**:
  - [`turbulence-iddes-sst.md`](turbulence-iddes-sst.md) — DES/IDDES。本計画は「純 LES (`LESorRANS==1`) + 壁応力モデル」で、IDDES の WMLES モード (RANS 層が壁モデルを兼ねる) とは別経路。共存可
  - [`turbulence-node-sst-wallfunction.md`](../accepted/turbulence-node-sst-wallfunction.md) / [`turbulence-node-wall-function-coverage.md`](../accepted/turbulence-node-wall-function-coverage.md) — node 壁関数の既存構造 (Normal_Neighbor / `Tau_Wall` 再スケール)。本計画はこの機構を流用する
  - [`turbulence-enhanced-wall-treatment.md`](../accepted/turbulence-enhanced-wall-treatment.md) — Reichardt 壁法則の初出
  - [`convection-keep-diss-lowmach-precond.md`](convection-keep-diss-lowmach-precond.md) — LES 対流基盤 (KEEP + matrix ES 散逸 + 低マッハ前処理)。WMLES 実行の前提スキーム
  - [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) — node 壁境界のゴースト撤廃 (壁から着手)。壁面状態の取得層と干渉するため実装順を調整する
- **related_notes**: [`turbulence-des-wmles-survey.md`](../../notes/investigations/turbulence-des-wmles-survey.md) — 平衡壁モデルの限界 (高 Re 剥離予測失敗) を含むサーベイ
- **created**: `2026-07-19`
- **owner**: `CFD Dev`

---

## 1. 目的

forge に WMLES (壁モデル LES) 用の代数壁応力モデルを実装する。壁境界面ごとに、マッチング点 (壁から 1 つ内側の解点) の瞬時 LES 状態と壁面状態から、壁面せん断応力 $\tau_w$ と壁面熱流束 $q_w$ を代数壁法則 (Reichardt + Kader) で算出し、壁の粘性流束をモデル値で置き換える。対象は圧縮性 (〜Mach 2 級)、断熱壁・等温壁の両方、cell / node (median-dual) 両対応。node を一次対象とする (プロジェクト方針)。

## 2. スコープ

- **やる**:
  - 壁法則デバイス関数群の共通化 (`wallLaw_d.cuh` 新設): 既存 SST 壁関数の Reichardt $u^+(y^+)$ / $du^+/dy^+$ を昇格し、Kader $T^+$ を追加
  - WMLES 用 $u_\tau$ ソルバ (warm start + 収束判定付き Newton。既存 SST 壁関数の固定 5 回反復とは別実装、SST 経路はビット不変)
  - $\tau_w$ の流束適用: cell = 壁面粘性流束の置換 (既存 `wallTreatment==1` 分岐の流用)、node = `Tau_Wall` による W↔I 内部双対面の traction 再スケール (既存 AddTauWall 機構の流用)
  - $q_w$ の流束適用: cell = 壁面 heatflux 項の置換、node = AddTauWall と対称な W↔I 熱流束再スケール (AddQWall・§4.6)
  - 有効化 config (`wallModelLES`) と bcondConfig の壁単位指定
  - 単体検証 (`tools/verify_wall_law.py` + `tests/unit/`) と チャネル流検証 (Reτ≈2000 目標、まず Reτ≈550 で較正)
  - チャネル駆動用の一様体積力 (streamwise body force) config — 前提の小実装 (§6.2)
- **やらない (将来拡張 / 別 plan)**:
  - ODE 型 (1D 埋め込み格子)・TBLE 型非平衡モデル (インターフェースは差し替え前提で設計)
  - 乱流流入生成 (SEM / recycling)。**超音速平板 M2 検証はこれが前提のためゲート** (IDDES plan の「やらない」と同一依存)
  - 種輸送・化学種拡散の壁モデル整合 (現状 node 半割面 skip のまま)
  - RANS SST 壁関数 (`wallTreatmentSST`) の挙動変更 — 共通化はデバイス関数の昇格のみ、数値経路はビット不変

## 3. 現状コードの実態 (2026-07-19 調査。本計画はこの資産の上に組む)

初稿の仕様書は取得層・壁法則・流束置換をすべて新規実装と想定していたが、**実態は主要部品の大半が SST automatic wall treatment として実装済み**である。差分だけ作る。

| 部品 | 既存実装 | 状態 |
| --- | --- | --- |
| Reichardt $u^+(y^+)$, $du^+/dy^+$ | `cuda_forge/ransWallFunction_d.cu:18-39` (`reichardt_uplus` / `reichardt_duplus_dyp`, κ=0.41) | **あり**。ただし無名 namespace でファイルローカル → ヘッダ昇格が必要 |
| $u_\tau$ Newton 反復 | 同 `:154-164`。残差形 $f=U_t/u_\tau-u^+$、初期値は粘性則、固定 5 回・収束判定なし・**warm start なし** | あり (要改良: §4.4) |
| 壁平行速度射影・よどみガード | 同 `:129-152` | あり |
| node マッチング点選択 (壁面→第一内点ノード) | 同 `:94-124`。SU2 Normal_Neighbor 方式: 壁ノードの双対面 CSR を走査し内向き法線 cos 最大の内点 `irep` と法線距離 `y=(x_I-x_W)\cdot(-n)` を取る | **あり** (仕様書の「取得層」がほぼこれ) |
| cell マッチング点 | `irep=ic`, `y=wall_dist[ic]` (同 `:92-93`) | あり |
| 壁面ごとの永続配列 | `bc.bvar_d["utau"/"ypls"/"twall_*"]` (`mesh.hpp:65`, 確保 `mesh.cpp:44-79`)。全 step 保持 | **あり** (warm start の置き場も既にある) |
| cell 壁面粘性流束の置換 | `cuda_forge/viscousFlux_d.cu:345-365`: `wallTreatment==1` で接線せん断を $\rho u_\tau^2$ (向き $-\hat e_t$) に置換、法線粘性・体積項は落とす | **あり** (WMLES はこの分岐にゲートを追加するだけ) |
| node 壁応力の伝達 | `viscousFlux_d.cu:148-165` AddTauWall: 片端のみ壁ノードの W↔I 内部双対面で解像 traction の接線成分を `Tau_Wall[icW]`$=\rho u_\tau^2$ に再スケール。**再スケール後の τ が粘性仕事項 (`:182`) にも入るためエネルギー整合** | **あり** |
| 壁 BC 種別 | `wall` (断熱: `T[ig]=T[ic]` ミラー) / `wall_isothermal` (`Ts` を YAML uniform で読み `T[ig]=Tsb`) — `boundaryCond.hpp:30-247`, `boundaryCond_d.cu` | あり |
| $q_w$ の流束フック | **なし**。壁熱流束は ghost 温度差 (cell) / 弱形式 $\nabla T\cdot S$ (node 半割面) で暗黙に決まるのみ | **新規** |
| Kader $T^+$ | なし | **新規** |
| thermo | `thermo_d.cuh`: `thermo_cp_mix`/`thermo_R_mix`/`thermo_mu_mix`/`thermo_lambda_mix` 等 (内部 double)。CPG は `cfg.cp`/`cfg.gamma` 定数。μ(T) は `viscMethod` 0=定数/1=Sutherland/2=分子論 (`gasProperties_d.cu:50-79`) | あり (可変比熱対応可) |
| 層流 Pr | **config キーは存在しない**。λ は μ と独立に決まる (定数 or Eucken) → Pr は $\mu(T_w)c_p(T_w)/\lambda(T_w)$ をその場評価する | 注意 |
| 精度 | `flow_float` = float (FP32, `flowFormat.hpp:6`)。thermo 内部は double | 注意 |

**重要な構造的事実 (仕様書の想定と異なる点)**:

1. **node モードでは壁ノードの運動量残差は Dirichlet ゼロ化される** (`nodeWallDirichlet=1` 既定, `zeroWallMomentumResidual`)。したがって「壁境界半割面の粘性流束を置き換える」だけでは τ_w は内点の運動量方程式に届かない。壁応力は **W↔I 内部双対面の AddTauWall 再スケールで効かせるのが forge の検証済み構造**であり、本計画もこれに従う (§4.5)。
2. 既存 wall function ゲートは `LESorRANS==2 && RANSmodel==1 && wallTreatmentSST` に閉じているため、LES (`LESorRANS==1`) では壁関数経路は一切走らない。**コードパス分離は自然に成立**しており、新ゲートを足すだけで SST と干渉しない。
3. RANS/SST の k/ω 境界処理 (`ransBoundary`) は `LESorRANS==2` ゲート内 → WMLES 実行時に呼ばれない (仕様書の前提どおり)。
4. node の k/ω・化学種拡散は境界半割面 skip が正 ([diffusion-node-boundary-halfface-skip])。WMLES はこれに触れない。

## 4. 設計方針

### 4.1 インターフェース (仕様書 §1 を維持、取得層は既存機構に対応付け)

壁モデル本体は離散化を知らない `__device__` 純関数とする。

```
入力:  マッチング点: ρ_m, u_i,m, T_m, μ_m
       壁面状態:    T_w_in, p_w_in (等温なら T_wall 指定値を優先)
       幾何:        単位法線 n_i (壁→流体), 壁面垂直距離 d
       warm start:  前 step の u_τ
出力:  τ_w (大きさ), e_∥,i (向き), q_w (壁→流体 正), u_τ (保存用)
```

取得層 (離散化ごと):

- **node**: マッチング点 = Normal_Neighbor 内点 `irep`、$d$ = 内向き法線距離 `bestDn` (既存 `ransWallFunction_d.cu:94-124` を共通関数に切り出して流用)。壁面状態 = 壁ノードの解 ($T_w$ = 壁ノード T、$p_w$ = 壁ノード Ps)。
- **cell**: マッチング点 = 壁第 1 セル中心 (`ic`)、$d$ = `wall_dist[ic]`。壁面状態 = ghost との面平均 $T_w=\tfrac12(T_g+T_{ic})$ (断熱ミラーなら $=T_{ic}$)、$p_w$ 同様。

ODE/TBLE 型への将来差し替えは同一シグネチャの別実装として追加する。

### 4.2 定式化 (仕様書 §2 のまま。変更点のみ注記)

壁平行速度: $u_{\parallel,i}=u_{i,m}-(u_{j,m}n_j)n_i$, $\hat e_{\parallel,i}=u_{\parallel,i}/u_\parallel$。

壁面物性: $\rho_w=p_w/(R_w T_w)$ ($R_w$: TP は `thermo_R_mix`(壁組成)、CPG は $(\gamma-1)c_p/\gamma$)、$\mu_w=\mu(T_w)$、$\lambda_w=\lambda(T_w)$、$c_{p,w}=c_p(T_w)$。μ/λ は `viscMethod` に整合するデバイス helper を新設 (0: 定数 / 1: Sutherland 式 / 2: `thermo_mu_mix`・`thermo_lambda_mix`) — `gasProperties_d.cu` のセル評価と同式。$\mathrm{Pr}=\mu_w c_{p,w}/\lambda_w$ をその場評価 (層流 Pr の config キーは無い・作らない)。γ 一定・cp 一定は式に焼き込まない。

Reichardt 則 (既存実装と同一定数):

$$u^+ = \frac{1}{\kappa}\ln(1+\kappa y^+) + 7.8\left(1 - e^{-y^+/11} - \frac{y^+}{11}e^{-y^+/3}\right),\quad \kappa=0.41$$

Kader 温度壁法則 (新規、等温壁のみ):

$$T^+ = \mathrm{Pr}\,y^+ e^{-\Gamma} + \left[\mathrm{Pr}_t\,(u^+ + P(\mathrm{Pr}))\right]e^{-1/\Gamma},\quad
\Gamma=\frac{0.01(\mathrm{Pr}\,y^+)^4}{1+5\,\mathrm{Pr}^3 y^+}$$
$$P(\mathrm{Pr}) = (3.85\,\mathrm{Pr}^{1/3}-1.3)^2 + 2.12\ln\mathrm{Pr},\quad \mathrm{Pr}_t=0.9$$

駆動温度差は回復温度 $T_r = T_m + r\,u_\parallel^2/(2 c_{p,w})$, $r=\mathrm{Pr}^{1/3}$ を用い

$$q_w = \frac{\rho_w\,c_{p,w}\,u_\tau\,(T_r - T_w)}{T^+}$$

断熱壁は $q_w=0$ で閉じる (壁面温度は解に含まれる摩擦加熱込みの値を物性評価に使う。回復温度による壁温推定はしない)。

### 4.3 Pr_t の整合 (仕様書からの追加注記)

内部場の乱流熱流束は `turbulentPrandtl` (既定 0.85, `viscousFlux_d.cu:174`) を使うが、これは SGS 渦粘性に掛かる係数で WMLES では近壁を支配しない。Kader の $\mathrm{Pr}_t=0.9$ は壁法則側の定数として独立に持つ (config で変更可、§4.8)。

### 4.4 $u_\tau$ Newton (仕様書 §2.4 + 残差形の変更)

WMLES ソルバは残差を **$F(u_\tau)=u_\tau f(y^+)-u_\parallel$** で組む ($u_\tau\to0$ で正則。既存 SST 形 $U_t/u_\tau-u^+$ は 0 近傍特異で warm start 失敗時に脆い)。$F'=f+y^+f'$。

- 初期値: `bvar_d["utau"]` の前 step 値。初回 / 無効値 (≤0, NaN) は層流推定 $u_\tau^0=\sqrt{\mu_w u_\parallel/(\rho_w d)}$
- 収束判定 $|\Delta u_\tau|/\max(u_\tau,\varepsilon)<10^{-6}$、最大 20 回。反復中 $u_\tau\le0$ は層流推定に置き直し
- 不収束時: 層流応力 $\tau_w=\mu_w u_\parallel/d$ にフォールバックし、デバッグカウンタを atomicAdd (黙って握りつぶさない)
- $u_\parallel<\varepsilon_u$: Newton スキップ、$\tau_{w,i}=\mu_w u_{\parallel,i}/d$ (既存よどみガードと同型)
- `bvar_d["utau"]` は restart ファイルに含めない (初回数 step で回復するため。既存 bvar と同運用)

SST 壁関数の Newton (固定 5 回) は**変更しない**。共通化するのは $u^+/du^+$ 関数のみ。

### 4.5 $\tau_w$ の流束適用

- **cell**: `viscousFlux_wall_d` の既存 `wallTreatment==1` 分岐 (`viscousFlux_d.cu:345-365`) と同じ位置に WMLES 分岐を追加: 接線せん断を $\tau_w=\rho_w u_\tau^2$、向き $-\hat e_{\parallel}$ (マッチング点瞬時速度に沿う。逆流時も瞬時向きに従う) で置換。法線粘性・体積項は既存分岐と同様に落とす。壁面 no-slip なので粘性仕事は 0。
- **node**: 壁ノードに `Tau_Wall[icW]=\rho_w u_\tau^2` を書き、既存 AddTauWall 再スケール (`viscousFlux_d.cu:148-165`) をそのまま使う。W↔I 面ごとの解像 traction の**向きを保って大きさだけ** τ_w に合わせる方式で、仕様書の「$\hat e_\parallel$ に沿わせる」とは厳密には異なるが、面ごとの瞬時解像方向を使うためむしろ局所的で、SST node 壁関数で検証済みの構造 (逆流も面の解像方向が自動で追随)。粘性仕事項への伝播も既存コードで整合済み (`:182` が再スケール後の τ を使う)。壁ノード運動量は従来どおり Dirichlet ゼロ化 (`nodeWallDirichlet=1` 必須)。境界半割面の弱形式運動量流束は W ノード残差がゼロ化されるため触らない。

### 4.6 $q_w$ の流束適用 (新規フック)

- **cell**: `viscousFlux_wall_d` の heatflux 項 (`viscousFlux_d.cu:371-384`) を WMLES 分岐で $q_w\,S$ に置換 (断熱は 0)。ghost 温度差による解像伝導は使わない (粗格子で過小評価するため)。
- **node**: AddTauWall と対称の **AddQWall** を新設: 片端のみ壁ノードの W↔I 面で、解像伝導熱流束 (`heatflux`) を壁ノード側の $q_w$ 値 (per-cell 配列 `Qw_Wall`) に再スケール (または置換)。境界半割面の弱形式熱流束は壁ノード CV にしか入らないため、内点へ熱を正しく届けるには W↔I 面で効かせる必要がある (運動量と同じ理屈)。
- **等温壁 node の現状確認 (実装ステップに含める)**: node モードで `wall_isothermal` の壁ノード温度が現在どう強制されているか (Dirichlet ピンか、弱形式のみか) を確認し、WMLES 時は壁ノード T を `Ts` に Dirichlet ピン + エネルギー残差ゼロ化 (運動量の `nodeWallDirichlet` と同型) に揃える。[`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) が壁境界を先に触る予定のため、着手前に同 plan と順序を調整する。
- 断熱壁: cell はミラーで解像熱流束が元々 0、node は $q_w=0$ を書くだけ。両モードとも Phase 1 でここまで、等温壁 (AddQWall + Kader) は Phase 2 とする (チャネル検証で等温壁が要るため Phase 2 までを検証前に完了)。

### 4.7 対流・SGS・勾配との関係 (仕様書 §3 のまま)

- 対流流束: 従来通り (壁面貫通質量流束ゼロ・圧力項不変)。KEEP/ES 散逸レイヤの壁面処理には触れない
- SGS 渦粘性: 壁面流束の評価には使わない (壁モデルが全応力を与える)。内部はそのまま。node W↔I 面は解像 traction (SGS 込み) を再スケールするので自動整合
- 勾配計算: 既存 no-slip 壁と同一
- `ypls_b` は WMLES 側でも書く (既存の「壁関数時は更新しない」ゲート `viscousFlux_d.cu:404` に WMLES 分岐を追加)

### 4.8 有効化と config

- **壁単位の有効化**: 新 BC 種は作らず、既存 `wall` / `wall_isothermal` の bcond `ints:` に `wallModelLES: 0/1` (type 3) を追加。断熱/等温は従来どおり BC 種で決まる (仕様書の「同じ BC 名で動作」を満たす)
- **ゲート**: `viscousFlux_wall_d` へ渡す `wallTreatment` を 3 値化 (0: なし / 1: SST 壁関数 / 2: WMLES)。2 は `LESorRANS==1` (または 0 = ILES) かつ当該 bcond の `wallModelLES==1` のとき。SST 経路 (`LESorRANS==2`) は従来ロジックのまま**ビット不変**
- **solverConfig (`turbulence` セクション、`getOptionalValidatedValue` で既定値付き)**:

| キー | 既定値 | 内容 |
| --- | --- | --- |
| `wmlesNewtonTol` | 1e-6 | Newton 相対許容 |
| `wmlesNewtonMaxIt` | 20 | Newton 最大反復 |
| `wmlesPrt` | 0.9 | Kader Pr_t |

  κ, C1–C3 (Reichardt) は既存 SST 壁関数と共有の compile-time 定数のまま (config 化しない。較正対象ではなく、SST 側との二重管理を避ける)。マッチング点は「壁に最も近い内側の解点」固定 (パラメータ化しない)
- **stale build 注意**: `solverConfig.hpp` 変更後は full rebuild 必須 (既知の罠 [stale-build-struct-layout-trap])

### 4.9 実装配置

- `cuda_forge/wallLaw_d.cuh` (新設): `reichardt_uplus` / `reichardt_duplus_dyp` を `ransWallFunction_d.cu` から昇格 (`__device__ __forceinline__`、κ 等の定数も移動)、Kader $T^+$/$P(\mathrm{Pr})$/$\Gamma$、WMLES Newton ソルバ、μ/λ(T) helper を追加。`ransWallFunction_d.cu` は昇格した関数を include して使用 (数値経路ビット不変をビルド後に回帰で確認)
- `cuda_forge/wmlesWallModel_d.cu` (新設): 壁境界面並列カーネル。取得層 (cell/node) → 壁モデル純関数 → `bvar_d["utau"/"ypls"/"twall_*"]` と node 用 `Tau_Wall`/`Qw_Wall` 書き込み。呼び出しは `applyBconds` 後・`viscousFlux` 前 (SST 壁関数と同じ位相)
- `viscousFlux_d.cu`: cell 壁面分岐の 3 値ゲート化、AddQWall 追加
- `variables.hpp`: `Qw_Wall` を cell 変数リストに追加 (`Tau_Wall` と同様、毎 step −1 初期化)

## 5. 実装ステップ

1. **現状確認**: node `wall_isothermal` の壁ノード温度強制の実態調査 (§4.6)。ghostless plan との順序調整
2. **methods 更新**: `methods/turbulence/` に壁モデル節 (理論: Reichardt/Kader/回復温度、実装: 取得層と AddTauWall/AddQWall 経路)、`methods/boundary.md` に `wallModelLES` 指定方法を追記
3. **`wallLaw_d.cuh` 新設**: 既存関数昇格 + Kader + Newton + 物性 helper。full rebuild + 回帰 (`tests/regression --smoke` と SST ケースのビット不変確認)
4. **`tools/verify_wall_law.py`**: §6.1 の単体検証 (python 参照実装との突合)
5. **WMLES カーネル (断熱・Phase 1)**: 取得層 + Newton + `Tau_Wall` (node) / cell 置換分岐 + config/bcond 配管。cell/node 両方
6. **等温壁 (Phase 2)**: Kader $q_w$、cell heatflux 置換、node AddQWall + 壁ノード T Dirichlet
7. **体積力 config (チャネル前提)**: `bodyForce` (運動量ソース $f_x$ + エネルギー $f_x u_x$) を追加 (小実装、別 plan にしない)
8. **チャネル検証** (§6.2) → 結果に応じ SGS/σ 較正
9. plan を `done` 化し accepted へ移動、README/methods/index 同期

## 6. 検証

### 6.1 単体 (`tools/verify_wall_law.py` + `tests/unit/`。CI は未整備のためローカルハーネス運用)

仕様書 §5 を踏襲:

- 既知解回復: $u^+=f(y^+)$ を満たす組を $y^+=0.5,5,30,100,300,1000$ で生成し Newton が $u_\tau$ を相対 1e-5 で回復
- $f, f'$ の正値・連続・単調性 ($y^+\in[0.1,2000]$ 対数スキャン)、$y^+\to0$ で $f\to y^+$ (<1%@0.1)、$y^+=1000$ で log 則との差 <2%
- 層流極限 (NaN/Inf なし)・逆流反転・warm start 3 反復以内収束
- 断熱 $q_w=0$ 厳密、等温 Kader $T^+$ の文献値整合 (Pr=0.71, $y^+=100$)、温度依存 cp で NaN/Inf なし
- 取得層の離散化間一致: 線形速度分布 + 一様温度の解析場で cell/node 両取得層の (T_w, p_w, u_∥, d) と $\tau_w$ が離散化誤差内で一致
- 追加 (実態due): `wallLaw_d.cuh` 昇格後に SST 壁関数経路がビット不変であること (回帰ハーネス)

### 6.2 検証ケース

1. **チャネル流 (本命 Reτ≈2000, Hoyas & Jiménez; 較正は Reτ≈550 の小規模から)**
   - 前提: streamwise 体積力 (§5-7)、周期境界 (node Cartesian 周期・3D median-dual 周期は実装済み)、`unsteady:1`、KEEP + matrix ES 散逸 (`keepDissType:2, keepDissJump:2, keepDissPrecond:1`) + 低マッハ (M~0.2-0.3 で運転し圧縮性効果を抑える)。壁は等温 (散逸熱の逃げ場。断熱では bulk 温度が漸増する) — **したがって Phase 2 完了が検証前提**
   - 格子: 第 1 内点 $y^+\approx50\text{–}100$、壁法線緩等比。`check_mesh_quality.py` VERDICT 併記
   - 判定: 平均速度 log 層 $u^+$ 誤差 <5%、$\langle\tau_w\rangle$ が体積力と釣合うこと (大域運動量収支)、$u_\tau$ 統計の定常性 (`check_quasisteady.py` 系の統計確認)。SGS モデル (WALE/sigma) と ES σ の組はここで確定
   - 判定は瞬時場一致ではなく統計量の許容帯 (atomicAdd 非決定性のため。node は決定性が高いが方針は同じ)
2. **超音速断熱平板 M≈2**: van Driest 変換速度分布と Cf 相関比較。**乱流流入生成 plan の完了にゲート** (それまで実施しない)
3. **本命形状 (ノズル / ピントル)**: チャネル合格後

### 6.3 既知リスク (サーベイ反映)

- 平衡壁法則 + 標準 SGS 係数は高 Re の緩剥離を見逃す実例あり (Iyer & Malik 2022 Gaussian bump, [survey](../../notes/investigations/turbulence-des-wmles-survey.md) §4)。剥離を含む本命形状では壁モデルだけでなく SGS 係数 (Vreman/WALE 係数、ES σ) の感度確認を必ず行う
- FP32: Newton は flow_float (float) で実施 (log/exp は 1–3 反復 × 少数回で性能・精度とも問題なし)。thermo 呼び出しは内部 double
- Reτ2000 は FP32 + 壁法線高 AR 格子でメッシュ品質ゲート (AR≤1000) と衝突しうる → WMLES 格子は $y^+\approx50$ 起点なので AR は緩い見込みだが VERDICT 必須

## 7. 影響範囲

- 触るファイル: `cuda_forge/wallLaw_d.cuh` (新)、`cuda_forge/wmlesWallModel_d.cu` (新)、`cuda_forge/ransWallFunction_d.cu` (関数昇格のみ)、`cuda_forge/viscousFlux_d.cu` (ゲート 3 値化・AddQWall)、`cuda_forge/boundaryCond_d.cu` (node 等温壁ピン)、`variables.hpp` (`Qw_Wall`)、`input/solverConfig.{hpp,cpp}`、`boundaryCond.hpp` (wall 種の ints に `wallModelLES`)、`main.cpp` (呼び出し 1 箇所 + 体積力)
- 既存ケースへの影響: `wallModelLES` 未指定 (既定 0) で全経路ビット不変。SST 壁関数・DES は無変更
- docs: `methods/turbulence/theory.md` / `implementation.md`、`methods/boundary.md`、`methods/index.md`

## 8. 完了条件

- [ ] `methods/` 更新 (§5-2)
- [ ] 単体検証 §6.1 全通過 + SST 経路ビット不変
- [ ] チャネル検証 §6.2-1 合格 (run パスと VERDICT を変更ログに記録)
- [ ] status を `done` 化し `plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-07-19` — 初稿 (外部ドラフト仕様をコード実態調査に基づきリバイズ)。主な変更: (1) Reichardt/Newton/Normal_Neighbor/`bvar_d["utau"]`/cell 置換分岐/AddTauWall は既存実装の流用に変更、(2) node の τ_w 適用を境界半割面置換から W↔I AddTauWall 再スケールに修正 (壁ノード残差ゼロ化のため)、(3) q_w 用に AddQWall を新設方針、(4) Newton 残差形を $u_\tau f-u_\parallel$ に変更 (0 近傍正則) + warm start 追加、(5) 層流 Pr はその場評価 (config キー非存在)、(6) チャネル検証の前提として体積力 config と等温壁 Phase 2 を明記、(7) M2 平板は乱流流入生成 plan にゲート。
- `2026-07-20` — **実装ステップ 2–6 完了** (§5-2 methods 更新 / §5-3 `wallLaw_d.cuh` 新設 /
  §5-4 単体検証 / §5-5 WMLES カーネル断熱 / §5-6 等温壁 Kader+AddQWall+node 温度ピン)。
  - 新規: `cuda_forge/wallLaw_d.cuh` (Reichardt/Kader/warm-start Newton 昇格・共通化)、
    `cuda_forge/wmlesWallModel_d.{cu,cuh}` (取得層→壁モデル→`bvar_d` + node `Tau_Wall`/`Qw_Wall`)、
    `tools/verify_wall_law.py`。config `wmlesNewtonTol`/`wmlesNewtonMaxIt`/`wmlesPrt`、bcond
    `ints: wallModelLES`、bvar `qwall`、cell 変数 `Qw_Wall`。`viscousFlux_wall_d` の
    `wallTreatment` 3 値化 (2=WMLES) と node W-I AddQWall、`zeroWallMomentumResidual` の
    等温 WMLES `res_roe` ピン (Qw_Wall マーカ)。
  - **単体検証 §6.1: `verify_wall_law.py` 17/17 PASS** (u_τ 回復 rel≤1e-5、層流極限、逆流反転、
    warm start ≤3 反復、Kader log/粘性極限、断熱 q_w=0 厳密)。
  - **SST 経路ビット不変 (cell)**: `case/18.backstep` `run_0128`/`run_0130` の 3-way 比較
    (wallTreatmentSST=1 稼働・wallLaw 昇格後 vs HEAD) で全 12 フィールドの差が atomicAdd
    ノイズフロアと同水準 → 昇格は数値的に不変。node 経路の同等確認は未実施 (残タスク)。
  - ビルドは full rebuild で成功 (clean configure は `-DCMAKE_CUDA_ARCHITECTURES=86` 必須)。
  - **機能スモーク (cell/node 両方)**: `case/24.laminar_channel_bl` `run_wmles_smoke_cell`
    (RK3) / `run_wmles_smoke_node` (implicit)、各 2000 step。両方とも完走・NaN 無し・WMLES 経路
    稼働 (ログ `[WMLES]` カウンタ)。層流チャネル (y+≪1) のため u_τ Newton は設計どおり
    laminar fallback が主で、乱流域 (y+~50-100) での実効検証はチャネル LES (§6.2) 待ち。
  - **残タスク**: §5-7 体積力 config (チャネル駆動)、§5-8 チャネル検証 (Reτ550→2000)、
    node SST 壁関数経路の回帰。
- `2026-07-20` — **等温壁の事前検証 (§5-1 の実態調査を厳密解ケースで実施)**。チャネル検証
  (等温壁前提) に先立ち、純伝導ボックス (350K/300K, 厳密解 = 線形 T・$q_w=\lambda\Delta T/H$) を
  cell/node で検証 (`case/24` `run_isoT_cond*` 一連、詳細は同 README)。
  - **cell 等温壁のバグ発見・修正**: ghost 温度が仕様 ($T_R=2T_w-T_L$, methods/boundary.md) に
    反して $T_w$ 直置きで、壁熱流束係数が 1/2 (dcc=2y₁ のため)。`wall_isothermal_d` を鏡像外挿
    (負温度ガード 0.2Tw) に修正 → **線形厳密解が機械精度の不動点・壁 q_w +0.02%**。
    回帰スイート 4/4 PASS (等温壁は回帰ケース外・既存経路不変)。**注意: 過去の等温壁ケース
    (hifire/flare/cutler 等) の壁熱流束は旧バグの影響下にあり、再評価が必要**。
  - **node 等温壁**: 強制ピンは無いが機能する (weak 閉包)。特性: 壁ノード T は BC 値でなく
    壁 CV 平均に緩む (ny=64 で 0.106 K オフセット、細分化で半減)・第 1 スペーシングの勾配
    −24% 歪み・**中央部の熱流束は厳密 (−0.05%)**。WMLES 等温壁 (強制ピン + AddQWall) とは別経路。
  - **node slip バグ発見 (未修正・別件)**: slip 境界 + 接線密度勾配で市松状スプリアス接線流
    (0.47 m/s) が定在。[調査ノート](../../notes/investigations/node-slip-tangential-density-spurious-flow.md)。
    **チャネル検証 (周期+壁のみ) には影響しない**が、slip を使う node ケースは要注意。
- `2026-07-20` — **node 等温壁の壁ノード温度ピンを実装** (前項の「壁ノード T が CV 平均に緩む」への
  対策。methods/boundary.md「node 等温壁の壁ノード温度ピン」)。状態ピン (`applyNodeIsothermalWallPin`,
  WMLES pin カーネル流用) + `res_roe` ゼロ化 + **block-DPLUR エネルギー行 decouple (`iso_wall_flag_d`,
  SU2 DeleteValsRowi 相当)** の 3 点セット。decouple 無しの状態ピンは implicit で数 step 発散する
  (今回実測 — **既存の WMLES 等温 node + implicit も同リスクだったが本 decouple が共通で解消**)。
  検証: case/24 `run_isoT_condN_node` — 壁ノード T 厳密・q_mid +0.00%・explicit 健全・implicit
  cfl_pseudo≤5 安定 (上限 ~5, 20 発散)。回帰 4/4 PASS。残: W-I 閉包の第 1 スペーシング O(Δ) バイアス −15%。
