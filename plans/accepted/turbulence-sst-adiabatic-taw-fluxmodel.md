# SST 断熱壁 T_aw の流束モデル置換再設計 (SU2 式, node) — 検証の結果**注入は棄却**、最終形は出力専用

## メタ

- **area**: `turbulence / boundary`
- **status**: `done`  <!-- 2026-08-11: 流束注入案は実測発散 (§9) で棄却。最終形 = Taw は境界出力
  (Tsb/Taw_diag) 専用、W-I 内部粘性拡散は常に DOF 状態のみ、壁熱は境界半割面 q_w (断熱=厳密0)。
  この「棄却の記録+最終方針」= output-only fallback が本 plan の成果物。SU2 熱結合の再実装は
  後継 plans/active/turbulence-sst-su2-taw-coupling.md (experimental mode 2) -->
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5(f)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) 同節
- **related_plans**:
  - [`../active/turbulence-sst-su2-taw-coupling.md`](../active/turbulence-sst-su2-taw-coupling.md)
    (後継: SU2 熱結合の再実装, experimental `sstThermalWallFunction: 2`)
  - [`architecture-node-boundary-gradient-dof-only.md`](architecture-node-boundary-gradient-dof-only.md)
    (前提: 境界勾配の owner-state 化)
  - [`../archived/turbulence-sst-thermal-wall-function.md`](../archived/turbulence-sst-thermal-wall-function.md)
    (superseded — 初代「弱閉包」設計)
  - [`turbulence-sst-thermal-flux-model.md`](../active/turbulence-sst-thermal-flux-model.md) (等温壁 Kader $q_w$)
- **created**: `2026-08-11`
- **owner**: `sano`

## 0. 最終方針 (2026-08-11 確定 — ユーザレビューによる)

**GG/LSQ の owner-state-only 方針
([architecture-node-boundary-gradient-dof-only](../accepted/architecture-node-boundary-gradient-dof-only.md))
は正しく維持。撤回したのは $T_{aw}$ を W–I 内部辺の粘性拡散へ注入する処理だけ**。

正しい熱経路:

```
物理壁面
  ├─ モデル壁面温度 Tsb (境界出力)
  └─ 壁熱流束 q_w
        ↓ 境界半割面流束 (res_roe[W] へ直接、外部境界流束なので ±F ペアにしない)
壁ノード W
        ↓ 通常の内部粘性拡散 (compact 温度差 = T[I]−T[W]、DOF 状態勾配)
内部ノード I
```

1. **W–I 内部辺**: 粘性熱拡散は常に DOF 状態 (`T[I]−T[W]`, DOF 勾配) で評価。
   `Tsb`/`Taw_diag` は一切参照しない。±F 保存的加算は不変。
2. **GG/LSQ 勾配**: owner-state-only を維持 (bvar 不参照、LSQ は内部隣接 node のみ)。
3. **壁面→壁ノード W の熱流束**: 境界半割面流束 `q_w` として実装する。
   断熱壁: q_w=0 (厳密, `viscousFlux_wall_d` の `adiabaticWall` 分岐で実装済み) /
   指定熱流束壁: 指定 q_w (将来) / 熱壁関数: 壁関数の q_w (将来) /
   等温 Dirichlet 壁: T[W]=Tw 強制、q_w は結果診断 (Dirichlet と Neumann を同時強制しない)。
4. **T_aw の役割**: `Taw_diag[W]` と `Tsb` (壁関数が再構成した物理壁面温度のモデル値、境界出力)
   に限定。状態ピン・GG/LSQ 勾配・W–I 内部流束・res_roe ゼロ化には使わない。
   出力上 `T[W]` (壁半 CV の計算状態温度) と `Tsb`/`Taw_diag` (モデル壁面温度) を明確に区別する。
5. **壁エネルギー残差**: SST 断熱壁の `res_roe[W]` はゼロ化しない。壁半 CV の熱収支は
   W–I 通常内部粘性流束 + 内部辺の粘性仕事 + 物理壁面熱流束 (q_w=0) の和として通常どおり解く。

---

# 【以下 §1–§8 は棄却済み旧案の記録】

**注意 (2026-08-11 整理)**: §1〜§8 は「W–I 流束注入」旧案の起票時仕様であり、**現行仕様ではない**。
現行仕様は §0 (output-only baseline)。旧案は §9 の実測発散により棄却された。
SU2 と同じ熱的結合の再実装は後継 plan
[`turbulence-sst-su2-taw-coupling.md`](../active/turbulence-sst-su2-taw-coupling.md)
(experimental mode `sstThermalWallFunction: 2`) が別途進める — 本 plan は
**output-only fallback と失敗履歴の記録**として保持する。

## 1. 目的 【棄却済み旧案】

SST automatic wall treatment の断熱壁温出力欠陥 (~200 K 過小評価) に対する `sstThermalWallFunction` を、
旧「壁面出力 bvar にのみ書く弱閉包」から **SU2 式の内部粘性流束モデル置換**へ再設計する。
$T_{aw}$ (Crocco 型回復温度) を「壁面出力の飾り」ではなく「熱流束評価用の壁 primitive temperature」として
明示的に流束層へ注入し、場が実際に $T_{aw}$ 近傍へ応答する構造にする。

## 2. 背景 (旧設計の何が問題だったか) 【棄却済み旧案】

`architecture-node-boundary-gradient-dof-only.md` の前提となったレビューで、旧弱閉包 (`Tsb=Taw_diag` を
bvar `Ts` に書くだけ) には 2 つの問題が判明した:

1. **「出力専用」は不正確だった**: node の境界 Green-Gauss/LSQ 勾配が当時 bvar を読んでいたため、
   $T_{aw}$ は境界勾配の接線補正項を経由して実際に場へ触れていた (ただし極めて弱く抑制されていた)。
2. **効果が弱すぎて場はほぼ動いていなかった**: case/40 `run_0038` (node 5e-3) の実データを検証すると、
   壁ノードの**実状態** $T[W]$ (bell 平均 1195.3 K) は OFF 基準 (1196 K) とほぼ不変で、
   「壁温 1417.9 K = SU2 と 4 K 一致」という報告は**出力配列に上書きした $T_{aw}$ 診断値同士の一致**を
   見ていたに過ぎなかった。生産値 (壁温 1400±15 K) のうち wf+閉包系列 (node run_0038/0040) はこの点で
   汚染されている (詳細: `case/40.nozzle_design_tool/README.md` の訂正注記)。

## 3. スコープ 【棄却済み旧案】

- **やる**: node × `kind:wall` (断熱) × SST × `wallTreatmentSST==1` × `sstThermalWallFunction==1` の
  壁ノードについて、内部粘性流束 (`viscousFlux_d` 主ループ) の熱流束 compact 項の端点温度を
  $T[W]\to T_{aw,\mathrm{Wall\_Flux}}[W]$ に置換する。対象は Normal_Neighbor の 1 面だけでなく、
  対象壁ノードに接続する**全内部辺**。
- **やらない**: `wall_isothermal` (既に物理壁温がある。等温壁の熱流束モデル置換は
  `turbulence-sst-thermal-flux-model.md` の Kader $q_w$ が別途担う)。cell モード
  (§4 の理由により流束モデル置換の対象外、bvar 出力の意味は変わらず現状維持)。
  $T_{aw}$ を状態としてピンすること (却下済み、§4)。

## 4. 設計方針 【棄却済み旧案】

### 4.1 3 層責務

- **状態層**: 壁ノード $T[W]$ は**通常の DOF**。Dirichlet ピンしない (`res_roe` も生かす)。
- **モデル層**: `ransWallFunction_d.cu compute_wall_friction_sst_d` が毎ステップ $T_{aw}=T_{rep}+r U_{t,rep}^2/2c_p$
  を計算し `Taw_diag` に書く (既存)。`applySstThermalWallFunction` (`nodeWallDirichlet_d.cu`) が対象壁ノードの
  per-CV マーカ `Taw_Wall_Flux[W]` に `Taw_diag[W]` を書く (非対象は $-1$)。$T_{aw}$ は**前ステップ収束状態
  ベースの lagged wall-model 値**として扱う (他の壁関数量と同じ)。壁面出力 (`res_wall_*` の bvar `Ts`) にも
  引き続き `Taw_diag` を書く (`set_wall_taw_output_d` は残す — 今度は流束へ実際に効く値の表示として妥当)。
- **流束層**: `viscousFlux_d` 主ループ (node 内部双対面) で、端点 (ic0 または ic1) が対象壁ノードなら
  熱流束の compact/法線差分項 `(Ts[ic1]-Ts[ic0])/dcc` の**その端点の温度だけ** `Taw_Wall_Flux` に差し替える
  (両端が対象なら両端とも)。面平均 Green-Gauss 勾配の接線補正項 (`dTdxf*k_x+...`) は変更しない (実状態の
  勾配のまま)。熱伝導率・粘性係数は当該ステップの状態物性のまま (Taw で再評価しない)。

### 4.2 マーカー配列

`Tau_Wall`/`Qw_Wall` と同型の per-CV `flow_float* Taw_Wall_Flux` (`variables.hpp` に登録、既定 $-1$)。
初期化は `ransWallFunction_d.cu init_wf_pk_d` (Tau_Wall 等と同じタイミングで $-1$ 化) に追加する。

### 4.3 ~~なぜ暴走しないか~~ → 【撤回・訂正 (2026-08-11)】この安定性主張は誤りだった

本節に書いていた 3 つの主張 —「モデル流束が有界なので安定」「res_roe を生かせば壁状態が平衡を解く」
「正帰還を持たない」— は**すべて撤回する**。実測 (§9) で、Taw 端点置換は正帰還ではなく**逆方向の
不安定 (壁エネルギー DOF の異常冷却)** で発散した。

**実際に起きた機構**: compact 項の端点温度を $T[W]\to T_{aw}$ ($>T[W]$) に置換すると、この面の流束は
「高温の仮想壁面 → 内部」向きに増大し、±F の保存的加算により **$W$ の半 CV から実在しないエネルギーが
流出し続ける** (物理壁は断熱で境界からの流入ゼロ、$W$ には補償source がない)。しかも置換後の
compact 項は $T[W]$ に依存しないため、**$W$ 自身の温度が下がってもこの流出は減らない (復元項の喪失)**。
結果、$T[W]$ は EOS 床まで単調に冷却し (step1000 で min 27.1 K / 15.5 K, §9)、近壁の物性・SST が
崩壊して roOmega NaN に至る。「有界な流束」であっても**壁半 CV 単体の収支を閉じない外部吸熱**として
働けば発散する — 旧・状態ピン版が「壁 CV の収支を捨てて加熱し続けた」のと鏡像の失敗であり、
教訓は同じ: **W–I 内部辺に手を入れる方式は、壁半 CV 単体のエネルギー収支を必ず監査せよ**
(内部辺全体の ±F 保存だけでは不十分)。

最終方針 (§0) はこの教訓に基づき、モデル値の場への注入を諦め、熱の出入りは物理壁面の境界半割面
流束 q_w (断熱=0) だけに限定する。

### 4.4 エネルギー保存 (必須維持事項)

- 断熱壁半割面の熱流束は変更なく厳密 $q_{wall}=0$。
- W–I 内部面のエネルギー流束は両端へ $+F/-F$ で保存的に加算 (置換したのは compact 項の**値**であって
  加算則ではない)。
- 断熱壁ノードの `res_roe` はゼロ化しない。
- implicit のエネルギー行 decouple はしない (等温壁の `iso_wall_flag` 機構とは無関係)。
- `roe` や状態 $T[W]$ を $T_{aw}$ にピンしない。

## 5. 実装ステップ 【棄却済み旧案】

1. `variables.hpp` に `Taw_Wall_Flux` を登録。
2. `ransWallFunction_d.cu init_wf_pk_d` で $-1$ 初期化に追加。
3. `nodeWallDirichlet_d.cu`: `applySstThermalWallFunction` に隣接して、対象壁ノード判定
   (`bc.bcondKind=="wall" && sstThermalWallFunctionActive(cfg) && discretization=="node"`) で
   `Taw_Wall_Flux[W]=Taw_diag[W]` を書くカーネル/wrapper を追加 (既存 `set_wall_taw_output_d` の出力書き込みは
   残す)。
4. `viscousFlux_d.cu` 主ループ: `Taw_Wall_Flux` 引数を追加し、`heatflux` の compact 項計算直前に
   端点温度を条件付き差し替え。
5. `architecture-node-boundary-gradient-dof-only.md` の実装 (境界勾配の owner-state 化) が先に完了して
   いることを前提とする (順序依存: 先に境界勾配を bvar 非依存にしないと、旧弱閉包と流束置換が二重に
   場へ効いてしまう)。
6. full rebuild。

## 6. 検証 【棄却済み旧案】

- **case/40 A/B**: 断熱壁ケースで旧弱閉包 (OFF 相当を含む) と新流束置換を比較。壁ノードの**実状態**
  $T[W]$ が $T_{aw}$ 近傍へ実際に応答すること (旧版のような「OFF とほぼ不変」にならないこと) を確認。
- **エネルギー収支監査**: 全内部面エネルギー流束の $\pm$ 和が丸め誤差内ゼロ。断熱壁半割面 $q_{wall}=0$。
  断熱壁 `res_roe` が生きており (非ゼロ・非射影)、implicit 行が decouple されていないこと。
- **SU2 同一メッシュ比較**: `su2-cross-check.md` 手順で断熱壁温を同一メッシュ SU2 と比較。
- **OFF 回帰**: `sstThermalWallFunction=0` でビット不変。
- **cell 非退行**: cell モードは本計画の対象外なのでビット不変 (出力専用のまま)。
- NaN/Inf・メッシュ品質・`check_convergence.py`・`check_quasisteady.py` を実行し、未収束を「収束・一致」と
  報告しない。
- 新規 `run_NNNN_<slug>` で実行し、case README の run 一覧を同期する (既存 run を使い回さない)。

## 7. 影響範囲 【棄却済み旧案】

- `solver_density_cuda/variables.hpp`, `cuda_forge/ransWallFunction_d.cu`, `cuda_forge/nodeWallDirichlet_d.cu`,
  `cuda_forge/nodeWallDirichlet_d.cuh`, `cuda_forge/viscousFlux_d.cu`
- `methods/turbulence/theory.md` §6.5(f) (更新済み), `implementation.md` 同節 (更新済み)
- `sstThermalWallFunction=1` を使う全 node 断熱壁ケース (case/40 等)。既定 OFF のためそれ以外は無影響。

## 8. 完了条件 【棄却済み旧案】

- [x] 関連 `methods/turbulence/theory.md`・`implementation.md` を更新済み
- [ ] 実装・検証完了 (本計画の §6 を満たす)
- [ ] 本計画の `status` を `done` に変更し、§9 に変更ログを記載
- [ ] ファイルを `plans/active/` → `plans/accepted/` へ移動
- [ ] [`plans/README.md`](../README.md) の一覧を同期
- [ ] `case/40.nozzle_design_tool/README.md` の生産値 (壁温 1400±15 K) の wf+閉包系列を新設計の値で再確認・訂正

## 9. 変更ログ

- `2026-08-11` — 起票。旧 `turbulence-sst-thermal-wall-function.md` (accepted) を superseded として
  archived へ移動し、本計画を新設計として起票。docs (theory.md/implementation.md) 更新済み、実装は
  `architecture-node-boundary-gradient-dof-only.md` の完了待ち。
- `2026-08-11 (2)` — **実装完了・検証で発散を確認、root cause 未特定のまま中断**。
  `Taw_Wall_Flux` per-CV マーカ (`variables.hpp`) + `set_wall_taw_flux_marker_d`
  (`nodeWallDirichlet_d.cu`) + `viscousFlux_d` 主ループの compact 項端点置換を実装しビルド完了。
  - **A/B で単離**: `sstThermalWallFunction=0` の同一入力 (run_0055) は一般 GG/LSQ 再設計と合わせて
    `ALL PASS`。`=1` (run_0053, 全内部辺へ適用する当初仕様) は node y+30 wall-function ケースで
    **step 1686 付近で roOmega が NaN 発散**。発散ノードは壁ノードではなく第一内部行 (壁から
    節点間隔 1 個分, `bad&wall=0` を確認)。
  - **仮説① (棄却)**: 「全内部辺への適用が、代表点以外の接線方向隣接ノードとの温度差まで
    流束化し不安定化させている」と推測し、`Taw_Rep_Id` (壁ノードの代表点 CV id を記録する新配列) を
    追加して**代表点との辺のみ**に適用を制限 (run_0057)。**効果なし** — step 1480 (制限前より早い)
    で同型発散。この仮説は誤りだったか、または副次的要因に過ぎない。
  - **診断データ**: 発散直前ダンプ (`res_nan_1480.h5`) で `Taw_diag` 自体が既に発散済み
    (max 4.1e8 K)・`T` 状態が破綻 (min 1e-4 K) — 何が最初に崩れたか (Taw_diag の計算自体が
    先に不安定化したのか、流束置換が場を壊してから Taw_diag が異常値を拾ったのか) は
    中間ステップのダンプが無く未特定。
  - **中断判断 (当時)**: 効果のあった修正 (境界勾配 DOF-only 化,
    `architecture-node-boundary-gradient-dof-only.md`) は確定させ、本機能は既定 OFF のため
    他への影響なくコミット可能と判断し、いったん `in_progress` で commit した (946a98f6)。
- `2026-08-11 (3)` — **root cause 判明 (ユーザレビュー)・流束注入を全面撤回・最終方針確定 (§0)**。
  - **崩壊起点の特定**: step 1000 の壁ノード**実状態温度** T[W] を 3 run で比較:

    | run | 構成 | T[W] min | T[W] mean | 結果 |
    | --- | --- | --- | --- | --- |
    | `run_0055_node_yp30_dofonly_tawoff` | Taw OFF | 1100.7 K | 1270.3 K | **ALL PASS** |
    | `run_0053_node_yp30_dofonly_regr` | Taw 全辺置換 | **27.1 K** | 640.0 K | DIVERGED (step1686) |
    | `run_0057_node_yp30_taw_repedge` | 代表辺限定 | **15.5 K** | 627.5 K | DIVERGED (step1480) |

    最初の異常は roOmega ではなく **Taw 端点置換による壁エネルギー DOF の異常冷却** (§4.3 の
    訂正参照: 置換流束が T[W] 非依存になり壁半 CV の復元項を失い、補償のない外部吸熱として
    EOS 床まで排熱し続ける)。roOmega NaN はその下流の症状。
  - **撤回実装**: `viscousFlux_d` の `Taw_Wall_Flux`/`Taw_Rep_Id`/compact 端点置換/代表辺限定/
    ゲート関数、`set_wall_taw_flux_marker_d`、`variables.hpp` の 2 配列、`ransWallFunction_d` の
    `Taw_Rep_Id` 書き込みをすべて削除。W–I 内部熱拡散は常に DOF 状態 (`Ts[ic1]−Ts[ic0]` +
    DOF 勾配) に復帰。`set_wall_taw_output_d` (Tsb=Taw の境界出力) のみ残置。
    viscousFlux_d に禁止事項コメントを恒久記載。
  - **§4.3 の誤った安定性主張 (「有界なので安定」「res_roe を生かせば平衡」「正帰還なし」) を
    撤回・訂正** (同節参照)。
  - **最終方針**: §0 のとおり。Taw は境界出力 (Tsb/Taw_diag) 専用の「壁関数が再構成した物理壁面
    温度のモデル値」。壁熱は境界半割面 q_w (断熱=厳密 0, 既存 `adiabaticWall` 分岐)。SST 断熱壁の
    res_roe[W] はゼロ化しない (該当ゼロ化経路が無いことをコードで確認済み)。
  - **検証 (撤回後, `run_0058_node_yp30_taw_outputonly` = run_0055 と同一入力で
    `sstThermalWallFunction: 1`)**:
    1. `check_convergence.py` **ALL PASS** (全列 3.5–5.1 桁低下, OFF の run_0055 と同一性状)。
    2. 場の解は OFF と max rel 2.1e-6 (ro) 〜 5.0e-5 (T) で一致 = node run-to-run atomicAdd
       再現幅内 → **W–I 内部熱流束は Taw ON/OFF で不変**を実証。
    3. 壁ノード実状態 T[W]: min 1100.6 K / mean 1270.3 K = OFF と同一水準 (**低温床崩壊なし**)。
    4. 差分は壁面出力 Tsb のみ: OFF 1134.9 K (状態値) → ON 1412.7 K (Taw モデル値, bell 平均)
       — SU2 壁関数出力 1418.9 K (run_0042) とモデル値同士で 6 K 差。
    5. 断熱壁境界エネルギー流束は厳密 0 (`viscousFlux_wall_d` `adiabaticWall` 分岐)・
       res_roe[W] 有効 (ゼロ化経路なし)・内部辺 ±F 保存 — いずれもコード検証。
    6. 非ゼロ q_w の大域収支検証 (Σ q_w A = 流体エネルギー流入) は、指定熱流束壁/熱壁関数 q_w を
       実装する将来 plan の必須検証項目として繰り越す (現スコープは断熱 q_w=0 のため対象なし)。
