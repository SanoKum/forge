# SST 断熱壁 T_aw の流束モデル置換再設計 (SU2 式, node)

## メタ

- **area**: `turbulence / boundary`
- **status**: `in_progress`  <!-- 2026-08-11: 実装済みだが y+30 wall-function ケースで root cause 未特定の
  発散あり。既定 OFF のため安全に commit 可能だが本番投入不可。§9 の診断記録を参照 -->
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5(f) (2026-08-11 全面改訂)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) 同節
- **related_plans**:
  - [`architecture-node-boundary-gradient-dof-only.md`](architecture-node-boundary-gradient-dof-only.md)
    (前提: 境界勾配の owner-state 化。本計画の設計はこれが完了して初めて「bvar `Ts` は勾配に効かない」
    という前提が成立する)
  - [`../archived/turbulence-sst-thermal-wall-function.md`](../archived/turbulence-sst-thermal-wall-function.md)
    (superseded — 旧「弱閉包」設計。本計画がその後継)
  - [`turbulence-sst-thermal-flux-model.md`](turbulence-sst-thermal-flux-model.md) (等温壁 Kader $q_w$ 流束置換。
    同じ「モデル置換の流束層」思想を断熱壁に適用したのが本計画)
- **created**: `2026-08-11`
- **owner**: `sano`

## 1. 目的

SST automatic wall treatment の断熱壁温出力欠陥 (~200 K 過小評価) に対する `sstThermalWallFunction` を、
旧「壁面出力 bvar にのみ書く弱閉包」から **SU2 式の内部粘性流束モデル置換**へ再設計する。
$T_{aw}$ (Crocco 型回復温度) を「壁面出力の飾り」ではなく「熱流束評価用の壁 primitive temperature」として
明示的に流束層へ注入し、場が実際に $T_{aw}$ 近傍へ応答する構造にする。

## 2. 背景 (旧設計の何が問題だったか)

`architecture-node-boundary-gradient-dof-only.md` の前提となったレビューで、旧弱閉包 (`Tsb=Taw_diag` を
bvar `Ts` に書くだけ) には 2 つの問題が判明した:

1. **「出力専用」は不正確だった**: node の境界 Green-Gauss/LSQ 勾配が当時 bvar を読んでいたため、
   $T_{aw}$ は境界勾配の接線補正項を経由して実際に場へ触れていた (ただし極めて弱く抑制されていた)。
2. **効果が弱すぎて場はほぼ動いていなかった**: case/40 `run_0038` (node 5e-3) の実データを検証すると、
   壁ノードの**実状態** $T[W]$ (bell 平均 1195.3 K) は OFF 基準 (1196 K) とほぼ不変で、
   「壁温 1417.9 K = SU2 と 4 K 一致」という報告は**出力配列に上書きした $T_{aw}$ 診断値同士の一致**を
   見ていたに過ぎなかった。生産値 (壁温 1400±15 K) のうち wf+閉包系列 (node run_0038/0040) はこの点で
   汚染されている (詳細: `case/40.nozzle_design_tool/README.md` の訂正注記)。

## 3. スコープ

- **やる**: node × `kind:wall` (断熱) × SST × `wallTreatmentSST==1` × `sstThermalWallFunction==1` の
  壁ノードについて、内部粘性流束 (`viscousFlux_d` 主ループ) の熱流束 compact 項の端点温度を
  $T[W]\to T_{aw,\mathrm{Wall\_Flux}}[W]$ に置換する。対象は Normal_Neighbor の 1 面だけでなく、
  対象壁ノードに接続する**全内部辺**。
- **やらない**: `wall_isothermal` (既に物理壁温がある。等温壁の熱流束モデル置換は
  `turbulence-sst-thermal-flux-model.md` の Kader $q_w$ が別途担う)。cell モード
  (§4 の理由により流束モデル置換の対象外、bvar 出力の意味は変わらず現状維持)。
  $T_{aw}$ を状態としてピンすること (却下済み、§4)。

## 4. 設計方針

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

### 4.3 なぜ暴走しないか (却下版との違い)

流束層の置換後、compact 項は**壁ノード自身の状態 $T[W]$ を含まない**ため、数学的には
「$W$–$I$ 間に $k_{\mathrm{eff}}(T[I]-T_{aw})/d\cdot\delta$ という**モデル流束**を課す」ことと同値になる
($\tau_w$ の AddTauWall と同じ「流束のモデル置換」)。$T_{aw}=T_{rep}+\Delta$ の $T_{rep}$ が通常
その面の相手ノード ($I$) の状態なら、この項は $\propto T[I]-T_{aw}=-\Delta$ で**絶対温度レベルに依らない
有界な量**になる。$T[W]$ 自身はこの有界な寄与と、実状態から計算される GG 補正項・他の内部面流束の合計を
自身の `res_roe=0` から solve する — 他の任意の内部ノードと同じ安定な拡散平衡構造であり、正帰還ループを
持たない。

対照的に、却下済みの「$T_{aw}$ を状態ピン + 壁ノード `res_roe` をゼロ化」構成は、壁 CV 自身のエネルギー
方程式を捨てながら W–I 面の流束だけは高い $T_{aw}$ を使って計算される非保存構成になり、内部ノードを
加熱し続けて $T_{rep}$→$T_{aw}$→(強制コピー)→さらに高い $T_{aw}$ という抵抗のない正帰還で暴走する
(node 1832 K, 実測 run_0038/0039 旧版)。**本設計は状態ピンをせず壁ノード自身の `res_roe` も生かした
ままなので、この機構は生じない**。

### 4.4 エネルギー保存 (必須維持事項)

- 断熱壁半割面の熱流束は変更なく厳密 $q_{wall}=0$。
- W–I 内部面のエネルギー流束は両端へ $+F/-F$ で保存的に加算 (置換したのは compact 項の**値**であって
  加算則ではない)。
- 断熱壁ノードの `res_roe` はゼロ化しない。
- implicit のエネルギー行 decouple はしない (等温壁の `iso_wall_flag` 機構とは無関係)。
- `roe` や状態 $T[W]$ を $T_{aw}$ にピンしない。

## 5. 実装ステップ

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

## 6. 検証

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

## 7. 影響範囲

- `solver_density_cuda/variables.hpp`, `cuda_forge/ransWallFunction_d.cu`, `cuda_forge/nodeWallDirichlet_d.cu`,
  `cuda_forge/nodeWallDirichlet_d.cuh`, `cuda_forge/viscousFlux_d.cu`
- `methods/turbulence/theory.md` §6.5(f) (更新済み), `implementation.md` 同節 (更新済み)
- `sstThermalWallFunction=1` を使う全 node 断熱壁ケース (case/40 等)。既定 OFF のためそれ以外は無影響。

## 8. 完了条件

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
  - **中断判断**: 効果のあった修正 (境界勾配 DOF-only 化, `architecture-node-boundary-gradient-dof-only.md`)
    は確定させ、本機能は **既定 OFF (`sstThermalWallFunction: 0`) のため他への影響なくコミット可能**
    と判断し、本 plan は `in_progress` のまま次セッションへ引き継ぐ。
  - **次セッションへの申し送り**: (1) 中間ステップ (例 step 1400–1480 を 10 step 間隔) をダンプし
    Taw_diag と T の時間発展を追って崩壊の起点を特定する。(2) `compute_wall_friction_sst_d` の
    Newton 反復や `Ut`/`utau` 自体が先に発散していないか確認 (Taw_diag は `T_rep + r·Ut²/(2cp)` で
    Ut に 2 乗依存するため、Ut の発散が Taw_diag を非線形に増幅する経路を疑う)。(3) 疑わしいなら
    Taw_Wall_Flux の値に上限 (例 Tt の数倍) を設けるクランプで応急止血できるか確認 (根治ではないが
    安全弁として)。(4) 深追いする場合は effort の高いモデル (Fable 5) での再着手を推奨
    (2026-08-11 セッションは Sonnet で実施し、この機構の収束まで到達できなかった)。
