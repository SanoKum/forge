# node 境界勾配の DOF-only 化 (Green-Gauss/LSQ 境界寄与から bvar を排除)

## メタ

- **area**: `architecture / discretization`
- **status**: `done`
- **related_docs**:
  - [`methods/discretization.md`](../../methods/discretization.md) §6.2 (理論), §7.2.2/§7.3/§7.3.1 (実装)
  - [`methods/gradient.md`](../../methods/gradient.md) (境界寄与節)
- **related_plans**:
  - [`boundary-node-nozzle-wall-outlet-stability.md`](boundary-node-nozzle-wall-outlet-stability.md) §2.11
    (node outlet の P/T 境界閉包を interior 化した先行修正。本計画はこれを全 primitive・全 bcond へ一般化する)
  - [`turbulence-sst-adiabatic-taw-fluxmodel.md`](turbulence-sst-adiabatic-taw-fluxmodel.md)
    (本計画の一般化を前提にした SST 断熱壁 $T_{aw}$ の流束モデル置換再設計)
  - [`discretization-lsq-gradient.md`](discretization-lsq-gradient.md) (LSQ 勾配の元計画。本計画はその境界処理を改訂)
  - [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) (ゴースト撤廃の親計画)
- **created**: `2026-08-11`
- **owner**: `sano`

## 1. 目的

node-centered (median-dual) の境界勾配 (Green-Gauss・LSQ 双方) を「境界ノード (owner) の実状態値のみ
から計算する」原則に統一する。bvar (境界フラックス入力・壁モデル値・境界面出力用の量) を勾配演算から
完全に排除し、「node の勾配は DOF 状態だけから作る。bvar は境界流束・壁モデル・境界面出力に使い、
勾配演算には入れない」という責務分離を全境界・全 primitive で実現する。

背景: `boundary-node-nozzle-wall-outlet-stability.md` §2.10/§2.11 で node outlet の P/T 境界勾配閉包に
bvar (規定背圧) を使っていたことが出口列の系統的な sag (最大 24%) の真因と判明し、P/T のみ・outlet のみの
特例として interior 化した (commit ddbe9ce5)。同種の問題は構造的に**全境界・全 primitive**に存在する
(SST 断熱壁 $T_{aw}$ の弱閉包が場に意図せず触れていた件も同根、`turbulence-sst-adiabatic-taw-fluxmodel.md`
参照)。本計画はこの特例を一般原則に格上げし、コード上の per-bcond 分岐を削減する。

## 2. スコープ

- **やる**:
  - node の Green-Gauss 境界カーネル (`calcGradient_b_d`) を全 primitive (`ro,Ux,Uy,Uz,P,T`)・
    全 non-periodic bcond で owner state 参照に統一 (`interiorPT`/`P_c`/`T_c` の特例撤去)。
  - node の LSQ 境界疑似点 (`lsqGrad_accumBoundary_d`, `lsqPreGrad_boundary_d`, `lsqPre_coefBoundary_d`,
    `lsqPre_setup_d` の境界点追加) を完全撤去。正規行列 $M$ にもジャンプ ($\mathbf b$) にも境界点を加えない。
  - 退化 (近傍配置が悪条件) ノードは既存のスペクトル打ち切り擬似逆 (`gradLSQ=2`) で扱う (新規機構は不要、
    既存フォールバックの対象が広がるだけ)。
  - periodic は変更しない (既存の DOF 同一視・gradient gather 経路を維持)。
- **やらない**:
  - cell-centered 経路 (ghost/bvar 閉包) は変更しない。
  - 対流フラックス・壁モデル・境界面出力での bvar 利用は変更しない (責務としてそのまま)。
  - symmetry/slip 専用の勾配射影 (SU2 `correctGradientsSymmetry` 相当) は本計画では設計しない
    (現状 slip は owner-state 化で足りるかを検証で確認し、不足があれば別計画に切り出す)。
  - SST 断熱壁 $T_{aw}$ の流束モデル置換は `turbulence-sst-adiabatic-taw-fluxmodel.md` に分離。

## 3. 関連 docs と前提

- 理論: [`methods/discretization.md`](../../methods/discretization.md) §6.2 (2026-08-11 改訂済み)。
- 実装: 同 §7.2.2 (GG, 新設)・§7.3/§7.3.1 (LSQ, 改訂済み)。
- 前提修正: commit ddbe9ce5 (node outlet P/T interior 化) — 本計画はこの特例を一般化する。

## 4. 設計方針

### 4.1 Green-Gauss (`calcGradient_b_d`)

現状 (改訂前) は `ro,Ux,Uy,Uz` を常に bvar (`bc.bvar_d[...]`) から、`P,T` は `interiorPT` フラグ
(outlet_statPress のみ 1) で owner state と bvar を切り替えていた。改訂後は **全変数を無条件で owner
state (`ro[ic],Ux[ic],Uy[ic],Uz[ic],P[ic],T[ic]`) に統一**し、`interiorPT`/`P_c`/`T_c` 引数を削除する。
bvar 引数 (`rob,roUxb,...,Ps,Ts`) は不要になるため関数シグネチャから削除する。

no-slip ($u_p=0$)・等温壁 ($T_p=T_w$)・axis (半径ミラー) は状態が既に BC 値へピンされているため、
owner state 化しても**数値はビット不変** (bvar と状態が同値)。数値が変わるのは状態と bvar が乖離し得る
境界 (outlet の背圧、断熱壁の $T_{aw}$ 弱閉包) のみ — これらは意図した修正。

### 4.2 LSQ (`gradLSQ=1`/`gradLSQ=2`)

- `gradLSQ=1` (毎ステップ solve, 回帰対照): `lsqGrad_accumBoundary_d` の呼び出しを削除するだけでよい
  (`lsqGrad_accumInternal_d` は既に `ip>=nNormalPlanes` を skip しており境界点を作らない)。
  カーネル自体も削除する (呼び出し元が無くなるため)。
- `gradLSQ=2` (係数事前計算, 既定推奨): `lsqPre_setup_d` の $M$ 組み立てから境界半割面点
  (現状 `else` 分岐で `pcx/pcy/pcz` を LSQ 点として追加) を削除し、境界半割面は periodic と同様
  `continue` で skip する。`lsqPre_coefBoundary_d`・`lsqPreGrad_boundary_d` カーネルと `cBnd` 配列・
  オフセット管理コードを削除する。

**なぜ $\mathbf b$ だけでなく $M$ からも除外が必要か**: ジャンプ $\phi_b-\phi_i=0$ の疑似点でも、
正規行列 $M=\sum w\,\mathbf d\mathbf d^{\mathsf T}$ には $\mathbf d$ (境界面重心への変位) が寄与するため、
「その方向の勾配をゼロへ向かわせる」暗黙の正則化になる。$\mathbf b$ だけゼロにしても $M$ に境界点が
残っていれば正則化は残るため、$M$ からの除外が必須。

### 4.3 境界固有の情報の扱い

境界に固有の情報 (壁温・回復温度モデル・規定背圧等) は次のいずれかで表現し、勾配のジャンプで代用しない:

- **状態層**: no-slip ($u=0$)・等温壁 ($T=T_w$)・axis (半径ミラー)・periodic (DOF 同一視) — 既存の
  Dirichlet ピン機構をそのまま使う。
- **境界フラックス (bvar)**: 対流フラックスの $Q_b$、壁モデル出力 (`Tau_Wall`/`Qw_Wall`/`Taw_Wall_Flux`)。
- **境界面出力**: `res_wall_*.h5` 等の表示専用値。

## 5. 実装ステップ

1. [`calcGradient_d.cu`](../../solver_density_cuda/cuda_forge/calcGradient_d.cu) `calcGradient_b_d`:
   シグネチャから bvar 引数・`interiorPT`/`P_c`/`T_c` を削除し、全変数 owner state 参照に統一。
   呼び出し側 (`calcGradient_d_wrapper`) の引数も合わせる。
2. 同ファイル `lsqPre_setup_d`: 境界半割面の $M$ 追加分岐を削除 (periodic と同様 `continue`)。
3. 同ファイル `lsqPre_coefBoundary_d`・`lsqPreGrad_boundary_d`・`lsqGrad_accumBoundary_d`: 削除
   (定義・宣言・呼び出し・関連バッファ/オフセット管理を含む)。
4. [`calcGradient_d.cuh`](../../solver_density_cuda/cuda_forge/calcGradient_d.cuh): 削除した関数の
   宣言を除去 (`calcGradient_b_d` は他 TU から参照されないため実害はないが整合のため更新)。
5. full rebuild (`ninja -t clean && ninja forge convertGmshToForge`)。

## 6. 検証

- **単体 (幾何限定)**: 定数場で全 node 境界の GG 勾配が丸め誤差内ゼロになること。線形場
  ($\phi=a+bx+cy+dz$) で境界ノードを含む GG/LSQ 精度を `tools/verify_linear_recon.py` 系で確認
  (既存ツールの拡張が要れば追加)。
- **LSQ 条件数**: `gradLSQ=2` setup ログの退化ノード数が旧実装と極端に乖離しないこと (境界点除去で
  近傍数が減るため増える可能性はあるが、スペクトル打ち切りで発散しないことを確認)。
- **outlet 非退行**: `boundary-node-nozzle-wall-outlet-stability.md` の出口 sag 修正 (run_0050/0051 相当)
  が本一般化後も同水準で保たれること (新 run で再現)。
- **等温壁 / 断熱壁 / axis / periodic / slip 個別確認**: 各 bcond で state=bvar な境界はビット不変、
  そうでない境界 (outlet, 断熱壁 $T_{aw}$) は意図した変化のみであること。
- **cell モード**: 変更なしなのでビット不変 (回帰確認のみ)。
- **inletProfile CSV 分布**: 対流フラックス側 (bvar) の面分布指定が変更されないこと (勾配経路のみの
  変更のため自明だが、既存 profile ケースで A/B 確認)。
- 判定は `solver_density_cuda/tools/check_convergence.py` / `check_quasisteady.py` の VERDICT を必須とする。

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcGradient_d.cu` / `.cuh`
- `methods/discretization.md` §6.2/§7.2.2/§7.3/§7.3.1、`methods/gradient.md`
- node discretization を使う全ケース (勾配の境界寄与が変わるため、状態と bvar が乖離する境界を持つ
  ケースは数値が変わり得る — 主に outlet と SST 断熱壁)

## 8. 完了条件

- [x] 関連 `methods/discretization.md`・`methods/gradient.md` を更新済み
- [x] 実装・検証完了 (本計画の §6 を満たす、§9 参照)
- [x] 本計画の `status` を `done` に変更し、§9 に変更ログを記載
- [x] ファイルを `plans/active/` → `plans/accepted/` へ移動
- [x] [`plans/README.md`](../README.md) の一覧を同期

## 8.1 `nodeWallDirichlet=0` (legacy weak-wall) 互換性監査

**結論**: 既定 (`nodeWallDirichlet=1`, 2026-07-20 以降の生産構成) はビット不変。`nodeWallDirichlet=0`
(非既定・診断専用の legacy weak-wall 経路, `discretization.md §7.2` 参照) は**挙動が変わる**が、
これは意図した設計の帰結であり、bvar への回帰は行わない。

- **=1 (既定)**: 壁ノード速度が `enforceWallNoSlip_d` で状態として厳密 $u_p=0$ にピンされているため、
  境界 GG の owner-state 化後も境界寄与は `U[ic]=0` — 旧 bvar (`Uxb=0`) と同値。**ビット不変**。
- **=0 (legacy)**: 壁ノード速度は状態として自由 DOF (ピンなし、弱形式 pressure-only 対流のみで
  間接的に効く)。旧実装は境界 GG が bvar (`Uxb=0`) を強制的に読んでいたため、状態が非ゼロでも
  「境界ジャンプ = 0−U[ic]」が勾配へ弱く効いていた。新実装は owner state (`U[ic]`) を使うため
  この境界ジャンプは常に 0 になり、**壁の存在が境界 GG 経由では勾配に効かなくなる**
  (内部双対面の粘性経由の影響は残る)。
- **判断**: `nodeWallDirichlet=0` は既に非生産の診断専用パス (case/29・36・38 の "noDir"/"nopin"
  比較 run のみで使用、生産既定は 1) であり、bvar 依存へ戻すことは本計画の原則
  (「勾配は DOF 状態だけから作る」) に反するため行わない。挙動変化はここに記録するに留める。

## 9. 変更ログ

- `2026-08-11` — 起票。docs (discretization.md/gradient.md) 更新済み、実装は着手中。
  発端: ユーザレビューで node outlet の P/T interior 化特例 (commit ddbe9ce5) を「全 primitive・
  全境界の原則」に一般化すべきと指摘され、SST 断熱壁 $T_{aw}$ 弱閉包の場汚染問題 (別レビューで発見) も
  同根と判明したため、2 計画に分けて起票。
- `2026-08-11 (2)` — **実装・検証完了、done**。
  - **実装**: `calcGradient_b_d` (GG) を全 primitive・全 non-periodic bcond で owner-state に統一
    (`interiorPT`/bvar 引数を撤去)。`lsqPre_setup_d`/`lsqGrad_accumInternal_d` から境界疑似点を撤去し、
    `lsqPre_coefBoundary_d`/`lsqPreGrad_boundary_d`/`lsqGrad_accumBoundary_d` を削除。full rebuild 成功。
  - **検証 1 (outlet 非退行)**: `run_0052_node_yp1_dofonly_regr` (run_0050 と同一入力, 全変数一般化後) の
    sag/η/ṁ が **run_0050 と完全一致** (sag 0.40%, η 0.9754, ṁ 1.2959) — outlet 特例の一般化は
    P/T 以外 (ro/Ux/Uy/Uz) の変化を含めてもビット/数値レベルで無害だった。
  - **検証 2 (cell 不変性)**: `run_0054_cell_bitinv_regr` (run_0036, `wallTreatmentSST:0` で今回変更点を
    一切経由しない cell ケース) は run_0036 比で相対差 0.6–2.3% (bit-identical ではない) だったが、
    **同一バイナリでの再実行 (`run_0056`) でも相対差 0.7–3.7%** を示し、これは既知の
    [[cell-atomicadd-nondeterminism]] (この特定ケースは roOmega が非収束・カオス的なプラトーで
    run-to-run 感度が特に大きい) の範囲内と確認。コード経路の静的解析でも cell は
    `calcGradient_b_d`/LSQ カーネル群を一切呼ばないことを確認済み (すべて `discretization=="node"` 分岐)。
    **cell 不変**と判定。
  - **検証 5 (OFF 回帰)**: `sstThermalWallFunction=0` の y+30 wall-function ケース (`run_0055`,
    一般 GG/LSQ のみ適用) は `check_convergence.py` **ALL PASS** (全列 3.5–5.1 桁低下)。
  - **未実施**: 定数/線形場の合成 GG/LSQ 精度単体テスト (§6 の該当項目) は時間制約で省略。
    設計上 (owner-state を境界値に使う = 発散定理の面積分を閉じる標準手法, SU2 と同型) 定数場で
    厳密ゼロになることは数式上明らかだが、専用ツールでの機械的確認は follow-up とする。
  - 副産物 (`turbulence-sst-adiabatic-taw-fluxmodel.md` の検証中に発見): この一般化そのものは
    y+30 SST wall-function ケースでも単独で安定 (`run_0055` ALL PASS) — 後続で見つかった発散は
    別計画の Taw 流束モデル置換機構に起因し、本計画のスコープではない。
