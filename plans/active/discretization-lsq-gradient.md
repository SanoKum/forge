# node-centered 最小二乗 (LSQ) 勾配

## メタ

- **area**: `gradient / architecture`
- **status**: `in_progress`  <!-- 実装済・GPU 検証未 -->
- **related_docs**:
  - `methods/discretization.md` (§7.3 LSQ 勾配)
  - `methods/gradient.md` (既知の制約: LSQ は node 限定)
- **related_plans**:
  - `discretization-median-dual.md` (親: cell/node 両対応化)
  - `discretization-node-boundary-ghostless.md` (兄弟: node 境界弱形式)
- **created**: `2026-06-14`
- **owner**: `CFD Dev`

## 0. 次の一手 = 動機の再検証ゲート (2026-07-19 ユーザー承認・必読)

**修正候補 (§9 の ①double 蓄積 / ②QR・列スケーリング / ③GG フォールバック) に着手しないこと。**
LSQ は近壁で発散する負結果 (§9) で停止中であり、かつ導入動機だった「GG の近壁 checkerboard」は
本ケース (case/29 node viscous) で **GG が問題なく収束** (rms_ro 5e-8)・`nodeMidpointFx` (fx=0.5) で
checkerboard roughness −36% 低減済み、と動機自体が揺らいでいる。

- **ゲート (安価・これだけやる)**: 固定した GG 収束場の上で勾配 checkerboard 指標 (近壁 `dUxdy` の
  隣接符号反転率・roughness) を測り、**せん断応力・壁面量に実害が出ているケースを特定する**。
  対象は case/29 node viscous と、動機の出所とされる平板系 (case/26)。
- **実害が確認できなければ本 plan を `archived/` へ移す** (動機消滅・LSQ 実装は gated のまま残置)。
  実害があった場合のみ、そのケースを検証ケースに据えて修正候補①〜③を再評価する。

## 1. 目的

node-centered (median-dual) の近壁で Green-Gauss (GG) 面勾配が持つ checkerboard (奇偶デカップリング)
モードが薄い粘性境界層のせん断応力を振動させる問題に対し、重み付き最小二乗 (LSQ) 勾配を node モード限定の
opt-in (`gradLSQ`) として導入する。完了時、node 粘性ケース (例 case/29) で BL の勾配振動が解消し、
cell モードおよび `gradLSQ=0` の node モードは従来とビット不変である状態を得る。

## 2. スコープ

- **やる**:
  - 近傍 CV 中心の逆距離二乗重み付き LSQ 勾配を node モードで計算する GPU カーネル群。
  - 境界は bvar 境界値 (壁=速度0 等) を半割面重心の LSQ 点として加え勾配を閉じる。
  - `mesh.gradLSQ` フラグ (既定 0) で GG / LSQ を切替。`discretization=="node"` のみ有効。
- **やらない**:
  - cell モードへの LSQ 適用 (cell は GG 維持)。
  - limiter / 高次再構成。
  - GG 経路そのものの変更 (LSQ は別経路で早期 return)。

## 3. 関連 docs と前提

- 理論: GG 勾配は [`methods/gradient.md`](../../methods/gradient.md#理論)。LSQ の動機・式・実装は
  [`methods/discretization.md`](../../methods/discretization.md#実装) §7.3。
- 前提: node 弱形式境界 (§7.2, commit 97e190a 系列) と bvar 境界値が利用可能であること。
- 親計画 `discretization-median-dual.md` の node-centered パイプラインに乗る。

## 4. 設計方針

セル $i$ の勾配 $\mathbf g$ を正規方程式 $M\mathbf g=\mathbf b$ で解く。
$M=\sum_j w_{ij}\mathbf d_{ij}\mathbf d_{ij}^{\mathsf T}$、$\mathbf b=\sum_j w_{ij}\mathbf d_{ij}(\phi_j-\phi_i)$、
$w_{ij}=1/|\mathbf d_{ij}|^2$、$\mathbf d_{ij}=\mathbf x_j-\mathbf x_i$。

- $\mathbf b$ は既存の勾配配列を流用 (accum で RHS を貯め、solve で in-place に勾配へ上書き)。
- $M$ (対称 3×3) は static scratch (`Mxx..Mzz`、`nCells` 変化時のみ再確保)。
- 2D (1 セル厚) は $m_{zz}\le\text{tol}$ で検出し $xy$ の 2×2 を解き $g_z=0$。
- 組み立て・solve は倍精度。LSQ 経路は wrapper で早期 `return` し GG 正規化・弱形式境界加算をスキップ。

## 5. 実装ステップ

1. `solverConfig.hpp` / `.cpp`: `mesh.gradLSQ` (既定 0) を追加。— **実装済**。
2. `calcGradient_d.cu`: `lsq_solve_sym3` / `lsqGrad_accumInternal_d` / `lsqGrad_accumBoundary_d` /
   `lsqGrad_solve_d` を追加し、`calcGradient_d_wrapper` で `gradLSQ && discretization=="node"` 分岐。— **実装済**。
3. docs §7.3 と gradient/theory.md の制約更新。— **実装済**。
4. ビルド (コンパイル) 確認。— 未 (C セッションの GPU 稼働終了後)。
5. 検証 run (§6)。— 未。

## 6. 検証

- **単体 / ビルド**: native もしくは Docker で `tools/build.sh` がコンパイル成功すること。
  cell モードの既存ケースが `gradLSQ` 追加で**ビット不変**であること (既定 0 経路を踏む)。
- **検証ケース**: node 粘性で BL 勾配振動が出ていたケース (第一候補 `case/29.bell_vs_conical` の
  node viscous、補助で平板系)。新規 `run_*` を複製作成し `gradLSQ=1` で実行。
- **判定基準**:
  - `tools/check_convergence.py` の VERDICT を根拠に収束/未収束を判定 (全残差列)。
  - GG (`gradLSQ=0`) との比較で近壁せん断/速度勾配の checkerboard 振動が低減すること。
  - 壁法則・物理量 (P≤Pt, ro>0, T>0) が保たれること。
- **未確認リスク** (検証で潰す):
  - ghost セル (`ic>=nCells`) の勾配を LSQ が更新しない点が node 弱形式で問題ないか。
  - 2D z 境界面が $m_{zz}$ に混入し誤って 3D 解に落ちないか。

## 7. 影響範囲

- 触るモジュール: `calcGradient_d.cu`、`input/solverConfig.{hpp,cpp}`。
- docs: `methods/discretization.md` §7.3、`methods/gradient.md`、`methods/index.md`。
- 既存実行手順: `gradLSQ` 未指定時は無影響 (既定 GG)。

## 8. 完了条件

- [x] 関連 `methods/gradient.md` 更新済み
- [x] 関連 `methods/discretization.md` (§7.3) 更新済み
- [ ] 実装・検証完了 (本 plan の §6 を満たす) — **GPU 検証が未**
- [ ] `plans/README.md` の状態を `done` に更新
- [ ] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `2026-06-14` — 初稿。中断セッションの未コミット実装 (LSQ カーネル群 + `gradLSQ` フラグ) を
  正規化。docs §7.3 / gradient 制約を整備。コード実装済・GPU 検証は後続。
- `2026-06-14` — **GPU 検証 → 負の結果 (LSQ は現状発散)**。
  - **手順**: GG 安定ベース `run_dual_visc_conical_node_expl_long` (node, axisym, explicit RK3,
    cfl=0.1, SLAU, laminar) を複製し `gradLSQ` のみ変えた 2 run を 40000 step で実行
    (build: `solver_density_cuda/build-lsq/`、binary はフルソース文脈でコンパイル成功 exit 0)。
    - `case/29.bell_vs_conical/run_lsq_gg` (`gradLSQ=0`): **収束** (step39999 rms_ro=5.24e-8)。
      step0 rms_ro が元 GG run とほぼ一致 = GG 経路の回帰 OK。
    - `case/29.bell_vs_conical/run_lsq_on` (`gradLSQ=1`): **発散**。運動量 (rms_roUx/roUy) が
      step~150 から指数増大 (step200 roUx=7.0 → step400 1e5 → step1114 inf/nan)。
    - `case/29.bell_vs_conical/run_lsq_diag` (`gradLSQ=1`, 120 step・30 step 出力): step120 の場で
      異常勾配を局在化。
  - **局在化** (run_lsq_diag/res_120.h5): 異常勾配は**近壁第一層ノードに集中** (top0.1% の
    `|dUxdy|` ノードは 100% が `wall_dist<1e-4`、軸近傍ではない)。`max|dUxdy|=1.0e9` vs
    median ~7e2 = 壁で約 6 桁のスパイク → 粘性応力爆発 → 運動量発散。
  - **根因 (精緻化, データ裏付け)**: 「float だから一様に悪い」のではない。GG (`gradLSQ=0`) は
    **行列を解かず面の重み付き総和**なので float でも有界・安定 (近壁 roughness median≈1・max
    `|dUxdy|`=1.2e8 とそれなりにラフだが破綻しない)。LSQ は正規行列 `M` を**逆行列で解く**ため、
    **少数の近壁ノードで M が退化 (near-singular)** し、その解だけが破滅的外れ値になる:
    発散前 step30 で大半のノードは GG と同等以上に平滑だが、近壁帯の roughness 99pct=3637・
    `max|dUxdy|`=1.2e9 と**一握りのノードが 1e9 級**に飛び、これが種となり運動量発散。
    退化の主因は近壁の高アスペクト比 (法線~1e-5 / 流れ方向~1e-3、比~100) + 歪み + ノード当たり
    近傍数不足で `det=mxx·myy−mxy²` が桁落ち (`M` を float `flow_float` 蓄積も悪化要因。`lsq_solve_sym3`
    は double 演算だが入力 `M` が既に float 精度)。
  - **次の修正候補** (未実施): ① `M`・`b` を double scratch で蓄積し double で solve、
    ② 逆距離二乗重みを近壁高アスペクトでロバストな重み/列スケーリング (正規化正規方程式 or QR) へ、
    ③ 退化 (det≈0) 検出時に GG へフォールバック。
  - **前提の再検討事項**: 本ケースでは **GG が問題なく収束** (rms_ro 5e-8)。LSQ 導入の動機だった
    「近壁 checkerboard」が本ケースで実害として出ているか自体を、固定 GG 場上での勾配 checkerboard
    指標で再確認すべき (動機の出所ケースの特定)。
- `2026-07-22` — **式レベルの補強証拠** (`tools/verify_linear_recon.py`, §0 ゲートとは独立の机上検証):
  LSQ 正規方程式は **double では線形場を機械精度で厳密再現** (AR=100+ジッタでも 1.5e-12)、
  一方 **float32 格納 (現行 `flow_float` 配列) では ~1e-6 に 6 桁劣化**。GG は非一様メッシュで
  線形場非厳密 (30% ジッタで相対誤差 66%)。→ §9 の近壁発散は LSQ の数学でなく格納精度+近特異幾何の
  問題で、修正候補① (double scratch 蓄積) の妥当性を支持。着手判断は §0 ゲートに従う (変更なし)。
- `2026-07-22` — **実メッシュ診断 (`tools/check_lsq_gradient.py` 新設) — 上記エントリの重心を訂正**。
  実 h5 メッシュから近傍セット (内部双対面 + 非 periodic 境界半割面 pc) を再構成し、
  ①退化センサス λ_min/λ_max(M̂=Σd̂d̂ᵀ) と ②線形場勾配誤差を (a) 現行相当 float32 格納 /
  (b) 全 double / (c) 係数事前計算 (setup double solve → float32 係数適用) / (d) GG で比較。
  - **case/39 hill 3D production** (`case/39.periodic_hills/run_0013_prod_ddes/hill_xc_160x100x60.h5`,
    99.2 万ノード): 退化ゼロ (worst 比 0.20)。LSQ は (a)=(b)=(c) 完全一致で max 5e-3
    (float32 場の差分ノイズ床、壁ノード)。**GG は同じ線形場で near-wall max 59% / p99.9 46%**
    (fx 補間の整合性誤差) → このクラスのメッシュでは **LSQ 化で近壁勾配 ~2 桁改善**の余地
    (勾配は粘性応力・keepDissJump・SGS の共通入力)。現行 float 実装でも数値的には安全な条件。
  - **case/29 nozzle 軸対称** (`case/29.bell_vs_conical/run_lsq_on/nozzle.h5`, §9 発散の本家):
    **真の面内退化が実在**: 比<1e-2 が 2.58% (709 ノード)、worst 4.3e-5、全て近壁 (wd~1e-4)
    = 法線エッジが接線から ~0.5° しか立っていない高AR+シア壁。**(a)=(b)=(c) が同一** (max 10%,
    壁) — つまり **float32 格納は発散の主因ではなく** (前エントリの強調を訂正)、double でも
    線形場で 10% 誤る「情報が無い」退化。§9 の発散は退化ノードの勾配誤差×粘性フィードバックが
    有力で、**必須修正は候補③ (退化検出→フォールバック) ないし stencil 拡張/壁 virtual 点**。
    候補① (double 化) 単独では不足。(d) GG は軸対称 planar 幾何が h5 に無く N/A。
  - 示唆: 修正の本命は「**係数事前計算** (setup 時 double solve + 退化ノードの GG/縮退方向
    フォールバック焼き込み) + ランタイム float32 gather」。実行時 double 不要・現行 LSQ より速く、
    退化処理を静的に済ませられる。§0 ゲートの動機再確認には hill の GG 46-59% が新証拠。
