# 次セッション用プロンプト: 非平衡凝縮 (4 モーメント方程式) の forge 実装

このファイルは「窒素の非平衡凝縮 (homogeneous nucleation + droplet growth, method of moments)」を
forge に実装するための**着手用プロンプト兼コンテキスト引き継ぎ**。次セッションの冒頭でこの
ファイルを読み、`docs`→`plan`→実装の順 (AGENTS.md 開発フロー) で進めること。

---

## ゴール

Lin, Cheng, Luo, Qin, *On nitrogen condensation in hypersonic nozzle flows*, Shock Waves 24:179–189
(2014) の物理モデルを forge に実装し、極超音速ノズルでの**窒素の非平衡凝縮**を計算できるようにする。
最終検証は case/34 の Arthur ノズルで、dry (凝縮なし) に対し**凝縮による静圧・静温上昇**が出ることを
示し、論文の傾向 (中心線静圧が等エントロピーより上振れ) を再現する。

## 何を追加するか — 「4 方程式追加」の中身

気相 (forge 既存の `ro, roUx, roUy, roUz, roe`) に加え、**液相を表す 4 本の輸送方程式**を足す。
論文の保存ベクトル `U = (ρ, ρu, ρv, ρE, ρg, ρQ2, ρQ1, ρQ0)` の後半 4 成分:

- `ρQ0` : 単位質量あたり液滴数密度
- `ρQ1`, `ρQ2` : サイズ分布のモーメント (∫ r^n f dr)
- `ρg` : 液相質量分率 (g = 4πρ_l Q3/3、実質 Q3 モーメント)

各方程式は**移流項 + 凝縮ソース (核生成 + 成長)** を持つ。気相の質量・運動量・エネルギー方程式の
ソースはゼロ (混合気体全体は保存) で、結合は**熱力学 (状態方程式・温度関係)** 経由:
`e = (1-g) e_v + g e_l`、潜熱 `L` が静温を上げる。

### 物理モデル (詳細は下記まとめ md に既に整理済み)

- **核生成**: 古典核生成理論 CNT (Becker–Döring) × **Iland et al. (2009) の窒素経験補正**
  `J = J_CNT · exp(A + B/T)`, `A=-55, B=4270 K`。臨界半径 `r* = 2σ/(ρ_l R T ln(p_v/p_sat))`。
  **核生成は気相温度 T で評価** (臨界核は自己加熱なしとみなす)。
- **成長**: 修正 Gyarmathy (Goodheart 2004)。`dr/dt = (1/ρ_l)·f_FS(Kn)·k R T²/L² · ln(p/p_sat(T_d))`。
  **成長は液滴温度 T_d で評価** (Hill モデル, エネルギーバランス式、Kelvin 効果 `p_d` 込み)。
- **液滴温度 T_d**: Hill のエネルギーバランス (純窒素で β_c 項を落とした陰的式) を反復で解く。
- **N2 物性**: 飽和蒸気圧・液密度・潜熱・表面張力のフィット (論文 Appendix 1)。

> 物理式・係数・T_d の解き方・「核生成=T / 成長=T_d」の区分け・Kantrowitz との違いは
> [`papers/on nitrogen condensation in hypersonic nozzle flows_summary.md`](../../papers/on%20nitrogen%20condensation%20in%20hypersonic%20nozzle%20flows_summary.md)
> に詳細整理済み。まずこれを読むこと。

### ローカル資料パス (リポジトリルートからの相対)

- 主論文 PDF: `papers/on nitrogen condensation in hypersonic nozzle flows.pdf`
  (Lin, Cheng, Luo, Qin, Shock Waves 24:179–189, 2014)
- 上記まとめ (日本語): `papers/on nitrogen condensation in hypersonic nozzle flows_summary.md`
- 検証ノズルの原典 (形状の一次出典): `papers/Arthur_pd_1952.pdf`
  (P.D. Arthur PhD thesis, Caltech/GALCIT, 1952。スキャン PDF=テキスト層なし、画像描画で読む)
- Ref.5 (製作詳細, 入手不可): Nagamatsu & Willmarth, GALCIT Memo No.6 = J. Appl. Phys. 23(10) 1089 (1952),
  DOI 10.1063/1.1701991 (AIP 有料、NTRS/DTIC・ローカルに無し)。スロート曲率はここにも無い公算大。

## 既存実装の最良テンプレート: RANS 2 方程式 (roK/roOmega)

forge は既に「移流される追加保存スカラー + stiff ソース項」を RANS (k-ω SST) で実装済み。
4 モーメントは**これと同じ骨格**で追加するのが筋。次のファイルを読んで作法を踏襲する:

- 保存変数名リスト: `variables.cpp:426`, `input/setInitial.hpp:427`, `main.cpp:1020` 付近
  (`roK`,`roOmega` をどこで宣言・確保・H2D コピーしているか)
- スカラー移流: `cuda_forge/ransTransport_d.cu`
- ソース項: `cuda_forge/ransSource_d.cu` (ここに核生成・成長ソースを並列で追加)
- 更新ループ: `cuda_forge/update_d.cu` / `update.cpp` (保存量更新に 4 モーメントを追加)
- 従属変数: `cuda_forge/dependentVariables_d.cu` (g→液相, T_d などの算出)
- 残差: `cuda_forge/residualMonitor_d.cu` (rms_roQ0 等を追加)
- 境界: `boundaryCond.cpp` / `cuda_forge/boundaryCond_d.cu` の inlet は液相モーメント=0
  (入口は dry)。多成分 TP 入口の `Y0..Yn` 動的登録 (`boundaryCond.cpp:131`) が参考になる。
- 熱力学: `thermalMethod` (0=熱量的完全気体, 2=NASA 多項式 TP)。二相 EOS `e=(1-g)e_v+g e_l` は
  新しい thermalMethod か既存への分岐追加で扱う。`dependentVariables` の T/P 逆算に潜熱項 gL を入れる。

## 数値解法の方針 (論文準拠)

- **fractional-step**: 「ソース無しの均質部 (既存 forge の対流)」と「凝縮ソース項を持つ非均質部」を分離。
  forge の implicit/explicit 時間積分にソース積分を挟む。stiff な核生成・成長は**部分陰 or 部分時間刻み**
  を検討 (核生成率 J は ln(p_v/p_sat) に指数的で極めて stiff)。
- **float 精度の注意**: 出口 T~27K, P~200Pa まで落ちる (case/34 で確認済み)。`r*`,`J` は
  指数・対数で桁が飛ぶので、ソース評価は double で行う (forge は混合精度可)。

## 設計上の重要論点 (実装前に方針を決めること)

着手前に plan で次の 3 点の方針を必ず決める。安易にスカラー追加するだけでは陰解法・Roe で破綻する。

### (A) 多成分凝縮への一般化 (H2O 等を見据える)

- **今は凝縮する分子を 1 種類しか想定していない** (N2)。だが将来 **H2O 凝縮やその他の凝縮**に
  使う可能性があるため、**最初から「凝縮種ごとに 4 モーメント + 物性セット」を持てる構造**で設計する。
  - モーメント変数は `roQ0_<sp>, roQ1_<sp>, roQ2_<sp>, rog_<sp>` のように**凝縮種インデックス付き**で
    確保し、種数 `nCondSpecies` でループする (RANS の roK/roOmega 固定 2 本とは違い、可変本数)。
  - 物性 (飽和圧 `p_sat(T)`, 液密度 `ρ_l(T)`, 潜熱 `L(T)`, 表面張力 `σ(T)`, 経験補正 A,B) を
    **種ごとのテーブル/構造体**にまとめ、N2 は最初の実装、H2O は係数差し替えで足せるようにする。
  - 二相 EOS `e = e_gas + Σ_sp g_sp (e_l,sp - e_v,sp)` も種で和を取る形に。気相は当面 1 成分 (キャリア)
    でよいが、将来「キャリアガス + 複数凝縮種」へ拡張できるよう既存の多成分 TP (`thermalMethod 2`,
    `nSpecies`, `Y_s`) との関係を整理する (凝縮種 = 気相成分の一部が液化、という対応付け)。
  - **当面の実装は nCondSpecies=1 (N2) で動かす**が、配列・ループ・物性 API は多種前提で書く。

### (B) flux Jacobian (対流) が変わる — Roe で強く効く / SLAU は要調査

- 保存ベクトルに 4×nCondSpecies 本のスカラーが増えるので、**対流フラックスの Jacobian (∂F/∂U) が
  拡大**する。モーメントは基本「気相速度で移流される受動スカラー」だが、**二相 EOS で圧力 p が g に
  依存する**ため、純粋な受動スカラーではなく**気相系へ弱く逆結合**する (p の固有構造が変わる)。
  - **Roe**: 完全な Roe 行列 (固有値分解) を使うと、追加変数で**固有ベクトル/固有値構造が変わり影響大**。
    モーメントを「移流のみ (固有値 = u_n) の付加波」として扱い、p の g 依存を Roe 平均にどう入れるかを
    決める必要がある。素朴に分けると保存性・上流性が崩れるので注意。
  - **SLAU**: 圧力流束と質量流束を分離する形式なので、追加スカラーは質量流束に乗せた upwind で
    扱える見込みだが、**g 依存の圧力項が SLAU の圧力フラックスにどう入るかは未確認**。まず SLAU で
    「モーメント=質量流束 upwind の受動スカラー、圧力は二相 EOS で評価」を実装して挙動を見る。
  - 実装初手は **Roe より SLAU 優先** (case/34 も SLAU で回している)。Roe 対応は後段で固有構造を
    検討してから。

### (C) 陰解法 (block-DPLUR / DPLUR) の方針 — 圧力微分とソース stiff 化

- 陰解法の対角/ブロック Jacobian は **∂(flux)/∂U** と **∂p/∂U** に依存する。二相化で:
  - **圧力の微分が変わる**: `p = p(ρ, ρe, ρu_i, g)` となり (e=(1-g)e_v+g e_l)、`∂p/∂(roe)`,
    `∂p/∂ro` に加え **`∂p/∂(rog)` が新たに出る**。既存の implicit (scalar/block DPLUR) の圧力微分
    (`solverConfig`/`timeIntegration` 周り、`cuda_forge/implicitCorrection_d.cu`,
    `timeIntegration_d.cu`, axisym 源 Jacobian の既存実装が参考) を**二相対応に拡張**する作戦が要る。
  - **核生成・成長ソースが極めて stiff**: `J ∝ exp(-ΔG*/kT)`, `dr/dt ∝ ln(p/p_sat(T_d))` は
    過飽和に指数的・対数的で、陽な積分では時間刻みが潰れる。**ソースの point-implicit 化**
    (`∂S/∂U` を対角に組み込む) が事実上必須。RANS の SST ソース damping
    (`run_*`/`ransSource` の point-implicit、flat_plate SST の commit 履歴) が**直接の手本**。
  - **結合度の選択**: ①モーメントを気相と**疎結合 (スカラー DPLUR で別個に陰)** + ソース point-implicit、
    ②気相+モーメントを**密結合 block** で解く、の二択。まず ① (実装が軽く安定) で立ち上げ、収束が
    悪ければ ② を検討、と plan に明記する。
  - **fractional-step との整合**: 「対流 (既存陰解法) → 凝縮ソース (point-implicit サブステップ)」の
    分離なら、対流側 Jacobian は圧力の g 依存だけ直せばよく、ソースの stiff 性はサブステップ側に閉じる。
    これを既定方針の候補とする。

> まとめ: **(B) 対流 Jacobian (特に Roe)** と **(C) 陰解法の圧力微分 + ソース point-implicit** は
> 「スカラーを足すだけ」では済まない中核課題。plan の段階で具体的な数式 (∂p/∂rog, ∂S/∂U) と
> 結合戦略 (疎結合+split を初手とする) を書き下してから実装に入ること。

## 検証計画

1. まず **dry 一致**: 凝縮ソースを切った状態 (g≡0) で case/34 run_0003 と一致することを確認
   (回帰防止)。
2. **凝縮 ON**: Arthur ノズルで貯気を下げる/上げて凝縮量を振り、中心線静圧が dry の等エントロピー
   線より**上振れ**することを示す (論文 Fig.11 / 本文の「凝縮で静圧上昇」)。
3. 可能なら論文の parametric study (T0 を振って出口 P/T/Ma・凝縮量) の傾向を再現。

## AGENTS.md 開発フロー (実装前に必須)

新規物理スキーム追加なので、実装前に:

1. `docs/condensation/theory.md` を新規作成 (物理: CNT+Iland, Gyarmathy/Goodheart, moment method,
   T_d, N2 物性)。**新エリア `condensation/` を作る** (docs/index.md 目次も同期)。
2. `docs/condensation/implementation.md` を新規作成 (forge への対応: 追加保存変数, ソース項配置,
   fractional-step, 二相 EOS, 境界・初期, RANS 実装の流用方針)。
3. `.github/plans/condensation-nonequilibrium.md` を `_template.md` から作成し、`related_docs` で
   上記 2 つをリンク。`.github/plans/README.md` 一覧にも追記。
4. 以上が揃ってから実装着手。

## この前のセッションで完了済みのこと (前提)

- **dry 膨張の再現が完了** = case/34.arthur_n2_nozzle:
  - 形状: Arthur 1952 の 2D source-flow ノズル (スロート0.254mm→出口25.4mm, 全角11°, A/A*≈100) を
    上半分・発散部のみ・双曲線スロートでメッシュ化 (`mesh/arthur_nozzle.geo`)。
  - `run_0002_slau_dry_cfl1` (1次) / `run_0003_slau_2nd` (2次) で収束。出口 Mach 1次6.75→2次6.87
    (等エントロピー6.93)。中心線・壁が等エントロピー線に一致。詳細は case/34 README。
  - 初期場 `arthur_n2` を `input/setInitial.hpp` に追加済み (M=1.05 throat 一様)。
  - 後処理 `case/34/postproc.py` (Mach vs A/A*, P/P0 vs x を等エントロピーと比較)。
- forge の BC 作法を確認済み: 超音速入口は `inlet_uniformVelocity` が `Umag>=sonic` で ρ,P,U 全量固定。
  入口は throat 直下に置き choking を回避 (収束部不要)。`probe.yaml` が必須 (空でよい)。
- native build: `solver_density_cuda/build-native/{forge,convertGmshToForge}`、
  再ビルドは `cmake --build build-native -j`。`.geo`→`.msh`: `gmsh -3 x.geo -o x.msh -format msh4`、
  変換は run ディレクトリで `convertGmshToForge x.msh x.h5`。

## 着手プロンプト (次セッションにそのまま投げる想定)

> forge に窒素の非平衡凝縮 (4 モーメント方程式: ρQ0,ρQ1,ρQ2,ρg) を実装したい。まず
> `.github/plans/condensation-nonequilibrium-session-prompt.md` と
> `papers/on nitrogen condensation in hypersonic nozzle flows_summary.md` を読み、AGENTS.md の
> 開発フローに従って (1) docs/condensation/theory.md, (2) docs/condensation/implementation.md,
> (3) .github/plans/condensation-nonequilibrium.md を作成してから実装に入って。実装は既存の RANS
> 2 方程式 (roK/roOmega) の追加保存スカラー+ソースの骨格を流用する方針で。ただし「設計上の重要論点」
> (A) 多成分凝縮 (H2O 等) への一般化 = 凝縮種ごとに 4 モーメント+物性を持てる構造 (当面 N2 1 種で動かす)、
> (B) 対流 flux Jacobian の変化 (Roe で強く効く・SLAU は要調査、SLAU 優先)、(C) 陰解法の圧力微分
> ∂p/∂rog とソース point-implicit (疎結合+fractional-step を初手) — の 3 点の方針を plan で先に
> 数式込みで決めてから着手すること。検証は case/34 の Arthur ノズルで dry 一致 → 凝縮 ON で静圧
> 上振れ、の順。
