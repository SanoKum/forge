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

### モデルは凝縮種ごとに切替可能にする (必須要件)

核生成・成長・物性は**固定実装にせず、凝縮種ごとに config で選べるプラグイン構造**にする。
予定している組み合わせ:

| 種 | 核生成モデル | 成長モデル | 表面張力 | 備考 |
| --- | --- | --- | --- | --- |
| N2 | CNT × Iland 経験補正 (A,B) | 修正 Gyarmathy (Goodheart) | N2 σ(T) フィット | 上記 Lin 2014 |
| H2O | **CNT + Kantrowitz 非等温補正** | **Hertz–Knudsen モデル** | H2O σ(T) | Wyslouzil 検証 (下記) |

- **凝縮種ごとに以下を独立に切替**:① 核生成モデル (`CNT` / `CNT_Iland` / `CNT_Kantrowitz` / …) と
  その**係数** (補正 A,B 等)、② 成長モデル (`gyarmathy` / `goodheart` / `hertz_knudsen` / …)、
  ③ 表面張力モデル/係数、④ 飽和圧・液密度・潜熱のフィット。
- 実装イメージ: 各モデルを**enum + 関数ポインタ/switch (device 側)** で選び、係数は種ごとの
  物性構造体に持たせる。`bcondConfig`/`solverConfig` か専用 `condensationConfig.yaml` で
  `species: [{name: N2, nucleation: CNT_Iland, growth: goodheart, ...}, {name: H2O, nucleation: CNT_Kantrowitz, growth: hertz_knudsen, ...}]` のように記述。
- **Kantrowitz 非等温補正**: 生成中クラスタの自己加熱で J を割り引く補正 (J_noniso = J_iso/(1+b·q²) 系)。
  N2 では使わず H2O で使う、という**種ごとの選択**をそのまま表現できること。
- まとめ md の「Kantrowitz との違い」節も参照 (N2 は採用せず成長側 T_d で自己加熱を扱う方式)。

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
  - **モデル切替を種ごとに**: 核生成/成長/表面張力/飽和圧などを enum+switch (device 側) で選び、
    係数は種ごと構造体。N2=CNT_Iland+Goodheart、H2O=CNT_Kantrowitz+Hertz–Knudsen を**同じ枠組みで**
    表現できること (上の「モデルは凝縮種ごとに切替可能にする」表を参照)。新モデル追加は enum と
    関数を足すだけで済む形にする。

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

### (D) float 精度 — 液滴数密度が巨大化する問題と打ち手

- **懸念は妥当**: 核生成率 J は ~1e20–1e30 個/(m³·s) になり、数密度 n=ρQ0 が ~1e18–1e22 個/m³ に達する。
  一方 r は nm スケール (~1e-9 m) なので、モーメントは
  `Q0~1e20`, `Q3~Q0·<r³>~1e20·1e-27~1e-7` と**桁が ~27 オーダー跨ぐ**。float32 (有効 ~7 桁) では
  Q0 と Q3 を同じスケールで保持できず、平均半径 `r=(Q3/Q0)^{1/3}` や「巨大な累積数に微小増分を足す」
  操作で**精度が崩壊**する。これは magnitude (1e20 は 3.4e38 内) ではなく**有効桁数の問題**。
- **打ち手 (plan で方針決定)**:
  1. **モーメントと凝縮ソース評価を double で持つ** = 最有力。forge は既に
     `build-double` / `build-mixed` / `build-mixed2` を持ち混合精度運用の実績あり
     ([[mixed-precision-axisym-refuted]] 系)。**気相対流は float のまま、凝縮の 4 モーメントと
     核生成/成長の積分だけ double** にする混合精度が現実解。`flow_float` とは別に
     `cond_float=double` 型を導入する案。
  2. **無次元化/リスケール**: モーメントを基準量で割って O(1) 化する。
     `μ_n = Q_n/(N_ref · r_ref^n)` (例 N_ref=1e18, r_ref=1e-9 m=1nm)。これでクロスモーメント比
     (平均半径) の計算が安定。double と併用するとさらに堅い。
  3. **質量分率 g は別扱いでよい**: g∈[0,~0.1] は float で十分。エネルギー結合は g が担うので、
     **g は通常精度、数密度モーメント Q0/Q1/Q2 だけ double+無次元化**という分担が効率的。
  4. **累積の round-off 対策**: 核生成の時間積分 (巨大数への微小加算) は Kahan 加算や、
     「新規核生成分」と「既存成長分」を分けて扱うと round-off を抑えられる。
  5. r=(Q3/Q0)^{1/3} 等のモーメント間演算は必ず double で。
- 結論として「float では太刀打ちできない」→ **モーメントは double (混合精度) + 無次元化**で対処可能。
  plan に「cond_float=double, μ_n 無次元化, g は通常精度」を既定方針として明記する。

> まとめ: **(B) 対流 Jacobian (特に Roe)**、**(C) 陰解法の圧力微分 + ソース point-implicit**、
> **(D) モーメントの double+無次元化** は「スカラーを足すだけ」では済まない中核課題。plan の段階で
> 具体的な数式 (∂p/∂rog, ∂S/∂U) と戦略 (疎結合+split / cond_float=double / μ_n 無次元化) を
> 書き下してから実装に入ること。

## 検証計画

1. まず **dry 一致**: 凝縮ソースを切った状態 (g≡0) で case/34 run_0003 と一致することを確認
   (回帰防止)。
2. **N2 凝縮 ON**: Arthur ノズル (case/34) で貯気を下げる/上げて凝縮量を振り、中心線静圧が dry の
   等エントロピー線より**上振れ**することを示す (論文 Fig.11 / 本文の「凝縮で静圧上昇」)。
3. **H2O 凝縮の検証ケースを追加** (要対応):
   - 既存の **case/16.nozzle_wys** (Wyslouzil 2D ノズル) を水蒸気凝縮の検証に使う。検証データ論文は
     **`papers/condensation/wyslouzil2000.pdf`**。まずこの論文を読み、ノズル形状・キャリアガス・
     水蒸気モル分率・貯気条件・測定量 (オンセット位置, 凝縮による圧力/温度, 液滴サイズ/数密度) を
     summary 化する (`papers/condensation/` に md でまとめる)。
   - case/16 に水蒸気凝縮の `run_*` を追加し、**H2O = CNT+Kantrowitz + Hertz–Knudsen** モデルで計算、
     Wyslouzil の実測 (圧力トレース・液滴サイズ等) と比較する。
   - これにより「**N2 (Arthur) と H2O (Wyslouzil) の 2 種を同じ枠組みで切替**」できることを実証する。
4. 可能なら Lin 2014 の parametric study (T0 を振って出口 P/T/Ma・凝縮量) の傾向も再現。

## AGENTS.md 開発フロー (実装前に必須)

新規物理スキーム追加なので、実装前に:

1. `methods/condensation.md` を新規作成 (物理: CNT+Iland, Gyarmathy/Goodheart, moment method,
   T_d, N2 物性)。**新エリア `condensation/` を作る** (methods/index.md 目次も同期)。
2. `methods/condensation.md` を新規作成 (forge への対応: 追加保存変数, ソース項配置,
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
> 開発フローに従って (1) methods/condensation.md, (2) methods/condensation.md,
> (3) .github/plans/condensation-nonequilibrium.md を作成してから実装に入って。実装は既存の RANS
> 2 方程式 (roK/roOmega) の追加保存スカラー+ソースの骨格を流用する方針で。ただし「設計上の重要論点」
> (A) 多成分凝縮への一般化 = 凝縮種ごとに 4 モーメント+物性を持ち、**核生成/成長/表面張力モデルと
> 係数を種ごとに config 切替**できる構造 (N2=CNT_Iland+Goodheart, H2O=CNT+Kantrowitz+Hertz–Knudsen を
> 同枠組みで; 当面 N2 1 種で動かす)、(B) 対流 flux Jacobian の変化 (Roe で強く効く・SLAU は要調査、
> SLAU 優先)、(C) 陰解法の圧力微分 ∂p/∂rog とソース point-implicit (疎結合+fractional-step を初手)、
> (D) 液滴数密度の巨大化に対し**モーメントは double(混合精度)+無次元化**、の 4 点の方針を plan で
> 先に数式込みで決めてから着手すること。検証は case/34 Arthur で dry 一致 → N2 凝縮 ON で静圧上振れ、
> さらに **H2O は case/16.nozzle_wys + `papers/condensation/wyslouzil2000.pdf`** で水蒸気凝縮を検証。
