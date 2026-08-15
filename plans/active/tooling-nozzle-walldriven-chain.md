# ノズル設計 wall-driven チェーン (①風洞): 円弧廃止 + CFD 特性抽出 + wave-cancellation MOC

## メタ

- **area**: `tooling / optimization`
- **status**: `in_progress`  <!-- W0–W3 済・W4 は平面近似止まり。W5 以降は保留 (2026-08-15、主線は tooling-nozzle-axismach-chain.md へ) -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (現在仕様。「wall-driven チェーン」節)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md)
  - 兄弟: [`tooling-nozzle-phase3-windtunnel.md`](tooling-nozzle-phase3-windtunnel.md) (現行 B8 チェーン。**生産構成として維持し壊さない**)
  - 出典調査: [`notes/investigations/nozzle-throat-curvature-shape-representation-survey.md`](../../notes/investigations/nozzle-throat-curvature-shape-representation-survey.md)
  - 引き継ぎ: [`notes/sessions/2026-08-15-handover-wind-tunnel-nozzle-design.md`](../../notes/sessions/2026-08-15-handover-wind-tunnel-nozzle-design.md)
- **created**: `2026-08-15`
- **owner**: `sano`

## 1. 目的

ユーザ持ち込みの外部推奨方針 (ChatGPT 作成の「wall-driven 設計ツール試作仕様」、以下「原案」) を
引き継ぎ §2 の確定事実・§3 の文献結論に照らして評価し、採否と修正を確定して実装する。
完成時には「U→T→D 低自由度多項式壁 → CFD (遷音速+初期超音速) → 特性情報抽出 →
D 下流 wave-cancellation MOC → δ\* 補正 → 全域 B-spline → NS/RANS 検証」の
チェーンが現行 B8 チェーンと**並存**して動き、同一要求 (①風洞 $M_d=4$) で
出口一様性とこぶ振幅を直接比較できる状態を得る。

## 2. 原案の評価 (確定事実 §2 / 文献結論 §3 との照合)

### 2.1 判定サマリ: **骨格は採用** — 3 経路のうち (b)+(c) を同時に実装し、(a) を dv 化する構成

こぶ問題の確定帰結 (引き継ぎ §2) は「効く経路は (a) $R$ を上げる / (b) スロート下流円弧を
設計出力に置き換える (接合廃止) / (c) 遷音速アンカーの高次化、の 3 つだけ。J 下流の壁表現
変更は的外れ」。原案を分解すると:

| 原案の要素 | 3 経路との対応 | 判定 |
| --- | --- | --- |
| 凍結円弧を廃止し U→D を低自由度多項式で設計 ($R_t$ はスロート曲率**条件**として残る) | **(b)** の変種。超音速域内の「κ=const 円弧 → 設計壁」接合 J が消滅。T→D 多項式の κ は $1/R_t$→0 へ**連続に減衰** (κ' 跳びの接合が超音速域から消える) | 採用 |
| 遷音速+初期超音速を CFD で解き特性情報を抽出、設計モデルにしない | **(c) の完成形**。B6 (CFD アンカー化でこぶ 58% 減) の実測と整合。Sauer の「$M''$ が構造的に存在しない」問題 (調査 §9-2) が消滅 | 採用 |
| D 下流は wave cancellation で壁を**従属出力**として生成 | CONTUR の「壁は MOC の従属出力」思想 (調査 §2.1) そのもの。データ線が D からの C⁻ 特性線になるので「縦線 starting line は古典に対応物がない」問題 (調査 §9-5) も解消 | 採用 |
| $M_c(x)$ を設計入力から診断量へ降格、出口 M は $\theta_D$ の root find | B10 が暴いたジレンマ (始端 C2 整合 vs 出口一様性、自由 CP 3 個で両立不能) が**構造的に消滅** — 整合させるべき target 曲線が存在しない | 採用 |
| $R_t$ が明示的な設計変数 | **(a)**。outer loop / MOO で $R$ を直接動かせる | 採用 |

### 2.2 正直な限界 (原案が言っていないこと)

- **軸こぶが「消える」保証はない**。到達不能域 ($x<x_{reach}$) の軸 M は U→D の壁形状と
  $R_t$ だけで決まる物理であり、新チェーンでも下流壁からは届かない (確定事実は不変)。
  新チェーンの利点は (i) スロート域が作る波が**実 CFD 場として正確に**下流設計に入り、
  相殺壁が実場基準になるので試験部への漏れが構造的に減る、(ii) $R_t$ と T→D の κ 分布が
  直接 dv なので、こぶ振幅そのものを outer loop で下げに行ける、の 2 点。
- **仮説 (検証対象・2026-08-15 レビューで表現訂正)**: 新チェーンで κ' 跳びが完全に
  消えるわけではない。正確には「**超音速域の J にあった跳び (κ=const 円弧との接合) は
  消えるが、T (U→T / T→D の接合点 = スロート) の κ' 跳びは残る**」(実装は
  `kappa_prime_jump_at_throat` で常時定量)。仮説は「接合特異性が J (超音速域) から
  T (遷音速点) へ移り、かつ T→D 区間内では κ が連続減衰するため、同じ $R_t$ でも
  円弧+J 構成より波源が弱い」。ただし κ' 不連続 (G3 破れ) の効果を G2 破れと分離定量
  した文献は存在しない (調査 §2.4 末尾) ので、これは実測で判定する (Phase W9 の A/B)。
- **Korte の相場観** (調査 §2.5): 遷音速モデル差起因の 0.02 級うねりは古典チェーンの実測
  相場。新チェーンはまさにその誤差源を CFD 実場で置換するので、0.02 を下回れるかが
  成否の定量指標になる。

### 2.3 原案から修正する点

1. **新規ディレクトリ構成 (原案 §26) は不採用**。`design/forge_design` の既存構成に統合する
   (原案自身が「既存 repository に合わせて変更してよい」と許容)。再利用資産:
   `_hermite_quintic` (wall.py) / `moc_kernel` (単位プロセス) / `cminus_cfd` (C⁻ 追跡) /
   `cfd_anchor` (局所多項式フィット抽出) / `mesh2d`+`runner_wt` (メッシュ・node Euler 評価) /
   `deltastar` (δ\* 相関) / `metrics`。
2. **判別実験 W0 (現行チェーン R=5 A/B、調査 §8) を並走枠として残す**。新チェーンの既定
   $R_t$ と「こぶの R 依存 vs 構成依存の分離」という物理情報を最安 (config+1 パス) で得る。
   新チェーン実装と独立に実行可能。**実行するかはユーザ判断** (引き継ぎ §8 の指示通り)。
   run 番号は `run_0039` 以降。
3. **U→T は既存 Bell–Mehta 5 次収縮と同型** (両端 6 条件の quintic Hermite)。差分は
   「円弧に C2 接続」でなく「throat で直接 $r''=1/R_t$」を課すことだけ。原案 §7 の単調条件
   $\mu=L_U^2/(R_t\Delta r)\le20$ は検証済み (§4.1) だが大収縮比では緩い制約で、実用制約は
   $\max|d\kappa/ds|$ と入口境界層 (Bell–Mehta の趣旨) の側にある。
4. **原案 §10 (直管→一定傾斜) は geometry option として実装・保存** ($\lambda\in[5/3,5/2]$、
   $\lambda=2$ で 4 次退化)。ただし①風洞の主経路では使わない。
5. **実在気体 MOC の抽象化 (原案 §14) は初期実装では最小限**。`moc_kernel` は γ=const の
   PM 関数ベースなので、ThermoModel 分離は TP 拡張の必要が生じた時に行う (原案も許容)。
   γ を関数引数に保つ現行流儀は維持する。
6. **D 下流の CFD ドメイン延長壁は θ_D 直円錐**。D からの C⁻ より上流の場は延長壁に
   依存しない (超音速の依存域) ため、データ線抽出には十分。
7. **BL 補正の λ 指定より tolerance 制約優先 (原案 §21) に同意** (MVC 汎関数、調査 §3.2/3.4
   と整合)。physical exit slope=0 を課さない (§22) も同意。ただし Phase W6–W8 は後段。

### 2.4 原案の数式の独立検証 (2026-08-15 実施、実装前に完了)

乱数 200 試行 + 境界 ±0.1% すれすれ試験で数値検証済み。**同一条件 (200 試行・適応求積) を
`design/tests/run_walldriven_tests.py` に恒久化してある** (2026-08-15 レビュー P3 反映。
判定閾値はテストでは相対 1e-12 — 実測誤差は機械精度 ~3e-16 でこれを大きく下回る):

- U→T 勾配閉形式 $\frac{dr}{dx}=\frac{\Delta r}{2L_U}\xi^2(\xi-1)[5\mu\xi-3\mu-60\xi+60]$ は
  quintic Hermite の厳密解と一致 (相対 1e-9 以下)。単調収縮 ⇔ $\mu\le20$ は**必要十分**
  (μ=20×0.999 単調 / ×1.001 非単調)。
- T→D の $r_D=r_t+\frac{L_D^2}{12R_t}+\frac{L_D}{2}\tan\theta_D$ は適応求積 (`scipy quad`)
  と一致 (実測 3e-16)。$r''$ 閉形式一致。全区間 $r''\ge0$ ⇔ $L_D\le3R_t\tan\theta_D$ は
  **必要十分**。
- §10 の $\lambda$ 窓 $[5/3,5/2]$ と $\lambda=2$ の 4 次退化も検証済み。

## 3. スコープ

- **やる**: Phase W1–W9 (§6)。①軸対称風洞のみ。ideal gas MOC。既存 B8 チェーンとの比較評価。
- **やらない**: 実在気体 MOC (抽象化のみ準備) / ②–⑤機種への展開 / 現行 B8 チェーンの改廃
  (比較で新チェーンが勝つまで生産構成は B8) / ソルバ側変更。

## 4. 設計方針 (差分のみ — 数式は methods 側)

幾何の定義式・検証条件は [`methods/design/overview.md`](../../methods/design/overview.md) の
「wall-driven チェーン」節に記載。実装上の判断:

- **無次元系は既存踏襲** ($r_t=1$、スロート $x=0$、U は $x=-L_U$)。
- **`geometry/wall_walldriven.py`** に `UpstreamThroatPoly` (U→T) / `ThroatExpansionPoly`
  (T→D) / `WallDrivenThroatRegion` (合成 + validator + 診断) を新設。$r,r',r'',r'''$ は
  全て解析式で持ち、κ と $d\kappa/ds$ も解析評価 (数値微分しない)。
- **validator は「解析条件 + 数値サンプル検査」の二重化**: 解析条件 ($\mu\le20$、
  $L_D\le3R_t\tan\theta_D$、$r>0$) と、密サンプルでの単調性・変曲・$\max|d\kappa/ds|$ を
  両方返す (解析条件は表現を差し替えたら失効するため)。
- **エラー方針 (2026-08-15 レビュー P1 反映)**: **構築時 `ValueError` = 入力が不正**
  (非有限・非正の長さ/半径/曲率半径・$\theta\notin(0,\pi/2)$・出口半径非正) — 係数生成前に
  検査し、NaN 座標が W3 のメッシュ生成へ流れることを遮断する。**`validate()` = 幾何は
  定義できるが設計条件違反** (μ 超過・変曲・λ 窓外・$d\kappa/ds$ 上限超過)。この 2 層を
  全クラスで統一する。
- **$\max|d\kappa/ds|$ の上限 (レビュー P2 反映)**: W1 では `validate(dkds_max=...)` の
  **設定可能な上限**として実装する (既定 `None` = 診断のみ)。**数値既定値の決定は W1 の
  完了条件に含めない** — 上限は幾何単独では決められず CFD 感度 (壁圧・軸 M うねりとの
  相関) が要るため、W3 で感度を測り W5 の outer loop 制約として確定する。
- **T 点の κ' 跳びは診断値として常時出力する** (両区分とも $\kappa(T)=1/R_t$ で κ 連続だが
  κ' は一般に不連続。跳び量は幾何診断 CSV に載せ、遷音速域内の跳びとして許容可否を
  CFD で確認する — 調査 §2.4 の教訓により「滑らかだから無害」とは推定しない)。
- **特性抽出は node 場 + 局所多項式フィット** (cell の $M''$ 雑音支配・node 124× 改善の
  教訓を踏襲)。データ線上の $P_0$ 変動を等エントロピー仮定の妥当性診断として記録する。
- **W4 cancellation MOC の仕様 (2026-08-15 レビュー P1 反映 — 着手前確定分)**:
  - **公開 API**: 新モジュール `geometry/moc_cancel.py` に
    `extract_dataline_cfd(run_dir, ...) -> DataLine` (D の壁点から C⁻ を軸まで追跡し、
    線上の $(x, r, M, \theta)$ と $P_0$ を局所多項式フィットで返す。`cminus_cfd` の
    追跡器と `cfd_anchor` の抽出器を再利用) と
    `cancellation_contour(dataline, gamma, theta_tol, ...) -> 壁点列` を置く。
  - **数理と未知数**: 軸対称・非回転・等エントロピー MOC。単位過程は `KernelMOC` の
    interior (予測子修正子)・軸点 ($\theta=0$ + 軸対称源項) を再利用し、各新点の未知数は
    $(x, r, \theta, \nu)$。**march 順序は `moc_inverse` と同じ三角充填**: データ線の隣接
    点対から interior 過程でレベルを 1 段ずつ上げ、軸で反射させながら下流へ埋める。
  - **非反射 (cancellation) 条件 — 2 回目の訂正 (2026-08-15、外部レビュー
    [Codex] で確認)**: 1 回目の訂正で「$K^-_{\rm target}=\nu(M_d)$ を全壁点に
    定数として課す」方式に変えたが、**これは軸対称 wave-cancellation MOC の
    厳密解ではなく、平面 simple-wave を模した暫定的な壁設計則にすぎない**と
    判明した。理由: (i) 真に出口から逆走するなら、出口の $C^-$ を実際に引き、
    反対族 $C^+$ との交点で局所 $M,\theta$ を更新しながら**逐次**壁点まで
    到達させる必要があるが、現行実装は「出口の $J^-$ を位置によらず直接
    代入」しているだけで、真の逆 MOC トレースをしていない。(ii) 軸対称では
    $J^-$ は特性線上で保存されず $dJ^-=S^-(M,\theta,r)\,ds$ という源項を持つ
    ため、真に出口から逆走するなら
    $J^-_W=\nu(M_d)-\int_W^e S^-(M,\theta,r)\,ds$ が必要 (積分内が未知で
    単純代入では閉じない)。$J^-_W=\nu(M_d)$ 一定は、この積分をゼロと置く
    **平面流近似**である。(iii) 終端条件 (下記) も壁点自身の $M,\theta$ を
    見ているだけで、**出口断面全体の一様性は未検証**。
    **したがって W4 を「軸対称 wave-cancellation MOC 完成」として扱わない。**
    実装 (`_wall_cancel`) はこの近似としては数値的に正しく動作する (退化
    チェック機械精度・流線自己整合性 ~1e-8・格子収束 <0.5% — 詳細は
    `moc_cancel.py` docstring) が、モデルそのものの限界として明記する。
    **真に厳密化する経路 (未実装)**: (a) 出口終端特性線から源項込みで
    backward MOC し $D$ 側 forward MOC (`build_field`) と matching する
    (出口位置・断面・壁・$\theta_D$ を同時調整する自由境界問題)、(b) 仮の壁
    (現行近似解) で全場を実 CFD で解き、出口全断面の非一様性を目的関数として
    壁または $\theta_D$ を反復調整する (W5 の outer loop と統合可能)。
  - **終端条件 (要見直し)**: 現行は壁角 $|\theta_w|\le\theta_{tol}$ かつ壁点
    自身の $M$ が $M_d$ の要求帯に入った時点で打ち切る (θ=0 交差の安全策込み)。
    **これは壁の自己充足的な条件で、出口断面全体の一様性 $M(r)\approx M_d,
    \theta(r)\approx0$ を保証しない** — 経路 (b) で全場を解いた際に別途確認が
    必要。
  - **$P_0$ 診断の合否**: データ線上の $|\Delta P_0|/P_0 \le 0.1\%$ = PASS /
    0.1–0.5% = WARN / $>0.5\%$ = FAIL。実測 (run_0042/0043, `extract_dataline_cfd`
    の線形補間版): 0.18% (WARN)。**初回実装の最近傍探索版は 2.3% の偽値を
    出した (バグ、修正済み) — P0 診断は必ず補間ベースであることを確認する**。
  - **照合状況**: (i) 放射源流厳密解: `build_field` (内点充填) は個々の点が
    <0.1% で一致 (検証済み)。壁生成 (`_wall_cancel`) は退化ケース (平面単純波)
    で機械精度一致・流線自己整合性 ~1e-8・格子収束 <0.5% (検証済み、ただし
    「モデル近似としての」内的整合性であり厳密解一致ではない)。(ii) 格子収束:
    検証済み。(iii) 自己一貫 (B8 場との照合): **未実施**。実 CFD (run_0043,
    $M_D\approx2.26\to M_d=4$、$J^-$ ギャップ ~20°) では、目標が局所値から
    大きく離れているとはしごが発散するケースを確認 (iter5 で $r\approx-33$)。
    原因未特定 — モデル近似の限界 ($D$ 直後で一気に目標へ飛ぶ非物理な第一歩)
    が主因の可能性が高い。
  - 回転 MOC は $P_0$ 診断 FAIL 時のフォールバックとして温存 (本仕様では未実装)。

## 5. 関連 docs と前提

- 現在仕様: `methods/design/overview.md` (実装と同期して更新)
- 確定事実・撤回履歴: `tooling-nozzle-phase3-windtunnel.md` §9.2 (B5–B10)
- AGENTS 必須ルール (新連番 run / メッシュ VERDICT / NaN・収束・準定常チェック /
  README run 一覧同期 / 検証後 commit) を全 Phase に適用。

## 6. 実装ステップ

| Phase | 内容 | 状態 |
| --- | --- | --- |
| **W0** | 判別実験: 現行 B8 チェーンで R=5 A/B (調査 §8 の計画のまま)。こぶの R 依存を測る | **済 (2026-08-15, `run_0039_r5_cold`)** — **判定: 円弧起因が支配**。同一コード 3 次トレンドでこぶ振幅 R2 0.0217→R5 0.0106 (比 0.49)、到達不能域のみでは 0.0217→0.0085 (比 0.39)、いずれも $1/\sqrt R$ 比 0.632 より深い。壁 $dM/dx$ の J 跳びは比 0.65→0.98 でほぼ消滅。**帰結**: 調査 §8 に従い MOO/新チェーンの既定は $R_t\ge5$ 帯とし、接合廃止 (E) 単独への投資は保留 — 本 plan の wall-driven チェーンは (b)+(c) を含むため方針変更なし、既定 $R_t=5$ で W3 へ。R5 では山 (x=1.32) が $x_{reach}$=1.10 の下流 = 設計到達域内へ移った (R2 では到達不能) 点も新チェーンに有利。留保: R5 は単一格子。詳細は [case README](../../case/41.wind_tunnel_design/README.md) の run_0039 行と `r5_ab_metrics.json` | 
| **W1** | 幾何生成器: `wall_walldriven.py` (U→T / T→D / 合成 / validator / 解析 κ・dκ/ds / §10 option) | **済 (2026-08-15)** |
| **W2** | 解析テスト + 可視化: 端点条件・境界すれすれ・違反ケースの検出、$r,\theta,\kappa,d\kappa/ds$ 図 | **済 (2026-08-15)** |
| W3 | メッシュ/CFD 接続 | **済 (2026-08-15, `run_0042_walldriven_w3c` → `run_0043_conecheck_staged` で訂正確定)**。`WallDrivenCFDWall` (直管+U→T→D+θ_D円錐+任意の戻しテーパ/円筒) + `runner_walldriven` (新 type `wind_tunnel_axisym_walldriven`、**段階起動必須** — cold 単段 conv1+cfl4 は step~9 発散: run_0040/0041)。**訂正**: 当初「円錐のみだと node 既知問題 [傾斜壁∩超音速流出コーナー] で発散する」と誤診断し戻しテーパ+円筒を導入した (run_0042) が、円錐のみ+段階起動だけの切り分け (run_0043、ユーザ指示) で NaN なし完走・うねり指標も同水準と判明。**真因は cold 単段の CFL が高すぎただけ** (`procedures/divergence-and-startup.md` の標準的教訓どおり)。**生産既定を円錐のみ (`L_turn=L_cyl=0`) に簡素化**。品質 PASS / NaN なし / ALL STEADY / 設計+データ線域は maxΔM ~4.5e-6 に凍結。**軸うねり +0.0009 (テーパ有) / +0.00088 (円錐のみ) — 円弧 R=5 の 1/11 — §2.2 仮説を強く支持**。W4 のデータ抽出は run_0042/0043 どちらでも有効 (幾何は同一、下流のみ差) | 
| W4 | 特性情報抽出 + cancellation MOC | **暫定近似が動作 (2026-08-15)、軸対称厳密化は未完了 (外部レビューで発見)**。`geometry/moc_cancel.py` 新設。**検証済み**: `extract_dataline_cfd` (D 発 C⁻ + $P_0$ 診断。線形補間で実測 WARN 0.18%)、`build_field` (内点充填、放射源流厳密解と個々の点が <0.1% 一致、格子収束確認済み)、`_wall_cancel`/`march_wall_from_dataline` (退化チェック機械精度・流線自己整合性 ~1e-8・格子収束 <0.5%、θ=0 交差の振動回避策込み — 68 件 ALL PASS)。**モデルとしての限界 (外部レビュー [Codex] で確認、2026-08-15)**: 壁の閉じ方 $J^-_W=\nu(M_d)$ (全壁点一定) は「出口から真に逆走した結果」ではなく、局所 $M,\theta$ を更新しながら特性線を逐次追跡していない。**軸対称では厳密に誤り** — $J^-$ は $dJ^-=S^-(M,\theta,r)\,ds$ という源項を持ち保存量でないため、真の逆走には $J^-_W=\nu(M_d)-\int_W^e S^-ds$ が要る (積分内が未知で単純代入では閉じない)。現行は積分をゼロと置く平面流近似。終端条件も壁点自身の状態のみを見ており出口断面全体の一様性は未検証。実 CFD ($M_D\approx2.26\to M_d=4$、大きな $J^-$ ギャップ) では march が発散するケースも確認 (原因: 近似の非物理な第一歩の疑い)。**厳密化の 2 経路 (未実装)**: (a) 出口終端特性線から源項込み backward MOC + $D$ 側 forward MOC との matching (自由境界問題)、(b) 仮の壁で実 CFD を解き出口非一様性を目的関数に壁/$\theta_D$ を反復調整 (W5 outer loop と統合)。次段階: 経路 (b) から着手が現実的 (既存 W5 outer loop 構想と親和性が高い) |
| W5 | outer loop: $M_e(\theta_D)$ root find (必要なら $R_t, L_D$ も)。$M_c(x)$ ・こぶ振幅を診断出力 | 未着手 |
| W6 | δ\* 補正: `deltastar` 相関 + smoothing spline + 法線オフセット | 未着手 |
| W7 | 全域 B-spline fitting: tolerance 制約 ($\max d\le\varepsilon_{geom}$) 下で $\int(d\kappa/ds)^2ds$ 最小化。throat 等の hard constraint | 未着手 |
| W8 | STEP/CSV 出力 (親計画 5a と合流) | 未着手 |
| W9 | NS/RANS 検証 + **B8 チェーンとの同条件比較** (ε_M / ε_θ / こぶ振幅 / 同一 R での円弧 vs 多項式 A/B) | 未着手 |

## 7. 検証

- W1–W2: 解析 unit test (`design/tests/run_walldriven_tests.py`) — 端点条件 (1e-12)、
  閉形式一致、境界 ±0.1% の判別、違反ケースで非単調/変曲が実際に発生すること、
  解析 κ・dκ/ds と数値微分の一致。
- W3 以降: 各 Phase で AGENTS 必須ゲート (メッシュ VERDICT / `check_convergence.py` /
  `check_quasisteady.py`) を通し、run は `case/41.wind_tunnel_design` の README 表に登録。
- 最終判定 (W9): 同一要求で B8 比較 — 出口 ε_M/ε_θ が要求内、こぶ振幅が B8 (0.019–0.023)
  以下、Korte 相場 0.02 を下回るかを報告。

## 8. 影響範囲

- 追加のみ: `design/forge_design/geometry/wall_walldriven.py`、`design/tests/run_walldriven_tests.py`。
- 既存チェーン (B8 生産構成・`wall_modef`・`moc_inverse`・帰還ループ) は無変更。
- W3 で `mesh2d`/`runner_wt` に option 追加 (既存経路は不変に保つ)。

## 9. 完了条件

1. W9 の同条件比較で新チェーンの採否が定量判定されている (勝てば生産構成の移行 plan を
   別途起票、負ければ本 plan を archived へ移し敗因を記録)。
2. `methods/design/overview.md` が最終仕様と同期している。

## 変更ログ

- 2026-08-15: 起票。原案 (ChatGPT 仕様) を引き継ぎ §2/§3 と照合して採否確定 (§2)。
  原案の設計式を独立数値検証 (§2.4、全て正)。W1–W2 実装・テスト 34 件 PASS・可視化完了。
- 2026-08-15 (レビュー反映, commit 2913d85f への指摘 5 件): [P1] 構築時入力検査を全クラスに
  統一 (非有限・非正・角度域外・出口半径非正 → `ValueError`。validate() との 2 層方針を
  §4 に明文化)。[P1] W4 cancellation MOC の仕様 (公開 API・非反射条件・march 順序・終端・
  $P_0$ 閾値・受け入れ照合 3 件) を §4 に追記。[P2] §2.2 の κ' 記述を「J の跳びは消えるが
  T の跳びは残る」に訂正。[P2] `validate(dkds_max=...)` を実装し数値既定値の決定は W3/W5 へ
  明示的に後送。[P3] 恒久テストを §2.4 の記録と同一条件 (200 試行・適応求積) に強化。
  不正入力拒否 12 件を含めテスト 48 件 ALL PASS。
- 2026-08-15 (W0 実施, ユーザ指示): `case/41.wind_tunnel_design/run_0039_r5_cold` で
  R=5 単発 A/B。品質 PASS / 収束 ALL PASS (全列 4.9–5.2 桁) / ALL STEADY / NaN なし。
  **円弧起因支配と判定** (振幅比 0.49、到達不能域のみ 0.39 — ともに $1/\sqrt R$=0.632 より
  深い)。新チェーンの既定 $R_t$ は 5 とし W3 へ進む。§6 W0 行に要約、一次データは
  case README run_0039 行・`r5_ab_metrics.json`・`analyze_r5_ab.py`。
- 2026-08-15 (W3 実施, ユーザ指示): `WallDrivenCFDWall` / `runner_walldriven` /
  probdef 新 type を実装、テスト 68 件 ALL PASS。**段階起動必須** (soft 1次+cfl0.5
  → 本段。cold 単段 conv1+cfl4 は D 発膨張波の軸反射が壁へ戻る帯 x≈8.6 で準 1D
  IC の過渡が破綻し step~9 発散 — run_0040/0041)。当初これを node 既知問題
  「傾斜壁∩超音速流出コーナー」と誤診断し、出口で壁を軸平行へ戻す戻しテーパ+
  円筒を導入 (`run_0042_walldriven_w3c`)。D 発 C⁻ データ線・P0 診断 (WARN
  0.196%)・軸うねり +0.0009 (円弧 R=5 の 1/11、こぶ実質消滅) を確認。
- 2026-08-15 (テーパ要否の切り分け, ユーザ指示): 円錐のみ (テーパなし) +
  段階起動を単独試験 (`run_0043_conecheck_staged`) → **NaN なし・ALL STEADY で
  完走、うねりも同水準 (+0.00088)**。**戻しテーパは不要だったと確定** — 真因は
  cold 単段の CFL (`procedures/divergence-and-startup.md` の標準的教訓どおり)。
  `WallDrivenCFDWall` の既定を `L_turn=L_cyl=0` に変更し生産構成を簡素化。
  W4 (特性抽出 + cancellation MOC) へ。
- 2026-08-15 (W4 着手, ユーザ指示「自走してよい」): `geometry/moc_cancel.py` 新設。
  `extract_dataline_cfd` (`cminus_cfd.load_field` に `with_pressure` オプション追加)
  と `build_field` は放射源流厳密解照合・格子収束を確認し**検証済み**。壁生成
  単位過程 `_wall_cancel` は WebSearch で古典 MLN 型公式 ($J^-$ 直前壁点継承 +
  $J^+$ 隣接内点) を確認して実装したが、**3 通りの march 構成すべてが厳密解照合で
  有意に乖離し未検証のまま** (詳細は §6 W4 行・`moc_cancel.py` docstring)。
  正直な状態管理として: 壁生成関数に `warnings.warn` を追加、恒久テストは
  検証済み部分のみに限定 (8 件 ALL PASS)、production では壁生成を使わない旨を
  明記。W4 は完了させず、次セッションでの継続調査対象として引き継ぐ。
- 2026-08-15 (ユーザ指摘で W4 の閉じ方を再設計): ユーザが「放射源流厳密解に
  合わなくてよいと思った」と指摘 — 検証方針 (manufactured solution) 自体は
  正しいが、閉じるための第二の情報 (単一スカラー $M_d$) が欠けていたと判明。
  $J^-_W=\nu(M_d)$ (壁で目標出口一様状態の反射側不変量を課す) を実装。
  デバッグ過程で 2 件の実バグを発見・修正: (1) `_wall_cancel` の軸対称源項
  補正 `Sm` が `sin(0.0)` で恒等的にゼロになっていた、(2) `front[-2]` を
  `_wall_cancel` の A とはしご最初の段の両方に使う二重消費で march が 1 反復
  で固定点に収束し停止していた。さらに $\theta_{tol}$ が粗いと march が
  収束点を通り過ぎ振動域で破綻することを発見し、θ=0 交差検出による安全策を
  実装。平面単純波退化チェック・流線自己整合性・格子収束のテスト 16 件
  ALL PASS。
- 2026-08-15 (外部レビュー [Codex] で W4 の位置づけを再訂正): 上記の
  $J^-_W=\nu(M_d)$ 一定方式は「軸対称からの真の逆 MOC」ではなく、平面
  simple-wave を模した暫定近似にすぎないと判明 (§4 の非反射条件節に理由を
  記載: 特性線の逐次追跡をしていない・軸対称源項 $S^-$ を無視している)。
  実 CFD (run_0043, $J^-$ ギャップ ~20°) では march が発散するケースも確認。
  W4 を「完成」と扱わず、モデルの限界として明記。厳密化の 2 経路 (backward
  MOC matching / 実 CFD 反復調整) を記録し次段階へ引き継ぐ。
- 2026-08-15 (**W5 以降を保留**): ユーザ判断により軸中心 Mach 分布指定を主線とする
  方針に転換。後継は [`tooling-nozzle-axismach-chain.md`](tooling-nozzle-axismach-chain.md)
  (Hall 遷音速 + 5次 Hermite 軸 Mach law + 逆 MOC + CFD-in-the-loop)。本計画は
  W4 厳密化 2 経路の記録と W1–W3 成果物 (段階起動 runner・cone 検証・
  `WallDrivenCFDWall`)・W0 計測手法の保管場所として active に残す。
