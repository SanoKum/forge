# [引き継ぎ] ①軸対称風洞ノズル設計ツール — 現況と次の一手

次セッションへの引き継ぎ文書。**ユーザは ChatGPT 推奨の方針を持ち込む予定**なので、
その方針を評価・実行できるだけの前提をここに集約する。

対象: `case/41.wind_tunnel_design`、$M_d=4$ / $P_t=1$ MPa / $T_t=800$ K / $r_t=10$ mm /
円弧 $R_u=R_d=2r_t$。座標は $r_t$ で無次元化。

---

## 1. 30 秒で分かる現況

- **設計チェーンは一通り動く**。目標軸マッハ $M_c(x)$ → 逆 MOC → 壁 → メッシュ →
  forge (Euler 軸対称) → 軸 $M$ 実測 → 特性線マップ経由で壁へ帰還 → 反復。
- **正式収束を達成済み** (`run_0026_b8`: masked $\|\Delta M\|_\infty$ = 0.0198 = 0.49% $M_d$ ≤ tol 0.5%)。
- **出口性能は風洞要求を大きく満たす**: $\varepsilon_M$ 0.023–0.032% (要求 <0.5%)、
  $|\theta|_{max}$ 0.015–0.018° (要求 <0.1–0.25°)。
- **CFD は node (median-dual) がベース** (2026-08-17 ユーザ決定)。cell は回帰対照のみ。
- **未解決の 1 点**: スロート円弧–設計壁の接合部が作る**軸マッハのこぶ**
  ($x/r_t\approx1.43$、振幅 0.019–0.023)。**設計壁のどの自由度でも届かない位置にある**。
- **文献調査は完了済**
  ([`notes/investigations/nozzle-throat-curvature-shape-representation-survey.md`](../investigations/nozzle-throat-curvature-shape-representation-survey.md))。
  推奨は **$R=5$ の単発 A/B を判別実験として先行**。

---

## 2. こぶ問題の確定事実 (ここは動かない)

| 事実 | 根拠 |
| --- | --- |
| こぶは $x/r_t\approx1.43$、振幅 0.019–0.023、谷は $x\approx1.86$ | 意図非依存指標 (軸 $M$ の 5 次トレンド差) |
| **格子収束する** (1x/2x/4x で 0.0218→0.0222→0.0227) が**鋭くなる** (最大勾配 +19%) | 接合部局所細分の 3 格子 |
| **cell/node 両方に残る** (node 0.0187 @1.44) → 数値雑音だけではない | run_0031–0036 |
| **波源は凍結円弧上** ($x_w=0.154 < x_j=0.2465$) | 実 CFD 場の C⁻ 後退トレース |
| $x<x_{reach,CFD}=1.704$ の軸 $M$ は**壁を変えても 4 桁目まで不変** | 3 種類の壁で実測 |
| 遷音速 Cauchy データを Sauer→CFD 実測に替えると**こぶ 58% 減** (0.118→0.050) | B6 (run_0018) |

**帰結**: こぶを消すには (a) $R$ を上げる、(b) スロート下流の円弧自体を設計出力に
置き換える (CONTUR 型・接合廃止)、(c) 遷音速アンカーの高次化 — の 3 経路しかない。
**J 下流の壁表現をどう変えても効かない** (G3 化・CST 化・κ(s) 化を J 下流だけに
適用するのは的外れ)。

---

## 3. 文献調査の要点 (次の判断に直結)

1. **Sivells/CONTUR は壁側で曲率を接合しない**。軸中心線分布を遷音速解・radial flow・
   出口一様条件に 1–2 階微分まで整合させ、壁は MOC の**従属出力**として一体生成する。
   本ケースの「円弧 + J で C2 接合」は CONTUR が批判した Foelsch 型の系譜に近い。
2. **接合部が軸に波を作る現象は 1963 年から文献化されている** (Darwell & Badham、
   Back & Cuffel 1966、Cuffel 1969)。軸対称では壁擾乱が軸に集束する。
3. **$R=2$ は遷音速理論の適用下限ちょうど** (Cuffel: "not less than about 2")。
   CONTUR の実設計例は $R=5.5$–$6$。
4. **Sauer の $M''$ は「精度が悪い」のではなく「存在しない」** — 軸上 $u$ が厳密に
   $x$ の 1 次なので $M''$ に幾何情報が構造的に無い。
5. **Korte の実在気体 Sivells 設計でも CFD 検証で ~0.02 Mach の残留変動**があり、
   「設計 MOC の遷音速モデルと実流れの差」に帰属されている。本ケースの 0.019–0.023 は
   **古典チェーンの実測相場と同水準**。ただし「不可避」の意味ではない。

**推奨実験** (調査 §8): **$R=5$ の単発 A/B**。表現・チェーン・dv 構成を一切変えず、
こぶ振幅の $R$ 依存を測る最安の分岐点。判定基準は「振幅が ~1/2 以下 ($\lesssim0.010$)
なら円弧起因が支配 → MOO に $R\ge5$ を課す」「ほぼ不変 ($\ge0.015$) なら構成依存が支配
→ 接合廃止 (CONTUR 型) を計画化」。

---

## 4. B10 (直近の試み) の結末 — 機構は正しいが自由度が足りない

**B10 = target と自由 dv の起点を $x_0$ から $x_{reach,CFD}$ へ移す方式**。
「自由壁が軸へ影響できる最初の位置」から設計すべき、という発想。

- **必須要件は全て満たした**: $x_{reach}$ 固定 / target は $x\ge x_{reach}$ のみ
  (上流呼び出しは `ValueError` で停止) / 始端で $M,M',M''$ が指定値と一致 /
  自由 CP 単調 (旧 B8 の 3.67→2.00 ジグザグを解消)。
- **B10 が暴いた本質**: $x_{reach}$ で旧 target は $M'=+0.511$, $M''=-0.196$ を
  要求していたが、CFD 実測は $M'=+0.290$, $M''=+0.190$。**$M'$ が 76% 違い $M''$ は
  符号が逆** — 旧 target は始端の整合性を犠牲にして出口を良くしていた。
- **代償**: CFD 値を C2 継承すると target 形状が最大 0.0758 変わり、**自由 CP 3 個では
  吸収できない**。接合帯は 0.0757→0.0082 と改善する一方、**出口が
  $\varepsilon_M$ 0.023%→0.128%、$|\theta|_{max}$ 0.015°→0.183°** に悪化 (要求の下限)。
- **結論**: 「始端の曲率整合」と「出口一様性」が現在の自由度で両立しない。
  未決の 3 択 = (i) 自由 CP 3→4/5 (ii) $x_{reach}$ で C1 に緩める
  (iii) CP を出口指標で最適化 (MOO へ委譲)。

**注意**: B10 の実装は `geometry.x_reach_cfd` を書かなければ**旧 B8 経路にフォールバック**
する。現行の生産構成は B8 (`problem_m4_b8.yaml`) のままで壊れていない。

---

## 5. 使えるインフラ (このセッションで整備したもの)

| 資産 | 場所 | 用途 |
| --- | --- | --- |
| `cminus_cfd.py` | `design/forge_design/geometry/` | **壁面状態を起点**にした C⁻ 追跡で $x_{reach,CFD}$ 抽出。感度 (出発オフセット/刻み/スナップ数) を自動記録 |
| `onesided_axis_anchor` / `axis_curve_true` / `anchor_sensitivity` | `feedback/cfd_anchor.py` | 上流片側の局所回帰で $M,M',M''$ (差分を使わない)。cell は偶関数フィットで真の軸値を復元 |
| `MachBezier.fit_free_cp` | `geometry/bezier.py` | 端点拘束下で自由 CP を参照曲線へ**拘束付き最小二乗射影** (手調整不要、CP 復元誤差 3.6e-15) |
| `_config_euler_node` | `evaluate/runner_wt.py` | node Euler 設定。**slip 壁では `nodeWallDirichlet` を書かない** |
| `Mesh2DParams.local_center/refine/width` | `meshing/mesh2d.py` | 接合部近傍のみ軸方向細分 (格子収束確認用) |
| `geometry.initial_wall_csv` | `feedback/euler_loop.py` | 帰還ループの初期壁を既存壁で指定 |
| `check_convergence.py` (改修) | `solver_density_cuda/tools/` | 低下桁数を**ピーク基準**に (node の IC は Uy≈0 で始まるため step0 基準だと誤判定)。プラトー呼吸の rising 誤判定ガードも追加 |

**node 化で見つけて直した cell 前提** (同種の監査を今後も要する):
`axis_mach` (双対 CV 重心を読んで非軸 DOF を 246 個混入させていた) /
`core_radius_traced` / `build_cminus_map_cfd`。`exit_uniformity` は node で健全と確認済。

---

## 6. 主要 run (詳細は [`case/41.wind_tunnel_design/README.md`](../../case/41.wind_tunnel_design/README.md))

| run | 位置づけ |
| --- | --- |
| `run_0026_b8` | **初の正式収束** (0.49% $M_d$)。生産構成 |
| `run_0027_c3` | 評価窓をリップまで延長。**仮壁の供給元** (pass_11/wall_next.csv = pass_12 が評価した壁) |
| `run_0030_ns_v1` | NS v1 (δ\* 相関壁 + y+1.5 node 低 Re SST)。軸 M が無帰還で 1.2% $M_d$ 内、δ\*/相関比 1.04–1.09 |
| `run_0031_baseline_cold` | cell の正式ゲート通過基準。`wall_provisional.csv` の置き場 |
| `run_0034/0035/0036_node_*` | **node の格子収束系列**。$x_{reach,CFD}=1.7042$、$M=1.9458$, $M'=+0.2896$, $M''=+0.19\pm0.03$ の出典 |
| `run_0037/0038_b10_*` | B10 の検証 (機構 OK・自由度不足を実証) |

**次に作る run は `run_0039` 以降** (0038 まで使用済み。調査 §8 が「run_0037_r5_cold」と
書いているのは番号衝突するので読み替えること)。

---

## 7. 繰り返してはいけない誤り (私が実際にやった)

1. **「支配方程式が切り替わるから不可避」は誤り**。CFD は全域で同じ Euler を解いており、
   切り替わるのは設計モデル。接合両側とも局所超音速。
2. **壁 $dM/dx$ の折れを「波源」と呼ばない**。壁 $M$ は境界条件でなく観測量で、因果が逆。
3. **Euler で「洗い流される」と言わない**。散逸はなく、特性線の広がり・波の相殺・数値散逸の結果。
4. **cell で $M''$ を測ろうとしない**。今回の抽出・再構成では 2 階微分が雑音支配
   (stencil 間 6.09、符号すら不一致)。node なら 0.049 に収まる。
5. **δ\* 抽出は探索窓をフリーストリームまで届かせる**。0.35 $r_t$ 窓では下流で
   打ち切られ「δ\* が下流ほど減る」という非物理な結論を出した (正しくは相関比 1.04–1.09)。
6. **帰還ループの初期壁を勝手に作り直さない**。仮壁 (既存の到達点) を使うこと。
   逆設計壁で始めて 170 µm 分の改良を捨てた。
7. **cell で成立した処理が node でも成立するとは限らない**。座標・重み・境界 DOF・
   補間処理は毎回監査する (§5 の実例参照)。
8. **帳簿の改善を物理の改善と混同しない**。評価参照を変えれば数値は良くなるが、
   意図非依存の物理指標 (うねり振幅・出口 ε) を必ず併記して判断する。

---

## 8. 次セッションでの進め方 (提案)

ユーザが ChatGPT 推奨方針を持ち込む。**それを §2 の確定事実と §3 の文献結論に
照らして評価する**のが最初の仕事。特に:

- その方針が **J 下流の壁表現だけを変えるもの**なら、§2 の「波源は円弧上」と衝突する
  ので、こぶには効かない見込みである旨を先に伝える (別の効能はあり得る)。
- その方針が **$R$・遷音速アンカー・接合構成のいずれかに触れる**なら、
  §3 の 3 経路と整合するので有望。
- どちらであれ、**先に $R=5$ の A/B (調査推奨) を走らせるかを確認する** — 表現を
  変えずこぶの $R$ 依存を測れるので、以後の判断が安くなる。

進める際は AGENTS の必須ルール (新連番 run・メッシュ品質 VERDICT・NaN チェック・
`check_convergence.py` / `check_quasisteady.py` の VERDICT 貼付・README run 一覧同期・
検証後 commit/push) を守ること。

---

## 9. 関連文書

- 経緯と全実測値: [`plans/active/tooling-nozzle-phase3-windtunnel.md`](../../plans/active/tooling-nozzle-phase3-windtunnel.md) §9.2
  (**B5/B6/B9/B10/B10-a/B10-b/B10-c。撤回済みの主張も明記してある**)
- 現在仕様: [`methods/design/overview.md`](../../methods/design/overview.md)
- 文献調査: [`notes/investigations/nozzle-throat-curvature-shape-representation-survey.md`](../investigations/nozzle-throat-curvature-shape-representation-survey.md)
- 調査依頼プロンプト: [`notes/sessions/2026-08-17-throat-curvature-shape-representation-survey.md`](2026-08-17-throat-curvature-shape-representation-survey.md)
- 論文 PDF: [`papers/nozzle_design/`](../../papers/nozzle_design/) (CONTUR 本体・Sauer 原典・Cuffel・Korte 2000 等)
