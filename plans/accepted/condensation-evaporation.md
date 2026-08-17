# 非平衡凝縮に蒸発 (負成長・液滴消滅) を追加する

## メタ

- **area**: `condensation`
- **status**: `done`  <!-- 2026-08-18 実装・検証済 (Phase A: 一温度, S<=1 蒸発)。Phase B (二温度 T_d<=T, S>1 亜臨界) は未着手 -->
- **related_docs**:
  - `methods/condensation.md` (現在仕様。「安定化」節に「dr/dt<0 は 0 (蒸発なし)」と明記されている)
- **related_plans**:
  - [condensation-nonequilibrium.md](condensation-nonequilibrium.md) (親: 4 モーメント凝縮の実装)
  - [tooling-nozzle-tp-split-h2o-condensation.md](tooling-nozzle-tp-split-h2o-condensation.md) (燃焼ガス風洞の H₂O 凝縮; 本 plan の動機となった run_0065)
- **created**: `2026-08-18`
- **owner**: `CFD Dev`

## 1. 目的

現行の凝縮モデルは **一方向 (凝縮のみ)** で、成長率 $dr/dt<0$ とソース $S_g<0$ を 0 にクランプし、
過飽和 $S<1$ では臨界半径 $r_*$ が定義されないため成長項自体を評価しない。その結果、いったん生成した
液相は**過熱域 (圧縮波・衝撃波・高温境界層) に入っても受動スカラーとして凍結**され、蒸発潜熱の吸収も
蒸気の復帰も起きない。本 plan では **蒸発 (負成長) と液滴の完全消滅**を 4 モーメント法に追加し、
$S<1$ の領域で液相が物理的な速度で蒸気に戻り、温度が対応して下がる状態を実現する。

## 2. 現状の実測 (2026-08-18)

対象: `case/42.isobutane_wt/run_0065_ib_R2_LU6_Lc8_h2o5_split_cond` (node Euler TP, slip 壁, H₂O 5 %,
Kw+HK, `res_12000.h5`)。

- 壁は **slip (Euler)** なので壁面加熱は起きていない。壁列 T≈288 K, S=1.00–1.09 (飽和線上)。
  「高温壁で蒸発しない」現象は**この run では未発生** (壁が熱くない)。粘性 (no-slip) 版に進めば
  回復温度 ~900 K の境界層に液相が凍結する形で顕在化する。
- 一方、試験部内の圧縮帯 (x≈2.6–3.07, y≈0.27–0.40; 接合波・内部衝撃の反射) で
  **S=0.71–0.95 かつ g≈1.2 % が凍結**した領域が 514 ノード (液相質量の 9.2 %) 存在する。
  ここが「蒸発すべきなのに蒸発しない」の実例。
- 液滴サイズ $\bar r=Q_1/Q_0$ ≈ 0.28 µm (中央値), 0.39 µm (圧縮帯), 数密度 $Q_0$≈5×10¹³ /kg。
- Hertz–Knudsen 蒸発率の見積り: S=0.7/290 K で $dr/dt$≈−6×10⁻⁴ m/s (0.39 µm 液滴が 0.6 ms で消滅、
  帯の滞在 ~0.4 ms → 帯内で g が数割減), S=0.3/400 K で 50 nm が 0.3 µs, 900 K では ns オーダー
  (**高温境界層では事実上瞬時**)。潜熱吸収は g=1.25 % 全量で ΔT≈−34 K。
- 圧縮帯では蒸発→冷却→$p_{sat}$ 低下で S が 1 に戻る**自己制限 (平衡蒸発)** になるため、Euler run
  への影響は「帯内で g が 1 割程度減り T が数 K 下がる」程度と見込む。NS 境界層では液相が消えて
  T が最大 ~30 K 下がるため、壁温・$q_w$・$\delta^*$ に効く。

## 3. スコープ

- **やる**
  - $S<1$ (過熱・不飽和) セルでの**負成長** (蒸発) を全成長則 (HK / Goodheart / Gyarmathy) に追加。
  - 蒸発時のモーメント更新の**実現可能性** (同一 θ 縮小、$|\Delta g|\le g$、$\Delta T$ 律速)。
  - **完全蒸発カットオフ** ($\bar r<r_{\min}$ or $g<g_{\min}$ で 4 モーメントを 0 化 = 液滴消滅)。
  - 点陰的ヤコビアン (`src_jac`) の蒸発分岐 (蒸気復帰による自己抑制)。本体 kernel と
    `cond_source_vector` の分岐を同一化。
  - config フラグ `condEvaporation` (既定 0 で回帰同一、検証後に既定 1 へ)。
  - 診断出力 (過飽和 $S$、$dr/dt$)。
- **やらない (後続)**
  - $S>1$ かつ $\bar r<r_*$ の **Kelvin 亜臨界蒸発** (核生成域で $\bar r$ が $r_*$ に張り付く
    問題の再燃を避けるため Phase A ではクランプ維持。`condEvapKelvin=1` で試験のみ)。
  - 二温度 (`condTwoTemp=1`) の蒸発 ($T_d\le T$ 側の Newton)。
  - 壁面付着 (液膜) や液滴の壁衝突・慣性 (モーメントは気相に追従、拡散なし)。
  - 固相 (氷) への切替。

## 4. 設計

### 4.1 統一駆動力 (核生成 $r_*$ に依存しない成長率)

現在の N2 Goodheart / Gyarmathy 成長率は $\ln S\,(1-r_*/\bar r)$ を駆動力にしているが、
$r_*\ln S = 2\sigma/(\rho_l R T)\equiv K_e$ なので

$$
\ln S\left(1-\frac{r_*}{\bar r}\right) \;=\; \ln S - \frac{K_e}{\bar r} \;=\; \ln\frac{p_v}{p_d(\bar r)},\qquad
p_d(\bar r)=p_{sat}(T)\exp\!\left(\frac{K_e}{\bar r}\right)
$$

と**恒等**。右辺は $S<1$ でも定義され (このとき常に負)、HK の駆動力 $p_v-p_d(\bar r)$ と同じ
Kelvin 平衡蒸気圧 $p_d$ を使う。実装は駆動力を $\ln S - K_e/\bar r$ (または $p_v-p_d$) で直接計算し、
`cond_nucleation` が返す $r_*$ に成長率を依存させない (現状 $S\le1$ で $r_*=0$ が返り分岐が壊れる)。

- **HK (H2O 既定)**: 既存式 $\frac{\alpha}{\rho_l}\frac{p_v-p_d}{\sqrt{2\pi R T}}$ をそのまま (符号付き)。
- **Goodheart / Gyarmathy**: 前因子はそのまま、駆動力を $\ln S - K_e/\bar r$ に置換。
  Gyarmathy の `if (p_v <= psat) return 0` を撤去。

### 4.2 分岐 (Phase A)

| 状態 | 核生成 $J$ | 成長 $dr/dt$ |
| --- | --- | --- |
| $S>1$, $\bar r>r_*$ | あり | 正 (現状どおり) |
| $S>1$, $\bar r\le r_*$ | あり | **0 (現状維持; Kelvin 亜臨界蒸発は `condEvapKelvin` で後段)** |
| $S\le1$, $Q_0>0$ | 0 | **負 (蒸発)** — `condEvaporation=1` のとき |

$S\le1$ では $r_*$ が不要なので `rstar>0` ガードに依らず `S`, `Q0` で分岐する。

### 4.3 モーメント更新と実現可能性

負成長でも同じモーメント式 $S_{Q_1}=Q_0\,dr/dt$, $S_{Q_2}=2Q_1\,dr/dt$,
$S_g=4\pi\rho_l Q_2\,dr/dt$ を使う ($Q_0$ は部分蒸発では不変 = 全液滴が同じだけ縮む monodisperse 仮定)。
1 step の過剰除去を防ぐため θ 律速を蒸発側に拡張する:

- $|\Delta g| \le g$ (存在する液相以上は蒸発できない)、$|\Delta g|\le$ `dg_max`、
  潜熱冷却 $|\Delta T| = |\Delta g|L/c_v \le$ `dT_max` (現状 1 K/step)。
- 全モーメントに同一 θ を掛け $g=Q_3$ 整合を保つ (現状の凝縮側と同じ)。
- `Sg<0→0` クランプは `condEvaporation=1` で解除。

### 4.4 完全蒸発 (液滴消滅)

monodisperse 更新では $\bar r\to0$ まで $Q_0$ が残り、$1/\bar r$ を含む Kn 項が発散する。また
再飽和時に微小残滓と新核が混ざり $\bar r$ を歪める。そこで **realizability クランプ kernel
(`cond_realizability_clamp_d`) に消滅条件を追加**する:

- $S\le1$ かつ ($\bar r<r_{\min}$ [既定 1 nm, `condEvapRmin`] または $g<g_{\min}$ [既定 $10^{-9}Y_w$])
  → `rog, roQ0, roQ1, roQ2` を 0。
- 条件を $S\le1$ に限定するのは、核生成域 ($S\gg1$ では $r_*$≈0.5 nm、$r_{nuc}=1.01r_*$) の
  新核を殺さないため。
- 質量は $Y_w$ (総水) が不変なので自動保存、エネルギーは二相 EOS $e=e_v-gL$ が $g\to0$ で
  潜熱吸収を自動処理 (残滓 $g<g_{\min}$ の温度誤差は無視できる)。

### 4.5 点陰的ヤコビアン

- $S_g$ 行: 蒸発では $p_v=\rho(Y_w-g)R_wT$ より $\partial p_v/\partial(\rho g)=-R_wT$、
  HK で $\partial S_g/\partial(\rho g) = -4\pi Q_2\,\alpha R_w T/\sqrt{2\pi RT}<0$ → `sj_g` $=+4\pi Q_2\alpha R_wT/\sqrt{2\pi RT}\ge0$
  (蒸気復帰による自己抑制)。既存の潜熱経由項 ($\partial S_g/\partial T\cdot\partial T/\partial(\rho g)$) は
  蒸発側でも同符号 (冷却→$p_{sat}$↓→蒸発↓) なので数値微分をそのまま流用。
- $S_{Q_1}$ 行: $\partial S_{Q_1}/\partial Q_1 = \partial(dr/dt)/\partial\bar r>0$ (Kelvin: 小さいほど速く蒸発)
  = 正帰還なので `sj_Q1=0` (現状のクランプ ≥0 で自然にそうなる)。安定化は θ 律速と 4.4 の消滅で担う。
- `cond_source_vector` (ヤコビアン用) に本体 kernel と**同一の分岐**を入れる (片側だけだと
  過大な $\partial S/\partial Q_1$ で roQ1 が潰れる既知の罠)。

### 4.6 config / 出力

| キー | 既定 | 意味 |
| --- | --- | --- |
| `condEvaporation` | **1** (検証後に既定化) | 1 で $S\le1$ の負成長・液滴消滅を有効化。0 で旧挙動 |
| `condEvapRmin` | 1e-9 | 完全蒸発とみなす $\bar r$ [m] |
| `condEvapKelvin` | 0 | 1 で蒸発駆動力に Kelvin 項 $p_d=p_{sat}e^{K_e/r_{30}}$ を含める (既定は平面 $p_{sat}$; §5.1-1)。$S>1$ 亜臨界蒸発は本 plan では扱わない |

診断出力: `condS_<s>` (過飽和 $p_v/p_{sat}$)、`condDrdt_<s>` [m/s] (負=蒸発)、`condR30_<s>` [m]
(蒸発分岐のみ)。現状は S を後処理で再計算しないと蒸発域が見えないため。

### 4.7 二温度 (Phase B)

`condTwoTemp=1` の HK 経路は Newton 内で `Td<T → Td=T` にクランプしている (凝縮中 $T_d\ge T$)。
蒸発では $L\,j = h(T_d-T)$ の $j<0$ で $T_d<T$ (蒸発冷却) になるので、クランプを
「$j$ の符号に応じ $T_d\ge T$ / $T_d\le T$」に置換する。Phase A は一温度のみ。

## 5. 安定性の見立て

- 蒸発は **負帰還** (蒸発→潜熱吸収で T↓→$p_{sat}$↓→S↑→蒸発↓、蒸気復帰でも S↑) で、
  凝縮の潜熱自己抑制と鏡像。θ 律速 (dT_max) と `sj_g` で 1 step 内の暴走はない。
- 高温境界層では時定数 ns で剛だが、θ 律速下では「dT_max/step で減っていく」形になり、
  擬似時間の定常収束には問題ない (非定常精度が要るなら dT_max を緩めるか sub-step)。
- 核生成域 ($S>1$) の挙動は Phase A では**一切変えない** (分岐が $S\le1$ 限定) ので、
  case/16 Wyslouzil / case/42 の凝縮 onset は回帰同一のはず。

### 5.1 Kelvin 正帰還 (小さいほど速く蒸発) への打ち手 (2026-08-18 ユーザ指摘で追加)

$p_d(\bar r)=p_{sat}e^{K_e/\bar r}$ は $\bar r$ が小さいほど蒸発を速める**正帰還**で、離散的には
「$Q_1$ が先に潰れる → $\bar r=Q_1/Q_0\to0$ → $e^{K_e/\bar r}$, $1/\bar r$ (Kn) 発散」の経路
(過去の roQ1 暴走崩壊と同じ) が開く。Phase A の既定は次で正帰還ループ自体を消す:

1. **蒸発分岐では Kelvin 項を落とす** (既定): 駆動力を平面飽和圧 $p_v-p_{sat}$ (Goodheart/Gyarmathy は
   $\ln S$) にする。$S\le1$ では Kelvin なしでも常に蒸発する。Kelvin 長さは水 290 K で
   $K_e=2\sigma/(\rho_lR_wT)\approx1.1$ nm、$p_d/p_{sat}=e^{K_e/r}$ は 100 nm で 1 %, 30 nm で 4 %,
   10 nm で 12 %。駆動力への寄与 $(e^{K_e/r}-1)/(1-S)$ は $S\to1$ では 100 nm でも無視できない
   (率は変わる) が、**質量収支には効かない**: $r_0$ から $r$ まで縮んだ時点の残量は $(r/r_0)^3$ で、
   $r_0$≈300 nm (run_0065) なら $r<10$ nm 段階の残量は $3.7\times10^{-5}$, $r<30$ nm でも $10^{-3}$。
   落として変わるのは液相質量の最後の $10^{-4}\sim10^{-3}$ が消える時間だけで、蒸発総量と潜熱
   ($\Delta T=\Delta g L/c_v$)・平衡蒸発の到達点 ($S=1$) は変わらない。Kelvin 込み蒸発は
   `condEvapKelvin=1` に退避。
2. **蒸発率の評価半径を質量ベースに**: $r_{30}=(g/(\tfrac43\pi\rho_lQ_0))^{1/3}$ を使い、
   帰還を θ 律速される $g$ 経由にする ($Q_1$ 崩壊で率が発散しない)。$r_{10}/r_{30}$ の乖離は診断出力。
3. **ステップ制限 + 消滅**: 縮小比 $\lambda=r_{\rm new}/\bar r$ を $\lambda\ge\lambda_{\min}=\tfrac12$
   でクランプ (1 step の質量減は高々 $1-\lambda_{\min}^3=7/8$、別途 dT_max の θ)。クランプ後でも
   $\lambda_{\min}\bar r<r_{\min}$、すなわち **$\bar r<r_{\min}/\lambda_{\min}=2r_{\min}$** なら中間状態を
   作らず 4 モーメントを一括 0 (=§4.4)。$r_0$≈300 nm から毎 step 半減しても ~7 step で到達し、
   消える瞬間の残量は $(2/300)^3\approx3\times10^{-7}$ で質量・エネルギーの飛びは無視できる。
4. **比スケール更新**: $\lambda=r_{\rm new}/\bar r$ で $Q_1\to\lambda Q_1,\ Q_2\to\lambda^2Q_2,\ g\to\lambda^3g$
   と一括スケール ($S_\phi=(\phi_{\rm new}-\phi)/\Delta t$ で残差経路へ)。$Q_1^2\le Q_0Q_2$ 等の
   実現可能性が自動で保たれ「$Q_1$ だけ 0 で $g$ が残る」不整合を防ぐ。
5. **数値ガード**: $1/\bar r$, $e^{K_e/\bar r}$ 評価前に $\bar r\ge r_{\rm floor}$; `sj_Q1=0`
   (正帰還項は陰的に入れない); 消滅発火セル数・最小 $\bar r$ をログ出力。
6. (Phase B, 必要時のみ) セル局所 0-D ODE の適応サブステップ積分。

$S>1$ 亜臨界域 ($\bar r<r_*$) は現状クランプ + `r_nuc=1.01r_*` のまま触らない。

## 6. 検証計画

1. **0-D 検証** (`solver_density_cuda/tests/` に host 側テスト or Python レプリカ):
   S=0.5 の空気中に $\bar r$=100 nm, g=1 % の液滴雲を置き、HK 解析解と $g(t)$, $\bar r(t)$ が一致・
   $r_{\min}$ で全モーメント 0・エネルギー $\Delta T=-\Delta g L/c_v$ を確認。
2. **回帰**: case/16 Wyslouzil (Fig.3 の $p/p_0$ 一致が不変)、case/42 `run_0064` dry (無影響) と
   `run_0065` → 新 run (`run_0066_..._evap`): 上流 ($S>1$ 域) はビット同一、圧縮帯 (x>2.6) で
   g 減少・T 低下・S→1 に収束することを確認。`check_convergence.py` / `check_quasisteady.py` を貼る。
3. **高温壁デモ (NS)**: case/42 の NS 版 (`prepare_ns` の no-slip 断熱壁) or case/16 NS で、
   境界層内 (T > 350 K) の g が 0、壁温が回復温度になること。凍結版との壁温差 (~30 K 上限) を報告。
4. NaN/発散チェック (AGENTS.md 規約)。

## 7. 変更ログ

- 2026-08-18: 起票。run_0065 の実測 (§2) を記録。統一駆動力の説明とユーザ指摘 (Kelvin 正帰還) を受け §5.1 追加。
- 2026-08-18: **実装・検証完了 (Phase A)**。
  - 実装: `cond_evap_rate` / `cond_evap_source` (統一駆動力, Kelvin off 既定, $r_{30}$, λ スケール, λ_min=0.5,
    dg_max/dT_max 律速, $r_{30}<2r_{\min}$ で一括消滅) と本体 kernel の蒸発分岐 + 数値 `sj_g`
    (`condensationSource_d.{cu,cuh}`); 消滅硬クランプ + g=0 モーメント塵掃除 (`cond_realizability_clamp_d`);
    config `condEvaporation`(既定 1)/`condEvapRmin`/`condEvapKelvin`; 診断 `condS_/condDrdt_/condR30_`;
    `interp_field.py` が凝縮モーメント dataset を新規作成 (restart で液相を運べるように)。
  - 0-D 単体 `tests/unit/test_cond_evaporation.cpp`: HK 解析解一致 (r 1e-15, g 1e-15)、$r_{10}=r_{30}$ 保持、
    実現可能性、λ_min 律速→8 半減で一括消滅、dT_max 律速 1.000 K、Kelvin on/off で $r<10$ nm の残量
    2.7e-5/2.4e-5 (質量収支不変, 蒸発時間 4.78e-4 vs 4.56e-4 s)、S>1 不発火、N2 経路負・有限。ALL PASS。
  - Euler 圧縮帯 (case/42 run_0066 vs 0067 対照, run_0065 から restart): S_min 0.67→0.87、帯内 T −4.4 K、
    g −10 %、S→1 平衡蒸発。上流不変 (|ΔT| 4e-3 K)。`condDrdt_0` は HK と 0.7 % 一致、λ≥0.86 (律速非発火)。
    step 4000 で定常 (4000/8000/12000 同一)。
  - 回帰 (case/16 run_0187/0188 vs 0186): 蒸発 OFF/ON とも run 間ノイズ内 (|ΔT| ≤8e-3 K)。Wyslouzil は液相域
    が S≥2 で不発火。case/16 SST no-slip (run_0172/0173, Tt=287 K の低温壁) は壁近傍 S 0.4–0.99 に有限速度の
    残液 (物理どおり) と restart ドリフト (既知)。
  - **NS 高温壁 (case/42 run_0068 dry NS → 0069 蒸発 OFF / 0070 ON / 0071 既定 ON)**: 壁 T≈920 K。OFF では
    T>400 K のノードに液相 43 個 (g≤2.5e-6) と液相の 46 % が S<1 に凍結; ON では T>400 K の液相 **0**、
    S<1 の液相質量 35 % に減り S 中央値 0.89→0.997 (平衡蒸発)、出口で最大 −11 K。壁列 (wall_dist<0.1 mm) は
    OFF でも液相 1.9e-24 = **液滴は境界層に入らない** (モーメントに拡散なし; 境界層流体は上流でも高温)。
    「壁面で蒸発しない」の実体は「壁近傍の過熱域で液相が凍結」であり、これは解消。
  - 残課題 (Phase B, 別 plan): 二温度 $T_d\le T$; $S>1$ 亜臨界蒸発; 液滴の乱流拡散 (RANS で境界層へ液滴が
    入らない); 消滅の硬クランプが前ステップ T で判定される点。
