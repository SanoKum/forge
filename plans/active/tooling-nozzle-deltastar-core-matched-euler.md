# ノズル排除厚さ更新計画: 積分法初期壁 + 固定 Euler 基準 CFD 反復

## メタ

- **area**: `tooling / design / boundary layer`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「粘性 δ\* 補正」節 (本計画で書き換え)
- **related_plans**:
  - 旧: [tooling-nozzle-axismach-viscous-deltastar.md](../accepted/tooling-nozzle-axismach-viscous-deltastar.md) (A12: 相関 δ\* + ρU_x 最大縁抽出 + x_lo 継ぎはぎ — 本計画で生産経路から置換)
  - 旧: [tooling-nozzle-axismach-physical-throat.md](../accepted/tooling-nozzle-axismach-physical-throat.md) (A13: 真のスロート探索・上流 Hermite 再生成 — 機構は維持)
  - 根拠: [notes/investigations/nozzle-deltastar-throat-review.md](../../notes/investigations/nozzle-deltastar-throat-review.md) (スロート δ\* 3〜12 倍過大の実測)
  - 影響: [design-isobutane-m6-d155.md](../accepted/design-isobutane-m6-d155.md) (Md トリム 6.0144 は本計画完了後に撤回)
- **created**: `2026-09-04`
- **owner**: `sano` (実装: Claude 自走)

## 1. 目的

排除厚さ補正を次の 2 段構成に置き換える。

1. **初回 NS 用の壁**は、入口から前進積分する境界層積分法で推定する (初期値生成専用)。
2. **NS 計算後**は、反復中固定した設計 Euler 場と NS 場の差から、**スロートを含む全域**の半径方向等価排除厚さ
   $\delta_r(x)$ を単一の抽出法で取り、固定点反復する。

温度縁・質量流束最大縁・$x_{lo}$・上下流の方式切替・相関とのブレンド・Md トリム・law 側 Mach 帰還は生産経路から廃止する。

到達目標 (case/45): NS/Euler 質量流量比 $1.000\pm0.003$、出口コア M $6.00\pm0.006$ を**トリムなし**で満たし、
初期壁の与え方 (相関 / 積分法) に依らず同じ固定点に収束する。

## 2. スコープ

- **やる**
  - 積分法初期推定器 `feedback/deltastar_integral.py` (断熱壁 / 指定壁温、CPG / semi-perfect)
  - 固定 Euler 基準抽出器 `metrics/deltastar.py::deltastar_from_core_matched_euler`
  - `PhysicalNozzleWall` の半径方向オフセット化、`prepare_ns` / `collect` の統合 (固定 Euler 参照・緩和・帳簿)
  - 単体試験 `design/tests/run_deltastar_tests.py`
  - V0 (既存 run で抽出器成立性) → V1 (case/45 相関初期壁から反復) → V2 (積分法初期壁) → V3 (等温) → V4 (移植性)
- **やらない (非スコープ・生産経路から廃止)**
  - 温度縁 / 質量流束最大縁の生産抽出 (旧 `deltastar_from_run` は比較用に残置)
  - $x_{lo}$、上流相関補完、積分法と CFD のブレンド
  - Md トリム、law 側 Mach 帰還
  - Sasman–Cresci 閉包そのもの (§4.1 参照。閉包差し替え口は残す)

## 3. 関連 docs と前提

- 診断ノートの事実: 与えたスロート δ\* は NS 実効値 ($0.0013$〜$0.0021\,r_t$) の 3〜12 倍。NS 質量流量が Euler 比 +0.8〜3.7 %。
  1D 感度 $d\ln M/d\ln(A/A^*)\approx 1/5$〜$1/6$ で試験部 M −0.2〜−0.7 %。観測残差と ±0.1 % で一致。
- 固定点では **NS のコア = 「設計壁 $r_{inv}$ を壁とする非粘性流」そのもの**。したがって同じ物理 $r$ で設計 Euler と一致する。
  この事実が抽出器の定義 (§4.3) の根拠。
- 設計 Euler run は全 case に存在する (case/45 `run_0001`, case/42 `run_0096`/`run_0055`, case/44 `run_0005`)。
- 構造格子 (node, index $= i\,n_j + j$) の列 $i$ は $x$ 一定線なので、断面積分は列の台形則で厳密に取れる (補間なし)。

## 4. 設計方針

### 4.1 初期推定器: CONTUR 型運動量積分法 (Sasman–Cresci の代替)

ユーザ計画は Sasman–Cresci (SC) を指定していたが、SC 閉包 (モーメント運動量方程式のせん断積分・形状係数関係) は原論文
(Sasman & Cresci, AIAA J. 1966) が手元に無く再現できない。同じ役割 (入口から前進積分、加速による薄化を表現、断熱/等温両対応)
を持ち、かつ**手元の一次資料で式が確定している** Sivells CONTUR (AEDC-TR-78-63 §5, Eq. 61–90,
[papers/nozzle_design/](../../papers/nozzle_design/)) の方法を実装する。CONTUR は超音速・極超音速風洞ノズルの粘性補正の
標準実装であり、目的に対して SC より劣る理由はない。閉包は `closure:` で差し替え可能にし、SC は原論文入手後の追加とする。

軸対称 von Kármán 運動量積分 (Eq. 61):

$$\frac{d\theta}{dx} + \theta\left[\frac{2 - M^2 + H}{M\,(1+\frac{\gamma-1}{2}M^2)}\frac{dM}{dx} + \frac{1}{r_w}\frac{dr_w}{dx}\right] = \frac{C_f}{2}\sec\phi_w$$

- プロファイル: $q/q_e = (z/\delta)^{1/N}$ (Eq. 67)、$\rho/\rho_e = T_e/T$、
  $T = T_w + a(T_{aw}-T_w)\,q/q_e + [T_e - a(T_{aw}-T_w) - T_w](q/q_e)^2$ (Eq. 69; $a=1$ Crocco / $a=0$ 放物)。
  断熱壁は $T_w = T_{aw} = T_e(1 + r_f\frac{\gamma-1}{2}M^2)$, $r_f = Pr^{1/3}$。
- $\theta,\ \delta^*$ は横曲率重み $(1 - z\cos\phi_w/r_w)$ 込み (Eq. 62/63)、$H=\delta^*/\theta$。
- 摩擦: $C_f = C_{f_i}/F_c$, $R_{\theta_i} = F_{R_\delta} R_{\theta_c}$ (Spalding–Chi 形, Eq. 70–73)、
  $F_c = [\int_0^1 (\rho/\rho_e)^{1/2} d(q/q_e)]^{-2}$、$F_{R_\delta} = \mu_e/\mu_w$、
  $C_{f_i} = 0.0773/[(\log R_{\theta_i}+4.561)(\log R_{\theta_i}-0.546)]$ (Eq. 75)。$R_{\theta_c}$ は平板形 $\theta_c$ (Eq. 74)。
- $N = N(R_\delta)$ は CONTUR Fig. 10 の曲線をテーブル化 (Eq. 77–82 の連鎖の結果)。
- 各ステーションで $\theta$ (状態量) から $\delta, N$ を内部反復で決め、$\delta^*_n$ (壁法線) を出す。
- 半径方向補正へ変換: $\delta_r = \delta^*_n/\cos\theta_w$。初回壁 $r_{phys} = r_{inv} + \delta_r$。
- 入口初期条件: YAML `theta0`/`H0` 直接指定、未指定時は仮想発達長 $x_v$ (既定 = 入口直管長) の乱流平板 $\theta_0 = 0.036\,x_v Re_{x_v}^{-0.2}$。
  初期条件と採用理由は出力 JSON へ記録し、感度 (θ0 ×0.5/×2) を診断として出す。
- 熱境界条件: `thermal_bc.mode: adiabatic | prescribed_temperature` (`Tw` 定数または `Tw_table: [[x, Tw], ...]`)。
  同じ積分器へ壁温・物性 (Sutherland、TP は `gas.T_of_M`・$\gamma(T)$) を供給する。
- **合格条件は CFD 一致率ではない**: 全域で有限・非負、壁が滑らか、スロート探索と曲率ゲートを通る、初回 NS が立つ、
  最終固定点が初期推定器に依存しない。

### 4.2 座標対応 (Euler 場を NS 断面へ写す)

固定点では NS コアは設計壁 $r_{inv}$ の非粘性流に一致するので、対応は**同じ物理 $r$** が正解である。反復途中は NS の有効壁
$r_{ref}(x) = r_{w,NS}(x) - \delta_{r,in}(x)$ (この pass に与えた補正を差し引いた壁) が設計壁に対応する。よって

$$\eta_{NS} = \frac{r}{r_{w,NS} - \delta_{r,in}},\qquad \eta_E = \frac{r}{r_{inv}}$$

とし、同じ $\eta$ で対応させる (固定点で $\eta_{NS}=\eta_E=r/r_{inv}$)。
物理壁 $r_{w,NS}$ で正規化すると固定点でも $r_{w,NS}/r_{inv}$ (M6 で 1.5〜7 %) の半径引き伸ばしが残り、
コア勾配の強い膨張部で $\delta_r$ の 10〜20 % 級の偽欠損を生むため採らない (ユーザ計画からの修正点 1)。

- 走査は $x$ 一定断面の半径方向、流束は $q=\rho U_x$。
- NS 各列の不等間隔半径格子を積分格子にし、Euler 場を同じ $x$ (列間線形)・同じ $\eta$ (列内線形) で補間する。
  $\eta_E > 1$ (Euler 壁の外側) は Euler 壁値で定数外挿。
- Euler メッシュ・解は反復中固定、HDF5 は書き換えない。$x$ の無次元化は設計スロート半径 `scale_m` で固定。

### 4.3 コア整合と等価半径

コア範囲 $0\le\eta\le\eta_c$ ($\eta_c=0.30$ 既定、感度 0.25/0.35)。軸対称面積重みでコア流量を厳密一致させる倍率

$$\alpha = \frac{\int_0^{\eta_c} q_{NS}\,\eta\,d\eta}{\int_0^{\eta_c} q_E\,\eta\,d\eta},\qquad q_{ref} = \alpha\,q_E$$

コア形状の一致度はコア相対 RMS $\sqrt{\langle(q_{NS}/q_{ref}-1)^2\rangle}$ (軸直近 3 点除外版も併記)。

符号付き質量欠損 (正値クリップしない) と等価半径:

$$D = 2\pi\int_0^{r_w}(q_{ref}-q_{NS})\,r\,dr,\qquad
2\pi\int_{r_{eff}}^{r_w} q_{ref}\,r\,dr = D\ \Rightarrow\ \delta_r = r_w - r_{eff}$$

内部名は `delta_r_equiv`。

### 4.4 抽出ゲート

断面ごとに CSV/JSON へ: `x, r_wall_euler, r_wall_ns, r_ref, alpha_core, core_rms, core_rms_noaxis, core_maxdev,
mass_deficit, delta_r_raw, delta_r_smooth, delta_r_c25, delta_r_c35, ok, reason`。

不合格 (`ok=false`) 条件: コア形状 RMS > 閾値 / コア範囲感度 > 閾値 / $D<0$ / 求根解なし / 補間範囲外 /
欠損がコア全体に分布 (累積欠損の 80 % が $\eta<0.8$ で発生) / 壁更新後の幾何ゲート違反。
暫定閾値: コア RMS 1 %、コア範囲感度 5 % (V0 で確定)。

**pass 0 では診断扱い、固定点で本ゲート** (ユーザ計画からの修正点 2): pass 0 は有効輪郭が 2 % 太い状態で
コア形状 RMS 1 % を広域で破り得る。破った断面も値は使い `ok=false` で記録し、合否判定は最終固定点で行う。

### 4.5 平滑化と壁更新

全域を同じ抽出法で処理 ($x_{lo}$ なし、相関補完なし、積分法混入なし、下流切替なし、端点クリップ外挿なし)。
平滑化は壁の波打ちを防ぐ最低限 (平滑化スプライン、生値との差を必ず出力)。

$$\delta_{in}^{k+1} = (1-\omega)\,\delta_{in}^{k} + \omega\,\delta_{ext}^{k},\qquad r_{phys}^{k+1} = r_{inv} + \delta_{in}^{k+1}$$

$\omega=0.5$ 初期値、安定なら 1.0。壁更新は**半径方向** (法線オフセット + $x$ シフトを廃止)。真の物理スロート探索と
上流 Hermite 再生成は維持する。

### 4.6 runner 統合と帳簿

- `prepare_ns(problem, run_dir, ic_from, euler_ref=, delta_r_csv=, omega=, initializer=)`
  - 初期壁: `deltastar_initializer` (YAML) または CLI。`dstar_blend` は削除。
  - `delta_r_csv` は全域 `delta_r_equiv.csv` を直接使用。
  - `prepare_info.json`: 初期推定器設定・熱境界条件・固定 Euler 参照・コア範囲・緩和係数・**実際に使ったスロート補正量**。
- `collect`: Euler/NS 質量流量比、質量流量由来の等価スロート補正量 $r_{t,W} - \sqrt{\dot m_{NS}/\dot m_E}$、
  抽出ゲート集計、固定点差 (max/median of $|\delta^{k+1}-\delta^k|/\delta^k$)、初期推定値と最終値の差。

## 5. 実装ステップ

1. 文書: 本 plan、`methods/design/overview.md` δ\* 節、`plans/README.md`、調査ノート §6–8。
2. `metrics/deltastar.py`: `deltastar_from_core_matched_euler` + 断面 CSV/JSON 出力 + ゲート。
3. V0: 既存 run (case/42 M5/M6・case/44・case/45) で抽出器を実行し成立性を判定 (追加計算なし)。
4. `geometry/wall_axismach.py`: `PhysicalNozzleWall(offset="radial", delta_r_x=...)`。
   `evaluate/runner_axismach.py`: `prepare_ns` 統合・`collect` 帳簿・反復ドライバ `iterate_deltastar`。
5. V1: case/45 `run_0019_ns_cm_pass0` (相関初期壁 = run_0004 の壁と同一入力, Md 6.0) → 抽出 → `run_0020_ns_cm_pass1` → 固定点。
6. `feedback/deltastar_integral.py` (CONTUR 法) + 単体試験。
7. V2: case/45 積分法初期壁から pass0/pass1、同じ固定点への収束を確認。
8. V3: 等温 (単体 + 初回 NS 開始性)。V4: case/42 M5 / case/44 へ設定不変で適用。
9. 文書同期、Md トリムの撤回、commit/push。

## 6. 検証

- **単体** (`design/tests/run_deltastar_tests.py`): $q_{NS}=c\,q_E$ で $\delta_r=0$ / 既知の壁側欠損から円環厚さ復元 /
  非一様 Euler プロファイルでも復元 / 壁半径が異なっても $\eta$ 写像で復元 / 点数・伸長率収束 / $\alpha$ で振幅差除去 /
  コア形状差をゲートで検出 / 25–35 % 感度 / 積分法 CPG 断熱・等温 / 等温 → 断熱の連続性 / 初期条件感度 / 法線→半径変換。
- **V0**: 全域で有限・非負、$\alpha(x)$ 滑らか、コア形状差小、コア範囲感度小、累積欠損が壁側で増える、下流で旧抽出器と整合、
  スロートで質量流量由来の等価値と整合。
- **V1 主ゲート**: NS/Euler 質量流量比 $1.000\pm0.003$、出口コア M $6.00\pm0.006$。品質・NaN・収束・準定常は AGENTS どおり。
- **V2**: 相関初期壁と積分法初期壁が同じ固定点 ($\delta_r$・質量流量・出口 M)。
- **V4**: case/42 M5 と case/44 にコア 30 %・平滑化・緩和係数を変えず適用。

## 7. 影響範囲

- `design/forge_design/{metrics/deltastar.py, feedback/deltastar_integral.py, geometry/wall_axismach.py, evaluate/runner_axismach.py}`
- `design/tests/run_deltastar_tests.py`
- `methods/design/overview.md`、`plans/README.md`、case/45・42・44 README run 表
- 旧 `feedback/deltastar.py` (平板相関) は比較・回帰用に残置。case 配下の `build_dstar_v3_*.py` は非推奨 (履歴として残す)。

## 8. 完了条件

- [ ] 全域で同じ CFD 抽出法を使用 / 積分法は初回壁のみ / 断熱・等温を同じ実装で扱える
- [ ] Euler は NS 反復中固定 / NS メッシュ変更へ補間で対応 / 半径方向補正へ統一
- [ ] case/45 で質量流量と出口 M のゲートを満たす / 初期推定方法によらず同じ固定点
- [ ] case/42 または case/44 で移植性を確認
- [ ] methods・plan・case README を同期、検証済みマイルストーンごとに commit/push

## 9. 変更ログ

- `2026-09-04` — 起票 (ユーザ計画を採択。修正 3 点: η 正規化を有効壁で / pass 0 ゲートは診断扱い / SC → CONTUR 運動量積分法)。
