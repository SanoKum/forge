# SST-DDES / SST-IDDES 実装計画

## メタ

- **area**: `turbulence`
- **status**: `in_progress`
- **related_docs**:
  - `methods/turbulence/theory.md`
  - `methods/turbulence/implementation.md`
- **related_plans**:
  - [`turbulence-des-wmles-survey.md`](../../notes/investigations/turbulence-des-wmles-survey.md) — 本計画の背景サーベイ
  - [`discretization-median-dual.md`](discretization-median-dual.md) — node-centered 両対応（§5.6 で整合を記述）
- **created**: `2026-06-21`
- **owner**: `CFD Dev`

---

## 1. 目的

Menter SST（`LESorRANS=2`）に DDES / IDDES の length scale 修正を追加し、
**衝撃波-境界層干渉（SBLI）・ピントルバルブ内部流れ**といった剥離支配の非定常解析を
RANS より高精度に、wall-resolved LES より低コストで行えるようにする。

完了時の状態：
- `DESmode: 1`（DDES）で剥離域のみ LES、付着 BL は RANS のまま
- `DESmode: 2`（IDDES）で上記 + 解像乱流入力がある場合に WMLES モードへ自動切替
- Forge の既存 SST ケース（平板・ノズル）はビット不変（`DESmode: 0` 既定）

---

## 2. スコープ

> **重要な前提（外部レビュー 2026-06-21 を反映）**: 「DDES length scale を動かす」ことと
> 「SBLI で実験と定量一致する」ことは別問題。**最大の失敗要因は DDES 式ではなく、
> ①対流スキームの数値散逸（風上散逸が LES 領域の解像乱流を殺す）②流入乱流の不足
> ③3D スパン解像 ④統計サンプル不足**。Phase 1 のまま Settles 24° に行くと
> 「剥離せん断層が育たず壁圧プラトー・再付着が RANS 寄りに残る」可能性が高い（§5.7 失敗予測）。
> このため検証フェーズを Phase 1（機能）と Phase 1.5（低散逸 flux = 物理検証の前提）に分ける。

**やる（Phase 1: DDES 機能実装）**:
- `DESmode: 1` の config 追加
- per-cell グリッドスケール $\Delta$ (Δmax) の事前計算
- $r_d$, $f_d$ シールド関数の GPU カーネル実装（**float32 ガード必須・§4.4**）
- $l_\mathrm{DDES}$ の計算・保存
- `ransSource_d.cu` の k 消滅項を $\beta^* \rho k \omega \to \rho k^{3/2}/l_\mathrm{DDES}$ に切替
- **診断フィールド一式の出力**（`f_d`, `r_d`, `delta_les`, `l_des`, `l_des/l_rans`・§4.2）
- `case/18.backstep`（3D）で**機能**動作確認（f_d 分布・NaN なし）

**やる（Phase 1.5: DES 用低散逸 flux — SBLI 物理検証の前提）**:
- 「LES 域は KE 保存中心差分、RANS と衝撃波は風上」の flux blending を導入（§4.8）。
- **短期（本命）**: 既存 `solver: KEEP_SLAU` の中心↔風上 blend 係数 `duc` を **`1-f_d` 主導**にする
  （`duc=max(0.05, 1-f_d, 1-ψ_limiter)`、SU2 `ROE_LOW_DISSIPATION=FD` と一致）。**ψ_limiter は衝撃
  センサにせず positivity フロア扱い**（limiter をセンサ流用する一次 DES 文献は皆無・[survey](../../notes/investigations/turbulence-des-flux-survey.md) §4.1）。
- **中期（SBLI/ノズル/ピントル）**: 実衝撃と剥離せん断層が共存するため、**改良 Ducros センサ**を
  SU2 `NTS_DUCROS` 相当の OR ブレンド `σ=σ_Ducros+(1-f_d)-σ_Ducros(1-f_d)` で層化（§4.8）。
- これなしでは backstep は見えても SBLI でせん断層 KH 成長が遅れる（レビュー指摘）。
- 設計の一次調査根拠は [`turbulence-des-flux-survey.md`](../../notes/investigations/turbulence-des-flux-survey.md)。

**やる（Phase 2: IDDES）**:
- `DESmode: 2` の config 追加
- $f_B$, $f_e$, $f_\mathrm{dt}$ の**完全な関数群**（§4.6 の簡略版は機能確認用。定量評価には不可）
- $l_\mathrm{IDDES}$ = DDES/WMLES 自動切替
- **WMLES モードの定量評価には流入乱流生成が前提**（下記「やらない」と矛盾するため、
  IDDES 定量評価は流入乱流 plan 完了後にゲートする）

**やらない（別 plan）**:
- 乱流流入生成（Synthetic Eddy Method / Recycling）—— WMLES モード定量評価の前提だが別タスク。
  これが無い間は IDDES を「定量一致」目的で使わない（DDES に留める）。
- ZDES ゾーン指定
- k-ε / SA ベースの DES（forge は SST ベース固定）

---

## 3. 理論

### 3.1 SST 内の RANS length scale

$$
l_\mathrm{RANS} = \frac{\sqrt{k}}{\beta^{*} \omega}
$$

$\beta^* = 0.09$。SST k-方程式の消滅項：
$$
D_k = \beta^* \rho k \omega = \rho k^{3/2} / l_\mathrm{RANS}
$$

> **訂正 (実装時に判明・2026-06-21)**: 初稿は $l_\mathrm{RANS}=\sqrt{k}/(\beta^{*1/4}\omega)$ と
> 書いていたが、これは上の恒等式 $D_k=\beta^{*}\rho k\omega=\rho k^{3/2}/l_\mathrm{RANS}$ を満たさず
> 誤り。整合する唯一の定義は $l_\mathrm{RANS}=\sqrt{k}/(\beta^{*}\omega)$ (Strelets 2001 /
> Gritskevich 2012)。$\beta^{*1/4}$ を使うと $D_k$ が標準 SST の $\beta^{*-3/4}\approx6.1$ 倍になり、
> $f_d=0$ の付着 BL でも mode1 が標準 SST に縮退せず **MSD で k が崩壊する** (検証で `roK` 100%
> 変化を観測 → $\beta^{*}$ 修正で `roK` 差 $1.8\times10^{-5}$ に解消)。§3.3 の $l_\mathrm{DDES}$ も同様。

### 3.2 LES グリッドスケール

$$
\Delta = \Delta_\mathrm{max} = \max_{\text{隣接面}} |\mathbf{cc}_{ic} - \mathbf{cc}_{jc}|
$$

**Δmax（隣接セル重心距離の最大）を主採用とする。** $\Delta = V^{1/3}$ は等方メッシュなら簡便だが、
**BL の高アスペクト比セル**（壁法線が極薄）では $V^{1/3}$ が接線スケールを大幅に過小評価し、
$l_\mathrm{LES} = C_\mathrm{DES}\Delta$ が小さくなりすぎて、シールド $f_d$ があっても
DES リミッタが BL 内で誤発火するリスクがある。Δmax は Spalart 系で標準。
（実装簡便のため最初に $V^{1/3}$ で動作確認するのは可。ただし最終的に Δmax に置換すること。）

> **逆向きの落とし穴（SBLI せん断層・レビュー指摘）**: Δmax は BL 保護には正しいが、
> **剥離せん断層でスパン方向の粗さに支配される**と $l_\mathrm{LES}$ が大きくなりすぎ、
> $f_d \approx 1$（LES 切替済）でも実効フィルタが粗く渦が出ず、時間平均が URANS 的になる。
> → 失敗予測: 壁近傍は守れるが剥離泡内で渦が立たない。
> 対策: **`delta_les` と $C_\mathrm{DES}\Delta$ vs $l_\mathrm{RANS}$ を必ず出力**し、剥離せん断層で
> $\Delta_x, \Delta_y, \Delta_z$ のどれが律速かを可視化。スパンが律速なら**スパン解像を上げる**
> （メッシュ側の対処）。IDDES では単純 Δmax でなく wall-distance aware な Δ 定義（Shur 2015 の
> $\Delta_\omega$ 等）を比較対象にすると改善余地がある（Phase 2 で検討）。

$$
C_\mathrm{DES} = F_1 \cdot C_{\mathrm{DES},k\omega} + (1 - F_1) \cdot C_{\mathrm{DES},k\varepsilon}
$$

定数：$C_{\mathrm{DES},k\omega} = 0.78$, $C_{\mathrm{DES},k\varepsilon} = 0.61$（Strelets 2001 調整値）

> **F1 の入手について**: `ransSource_d.cu` は F1 を内部で計算しているが（L109–115）、
> `turbulent_viscosity_d.cu`（`sst_eddy_viscosity_d`）は F2 のみで F1 を持たない。
> l_des を `sst_eddy_viscosity_d` 内で計算する場合は F1 を同カーネルで再計算する必要がある
> （k, ω, wall_dist, vis_lam, CD_kw から。CD_kw は ∇k·∇ω が要るので dKd*/dOmegad* も渡す）。
> 簡略化したい場合、初版は $C_\mathrm{DES}=0.78$ 固定（F1 ブレンドなし）でも実用上問題ない。

### 3.3 DDES シールド関数（Spalart et al. 2006）

$$
r_d = \frac{\nu_t + \nu}{\kappa^2 d^2 \sqrt{\partial U_i/\partial x_j \cdot \partial U_i/\partial x_j}}
\qquad (\kappa = 0.41)
$$

$$
f_d = 1 - \tanh\!\bigl([8\, r_d]^3\bigr)
$$

| $r_d$ | 物理的意味 | $f_d$ | 動作 |
|-------|-----------|-------|------|
| $\approx 1$ | BL 内部（渦粘性大） | $\approx 0$ | 純 RANS |
| $\to 0$ | 剥離域・自由せん断層 | $\to 1$ | DES limiter 有効 |

$$
l_\mathrm{DDES} = l_\mathrm{RANS} - f_d \cdot \max(0,\; l_\mathrm{RANS} - C_\mathrm{DES}\Delta)
$$

### 3.4 IDDES length scale（Shur et al. 2008 / SST 版 Gritskevich et al. 2012）

> **確定版 (2026-07-19 Phase 2 着手時に訂正)**: 初稿の
> $l_\mathrm{IDDES}=\max(\tilde f_d(1+f_e)l_\mathrm{RANS},\,l_\mathrm{WMLES})$ /
> $f_e=\max(f_{e1}f_{e2}-1,0)$ は Shur 2008 の正準形と異なっていた（外部要約由来の誤記）。
> 以下は Gritskevich, Garbaruk, Schütze & Menter (2012, Flow Turb. Combust. 88:431) の
> SST-IDDES を正とする。全関数形状は `tools/verify_iddes_functions.py` で
> 検証済（RANS 縮退・WMLES 分岐・float32 ガードを含め VERDICT: PASS）。

**グリッドスケール（壁距離依存）**: Shur 原型の $h_{wn}$（壁法線刻み）項を落とした
Gritskevich 簡略形を採用する。非構造格子で壁法線方向の格子刻みを定義する必要がなく、
`wall_dist` と既存 `delta_les`（$=h_\mathrm{max}=\Delta_\mathrm{max}$）だけで組める:

$$
\Delta_\mathrm{IDDES} = \min\bigl(C_w \max(d_w,\; h_\mathrm{max}),\; h_\mathrm{max}\bigr),
\qquad C_w = 0.15
$$

壁で $0.15\,h_\mathrm{max}$、壁から離れると $h_\mathrm{max}$ に漸近（$[0.15h, h]$ に有界・単調）。
近壁で $\Delta$ が縮むことで $l_\mathrm{LES}$ が下がり、WMLES 分岐が log 層まで LES を入れられる。

**近壁ブレンド関数**（$h_\mathrm{max}=$ `delta_les`）:

$$
\alpha = 0.25 - d_w/h_\mathrm{max},\qquad
f_B = \min\bigl(2\exp(-9\alpha^2),\; 1\bigr)
$$

$f_B=1$（RANS 強制）は $d_w \le 0.53\,h_\mathrm{max}$、$d_w \gtrsim h_\mathrm{max}$ で $\approx 0$。

**log 層ミスマッチ補正** ($r_{dt}$: 乱流版、$r_{dl}$: 層流版の $r_d$):

$$
f_{e1} = \begin{cases} 2\exp(-11.09\alpha^2) & \alpha \ge 0 \\ 2\exp(-9\alpha^2) & \alpha < 0 \end{cases}
\qquad
f_{e2} = 1 - \max(f_t, f_l)
$$

$$
f_t = \tanh\bigl[(C_t^2 r_{dt})^3\bigr],\quad
f_l = \tanh\bigl[(C_l^2 r_{dl})^{10}\bigr],\qquad
C_t = 1.87,\; C_l = 5.0\ (\text{SST 較正値})
$$

$$
f_e = f_{e2}\,\max(f_{e1} - 1,\; 0)
$$

（SA 版にある低 Re 補正 $\Psi$ は SST 版では 1。付着 log 層 $r_{dt}\approx1$ では
$f_t\to1\Rightarrow f_e=0$ で補正は自動停止する。）

**WMLES/DDES 自動切替**:

$$
f_{dt} = 1 - \tanh\bigl[(C_{dt1}\, r_{dt})^3\bigr],\qquad C_{dt1} = 20\ (\text{SST 較正; SA/DDES の 8 と異なる})
$$

$$
\tilde f_d = \max\bigl(1 - f_{dt},\; f_B\bigr)
$$

$$
l_\mathrm{IDDES} = \tilde f_d\,(1+f_e)\,l_\mathrm{RANS} + (1-\tilde f_d)\,l_\mathrm{LES},
\qquad l_\mathrm{LES} = C_\mathrm{DES}\,\Delta_\mathrm{IDDES}
$$

挙動（verify スクリプトで数値確認済）:
- **付着 BL（モデル乱流 $r_{dt}\approx1$）**: $f_{dt}\to0 \Rightarrow \tilde f_d\to1$、かつ $f_e=0$
  → $l_\mathrm{IDDES}=l_\mathrm{RANS}$ 厳密（DDES 相当のシールド）。
- **解像乱流あり（瞬時 $\nu_t$ が SGS 級 → $r_{dt}\ll0.02 \Rightarrow f_{dt}\to1$）**:
  $\tilde f_d \to f_B$ →
  $l_\mathrm{IDDES}=f_B(1+f_e)l_\mathrm{RANS}+(1-f_B)l_\mathrm{LES}$ = WMLES 分岐
  （壁直近 $f_B=1$ は常に RANS、$d_w\gtrsim h_\mathrm{max}$ で LES）。
- 乱流流入がない場合は $f_{dt}\to0$ 側に留まり IDDES は DDES 同等に縮退する。

---

## 4. 実装設計

### 4.1 Config 拡張

`solverConfig.hpp`：既存 `wallTreatmentSST`（L120）の隣に追加。

```cpp
int DESmode = 0;          // 0=RANS(既定, ビット不変) / 1=DDES / 2=IDDES
flow_float C_DES_kw = 0.78;
flow_float C_DES_ke = 0.61;
```

`solverConfig.cpp`：既存 `wallTreatmentSST` のパース（L320–323）に倣い、`turbulence` セクションで
`getOptionalValidatedValue<int>(turb, "DESmode", 0, "turbulence")` で読む。範囲 `0..2` をバリデート。
C_DES_kw/ke も `getOptionalValidatedValue<flow_float>` で既定値付きで読む。

`LESorRANS` は引き続き `2`（SST）のまま。`DESmode` が length scale 修正の on/off を制御。
**`DESmode==0` のとき既存経路と完全ビット一致**であること（分岐で従来式を保持）。

### 4.2 新フィールド

`variables.hpp` の `cellValNames` リスト（L30–）に 2 名を追加するだけ。
`variables()` コンストラクタ（variables.cpp L18–22）と `allocVariables` がこのリストを走査して
host (`c`) / device (`c_d`) を自動確保する。**個別の emplace やコンストラクタ変更は不要。**

```cpp
// variables.hpp cellValNames に追記（src_jac_omega の近く）
"delta_les",   // per-cell グリッドスケール Δmax（幾何セットアップ時に1回計算）
"l_des",       // per-step l_DDES / l_IDDES
"fd_shield",   // per-step シールド関数 f_d（0=RANS, 1=LES limiter 有効）
"rd_des",      // per-step r_d（シールドの素。飽和診断に必須）
```

**診断出力は必須**（レビュー指摘: `l_des` だけでは f_d の誤作動を見逃す。float32 で
`f_d` が 0/1 に張り付くのは NaN より危険）。`output_cellValNames` に
`"l_des"`, `"fd_shield"`, `"rd_des"`, `"delta_les"` を追加し、`res_*.h5` で:

- `f_d` 分布: BL 内 ≈0・剥離域 ≈1 になっているか
- `r_d` 分布: clamp 上限に張り付いていないか（飽和＝シールド誤作動）
- `l_des / l_rans` 比（後処理で算出）: 剥離せん断層でどこが LES 化しているか
- `C_DES·Δ` と `l_rans` の大小関係: LES 化の律速がどの方向の Δ か

を必ず可視化する。これらは backstep / SBLI の合否判定の一次情報。

### 4.3 グリッドスケール Δ の事前計算

**幾何セットアップ時に device で 1 回計算**する。場所は `variables.cpp` の幾何セットアップ
（`setGeometricVariables` 相当、volume/ccx を device へ転送する L256–424 付近）の直後が自然。

Δmax は面ループで隣接重心距離の最大を取る:

```cpp
// 擬似コード（atomicMax か、面→セルの既存ループ機構を流用）
// 初期化: delta_les[ic] = 0
// 内部面 ip ごとに:
//   d = |cc[ic] - cc[jc]|
//   atomicMax(&delta_les[ic], d);  atomicMax(&delta_les[jc], d);
// 境界面: d = 2*|cc[ic] - face_center| で補完
```

forge には面→セルの集約カーネルが既にある（勾配計算等）ので、その実装パターンを流用する。
**初版は簡便に `delta_les[ic] = cbrt(volume[ic])` で動作確認 → Δmax に置換**でも可（§3.2 の注意参照）。

> atomicMax が float に直接効かない場合は、`__int_as_float`/`__float_as_int` トリック、
> または CPU 側（mesh 構築時に host で計算してから H2D 転送）でも良い。
> 静的量なので毎 step 計算は不要、1 回で十分。

### 4.4 DDES カーネル（Phase 1）

`turbulent_viscosity_d.cu` に新 device 関数を追加：

```cpp
__device__ flow_float compute_rd(
    flow_float nu_t, flow_float nu, flow_float y,
    flow_float dUxdx, ...) {
    // |∂Ui/∂xj|^2 = S^2 + Ω^2 (9 成分の二乗和 = velocity gradient magnitude)
    flow_float grad2 = dUxdx*dUxdx + dUxdy*dUxdy + ... ; // 全 9 成分
    // ★ float32 ガード（レビュー指摘・最重要）:
    //   分母は y² で効くので近壁・よどみ・一様流で極端値になりやすい。
    //   分母に「物理スケール付き floor」を入れ、grad2 にも floor を入れる。
    const flow_float gradMag = sqrt(max(grad2, GRAD_FLOOR)); // GRAD_FLOOR: ケールスケール由来
    flow_float denom = kKarman * kKarman * y * y * gradMag;
    denom = max(denom, DENOM_FLOOR);
    flow_float rd = (nu_t + nu) / denom;
    // ★ rd を 0..10 程度に clamp（飽和で f_d が 0/1 に張り付くのを防ぐ）。
    return min(rd, static_cast<flow_float>(10.0));
}

__device__ flow_float compute_fd(flow_float rd) {
    flow_float arg = static_cast<flow_float>(8.0) * rd;
    return static_cast<flow_float>(1.0) - tanh(arg * arg * arg);
}
```

> **float32 ガードの設計指針（レビュー反映）**:
> forge は近軸 SST で float32 起因の問題を経験済み（mixed-precision は棄却）。r_d は d² で
> 効くため cell モードでも飽和しやすく、node モード（壁ノード d→0）ではさらに顕著。
> - `grad2` の floor: 一様流（∇u→0）で r_d→∞→f_d=0（RANS）になる。これは**正しい極限**だが、
>   0 割回避に floor が要る。floor は流れの代表ひずみ率スケール（例 `(U∞/L)²`）で無次元的に置く。
> - `rd` clamp [0, 10]: f_d=1-tanh((8rd)³) は rd≳0.5 で実質 0 に飽和するので、上限 10 で十分。
> - **`f_d` と `r_d` を必ず出力**（§4.2）。NaN がなくても「どこでも f_d=0/1 張り付き」は
>   DES が機能していない兆候で、l_des だけでは見抜けない。

`sst_eddy_viscosity_d` または新カーネル `compute_des_length_d` で：

```cpp
float l_rans = sqrt(k_c) / (pow(kBetaStar, 0.25f) * w_c);
float l_les  = C_DES * delta_les[ic];
float rd     = compute_rd(...);
float fd     = compute_fd(rd);
float l_ddes = l_rans - fd * max(0.0f, l_rans - l_les);
l_des[ic]    = max(l_ddes, 1e-20f);
```

### 4.5 k 消滅項の切替（ransSource_d.cu）

`rans_sst_source_d` カーネル（既存 Dk は L140）に引数 `int DESmode` と `flow_float* l_des` を追加。
`ransSource_d_wrapper`（L183–206）の呼び出しに `cfg.DESmode, var.c_d["l_des"]` を追加。

```cpp
// 既存（L140）: const flow_float Dk = kBetaStar * rho * k_c * w_c;
// → 分岐に置換:
flow_float Dk;
flow_float jac_k;
if (DESmode > 0) {
    const flow_float l_hyb = max(l_des[ic], static_cast<flow_float>(1e-20));
    Dk    = rho * pow(k_c, static_cast<flow_float>(1.5)) / l_hyb;
    // 陰解法ヤコビアン: ∂Dk/∂(ρk) = 1.5 √k / l_hyb
    jac_k = static_cast<flow_float>(1.5) * sqrt(k_c) / l_hyb;
} else {
    Dk    = kBetaStar * rho * k_c * w_c;       // 従来式（ビット不変）
    jac_k = kBetaStar * w_c;
}
```

**ω 方程式は変更しない**（標準 SST-DES は k 消滅項のみ修正）。`src_jac_k` には上記 `jac_k` を入れる
（既存 L162 `src_jac_k[ic] = kBetaStar * w_c;` を `src_jac_k[ic] = jac_k;` に）。
automatic wall treatment の wf_pk 置換（L135–137）と ω ピン留め（L169–172）は DES と独立に
そのまま残す（wall-adjacent セルは RANS 扱いで問題ない）。

### 4.6 IDDES 追加関数（Phase 2・確定設計 2026-07-19）

**§3.4 確定版（Gritskevich 2012 完全関数群）をそのまま実装する。簡略版（初稿の
`max(f̃_d(1+f_e)l_RANS, l_WMLES)` 形）は廃止**（式自体が誤記だった。§3.4 冒頭注）。

実装は既存 `compute_des_length_d`（`turbulent_viscosity_d.cu`）に `DESmode==2` 分岐を追加する:

- 入力は既存 DDES と同一（`vis_turb`/`vis_lam`/`wall_dist`/`delta_les`/∇u/k/ω）。
  **新規の幾何量・事前計算は不要**（$\Delta_\mathrm{IDDES}$ が Gritskevich 簡略形のため。§3.4）。
- $r_{dt}=\nu_t/(\kappa^2 d^2|\nabla u|)$、$r_{dl}=\nu/(\kappa^2 d^2|\nabla u|)$ は既存
  `compute_rd_ddes` を流用（floor + clamp[0,10] の float32 ガード込み。$f_t/f_l/f_{dt}$ の
  tanh 引数は clamp 上限 10 で完全飽和するので clamp の影響なし）。
- 定数: $C_{dt1}=20$, $C_t=1.87$, $C_l=5.0$, $C_w=0.15$（SST 較正値・compile-time 定数。
  config 化しない）。$C_\mathrm{DES}$ は DDES と同じ `C_DES_kw` 固定（F1 ブレンドは将来課題）。
- **DESmode:1 (DDES) 経路はビット不変**: `f_d` の $C_{d1}=8$（Spalart 2006 較正）を含め一切触らない。
  IDDES の $f_{dt}$ は $C_{dt1}=20$ の別関数として実装する（Gritskevich も DDES=20/IDDES=20 だが
  forge の DDES は 8 で検証済のため変更しない。挙動差は backstep 回帰で観察）。
- **診断出力**: `fd_shield` には $\tilde f_d$、`rd_des` には $r_{dt}$ を書く（DDES と同じ配列を
  モードで意味替え）。$f_e$ は新フィールド `fe_iddes` に出力（WMLES 分岐と log-layer mismatch
  補正の活性を見る一次情報）。$f_B$・$\Delta_\mathrm{IDDES}$ は `wall_dist`/`delta_les` から
  後処理で厳密再計算できるためフィールドは追加しない。
- 事前検証: `tools/verify_iddes_functions.py`（関数形状・RANS 縮退・WMLES 分岐・float32 整合、
  VERDICT: PASS 済 2026-07-19）。

> **⚠ 定量評価の前提（レビュー指摘・維持）**: WMLES モードの定量評価
> （チャネル log 層 / mismatch）には解像乱流の存在が前提。周期チャネルは流入不要で
> 自己維持するが **streamwise 体積力 config が未実装**（[turbulence-wmles-wall-stress plan](turbulence-wmles-wall-stress.md)
> §5-7 と共通の前提部品。同 plan は純 LES + 代数壁応力の別経路で IDDES と共存可）。
> 非周期ケース（SBLI 等）の IDDES 定量評価は **流入乱流生成（別 plan）にゲート**。
> それまで定量目的は DDES に留める（機能確認 = f̃_d/f_e/l_des 分布の妥当性までを Phase 2 の
> 完了範囲とする）。

### 4.7 呼び出しフロー（main loop）— ★重要：依存順

**`r_d` は $\nu_t = $ `vis_turb`/ρ を必要とする。`vis_turb` は `turbulent_viscosity_d_wrapper`
（main.cpp:937, `sst_eddy_viscosity_d`）で初めて計算される。** したがって l_des の計算は
**`turbulent_viscosity` の後・`ransSource` の前**でなければならない。

実際の main.cpp 呼び出し順（確認済）:

```
main.cpp:920  calcGradient_d_wrapper        … ∂Ui/∂xj 更新（gradients ready）
main.cpp:937  turbulent_viscosity_d_wrapper … vis_turb 計算（nu_t ready）
main.cpp:943  ransTransport_d_wrapper        … k/ω 移流・拡散
main.cpp:954  ransSource_d_wrapper           … Dk 計算（ここで l_des を読む）
```

**推奨実装**: `l_des` の計算を `turbulent_viscosity_d_wrapper` の**末尾に追加**する
（vis_turb・wall_dist・gradients・delta_les が全て揃っており、追加カーネル呼び出しが
1 wrapper 内で完結する）。具体的には:

- 案A（最小）: `sst_eddy_viscosity_d` カーネル内で vis_turb を計算した直後に、同じローカル
  `mu_t` を使って l_des もその場で書き込む（追加の global read 不要・最速）。
  ただし F1 が要るので k/ω 勾配ポインタを追加引数で渡す（§3.2 の注記）。
- 案B（分離）: `turbulent_viscosity_d_wrapper` 内で `sst_eddy_viscosity_d` の後に
  新カーネル `compute_des_length_d` を呼ぶ（vis_turb を global から読む）。可読性優先。

`DESmode==0` のときは l_des カーネルを呼ばない（または l_des=l_rans を書く）。`ransSource` 側で
`DESmode==0` なら従来 Dk、`>0` なら l_des 経由 Dk に分岐する。

> **注意（軸対称はスコープ外）**: 既存 SST は軸対称フープひずみ `axisym_uy_over_r` を S² に
> 加えるが、DES 検証ケース（backstep・SBLI）は全て 3D 直交メッシュなので、初版の r_d /
> gradient magnitude に軸対称項は含めない。軸対称 DES は将来の別 plan とする。

### 4.8 DES 用低散逸 flux（Phase 1.5 — SBLI 物理検証の前提）★

> **★ 設計更新 (2026-07-19 後段, ユーザー承認済み・現行の本命)**: Phase 1.5 は flux ブレンドを
> やめ、**全域単一スキーム `KEEP` + matrix ES 散逸 (`keepDissType:2, keepDissJump:2,
> keepDissPrecond:1, σ0.05`) を第一候補**とする。根拠:
> (1) RANS が upwind を要求するわけではなく、中心+人工散逸の RANS は JST (SU2 既定の一つ) として
> 主流実績がある — KEEP+matrix ES は構造的に JST 同族 (σ≈k₄)、
> (2) `keepDissPrecond` (2026-07-19 実装, [convection-keep-diss-lowmach-precond](../active/convection-keep-diss-lowmach-precond.md))
> により低マッハ市松は行列構造レベルで解決済み (真の市松 142×減衰・物理コスト c' 以下・ES 維持)、
> (3) 単一スキームなら f_d ブレンドの実装・グレーゾーン不整合・LES 域散逸の f_d 依存が全て消える。
> **前提検証 (この順)**: ① flat plate SST を KEEP+ES(+implicit) で回し BL プロファイル/Cf が
> SLAU・SU2 と一致すること (KEEP+implicit は未検証領域)、② backstep DDES の再実行で NaN 無し・
> f_d 分布維持、③ SBLI 系へ伸ばす時点で**衝撃センサ駆動の σ ブースト (JST k₂ の ES 版, Ducros×σ
> ランプ) を別 plan で追加** (現 ES レイヤは shock capturing を持たない)。
> ①② が躓いた場合のフォールバックが下記の f_d 主導ブレンド (旧本命・設計は有効なまま残す)。
>
> **前提検証の結果 (2026-07-19 後段, ①②とも合格)**:
> ① flat plate SST (case/26 `run_0012`/`run_0013`): KEEP+ES+implicit cfl20 で安定・Cf 完全収束
> (ドリフト ≤0.04%/20k)・**Cf/Schl 0.88/0.91/0.94 = SLAU 基準比 −3.6〜−3.9%** (本 case の
> スキーム間ばらつき node vs cell MUSCL −3.4〜−3.9% と同帯) → **KEEP+ES の SST-RANS は成立**。
> 注: 最初の試行 (run_0004/0005) は Cf 未検証の run_0003 土俵で不適だった。
> ② backstep DDES (case/18 `run_0125`): 3D 857k cells + implicit で 500 step NaN 無し、
> f_d = せん断層 0.974 (SLAU 0.976)・泡内 0.000 と zoning 正常。付着 BL 帯の f_d が SLAU より
> やや高い (0.43 vs 0.28, 遮蔽は保持) — Phase 2 で νt 依存性として要観察。
> → **単一スキーム案で Phase 1.5 は実質完了。次は Phase 2 (IDDES 関数群) に進める。**

> **⚠ 整合メモ (2026-07-19・Phase 1.5 着手前に必読)**: 本節の本命設計「`solver: KEEP_SLAU` の
> `duc` に f_d 注入」は**現行コードでは前提が消失している**。`KEEP_SLAU_d` は
> `cuda_forge/convection/legacy/` に退避され dispatch から到達不能 (現行 `cfg.solver` は
> SLAU/SLAU2/HLLE/ROE/KEEP のみ)。KEEP 系はその後「純 KEEP 中心流束 + 独立な ES 散逸レイヤ
> (`keepDissType` / `keepDissJump` / σ)」に再編された
> ([`convection-keep-es-dissipation`](../accepted/convection-keep-es-dissipation.md) /
> [`convection-keep-diss-recon-jump`](../accepted/convection-keep-diss-recon-jump.md)、σ 較正は
> memory [[keep-es-dissipation-status]])。Phase 1.5 の実装時は KEEP_SLAU 復活ではなく、
> **ES 散逸レイヤの散逸強度 (σ または blend) を f_d で駆動する**設計に読み替えて再設計し、
> 先に本節を更新してから実装すること (LES 域 σ≈0.02 ↔ RANS/衝撃域フル散逸、という対応は
> 本節の設計原理 (§survey の Travin 型連続ブレンド・f_d 主導) をそのまま踏襲できる)。
> 以降の本文は 2026-06-21 時点の設計記録として残す。

> **✅ 再設計・実装 (2026-07-22, 上記メモの指示に従う現行設計)**: ES 散逸レイヤの σ を f̃_d で駆動する。
>
> $$\sigma_f = \max\!\big(\sigma_\mathrm{min},\ (1 - \tilde f_{d,f})\,\sigma_\mathrm{max}\big),\qquad \tilde f_{d,f} = \tfrac12\big(f_d[ic_0] + f_d[ic_1]\big)$$
>
> - **config**: `keepDissFdBlend` (0: off 既定=ビット不変 / 1: on)、`keepDissCoeffMax` (既定 1.0)。
>   σ_min は既存 `keepDissCoeff` (case/39 A/B で 0.05 採用 — 一様 σ0.05 は解像乱流無傷で格子スケール
>   滑らかさ 10-25% 改善、`run_0009_coarse_s005`)。SU2 `FD` (σ=max(0.05,1−f_d)) と同形。
> - **実装箇所**: `convectiveFlux_keep_d.inc.cuh` の散逸 3 分岐で per-face σ を計算 (`fd_shield` は
>   `turbulent_viscosity_d` が convectiveFlux より先に毎サブ反復更新するため鮮度整合)。
>   ゲート: `DESmode>0 && keepDissFdBlend==1` のときのみ `fd_shield` を渡す (他は nullptr=従来)。
> - **効果の狙い**: RANS 帯 (f̃_d≈0, 付着 BL・薄セル) で σ→1 のフル ES 散逸 = 他コードの風上相当、
>   LES 域 (f̃_d≈1) はフロア σ_min のみ。UZEN/Probst (arXiv:2511.08257, KE 保存 LD2 + f̃_d Travin σ)
>   の出版済みパターンと同型 ([survey §7](../../notes/investigations/turbulence-des-flux-survey.md))。
> - 旧設計の `desLowDissipation` config 名は廃止 (KEEP_SLAU 前提だったため)。中期の改良 Ducros OR
>   ブレンドは本形の σ_f に対する OR (σ_f ← σ_Ducros + σ_f − σ_Ducros·σ_f) として据え置き互換。
> - **検証済 (2026-07-22, case/39 粗格子 3 腕比較・同一 IC/同一窓/5 snap 平均)**: fdblend は
>   風下壁帯 (RANS/grey 帯) の市松度 r −0.50→−0.79・d2 −37% (vs σ0.02 一様)、コア Uz rms
>   6.63±1.2 vs 6.89±1.1 = 解像乱流無傷。**採用** (`run_0010_coarse_fdblend`)。
> - **制約 (実測で確定)**: `keepDissPrecond` と非互換 — Δp 散逸増強 (∝c²/Ur) は壁近傍で
>   Ur=εc フロアに張り付き、σ→1 との積が即発散 (σmax0.3 でも step2 NaN / precond off は
>   σmax1.0 安定, `run_diag_fdb_*`)。config 検証で禁止し `keepDissCprime` へ退避。
>   LES 側 precond + RANS 側標準 Roe の**演算子ブレンド**は将来課題。
> - **fd_shield の向きの罠**: DDES (f_d, 1=LES) と IDDES (f̃_d, 1=RANS) で逆
>   (variables.hpp)。カーネルは `fdLesOne` フラグで両対応済み。

**外部レビューが最大の失敗要因と指摘した箇所。** DDES length scale が正しくても、
LES 領域（f_d≈1）で対流スキームの数値散逸が強いと解像乱流が育たず、SBLI で
「剥離せん断層の KH 成長が遅れ、剥離泡が短く、非定常圧力変動が弱い」結果になる。

> **設計根拠の一次調査**: 本節の設計は [`turbulence-des-flux-survey.md`](../../notes/investigations/turbulence-des-flux-survey.md)
> （deep-research 98 エージェント / SU2 master ソース + Travin/DEShybrid/Ducros 一次文献）に基づく。
> 要点: (1) 正準形は `weight=(1-σ)w₁+σw₂` の中心↔風上連続ブレンド（Travin 2000 / OpenFOAM
> `DEShybrid` / SU2 `ROE_LOW_DISSIPATION`）、(2) **ブレンド係数の駆動量は乱流解像物理（f_d/NTS）で、
> 衝撃センサとは直交**、(3) SU2 の DES 既定は `ROE_LOW_DISSIPATION=FD`（σ_FD=max(0.05,1-f_d)）で
> **forge の `1-f_d` 注入と一致**。**ただし MUSCL limiter を衝撃センサに流用する設計は一次 DES/LES
> 文献に一例も無い**（SU2 は実衝撃に専用 Ducros を `FD_DUCROS`/`NTS_DUCROS` で層化）。

**現状の forge（2026-06-22 前提更新）**:
- `solver: SLAU/SLAU2/ROE/HLLE` は全面風上（散逸はリーマン解に内在）。
- **Ducros センサは現状排除**。`apply_ducros_limiter` / `ducrosSensor_d.*` は残置するが、SLAU では
  1 次化 ON/OFF で場の差 <0.1%（SLAU 自身の散逸が支配）で**センサとして当てにならない**ことが
  判明し `ducrosLimiter` 既定 0。**この観測は「limiter ψ も衝撃センサに使えない」ことの forge 内部証拠**
  でもある（ψ も SLAU 散逸下では不感）→ 短期は f_d 主導、中期で**改良 Ducros を別途復活**（下記）。
- **重要: forge は既に「中心基本＋センサで風上」を実装している**。`solver: KEEP_SLAU`
  ([`KEEP_SLAU_d`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu), L2210–) は
  $\mathbf{F} = (1-\mathrm{duc})\mathbf{F}_\mathrm{KEEP} + \mathrm{duc}\,\mathbf{F}_\mathrm{SLAU}$,
  $\mathrm{duc}=\max(\text{ducros},1-\psi_\mathrm{lim})$ の blend（L2452–2456）。$\mathrm{duc}=0$ で
  **KEEP（運動エネルギー保存中心差分・無散逸）**、$\mathrm{duc}=1$ で SLAU（風上）。圧縮性 LES/DES の
  定石（KE 保存中心を基本に、不連続のみ風上へセンサ切替）そのもの。

**短期設計（Phase 1.5 本命: KEEP_SLAU の blend に f_d を主導注入）**:

LES は KE 保存中心差分を基本にして渦を保存し、RANS でのみ風上にするのが正しい。
σ で SLAU の散逸を薄める（旧案＝下記案B）より、**既存 `KEEP_SLAU` の blend 係数 `duc` を `1-f_d`
主導にする**方が KE 保存性で優れ新規コードもほぼ不要。**SU2 `FD`（σ_FD=max(0.05,1-f_d)）と一致**：

```cpp
// 旧: duc = max(ducros, 1-lim)
// 新 (DESmode>0 かつ desLowDissipation=1):
flow_float fd_face = 0.5*(fd_shield[ic0] + fd_shield[ic1]);  // f_d を面へ補間
flow_float duc_des = fmaxf(0.05f, 1.0f - fd_face);           // ★ DES 主制御。SU2 FD と同形・フロア 0.05
duc = fmaxf(duc_des, 1.0f - lim);   // 1-lim は衝撃センサではなく positivity/robustness フロア扱い
```

- **`1-f_d` が DES 主制御**（SU2 `FD` 同形, フロア 0.05 で LES 域の残留散逸を最小化）。
- **`1-ψ_limiter` は衝撃センサにしない**。limiter は急勾配（せん断層・乱流極値・接触面）全部で発火し
  衝撃と乱流を区別しないため、衝撃センサに使うと**守りたい解像乱流を過散逸させる**（一次文献は専用
  Ducros 等を使い、limiter のセンサ流用例は皆無）。`max(…,1-lim)` 項は **positivity/robustness の
  フロア**としてのみ残す（実際 SLAU 下では <0.1% 不感なので実質 `duc≈1-f_d` に縮退する公算大）。

| 領域 | f_d | duc（短期） | 実効スキーム |
|------|-----|-------------|-------------|
| 付着 BL（RANS） | ≈0 | **≈1** | SLAU（風上・頑健、モデルが効く） |
| 剥離 LES（滑らか） | ≈1 | **0.05** | **KEEP（中心・低散逸 → 渦が生きる）** |

backstep は**亜音速・衝撃波なし**なので短期はこの f_d 単独 blend で十分（衝撃センサ不要）。

**中期設計（SBLI/ノズル/ピントル: 改良 Ducros を OR ブレンドで追加＝SU2 `NTS_DUCROS` 相当）**:

実衝撃と剥離せん断層が**共存**する forge の本命ターゲット（SBLI・超音速ノズル・ピントル）では、
f_d 単独では「LES 域に入った衝撃で散逸が落ちて捕獲できない」。SU2 はこのため BSCW 遷音速バフェット
DDES で `NTS_DUCROS` を採用し、**壁と衝撃で σ=1、せん断層/wake で σ→0.05** にしている。forge も
同じ **確率OR ブレンド**で改良 Ducros を層化する：

$$
\sigma_\mathrm{Ducros} = \frac{(\nabla\cdot\mathbf{u})^2}{(\nabla\cdot\mathbf{u})^2+\omega^2+\epsilon},\qquad
\mathrm{duc} = \sigma_\mathrm{Ducros} + (1-f_d) - \sigma_\mathrm{Ducros}\,(1-f_d)
$$

- **Ducros は原型でなく改良版**（Hendrickson-Kartha-Candler 2018）を使う。原 Ducros は低渦度域で
  発散に過敏で解を改変する失敗モードがある（[survey §2.2](../../notes/investigations/turbulence-des-flux-survey.md)）。
- これにより付着 BL（f_d→0）でフル散逸・衝撃でフル散逸・剥離せん断層でのみ低散逸、を同時成立。
- 既存の `apply_ducros_limiter` を流用せず、**DES 用 `σ_Ducros`（連続値）を別途実装**（既存は離散リミタ）。

**代替案（本命が不調なら）**:
- 案B: SLAU 散逸項に係数 $\sigma\in[\sigma_\mathrm{min},1]$。KEEP を使わず最小改修だが KE 非保存で渦保存に劣る。
- 案D（WENO/TENO）/ 案E（LD2 勾配再構成, DLR J054956）は長期。非構造FV・GPU・float32 と非親和な D は
  見送り、E は LES 枝の精度向上策として長期検討。

**config**: `desLowDissipation`（0:off 既定 / 1:f_d のみ短期 / 2:f_d+改良Ducros 中期）を新設。
**前提: `solver: KEEP_SLAU`**（SLAU 単独では KEEP 経路が無いので無効）。Phase 1 機能検証は off で実施済、
Phase 1.5 で 1（f_d）→ T3 で 2（+Ducros）と段階的に上げる。

**判定**: backstep でこの機構 on/off の解像乱流（速度変動 RMS・スペクトル傾き）を比較し、
on で −5/3 慣性域が伸びること。off では「f_d は立つが乱流が死ぬ」灰色領域になりやすい。
SBLI（T3）では `desLowDissipation=2` で衝撃捕獲とせん断層低散逸の両立を SU2 クロスチェックで確認。

**float32/GPU リスク**: σ ブレンドは局所量から face ごとに計算する有界・平滑な乗数（フロア 0.05）で
GPU 並列・float32 安定。**真のリスクはブレンドでなく既知の非直交 free-stream 桁落ち**（中心化で散逸が
減るほど顕在化）。LES 域の中心比率を上げる前に free-stream preservation を確認すること。

> 触るファイル: `convectiveFlux_d.cu`（`KEEP_SLAU_d` の `duc` に f_d 注入・`fd_shield` 引数追加、
> 中期で連続 `σ_Ducros` 追加）、config（`desLowDissipation`）。**§7 影響範囲の「影響なし: convectiveFlux」
> は Phase 1.5 で覆る**。既存 `ducrosSensor_d.*`（離散リミタ）は流用せず、DES 用連続 σ_Ducros を別実装。

---

## 5. 実装ステップ

1. **Config 追加** (`solverConfig.hpp/.cpp`)  
   `DESmode`, `C_DES_kw`, `C_DES_ke` を追加・解析。ビルド確認。

2. **フィールド追加** (`variables.hpp/.cpp`)  
   `delta_les`, `l_des` を host/device 両方に追加。

3. **Δ の事前計算**（`mesh/mesh.hpp` または `setInitial.hpp`）  
   全面ループで per-cell max 距離を計算、`delta_les` に書込み。

4. **DDES カーネル実装** (`turbulent_viscosity_d.cu`)  
   `compute_des_length_d` カーネルを新規追加（device 関数 `compute_rd`, `compute_fd`, `compute_l_ddes` を用いる）。  
   `turbulent_viscosity_d_wrapper` に `DESmode > 0` 分岐で呼び出し。

5. **k 消滅項切替** (`ransSource_d.cu`)  
   `rans_sst_source_d` に `l_des` ポインタと `DESmode` を追加引数として渡す。  
   `Dk` の計算を分岐。`src_jac_k` も対応するヤコビアンに変更。

6. **ビルド＋回帰確認**  
   `DESmode: 0` で `case/05.sod_shock_tube` と `case/26.flat_plate_sst` がビット不変であることを確認。

7. **DDES 動作確認** (`case/18.backstep`)  
   非定常計算（`timeIntegration: 3`）で剥離域に $f_d \approx 1$ (LES 活性) ・
   付着域に $f_d \approx 0$ (RANS 保持) を `l_des` フィールドの可視化で確認。

8. **IDDES Phase 2**（Phase 1 完了後）  
   `f_B`, `f_e` 関数追加、`l_IDDES` 切替。DDES との l_des 差分を同ケースで比較。

---

## 5.5 実装上の落とし穴（精査で判明・必読）

1. **依存順（最重要）**: `r_d` は `vis_turb`（=ρν_t）を読む。`vis_turb` は main.cpp:937 の
   `turbulent_viscosity_d_wrapper` で計算される。l_des はその**後**で計算すること（§4.7）。
   うっかり calcGradient 直後に置くと vis_turb が前 step の古い値になる（致命的でないが不整合）。

2. **ビット不変の死守**: `DESmode==0` は既存 RANS と完全一致が要件。Dk・src_jac_k の分岐は
   `DESmode==0` 側で**従来式そのまま**（係数・順序・キャスト一致）。検証は
   `case/26.flat_plate_sst` を新旧バイナリで回し `res_*.h5` の差が 0 であることで確認。

3. **stale build に注意**（既知の罠）: `solverConfig.hpp` に新メンバを足すと差分ビルドが
   CUDA obj を取りこぼし step0 で dt=0/NaN 凍結することがある。**config 変更後は full rebuild** すること
   （プロジェクトの既知事例。`make clean` 相当を挟む）。

4. **Δ の選択**: BL の高 AR セルで `V^(1/3)` は接線スケールを過小評価 → DES リミッタ誤発火。
   最終的に Δmax を使う（§3.2）。初版 V^(1/3) で動作確認したら必ず Δmax に置換し、
   `case/26.flat_plate_sst` で BL 内 f_d<0.05 を再確認。

5. **F1 ブレンド**: `sst_eddy_viscosity_d` は F1 を持たない（F2 のみ）。l_des を同カーネルで計算するなら
   F1 を再計算するか、初版は `C_DES=0.78` 固定で割り切る（§3.2）。

6. **軸対称はスコープ外**: 検証ケースは 3D 直交。r_d の gradient magnitude に
   `axisym_uy_over_r` は含めない。

7. **3D メッシュ確認**: `case/18.backstep` の既存メッシュが 2D か 3D（スパン方向セルあり）かを
   実装前に確認。DDES はスパン方向解像がないと LES モードが立たず RANS 同然になる。
   2D なら 3D メッシュを別途用意（スパン幅・周期 BC）。

---

## 5.6 node-centered（median-dual）両対応

forge は cell-centered 既定だが `discretization: cell | node` で node-centered（中点双対
median-dual）に切替できる両対応化が進行中（[`discretization-median-dual.md`](discretization-median-dual.md)）。
**DDES は本質的に per-CV 代数（volume→Δ, wall_dist→d, |∇u|→r_d, vis_turb→ν_t を読んで l_des を吐く）
なので、median-dual の設計とそのまま整合し、DES カーネルに cell/node 分岐は不要**。

### なぜ無変更で成立するか

median-dual は「GPU カーネルは CV+面+接続の抽象だけで書かれ、内部スキームは双対上で無変更で成立」
を前提に、前処理で `replacePrimalWithDual()` が primal を dual に差し替え **solver 本体は無変更**で
走る方式（median-dual plan M1 で確定済）。node モードでは `nCells==nNodes`、
`volume`/`ccx`/`wall_dist`/勾配/`vis_turb` が全て「双対 CV 上の量」になるだけ。
DES はこれらを読むだけなので双対 CV 上でも同じ式が成立する。

### 両対応のための実装規約（cell 実装時から守る）

- **Δ は実行時の `volume`/`ccx`（device 配列）から計算する**。gmsh の primal
  `msh.cells[].volume` を直接参照しない。これで node モードでは自動的に双対 CV の Δ になる
  （§4.3 の方針はこれを満たす）。
- DES カーネルは `cfg.discretization` を**参照しない**。CV 抽象（volume/wall_dist/勾配/vis_turb）
  のみに依存させる。

### node モード固有の注意（検証時）

1. **近壁勾配の checkerboard（最重要）**: node-centered median-dual は近壁 Green-Gauss 面勾配が
   checkerboard を持つ（`gradLSQ`・`dualEdgeMidInterp` フラグが存在する理由）。**r_d は ∇u に依存
   するため、勾配振動が f_d シールドを誤作動させ得る**（BL を守れず MSD で偽剥離、または過剰
   シールド）。→ **node モードで DES を回すときは `gradLSQ=1` を前提**とする。
2. **壁ノードで d→0**: node では壁ノードが壁面に乗り `wall_dist=0`。r_d 分母の d² → 壁で
   r_d→∞→f_d→0（純 RANS）= 物理的に正しい極限だが、float32 のゼロ割回避に `max(d, kSmall)`
   ガードが cell 以上に効く（§4.4 に既にある）。関連: [`diffusion-node-wall-viscous-distance.md`](../accepted/diffusion-node-wall-viscous-distance.md)。
3. **実効解像度の差**: 同一物理メッシュでも node 双対は DOF が異なる（tet で ≈1/5.5）。Δ が変わり
   DES 解像内容も変わる。cell vs node の DES 比較は「同一メッシュ」でなく「同一実効解像度」で揃える。

### 順序依存

**node モード DES は node モード SST の成熟度に乗る**。median-dual は現在 M1（Euler 1 次が回る段階）
で RANS/SST の node 検証はこれから。現実的な順序は:

```
本 plan: cell DDES 実装・検証（DES カーネルを CV 抽象に対して書く）
  ↓
node SST が検証通過（median-dual 側の別 plan）
  ↓
node DDES = 検証マトリクスに 1 行追加するだけ（コード流用・gradLSQ=1 で回す）
```

→ 今回は cell で実装すれば足りる。CV 抽象に対して書いておけば node 化は後から検証のみで済む。

---

## 5.7 失敗予測（外部レビュー 2026-06-21・実装者は必読）

**このまま Phase 1 のみで Settles 24° に行くとどうなるか**:
> 衝撃波は捕まるが、剥離せん断層の乱流化が遅れ、**壁圧プラトーが実験より平滑・
> 非定常圧力変動が弱い**結果になりそう。原因は DDES length scale ではなく、
> **SLAU/Roe 系の風上散逸と流入乱流不足**。

実装者が「壁圧がそれっぽく合った」で合格にしないための注意:
- 壁圧は衝撃捕獲＋平均剥離泡で**それっぽく合うことがある**（偽の合格）。
  Cf・剥離/再付着位置・圧力 RMS・スパン相関まで見ないと「乱流が立っている」とは言えない。
- `f_d > 0.9` だけでも不十分。**f_d が立っても風上散逸で乱流が死ぬ**から（§4.8）。
- 「解像乱流もモデル乱流も足りない灰色領域」が最悪。Dk だけ強まり Pk は RANS 上限
  （[ransSource_d.cu:121](../../solver_density_cuda/cuda_forge/ransSource_d.cu) の limiter）に縛られ、
  モデル k が急落する一方、数値散逸が解像乱流を消す。

**優先対処の結論（レビュー）**: Phase 1 DDES は進めてよい。ただし**物理検証の前**に、
①診断出力（f_d, r_d, delta_les, l_des/l_rans・§4.2）②3D backstep の乱流統計確認（§6 T1-B 強化）
③DES 用低散逸 flux（§4.8）を揃えること。SBLI 定量一致を目標にするなら低散逸ブレンドは実質必須寄り。

---

## 6. 検証

### 全体方針

DDES/IDDES の検証は **4 段階のティア**で行う。
ティアが上がるほど物理的に目標適用に近くなるが、計算コストも上がる。
**ティア 1–2 を先に完走させてから、ティア 3–4 に進む**こと。

---

### ティア 1：回帰・シールド関数の正常動作確認（既存ケース・安価）

#### T1-A: `case/26.flat_plate_sst`  — BL シールド確認

**目的**: 付着乱流 BL 全域で $f_d \approx 0$（DDES が RANS のまま）であることを確認。
DDES の最重要保証「BL を壊さない」をまず単純な流れで検証する。

**参照データ**: 既存 SST 収束済み run（`DESmode: 0`）

**手順**:
1. 既存収束済み場を restart して `DESmode: 1` に切替え
2. `l_des` フィールドを出力し $l_\mathrm{des} / l_\mathrm{RANS}$ マップを確認
3. BL 内 $l_\mathrm{des}/l_\mathrm{RANS} > 0.99$ であること（LES リミッタが効いていない）
4. Cf プロファイルが RANS と < 0.1% の差であること

**判定基準**:
- [ ] BL 全域で $f_d < 0.05$
- [ ] Cf の RANS との差 < 0.1%
- [ ] NaN/Inf なし

---

#### T1-B: `case/18.backstep` — 剥離域 LES 活性確認

**目的**: 後退段で剥離した直後に $f_d \to 1$（LES 活性）となり、
再付着長が RANS より実験値に近づくことを確認。

**Primary validation: NASA/TMR 2DBFS = Driver & Seegmiller (1985)**

TMR が validation case として整理し、**実験データファイルが落とせる**（ERCOFTAC Classic
Collection C.30 としても流通）。原典: Driver & Seegmiller, *AIAA Journal* 1985。
入手: NASA Turbulence Modeling Resource (tmbwg.github.io) の 2DBFS ページ。

| 量 | 値 |
|----|-----|
| $Re_H$（ステップ高さ基準） | ≈ 36,000 |
| upstream $Re_\theta$ | 5,000 |
| upstream BL 厚 $\delta$ | ≈ 1.5 H |
| マッハ数 $M$ | ≈ 0.128（**ほぼ非圧縮**。低マッハ前処理 `lowMachPrecond` 検討） |
| **再付着位置** $x_R/H$ | **6.26 ± 0.10** ← 主判定値 |
| 比較プロファイル位置 | $x/H = 1, 4, 6, 10$ |
| 公開比較量 | 壁面 $C_f$、壁面 $C_p$、平均速度 $U/U_\mathrm{ref}$、turbulent shear stress $\overline{u'v'}$ |

**補助参照（合否には混ぜない・解釈用）**:
- Le, Moin & Kim, *JFM* 1997: 低 Re_H の DNS（平均・RMS・Reynolds stress、$x_R/H \approx 6.3$ 系）。
  ただし Re・流入・上壁条件が違うので case/18 の合否基準に混ぜない。
- Jovic & Driver / Driver 系 NASA: 低〜中 Re backstep 実験、流入 BL・再付着詳細向き。
- Bin, Park, Lv & Yang 2023 (arXiv): backstep を含む canonical separated-flow DNS/LES/RANS benchmark。

> **注意（2D 実験 vs 3D 計算）**: 2DBFS は 2D 形状の実験だが実流れは 3D 乱流。3D DDES の
> validation target として使える。ただし公開比較量は**平均場・壁面量・乱流せん断応力が主**で、
> LES 的スペクトルや全 Reynolds stress tensor は無い。→ それらは forge 側で自前取得する（下記）。

**forge 既存状況**: RANS・LES の run あり（`run_0001` 〜）。  
DDES 用に新 run を追加する。3D スパン方向解像が必要（2D は DES でも RANS 的動作）。

**手順**:
1. 既存 RANS 収束場から DDES `DESmode: 1` で非定常（`timeIntegration: 3`）に切替え
2. **3D・十分なスパン幅・周期 BC**で回す（2D は DES でも RANS 的動作）
3. 時間平均を取り、$x_R/H$・$C_f$・$C_p$・$U/U_\mathrm{ref}$・$\overline{u'v'}$ を
   $x/H=1,4,6,10$ で TMR 2DBFS データと比較
4. `l_des`・`f_d` フィールドで「剥離せん断層で LES 化」を可視化
5. **自前で** $u_\mathrm{rms}, v_\mathrm{rms}, w_\mathrm{rms}, \overline{u'v'}$・できればスパン相関 or
   スペクトルを取り、解像乱流が立っているかを確認（TMR にはない量）

**判定基準（強化・レビュー反映: 壁圧・f_d>0.9 だけでは不十分）**:
- [ ] **再付着長 $x_R/H = 6.26 \pm 0.10$** に整合（主判定値。RANS SST は ≈5.5 で過少評価）
- [ ] **壁面 $C_f$・$C_p$** が TMR 2DBFS 実験データと整合（再付着点・回復領域）
- [ ] **平均速度 $U/U_\mathrm{ref}$ と $\overline{u'v'}$** が $x/H=1,4,6,10$ で実験と整合
- [ ] 剥離域の $f_d > 0.9$（LES への切替 — **必要条件にすぎない**）
- [ ] **自前 $u_\mathrm{rms}$ 等が有意**（せん断層で解像乱流が立つ。f_d 立っても風上散逸で死ぬので必須）
- [ ] **エネルギースペクトルに −5/3 慣性域**が見える（§4.8 低散逸 flux の on/off で比較）
- [ ] **スパン相関**がスパン幅内で減衰（スパンが狭すぎないことの確認）
- [ ] NaN/Inf なし、統計サンプルが十分（最低 5 流れ通過時間以上）

> **backstep は SBLI への門番**: ここで「f_d は立つが乱流が立たない」なら、SBLI に行っても
> 必ず失敗する。低散逸 flux（§4.8）の効果をこの安価な 3D ケースで先に確立してから T3 へ。

---

### ティア 2：（削除）

`case/06.mach3_wind_tunnel` は `LESorRANS: 0`（オイラー計算、粘性なし・乱流モデルなし）のため、
DDES 検証ケースとしては対象外。超音速付着 BL の RANS 保護確認は T3 ケースのコーナー上流平板部で兼ねる。

---

### ティア 3：SBLI 精度検証（新規ケース・高コスト・目標適用に直結）

#### T3: 超音速圧縮コーナー / 入射衝撃波-BL 干渉

超音速ノズル・ピントルバルブの核心物理（衝撃波起因の剥離と再付着）を
実験データと定量比較するための標準 SBLI ベンチマーク。
**IDDES-SST が DDES より優位かどうかをここで判断する。**

**推奨ケース比較**:

| 候補 | 条件 | 剥離点 | 参照データ | 難易度 |
|------|------|--------|-----------|--------|
| **Settles 圧縮コーナー** ← **推奨** | $M=2.84$, $\theta=24°$ | 幾何固定（コーナー） | Brown (2014) NASA TM-218353 より Tables 抽出済 | ★★☆ |
| Dupont 入射衝撃波 | $M=2.3$, 8° くさび | 非固定（圧力勾配） | Dupont et al. (2006) J. Fluid Mech. 559 | ★★★ |

**Settles 24° を推奨する理由**:
- 剥離点がコーナー幾何に固定 → DDES シールド関数が「剥離した」と確実に認識できる
- 付着 BL（コーナー上流）と剥離泡（コーナー後）が明確に分離 → T2 で落とした「超音速 BL 保護」も兼ねて確認できる
- シンプルな 2D 楔形状でメッシュ構築が容易

**forge 入力条件（NASA TM-2014-218353, Tables 1–7 より抽出・検証済）**:

| パラメータ | 値 | 用途 |
|-----------|-----|------|
| $M_\infty$ | 2.84 | 流入 BC |
| $p_\infty$ | 23,600 Pa | 流入 BC |
| $T_\infty$ | 100.3 K | 流入 BC |
| $\rho_\infty$ | 0.82014 kg/m³ | 確認用 |
| $U_\infty$ | 570 m/s | 確認用 |
| $\mu_\infty$ | 7.002×10⁻⁶ kg/(m·s) | 確認用 |
| 単位 Re | 67.16×10⁶ /m | メッシュ解像度設計 |
| 壁面条件 | 断熱（$T_\text{wall}=276\,\text{K}$） | BC |
| 楔角 | 24° | ジオメトリ |
| 衝撃波角 | 44.3° | 確認用 |
| $p_2/p_1$（コーナー後） | 4.42 | 出口 BC・確認 |
| $\delta_0$（上流 BL 厚） | 23.72 mm | 入流 BL 設定 |
| $\tau_{w,0}$（上流壁面せん断） | 145.0 Pa | 確認用 |

**壁面圧力の参照値**（自由干渉理論より、グラフ digitize の前に合否判定に使用可能）:

| 位置 | $p_w/p_\infty$ |
|------|----------------|
| 上流（基準） | 1.00 |
| 剥離点 $p_{w,s}$ | **1.68** |
| 剥離バブルプラトー $p_{w,b}$ | **1.97** |
| コーナー後（無粘性理論値） | **≈ 4.42** |

**参照ソース**:

| ソース | 種別 | 入手 |
|--------|------|------|
| Brown (2014) *NASA TM-218353* | CFD+実験比較（Tables 1–7 数値取得済） | NASA NTRS 無償公開 |
| Settles, Bogdonoff & Vas (1975/1976) | 実験原典 | AIAA Paper 1975-7 / AIAA J. 14(1) |
| Settles & Dodson (1991/1994) | 実験データベース | *NASA CR-177577/177638* |
| Wu & Martin (2008) *AIAA J.* 46(8) | DNS $M=3$, 16° | 条件近似・参照用 |
| Troshin & Bakhne (2024) | IDDES-SST 計算 | Springer（サーベイ既取得） |

**実践手順**:
1. 上記の流入条件を `bcondConfig.yaml` に直接入力（数値は本計画に記載済）
2. 平板部で BL を育成（Re 単位 67.16×10⁶/m × 1.98 m ≈ $Re_x=1.3\times10^8$、もしくは BL プロファイルを直接入流に指定）
3. $p_w/p_\infty$ を壁面に沿って抽出し、上記参照値と比較
4. 詳細な点群データが必要な場合：Brown (2014) の Fig. 27/28 を WebPlotDigitizer で digitize（20–30 点）

**新規 case ディレクトリ**: `case/37.sbli_compression_corner`

**フロー概略**:
```
  → M=2.84 一様流
  ┌──────────────────────┐
  │ 平板（付着 BL, RANS） │ ← f_d ≈ 0 を確認
  └──────────────────────┘
                          \← 24° コーナー
                           \
                            \ 剥離泡（LES 活性）← f_d → 1 を確認
                             \___________________________
```

**メッシュ要件**:
- 壁法線: automatic wall treatment で $y^+ \approx 30$–80（コスト削減）
- コーナー上流平板: RANS 相当の接線解像度で可
- コーナー後剥離域: スパン方向に 3D 解像 $\Delta z \approx 0.1 \delta$

**前提（レビュー反映・必須）**:
- T1-B（backstep）で**解像乱流が立つこと**（スペクトル・RMS）を先に確認済みであること。
- **§4.8 の低散逸 flux（Phase 1.5）を有効化**して臨むこと。風上散逸のままだと
  せん断層 KH 成長が遅れ、壁圧プラトーが平滑・非定常圧力が弱く出る（§5.7）。
- コーナー後剥離域は**スパン方向 3D 解像**が要る。スパンが Δmax を律速すると LES 化が遅れる（§3.2）。

**判定基準（強化・壁圧だけで合格にしない）**:
- [ ] 壁面圧力プロファイルの Settles 実験値との乖離 < 10%（$p_w/p_\infty$ 分布。**sanity check 扱い**）
- [ ] **Cf プロファイル**が実験/参照 CFD と整合（壁圧より厳しい・剥離再付着に敏感）
- [ ] **剥離点・再付着点の位置**が実験と一致（壁圧プラトーの平滑さで偽合格しないため）
- [ ] **壁面圧力 RMS**（非定常変動）が有意で、せん断層の非定常性が出ている
- [ ] **スパン相関**がスパン幅内で減衰
- [ ] コーナー上流平板で $f_d < 0.1$（RANS 保護の超音速動作確認、T2 の代替）
- [ ] DDES の結果が信頼できてから、IDDES は §4.6 の前提（完全関数群＋流入乱流）を満たした上で比較

> **偽合格の罠（レビュー）**: 壁圧は衝撃捕獲＋平均剥離泡で「それっぽく」合うことがある。
> Cf・剥離/再付着・圧力 RMS・スパン相関を併せて見ること。

---

### ティア 4：超音速キャビティ（ピントルバルブ類似・将来）

#### T4: `case/01.cavity` または新規超音速キャビティ

**目的**: ピントルバルブの環状ギャップに類似した「キャビティ」形状での
Rossiter 振動モードを DDES/IDDES で再現する。
非定常・音響的フィードバックを含む最も複雑な検証。

**参照**: AIAA/NASA キャビティ DDES ベンチマーク  
$M=1.5$, $L/D = 4$（Colonius, Lele & Moin 1993; Rossiter 1964）

**現状**: `case/01.cavity` の既存ケースが転用できる可能性があるが、
超音速 SST + DDES での計算は新規になる見込み。
**本ティアは ティア 3 完了後に着手する。**

---

### 検証ケース全体まとめ

```
ティア  ケース                           Phase   コスト    何を確かめるか
─────────────────────────────────────────────────────────────────────────────
T1-A   case/26.flat_plate_sst            DDES    低        BL シールド（f_d≈0）
T1-B   case/18.backstep                  DDES    中        剥離域 LES 活性・再付着長 x_R/H=6.26±0.10
                                                           参照: NASA/TMR 2DBFS = Driver & Seegmiller (1985)
                                                           実験データDL可 (Cf,Cp,U,u'v' @ x/H=1,4,6,10)
T2     削除（case/06 はオイラー計算）    —       —         —
T3     case/37.sbli_compression_corner   DDES/   高        SBLI 剥離予測・実験比較
       M=2.84, θ=24°                     IDDES             超音速 BL 保護も兼ねる
                                                           参照: Settles & Dodson (1994) NASA CR-177577
                                                                 Wu & Martin (2008) DNS
T4     case/01.cavity or 新規            IDDES   高        Rossiter 振動・音響
       (ピントルバルブ類似)                                 参照: AIAA キャビティ標準ケース
```

### 共通判定基準（全ティア）

- `check_convergence.py` を使って残差の収束を定量確認（PASS 必須）
- 時間平均統計は最低 5 流れ通過時間以上取ること
- `l_des`・`f_d`・`r_d` フィールドを毎ケース出力・可視化し、RANS/LES 分布が物理的に合理的かを確認
- NaN/Inf: `detectNaN: 1` を使い序盤 100 ステップ以内に NaN がないことを確認

### 時間積分の方針（レビュー反映）

- 物理時間精度の**陽的 RK3（`timeIntegration: 3`）は CFL 制限で高コスト**。600k セル非定常は
  可能だが SBLI 統計まで取ると実用的に厳しい。
- **手順**: まず**小さめ 3D backstep で陽的 RK3 の基準解**を作る（位相・乱流エネルギーの真値）。
  次に **dual-time（陰的非定常, `dualTime: 1`）が同じ統計を再現するか確認**してから SBLI へ。
- dual-time は内反復（`nSubIterDualTime`）が甘いと非定常の位相・乱流エネルギーを壊すので、
  backstep で内反復数を校正してから本番に使う。

---

## 7. 影響範囲

**変更ファイル（Phase 1: DDES 機能）**:
- `solver_density_cuda/input/solverConfig.hpp` / `.cpp`（Config 追加）
- `solver_density_cuda/variables.hpp` / `.cpp`（フィールド追加: delta_les, l_des, fd_shield, rd_des）
- `solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu` / `.cuh`（Δ計算 + DES カーネル + 診断出力）
- `solver_density_cuda/cuda_forge/ransSource_d.cu` / `.cuh`（k 消滅項切替）
- `solver_density_cuda/variables.cpp`（Δ事前計算: 幾何セットアップ直後）

**変更ファイル（Phase 1.5: DES 用低散逸 flux・§4.8）**:
- `solver_density_cuda/cuda_forge/convectiveFlux_d.cu`（散逸係数 σ の f_d/Ducros 連動スケーリング）
- `solver_density_cuda/cuda_forge/ducrosSensor_d.cu`（センサを「散逸を下げる」用途に流用）
- config（`desLowDissipation` 等）

**影響なし（Phase 1 時点）**:
- `ransTransport_d.cu`（k/ω の対流・拡散は不変）
- `ransBoundary_d.cu`（BC は不変）
- `ransWallFunction_d.cu`（wall function は不変、DES と直交）
- ※ `convectiveFlux_d.cu` は **Phase 1 では不変だが Phase 1.5 で変更**（SBLI 物理検証の前提）

**docs 更新**:
- `methods/turbulence/theory.md`：§8 に DES/DDES/IDDES 理論節を追加
- `methods/turbulence/implementation.md`：§3.4 に `DESmode` 実装方針を追記
- `methods/index.md`：変更があれば同期

---

## 8. 完了条件

- [x] `methods/turbulence/theory.md` に DES/DDES/IDDES 理論節追加 (§8、2026-06-21)
- [x] `methods/turbulence/implementation.md` に DESmode 実装節追加 (§3.8、2026-06-21)
- [x] Phase 1 (DDES) 実装・ビルド成功（診断出力 f_d/r_d/delta_les/l_des 込み）
- [x] `DESmode: 0` で `case/26.flat_plate_sst` がビット不変（atomicAdd 非決定性の床内で
      新旧バイナリが一致＝1 step で base–base と base–des の差が同値・§9 検証ログ）
- [x] `case/18.backstep`（3D・857k cells・周期スパン）で DDES 機能動作確認
      （剥離せん断層 f_d≈0.67, 付着 BL f_d≈0.02, NaN なし・§9 検証ログ）
- [ ] **（SBLI 物理検証を目標にする場合）Phase 1.5 低散逸 flux 実装＋ backstep で解像乱流
      （RMS・スペクトル −5/3）確認** — これ未達のまま SBLI 定量一致を主張しない（次セッション）
- [x] `case/18.backstep`・`case/26.flat_plate_sst` の README.md に run 追記
- [x] `.github/plans/README.md` の status 更新（Phase 1 完了・status は in_progress 継続）
- [ ] 本 plan の `status` を `done` に更新（Phase 1.5/2 が残るため `in_progress` 継続）

---

## 9. 変更ログ

- `2026-06-21` — 初稿。DES/WMLES サーベイ（turbulence-des-wmles-survey.md）を受けて作成。
  Phase 1 (DDES) / Phase 2 (IDDES) の 2 段構成に設定。
- `2026-06-21` — 外部レビュー（GPT-5.5）を反映。主な追加:
  ① **Phase 1.5 = DES 用低散逸 flux**（§4.8）を新設＝SBLI 物理検証の前提（最大の失敗要因は
  風上散逸が解像乱流を殺すこと。既存 Ducros は limiter を落とす＝逆効果）。
  ② **診断出力 f_d/r_d/delta_les/l_des**を必須化（§4.2、float32 で f_d が 0/1 張り付くのは
  NaN より危険）。③ **r_d の float32 ガード**（clamp[0,10]＋物理スケール floor、§4.4）。
  ④ **Δmax の逆向き落とし穴**（SBLI せん断層でスパン粗さ律速→LES 化遅れ、§3.2）。
  ⑤ **backstep/SBLI の合否基準強化**（壁圧・f_d だけで合格にしない: RMS・スペクトル・Cf・
  剥離再付着・スパン相関、§6 T1-B/T3）。⑥ **IDDES 簡略版は定量評価不可**＝完全関数群＋
  流入乱流が前提（§4.6/§2）。⑦ **時間積分**は陽的 RK3 基準解→dual-time 再現確認の順（§6）。
  ⑧ **失敗予測**節を新設（§5.7）。
- `2026-06-21` — **§4.8 (Phase 1.5 設計) 改訂**: ① Ducros センサは SLAU で実質排除 (`ducrosLimiter`
  既定 0・場差<0.1%) のため、低散逸 flux を Ducros に依存させない方針に変更。② forge は既に
  `solver: KEEP_SLAU` で「KE 保存中心 (KEEP) を基本にセンサ `duc` で SLAU 風上へ切替」を実装済みと
  判明 ([`KEEP_SLAU_d`] L2452 の `(1-duc)·KEEP + duc·SLAU`)。本命設計を **σ-scaling SLAU (旧案1) から
  「KEEP_SLAU の `duc` に f_d を注入」(`duc=max(1-f_d, 1-ψ_limiter)`) に変更** — KE 保存性で優れ
  新規コードもほぼ不要。衝撃波保護は既存 MUSCL リミタ ψ が肩代わりし排除した Ducros は不要。
  `desLowDissipation` は `solver: KEEP_SLAU` 前提。σ-scaling は案 B (代替) へ降格。
- `2026-06-21` — **Phase 1 (DDES) 実装・検証完了** (status → `in_progress`、Phase 1.5/2 継続)。
  - **実装**: config (`DESmode`/`C_DES_kw`/`C_DES_ke`)、診断 4 変数 (`delta_les`/`l_des`/`fd_shield`/
    `rd_des`, HDF5 出力)、Δmax 事前計算 (`variables.cpp`, 実行時重心+面接続から host 1 回)、
    DDES カーネル (`turbulent_viscosity_d.cu`, `compute_des_length_d`+`compute_rd_ddes`/`compute_fd_ddes`,
    float32 ガード)、k 消滅項切替 (`ransSource_d.cu`, `DESmode>0` で $\rho k^{3/2}/l_\mathrm{des}$)。
    依存順は vis_turb 計算直後に l_des を計算 (§4.7)。clean rebuild は `FORGE_CUDA_ARCHITECTURES=86`
    必須 (デフォルト configure が sm_52 を選び `atomicAdd(double*)` で失敗する環境メモ)。
  - **バグ修正**: §3.1 の $l_\mathrm{RANS}$ は $\sqrt{k}/(\beta^{*}\omega)$ が正 ($\beta^{*1/4}$ は誤、
    上記 §3.1 訂正注)。初版 $\beta^{*1/4}$ では mode1 の付着 BL で MSD により `roK` が ~100% 変化
    したが、$\beta^{*}$ 修正で `roK` 差 $1.8\times10^{-5}$ (= mode0 と実質一致) に解消。
  - **T1-A 回帰 (`case/26.flat_plate_sst`)**: forge は atomicAdd flux scatter の順序により run-to-run
    非決定 (同一 baseline binary でも `roe` ~0.14, `omega` ~624 ばらつく) のため厳密ビット一致は
    検証不能。1 step で base–base と base–des の max|Δ| が **全保存量で同値** (ro 1.19e-7, roe 1.56e-2,
    roK 3.73e-9, roUz 1e-15) → DESmode:0 は baseline と数値的に同一経路 (atomicAdd 床内一致)。
    根拠: `run_0019_des_regr_mode0`。
  - **T1-A シールド (`case/26.flat_plate_sst`, 発達 BL から mode0/mode1 を 1000 step)**: 発達した
    乱流場 (nut/nu med 28–65) で mode1 vs mode0 の `roUx` relL2 = **5.7e-7** (wall-resolved) / 4.7e-5
    (y+30 wall-modeled) ≪ Cf 0.1% 基準、`roK` relL2 = 1.8e-5。DES リミッタは BL 外縁の 1.6–3% の
    セルのみで発火 (内側 BL はシールド $f_d{\approx}0$ + グリッド規準 $l_\mathrm{RANS}<C_\mathrm{DES}\Delta$
    の二重保護)。NaN なし。根拠: `run_0019`(mode0)/`run_0020`(mode1, wall-resolved),
    `run_0021`(mode0)/`run_0022`(mode1, y+30)。**注意**: `run_regr_cf` 等の nut/nu≤1 の未発達場では
    シールドが効かない (νt≈0 で rd 小) ため、検証は必ず発達乱流場 (`run_0007`/ewt) から restart する。
  - **T1-B 機能 (`case/18.backstep`, 3D 857k cells, 周期スパン 4H, DESmode:1, 500 step implicit)**:
    3D メッシュ変換 (mesh quality VERDICT **PASS**: AR max 2.3, skew 0.333) → 2D RANS 場から
    `interp_field` cross-mesh restart → 実行。**NaN なし**。f_d 分布: 付着近壁 BL `f_d≈0.02`
    (frac<0.1=0.96)、ステップ後方の剥離せん断層 (x∈[3,9)) `f_d≈0.67`=LES 活性、再循環泡内部
    (Ux<0, 高 νt) `f_d≈0`=RANS modeled、再付着 ~6H。Δmax≈0.13, rd は [5e-4,10] に分布 (飽和張り付き
    なし)。根拠: `run_0035_des3d_ddes`。**収束 VERDICT**: `check_convergence.py` は NOT CONVERGED
    (DDES は本質的に非定常で残差は下げ止まらない・rms は falling)。これは想定どおりで、Phase 1 の
    機能確認 (f_d 分布 + NaN なし) が合格基準。**解像乱流の定量化 (RMS・−5/3・x_R/H) は Phase 1.5**
    (低散逸 flux) 後の仕事で、本セッションでは未実施。
- `2026-06-21` — T1-B backstep の定量基準を確定。**Primary = NASA/TMR 2DBFS
  (Driver & Seegmiller 1985, 実験データDL可)**: $Re_H$≈36,000, $x_R/H$=6.26±0.10,
  $C_f$/$C_p$/$U$/$\overline{u'v'}$ @ $x/H$=1,4,6,10, $M$≈0.128(ほぼ非圧縮)。
  補助参照に Le-Moin-Kim 1997 DNS / Jovic-Driver / Bin+2023 benchmark（合否には混ぜない）。
  LES 的量（rms・スペクトル・スパン相関）は TMR にないため forge 側で自前取得する旨を明記。
- `2026-07-19` — §3.4 を Gritskevich 2012 正準形に訂正 (初稿の $\max$ 形 $l_\mathrm{IDDES}$ /
  $f_e$ 式は外部要約由来の誤記)、§4.6 を確定設計に更新。関数群の事前検証スクリプト
  `tools/verify_iddes_functions.py` を追加 (形状・RANS 縮退・WMLES 分岐・float32、VERDICT: PASS)。
  §4.8 は単一スキーム KEEP+ES を本命に更新し前提検証①② 合格 (同日、既存記載)。
- `2026-07-20` — **Phase 2 (IDDES) 実装・機能検証完了**。
  - **実装**: `compute_des_length_d` に `DESmode==2` 分岐 (`compute_l_iddes`: $f_B/f_{e}/f_{dt}/
    \tilde f_d$, $\Delta_\mathrm{IDDES}$, 定数 $C_t{=}1.87/C_l{=}5.0/C_{dt1}{=}20/C_w{=}0.15$
    compile-time)。$r_{dt}/r_{dl}$ は `compute_rd_ddes` 流用 (clamp[0,10] は tanh 飽和域で無影響)。
    診断: `fd_shield`←$\tilde f_d$ (**mode2 は 1=RANS で DDES と向きが逆**)、`rd_des`←$r_{dt}$、
    新フィールド `fe_iddes` (HDF5 出力)。新規幾何量・事前計算なし。
  - **DDES 回帰 (mode1 ビット不変)**: `case/18.backstep` `run_0128` (新バイナリ) vs `run_0125`
    (実装前) vs `run_0130` (HEAD 再実行=atomicAdd ノイズフロア実測)。同一 restart (md5 一致)
    500 step implicit で、全 12 フィールド (保存量+P/T+fd/l_des/rd) の max|Δ|/relL2 が
    **signal ≒ ノイズフロア** (例 roUx max|Δ| 4.5e-4 vs 4.5e-4) → **DESmode:1 経路は数値的に不変**。
  - **機能確認① (RANS 的場, `run_0129`)**: run_0125 と同一入力で DESmode:2 のみ差分・500 step。
    NaN 無し。付着 BL $\tilde f_d$=1.000・$l/l_\mathrm{RANS}$=1.000 厳密 (シールド)、せん断層
    $\tilde f_d$ med 0.036 (LES)、泡内 1.000 (RANS)。$f_e$ は付着 BL 近壁のみ活性 (24%, max 0.98)。
  - **機能確認② (発達乱流場, `run_0131`)**: run_0126 res_30000 (解像 3D 乱流飽和) から
    DESmode:2 で 3000 step unsteady RK3。NaN 無し。**WMLES 分岐が設計どおり動作**:
    遠壁 ($f_B{<}0.01$) で $f_{dt}$ med=0.992・frac>0.9=67% (解像乱流検知)、$f_B{>}0.999$ 帯で
    min $\tilde f_d$=1.0000 (壁直近 RANS 強制厳密)、せん断層 $l/l_\mathrm{RANS}$ 0.537≒DDES 0.524、
    解像乱流維持 (Uz rms 4.97→5.13 m/s)。$f_e$ は近壁帯のみ活性 (20-29%)。
  - **Phase 2 の完了範囲 (機能確認まで) を達成**。WMLES モードの定量評価 (チャネル log 層 /
    mismatch) は従来どおり streamwise 体積力 config +
    [turbulence-wmles-wall-stress plan](turbulence-wmles-wall-stress.md) のチャネル検証系、
    非周期ケースは流入乱流生成 plan にゲート。
  - 環境メモ: clean configure は `-DCMAKE_CUDA_ARCHITECTURES=86` (cache 変数名はこれが正。
    2026-06-21 記載の `FORGE_CUDA_ARCHITECTURES` は誤記) 必須。
