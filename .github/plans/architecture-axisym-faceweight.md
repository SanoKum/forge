# 軸対称 近軸 k/ω の数値特異点 (調査・修正)

## メタ

- **area**: `architecture`
- **status**: `draft`
- **related_docs**:
  - `docs/turbulence/theory.md` (§7.5)
- **related_plans**:
  - `architecture-axisymmetric.md` (B 流儀 r 重み実装)
  - `architecture-axisym-sst.md` (軸対称 SST 幾何項)
  - `turbulence-kato-launder.md` (症状緩和)
- **created**: `2026-06-12`
- **owner**: `CFD Dev`

## 1. 目的

軸対称 RANS でスロート**軸中心線の $k$/$\omega$ が非物理的に過大**になる (case 29 conical)。
**軸対称固有 (2D で出て 3D で消える) は確実**だが、機構は未確定。当初の「グリッド発散の数値特異点」は
uniform 格子の k 非収束により**未確証**。正しい検証で機構を特定し根治する。

## 2. スコープ

- **やる**: 軸第1セル $k$ の数値特異点の機構特定と修正。候補 (§4) を切り分け。
- **やらない**: Kato–Launder (症状緩和は別 plan で実装済)。

## 3. 確定している事実 (根拠)

[`docs/turbulence/theory.md`](../../docs/turbulence/theory.md) §7.5 /
[`case/29.bell_vs_conical/README.md`](../../case/29.bell_vs_conical/README.md):

- **グリッド発散 (確定的)**: 軸第1セル $k$ = 17 → 836 → 1956 (軸セル 0.50→0.10→0.05 mm)。
  中間半径背景 $k$=4.7・壁 BL $k\sim3.8\times10^4$ は**グリッド収束**。軸だけ発散 = 数値特異点。
- **3D で消える**: 90° セクタ butterfly (壁第1層 8µm 全周一様), 低背景で軸平坦 (第1セル/核 ~0.85)。
- **平均流は滑らか**: 軸近傍で $U_x,T,\rho$ 平坦、$k$ のみ尖る。

## 4. 機構候補 (当初仮説の訂正を含む)

> **訂正**: 当初「面積が planar で体積と不整合」と推定したが誤り。flux 用の面積・体積は
> 両方 $r$ 重み済 ([`variables.cpp`](../../solver_density_cuda/variables.cpp) `ss *= r_face`,
> `volume *= r_cell`, 軸面 `r_floor=1e-20`。B 流儀)。h5 `surfArea` は重み付け前の幾何値だった。

候補:
1. **勾配の面積非対称**: 軸対称で勾配 (Green-Gauss, 拡散が使用) は **planar 面積**
   ([`calcGradient_d.cu`](../../solver_density_cuda/cuda_forge/calcGradient_d.cu): `ss_planar`/`A_planar`)
   一方 flux は $r$ 重み。近軸でこの非対称が $k$/$\omega$ 拡散勾配を乱す可能性。
2. **近軸 stiffness + 生産正帰還**: $r$ 重み FV は近軸で $V,S_f\to0$ (+`r_floor`) と硬い。
   生産リミタ $P_k=\min(\mu_t S^2, 10\beta^*\rho k\omega)$ は $P\propto k$ で、近軸の微小不均衡を
   指数増幅し得る。
3. **拡散の $\theta\theta$ 寄与 / 軸 BC**: $k$/$\omega$ 拡散の近軸幾何項・軸面 (`r_floor`) の扱い。

## 5. 調査ステップ

0. **正しいグリッド試験 (最優先)**: 壁を保ち近軸だけ細分した 2 面クラスタ格子を **$k$ 収束まで**回し、
   軸 $k$ がグリッド収束か発散かを判定 (uniform 格子は $k$ 非収束で無効だった)。
1. `scalarDiffusion=0` で軸スパイクが消えるか (※uniform 格子では 836→114 と拡散主因を示唆したが
   未収束場のため要再確認)。
2. 候補1: 拡散勾配を $r$ 重み面積で計算したら改善するか。
3. 候補2: 生産リミタ無し/緩めで軸スパイクのグリッド依存が変わるか。
4. 機構特定後、最小修正を実装。

## 6. 検証

- **本症状**: `case/29` conical RANS で軸第1セル $k$ が**グリッド収束** (uniform ny=100/200 で一致)。
- **回帰 (最重要)**: 全軸対称ケース (`case/27` CEA, `case/23`, `case/29` 推力 mdot/λ) で悪化なし。
- **3D 整合**: 2D 軸と 3D 軸の $k(r)$ が一致方向。

## 7. 影響範囲

- 機構次第 (勾配/拡散/source/軸 BC)。非軸対称・3D は不変が条件。

## 8. 完了条件

- [ ] 機構特定 (§5)
- [ ] 最小修正で軸 $k$ のグリッド収束
- [ ] 全軸対称回帰で悪化なし
- [ ] docs §7.5 を根治済みに更新

## 9. 変更ログ

- `2026-06-12` — 初稿 (誤診の訂正含む)。case 29 軸 $k$ 過大は**軸対称固有 (2D 出/3D 消) は確実**だが
  機構は未確定。誤りだった推定: ①planar 面積 vs r 重み体積の不整合 (flux は両方 r 重み済)、
  ②拡散勾配が planar 面積を使うのが bug (勾配は planar 作用素で正しい)、③グリッド発散の数値特異点
  (根拠の uniform 格子は壁粗く k 非収束=無効)。正しい検証は壁保持の 2 面クラスタ格子を k 収束まで。
