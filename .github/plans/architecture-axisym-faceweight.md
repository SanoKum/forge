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
**グリッド発散により数値特異点と確定**したが、機構は未確定。原因を特定し根治する。

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

1. ✅ **済**: `scalarDiffusion=0` で軸第1セル $k$ が **836→114 (7×減)**、ピークも軸外 (r=0.15) へ。
   → **拡散が軸スパイクの主因** (候補1 を支持)。残り (114) は生産/移流成分。
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

- `2026-06-12` — 初稿。case 29 で軸 $k$ 過大を**グリッド発散により数値特異点と確定**
  (背景/BL は収束、3D 非再現、平均流平滑)。当初の planar 面積仮説は誤り (flux は $r$ 重み済) と訂正。
  機構候補を整理。**`scalarDiffusion=0` 試験で拡散が主因と判明** (軸 $k$ 836→114)、候補1 (近軸拡散勾配が planar 面積) が有力。
