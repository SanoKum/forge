# 09. Taylor-Green Vortex (3D)

3 次元 Taylor-Green 渦 (TGV)。全方向周期境界の立方体 $[0,2\pi]^3$、構造化 hex $32^3$ (節点 $33^3$)。
低散逸スキームの **運動エネルギー (KE)・エントロピー総和の保存性** を見るための検証ケース。

- メッシュ: [`mesh/Taylor-Green.geo`](mesh/Taylor-Green.geo) (`nx=ny=33`, z 方向 32 層 extrude)。
- 初期条件: `initial: "Taylor-Green"` (無次元, $M_0=0.4$, $\rho_0=1$, $P_0=1/\gamma\approx0.714$)。
  $u=M_0\sin x\cos y\cos z$, $v=-M_0\cos x\sin y\cos z$, $w=0$。
- 数値: pure KEEP (`solver: KEEP`, `keepDissipation: 0` = Roe 散逸無し・非散逸中心流束のみ), 陽解法 RK4 (`timeIntegration: 4`), 乱流モデル無し。

## 重要な設定上の注意

- **無次元・低圧ケースなので EOS フロアを下げること**: 既定の圧力フロア `pMin=1.0` Pa は本ケースの圧力場
  ($P\in[0.65,0.78]$) を全域クランプし、初期場を破壊する (step0 で $P\equiv1.0$)。`physProp.pMin/roMin/tMin`
  を小さく (例 `1e-6`) 指定する。floor の config 化は [`plans/accepted/thermophysics-eos-positivity-floor-config.md`](../../plans/accepted/thermophysics-eos-positivity-floor-config.md) 参照。
- **node モードのメッシュは node 設定で変換すること**: `discretization: "node"` のとき `convertGmshToForge` は
  median-dual を構築して primal を置換する。cell 設定で変換した h5 を node run に流用すると primal hex 上で
  解いてしまい不正 (双対 nCV=35937 ≠ primal hex 32768)。各 node run dir で自身の node config を使って変換する。

## 計算 run 一覧

| run_* | 離散化 | 設定差分 | 主要結果 | 状態 |
| --- | --- | --- | --- | --- |
| `run_0004_node_keep_pmin` | node (median-dual, nCV=35937) | pure KEEP, 粘性 visc=0.0025 (Re≈160), pMin=1e-6 | 2000 step 完走・NaN なし。KE 物理減衰 (K/K0→0.64)、E 保存 ~1e-4 | active |
| `run_0005_cell_keep_pmin` | cell (primal hex 32768) | pure KEEP, 粘性 visc=0.0025, pMin=1e-6 | **step382 で発散** (K/K0=1.15 とスプリアス KE 生成)。非散逸中心流束は collocated で不安定 | active (発散例) |
| `run_0006_node_keep_euler` | node (median-dual) | pure KEEP, **非粘性** visc=0/thermCond=0, pMin=1e-6 | 2000 step 完走。**KE ±1.6% / エントロピー ±0.6% 保存**、E ~1e-4。pure KEEP の保存性を実証 | active |
| `run_keep`, `run_keep_M0.1` | cell | 旧 nagare 時代の参照入力 (旧スキーマ config, 現バイナリでは要修正) | — | ref |

成果物: 保存量比較図 [`conservation_compare.png`](conservation_compare.png) (KE と総エントロピーの時間変化、node vs cell)、
解析スクリプト [`analyze_conservation.py`](analyze_conservation.py)。各 run の `residual_history.png`。

### 結論

- **pure KEEP (散逸ゼロ) は node (median-dual) で安定、cell (collocated) で不安定**。cell は非散逸中心流束の
  odd-even / aliasing でエネルギーがスプリアスに増殖し発散、node は KE 保存性により有界。
- node 非粘性で KE・エントロピー総和が数 % 以内で保存され、pure KEEP の設計どおりの保存性を確認。
- 旧 `run_keep` config (`solver: smac`, 旧 bcond, floor 既定) は現バイナリでは①スキーマ不整合②圧力フロア
  クランプで「昔は計算できたのに回らない」状態になる。原因は solver の退行ではなく config スキーマ進化 +
  フロア設定であり、modern config + `pMin` 引き下げで解消する。
