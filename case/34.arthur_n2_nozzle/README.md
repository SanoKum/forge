# case/34 Arthur (1952) 2D source-flow 極超音速ノズル — dry 超音速膨張

P.D. Arthur の博士論文 (Caltech/GALCIT, 1952) で用いられた **2 次元 source-flow
ノズル**を再現し、**乾いた (dry) 非粘性の超音速膨張**を forge で計算するケース。
窒素凝縮を扱った論文
[`papers/on nitrogen condensation in hypersonic nozzle flows_summary.md`](../../papers/on%20nitrogen%20condensation%20in%20hypersonic%20nozzle%20flows_summary.md)
の検証ノズル (Arthur 1952) に対応する。

> **注意**: forge には核生成・液滴成長など**相変化 (凝縮) モデルは無い**。本ケースが
> 再現するのは「凝縮なし」の等エントロピー的膨張 (Arthur Fig.4/11 の dry 1D 理論線) まで。
> 凝縮による静圧上昇そのものは対象外。

## ノズル形状 (出典: Arthur 1952, II.A Apparatus)

- 2 次元 source-flow ノズル: スロート全高 **0.010 in = 0.254 mm**、出口全高
  **~1 in = 25.4 mm**、幅 1 in 一定、**全開き角 11° (half-angle 5.5°)**、面積比 **A/A*≈100**。
- 壁は仮想ソース点から出る直線 (source flow)。**スロート曲率は原典に記載が無い**
  (製作詳細は Ref.5 = Nagamatsu & Willmarth, GALCIT Memo No.6 / JAP 23(10) 1089, 1952 に委ねられ、
  入手できず)。source-flow は本来直線壁なので曲率は未定義。
- 本ケースは **発散部のみ・上半分**をモデル化 (中心線 y=0 を対称 slip)。スロート直下から
  始め、**入口に超音速インレット (M=1.05) を置いて choking を回避**。
- 壁形状はメッシュ上 **双曲線 `y(x)=√(yt²+(x·tanα)²)`** で表現:スロートで水平接線
  (dy/dx=0) となり鋭い凸コーナーの膨張特異点を回避、下流は 5.5° 直線 (Arthur の source-flow)
  に漸近。実効スロート曲率 R_t = yt/tan²α ≈ 13.7 mm (滑らかな代表値)。

メッシュ生成: [`mesh/arthur_nozzle.geo`](mesh/arthur_nozzle.geo) (gmsh, 400×60 の 1 層押し出し,
24,000 hex)。境界 physID: 1=inlet, 2=outlet, 3=wall(上壁+中心線+z両面, 全て slip)。

## 計算条件

- 貯気: P0 = 844,037 Pa (8.33 atm), T0 = 290 K (Arthur Run 34-1 / Fig.4,11 相当)。窒素 γ=1.4。
- 非粘性 Euler (viscMethod 0)、SLAU、1 次精度、陰解法 block-DPLUR。
- 初期場 `arthur_n2` ([`solver_density_cuda/input/setInitial.hpp`](../../solver_density_cuda/input/setInitial.hpp)):
  スロート超音速状態 (M=1.05) の一様場。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_slau_dry` | 直線壁(スロート鋭頂点) + cfl_pseudo=5 | step 467 で発散 (スロート上端の凸コーナー膨張特異点; Ux~1.5e8)。`res_nan_467.h5` | 破棄 (コーナー特異点の記録) |
| `run_0002_slau_dry_cfl1` | 双曲線スロート(滑らか) + cfl_pseudo=1 + **1次精度** (convMethod 0) | **収束** (rms 6.7–6.9 桁低下)。出口 **M=6.75** (等エントロピー 6.93)、A/A*=99.6、P/P0=2.79e-4 (理論 2.58e-4)。`postproc_arthur.png`, `residual_history.png` | active (1次基準) |
| `run_0003_slau_2nd` | run_0002 の収束場を引き継ぎ **2次精度** (convMethod 1, limiter 2 Venkatakrishnan) | 出口 **M=6.87** (等エントロピー 6.93、誤差 2.5%→**0.8%**)、P/P0=2.71e-4。Mach が A/A*=100 まで等エントロピー線にほぼ完全一致。残差は **~3e-4 でプラトー** (リミタ起因のリミットサイクル; 場は準定常)。`postproc_arthur.png`, `residual_history.png` | active (高精度) |
| `run_0004_cond_skeleton_dry` | 非平衡凝縮 Phase 1 (`condensation:1, nCondSpecies:1`)。run_0003 と同条件で **4 モーメント (ρg,ρQ2,ρQ1,ρQ0) を受動スカラー (ソース=0) として輸送**する骨格の回帰確認 ([plan](../../.github/plans/condensation-nonequilibrium.md)) | **dry 回帰一致**: 出口 **M=6.874**, P/P0=2.71e-4 で run_0003 と同一。液相モーメントは全セル厳密 0、NaN/Inf なし。cond ON vs OFF の場差 (ro maxrel 9.5e-4) は cond OFF を 2 回回した run-to-run ノイズ (1.0e-3, リミットサイクル+GPU atomic 非決定) **以下**=凝縮スカラーは気相に不干渉。`postproc_arthur.png`, `residual_history.png` | active (Phase 1 骨格) |
| `run_0005_cond_off_ref` | 上の回帰用 cond OFF 参照 (現行バイナリ, `condensation:0`)。run_0004 とのビット差が plateau ノイズ由来であることの対照 | run_0004 との場差は run-to-run ノイズ内 (上記) | ref (対照, 成果物は破棄可) |

> `run_0002` の VERDICT は形式上 NOT CONVERGED と出るが、これは**面外 z 方向の `rms_roUz`
> (≈1e-14 のノイズ; 単層押し出しで物理が無い) のみ**による偽陽性。物理 4 残差
> (`rms_ro`,`rms_roUx`,`rms_roUy`,`rms_roe`) は 6.7–6.9 桁低下し平坦=実質完全収束。

## 結果まとめ

- 出力先: `case/34.arthur_n2_nozzle/run_0002_slau_dry_cfl1/`
- 主要図:
  - [`run_0002_slau_dry_cfl1/postproc_arthur.png`](run_0002_slau_dry_cfl1/) — 左: Mach vs A/A* (Arthur Fig.4 相当),
    右: P/P0 vs x (Arthur Fig.11 相当)。forge の中心線・壁が等エントロピー線に一致し、
    かつ中心線と壁が互いに一致 (Arthur Fig.11 の知見を再現)。
  - [`run_0002_slau_dry_cfl1/residual_history.png`](run_0002_slau_dry_cfl1/)
- 後処理: [`postproc.py`](postproc.py) (`python3 postproc.py run_0002_slau_dry_cfl1`)。
- forge 出口 Mach の等エントロピー (6.93) との差は **1次 6.75 (2.5%) → 2次 6.87 (0.8%)** と縮小。
  1次の不足は数値拡散によるもので、2次精度化 (`run_0003`) でほぼ理論線に一致する。
  → A/A*≈100 の Arthur 形状の dry 出口 Mach ~6.9 は妥当。
- **マッハ数の解釈** (重要): A/A*≈100 → 等エントロピー dry で M≈6.93 は正しい。これは「凝縮しなければ」
  の値で、実機・論文では**凝縮の潜熱で M はこれより低下**する (論文 parametric study で凝縮あり
  Ma 5.40–6.67)。論文の「設計 Ma=6.0」ノズルは A/A*≈53 の**別形状**で、本ケース (Arthur 検証ノズル,
  A/A*=100) とは異なる。
- 2次精度の残差プラトー (~3e-4) は Venkatakrishnan リミタのリミットサイクルで、forge で散見される
  挙動。出口諸量は準定常で安定しており Mach 比較には使える。

## 実行手順

```bash
cd case/34.arthur_n2_nozzle/mesh
gmsh -3 arthur_nozzle.geo -o arthur_nozzle.msh -format msh4
cd ../run_0002_slau_dry_cfl1
cp ../mesh/arthur_nozzle.msh .
/home/sano/work/forge/solver_density_cuda/build-native/convertGmshToForge arthur_nozzle.msh arthur_nozzle.h5
/home/sano/work/forge/solver_density_cuda/build-native/forge
python3 ../postproc.py .
```
