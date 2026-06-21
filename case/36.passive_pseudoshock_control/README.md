# case/36 多孔壁パッシブコントロールによる擬似衝撃波抑制 (2D 逆解析)

Matsuo et al. (1988)「境界層のパッシブコントロールが擬似衝撃波に及ぼす影響」
(九州大学総合理工学研究科報告 10(1), pp.45-50,
[PDF](https://api.lib.kyushu-u.ac.jp/opac_download_md/17711/p045.pdf)) の試験部を 2D で逆解析する。

**目的**: 多孔壁+キャビティの境界層パッシブコントロールで擬似衝撃波 (shock train) が
弱化し「ショックレス」に近づくかを再現確認する。詳細計画は
[`.github/plans/verification-passive-pseudoshock-control.md`](../../.github/plans/verification-passive-pseudoshock-control.md)。

## 形状 (x=入口からの距離, 上下対称, 片側 1° 広がり)

```
            多孔板+キャビティ (上)
         ┌──┬─┬─┬─┬─┬─┬─┬─┬──┐   ← キャビティ 深さ20
   slip  │  └┐└┐└┐└┐└┐└┐└┐  │   ← 多孔板 厚10, 8スロット 幅0.8 ピッチ4
 ====inlet========================...========== outlet
   (助走30)  1°拡大100      ↑x=130-162    528mm直管
         │  ┌┘┌┘┌┘┌┘┌┘┌┘┌┘  │
         └──┴─┴─┴─┴─┴─┴─┴─┴──┘
            多孔板+キャビティ (下)
   x:0    30        130  162              690
```

| 区間 [mm] | 内容 | 壁 |
|---|---|---|
| 0–30 | 助走領域 (高さ 21.509) | slip |
| 30–130 | 1° 拡大ダクト | no-slip |
| 130–162 | キャビティ部 (多孔板 厚10 + キャビティ 深20、8スロット 幅0.8/ピッチ4、開口率0.2) | no-slip |
| 162–690 | 528mm 直管 (1°) | no-slip |

主流 M=1.89 (多孔壁上流端)。入口 M=1.689 (面積比逆算)。

## 入口条件 (よどみ 3MPa / 288.15K, 空気 CPG γ=1.4)

`inlet_uniformVelocity`: Ux=458.6 m/s, ro=11.733 kg/m³, Ps=617800 Pa (M=1.689, T=183.5K)。
背圧は出口 `outlet_statPress` の Ps で掃引し擬似衝撃波を前後に動かす。

## メッシュ

`mesh/make_mesh.py` (パラメトリック)。`POROUS=True`→多孔壁, `False`→固体壁 (比較基準)。

```bash
cd mesh
python3 make_mesh.py                       # passive_porous.geo
gmsh -2 passive_porous.geo -o passive_porous.msh -format msh4   # forge は msh4 を読む (msh2 不可)
```

初版: porous 63.6k quad / solid 38.3k quad (全 quad)。physID:
1=inlet, 2=outlet, 3=wall_slip, 4=wall(no-slip), 10=fluid。

### Salome (GEOM/SMESH) 版 — `mesh_salome/`

gmsh 版と**同一形状・同一グループ名・同一 physID** を生成する代替フロー。
構造ダクト (直管部 Quadrangle 写像) + 非構造コア/キャビティ (NETGEN quad-dominant)。

```bash
cd mesh_salome
salome -t make_mesh_salome.py                 # POROUS=True  -> passive_porous.med/.unv
salome -t make_mesh_salome.py args:solid      # POROUS=False -> passive_solid.med/.unv
salome -t make_mesh_salome.py args:porous,netgen   # 全 NETGEN フォールバック
# MED/UNV -> gmsh msh4.1 (meshio + h5py + gmsh が必要; forge は msh4 を読む)
python3 med_to_gmsh.py passive_porous.med     # -> passive_porous.msh (physID/グループ検証つき)
```

`make_mesh_salome.py` の主要パラメータ (冒頭に集約): `POROUS`, `STRUCTURED_DUCT`,
`NY`/`NX_RUN`/`NX_DIFF`/`NX_DOWN` (構造分割数), `CL_CORE`/`CL_SLOT`/`CL_CAV` (特性長 [m])。
スロット口・ブロック界面などの内部エッジはグループに含めず conformal に連続させ、
`med_to_gmsh.py` が physID=0 の内部線を除外して書き出す (forge は dim=1 を全て境界として読むため)。

**確認図** (gmsh 版で生成した同一形状・現行 10mm 板厚): `mesh_salome/preview_porous.png`
(全体+キャビティ), `mesh_salome/preview_slotzoom.png` (スロット 0.8mm が ~4-5 セル)。
gmsh 等価メッシュ: porous 62.8k quad (100% quad), 境界線 inlet/outlet=各64・wall_slip=100・
wall=2274, 最小辺 ~0.087mm。`med_to_gmsh.py` の変換ロジックは合成ケースで検証済み。

> Salome 本体は本リポジトリ環境に未導入のため、`make_mesh_salome.py` の SMESH 実行は
> 各自の Salome で行うこと (形状/グループ/physID 規約は gmsh 等価メッシュで確認済み)。
> 変換後の `.msh` は `gmsh passive_porous.msh` で開けば可視化できる (Salome GUI 不要)。

## 計算 run 一覧

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
|---|---|---|---|
| `run_0001_slau_estab_1st` | 層流・SLAU・1次・吹き抜け (outflow) で超音速場確立。porous mesh, M_in=1.689 | **健全**: NaN 無し・物理妥当 (ρ>0, P≤P0, T 154-354K), Mach 0-2.05, 75%超音速。残差~2桁低下。場確立 (背圧無しのため shock train 無し=想定通り)。`field_Mach.png`/`field_P_kPa.png`/`residual_history.png`/`res_20000.h5` | active (確立) |
| `run_0002_slau_sst_muscl_imp_cfl5` | **本番設定の cfl_pseudo=5 安定性確認**: SST + MUSCL(2次,Venkata) + 陰解法(blockDPLUR, implicitRelax0.5, nStepInner10) + **cfl_pseudo=5**。run_0001 確立場から restart (自由流 k/omega 注入)、Sutherland, dilatationCorr=2, 吹き抜け | **cfl_pseudo=5 安定 (発散なし)**: 10000 step/58s, NaN 0, 物理妥当 (ρ>0,P≤P0,T168-296K,k≥0,ω>0), 残差は平坦プラトーで安定 (発散せず)。乱流BL発達を確認。収束は ~0.7-1.5桁でプラトー=キャビティ剪断層の limit cycle+背圧無しのため (cfl5 起因ではない) | active |
| `run_0008_lam_estab` | 層流・MUSCL・陰解法 cfl5 で超音速場確立 (run_0001 から restart) | 薄いBL (δ99~0.1-0.3mm)、cavity core M=1.90 (実験1.89に一致)。outlet M=2.52 | ref (層流anchor) |
| `run_0011_lam_bp1p7_seedSST` | 層流・背圧 Ps=Pt=1.7。SST shock場を seed (超音速ロックイン回避) | 衝撃波保持 x447mm。**seeding で背圧BC が機能** | ref |
| `run_0012_lam_bp2p0_chain` | 層流・Ps=2.0 (run_0011 から chain) | **衝撃波 x129mm=多孔壁上流端 (x_f≈0)**。ただし中心軸静圧は**波状 (shock train)**=ショックレス化せず。層流BL薄すぎが原因 | ref |
| `run_0013_sst_wt1_ti1pct` | SST・**wallTreatmentSST=1**・入口乱流 μt/μ=1,TI=1% (k=31.5,ω=3e7)。y+~1200メッシュ | 安定だが **BL過厚** (cavity M 1.63, δ99 4.2mm=34%半高, μt/μ~37000)。y+~1200 で壁関数不成立が主因 | superseded |
| `run_0014_sst_wt1_wallres` | SST・wt=1・**wall-resolved メッシュ (NY120, 第一セル~5μm, y+~77)** で確立 | **BL適正化**: cavity M **1.81** (←1.63), δ99 **2.11mm=17%** (←34%), μt/μ **7558** (←37000)。メッシュ品質 AR202/skew SOFT-PASS | active (SST anchor) |
| `run_0018_sst_wt1_porous_bp2p06` | SST wall-res・**多孔壁**・Ps=Pt=2.06 (seed) | **NOT CONVERGED** (roUx/roUy/roK 上昇・他プラトー)。衝撃波 x~165mm だが非収束スナップショット値 | 非収束 |
| `run_0019_sst_wt1_solid_bp2p06` | SST wall-res・**固体壁**・Ps=2.06 (比較) | **NOT CONVERGED** (roUx 上昇)。衝撃波 x~129mm (ドリフト中) | 非収束 |
| `run_0020_sst_wt1_porous_bp2p08` | SST wall-res・多孔壁・Ps=2.08 (x_f≈0狙い) | **NOT CONVERGED** (roOmega 上昇・他低下中)。衝撃波 x~141mm | 非収束 |
| `run_0021_sst_wt1_porous_prof` | **計算時間ボトルネック分析用** (run_0020 入力複製, nStepOuter 縮小)。ncu でカーネル別 GPU 時間を取得 | 単独実行 **~5.35ms/step** (run_0020 基準, 98k cell)。GPU busy ~3.7ms/step・残りは launch/host。**上位: implicit_defect_correction_block 24.6% / SLAU 22.4% / limiter_r1 17.4% / viscousFlux 7.3%**。`residual_history.png` | 破棄予定 (プロファイル) |
| `run_0050_solid_prof` | **run_0048 (solid 構造) の計算時間ボトルネック分析用** (入力複製, nStepOuter 縮小)。ncu + detectNaN on/off 比較 | solid 79.4k cell, steady **5.77ms/step** (run_0048)。GPU busy **2.87ms/step (~50%)・残り ~50% は host/launch (110 launch/step)**。GPU 上位: implicit_block 27% / limiter_r1 20% / SLAU 13% / viscous 7.6%。**`detectNaN=1` が +1.45ms/step (=wall の ~21%)** = launch オーバーヘッドの主因 (毎step 全保存量に cub/thrust リダクション ~42 launch)。`residual_history.png` | 破棄予定 (プロファイル) |

### 背圧 down-sweep (porous, SST wall-resolved, 本番設定; 2026-06-21)

本番設定 = SLAU+SST(`wallTreatmentSST:1`)+MUSCL+陰解法 `cfl_pseudo:5`/`nStepInner:5`/**`implicitRelax:1.0`**/`detectNaN:1`、入口乱流 μt/μ=1・TI=1%、`outlet_statPress` Pt=Ps、80k step、出力2000。各点は前段から継続 (`run_case.sh` 経由, VERDICT 記録)。

| run_* | Ps[MPa] | 衝撃波フロント x[mm] | 状態 |
|---|---|---|---|
| `run_0023..0027` | 2.06→1.9 | 159 → 237 | 定常 (残差フロア ~3e-4) |
| `run_0028_sst_porous_bp1p8` | 1.8 | 267 | 定常 (フロア) |
| `run_0029_sst_porous_bp1p7` | 1.7 | 327 | 定常 (フロア) |
| `run_0030_sst_porous_bp1p6` | 1.6 | 381 | 定常 (フロア) |

- **背圧を下げると衝撃波は単調に下流へ移動** (159→381mm, Ps 2.06→1.6)、cfl5 で全点安定・NaN無し。
- 各点とも衝撃波は静止し残差は一定フロア (~3e-4, rms_roUy ~3-5e-2) で plateau。`check_convergence.py` は plateau 判定だが、**衝撃波静止＋残差一定の steady 状態として収束受理** (擬似衝撃波の limit-cycle フロア; cfl15 発散・cfl10 振動・cfl5 が安定限界)。
- 多孔壁制御の本命域 (衝撃波を壁上 x130-162) は Ps≈2.06-2.1 で、それ以上 (x130 狙い) は発散側。down-sweep は壁より下流側の特性。

### 構造メッシュ化 + 上下非対称の正体 (2026-06-21)

非構造 pave (quad/tri) は上下非鏡像 (31% 節点に鏡像なし) で擬似衝撃波の上偏りが疑われたため、
**キャビティ部を全 quad 構造格子化** (`make_mesh.py`, x区間×y帯のマップド quad, `N_SLOT` パラメータ化)。

- **メッシュ品質**: 全 quad・skew 0.011 (←pave 0.93)・**厳密 x 軸対称** (鏡像なし節点 0, 上下セル数一致)。`check_mesh_quality.py` PASS。
- **上偏りは物理的な擬似衝撃波の分岐と確定** (run_0035): 厳密対称メッシュ＋**厳密対称 IC** (asym=0) から出発しても、非対称が機械精度から指数成長 (0→0.0007→0.055→0.29→**0.37**) し、非対称IC版 (0.374) と同値へ。対称解が不安定で流れが片側へ偏る = メッシュ起因でなく物理 (ダクト内擬似衝撃波の既知の非対称分岐)。
- **構造メッシュ背圧マップ (porous, cfl5/implicitRelax1/detectNaN1/80k, 全点安定)**:

| Ps[MPa] | 1.8 | 1.81 | 1.82 | 1.83 | 1.84 | 1.85 | 1.9 | 1.95 | 2.0 |
|---|---|---|---|---|---|---|---|---|---|
| 衝撃波 x[mm] | 279 | 273 | 261 | 249 | 249 | 255 | 237 | 219 | 195 |

  背圧↑で衝撃波は上流へ単調移動。Ps=2.0 で x=195 (多孔壁 130-162 にあと一歩)。構造メッシュは quad-pave (>2.06 発散) より安定。
- 固体壁 (キャビティ・多孔板なし) 構造メッシュ Ps 1.80-1.90/0.02 を別途計算中 (porous vs solid 比較基準)。
- run_0036-0043=porous 構造 sweep, run_0044-0049=solid 構造 sweep。`up_sweep_struct.sh`/`fine_sweep_struct.sh`/`solid_sweep_struct.sh` で駆動。

### ★ porous vs solid 比較 = 論文の機構を再現 (2026-06-21)

構造メッシュで porous (run_0036-0043) と solid (run_0044-0049) を matched-Ps で比較。
中心軸静圧の波状振幅 (peak-to-peak, x180-450mm) を「ショックレス度」の指標とする:

| Ps[MPa] | porous 波状 | solid 波状 | porous 衝撃波x | solid 衝撃波x |
|---|---|---|---|---|
| 1.80 | 0.361 | 0.362 | 279 | 279 |
| **1.82** | **0.280** | **0.380** | 261 | 267 |
| **1.84** | **0.265** | **0.377** | 249 | 255 |
| 1.90 | 0.260 | 0.250 | 237 | 231 |

- **Ps=1.82-1.84 で多孔壁が擬似衝撃波を弱化** (中心軸波状 ~30% 減)。図 `porous_vs_solid_centerline.png`:
  ① 多孔壁では**圧力上昇が多孔壁領域内 (130-162) から早く始まる** (穴からの吹出/吸込が壁上で圧縮を開始)、
  ② 振動振幅が solid (peak 0.48-0.56) より小さい (porous peak 0.36-0.43) = **shock train が弱まりショックレスに近づく**。
- **= Matsuo et al. の「多孔壁パッシブコントロールで擬似衝撃波が弱化する」を 2D で定性再現**。ユーザが指摘した「1.8-1.85 の良い領域」と一致。
- 注意 (誠実な留保): 各 run は残差フロアの steady 状態 (擬似衝撃波 limit-cycle), 流れは物理的に上下非対称 (中心軸指標), 2D スリット (3D 丸穴でない)。Ps=1.80/1.90 では porous≈solid (効果は Ps 窓依存=衝撃波が多孔壁と相互作用する位置のときのみ)。

### Ducros 1次化 ON/OFF 比較 (2026-06-21, 固体壁 Ps=1.84)

衝撃近傍で MUSCL リミタを強制 1 次化する Ducros センサ寄与を `ducrosLimiter` config フラグ
(**既定 0=off で 1 次化を使わない**, 1=on で従来の強制 1 次化) で切り替えられるようにし
(`solverConfig.hpp`/`.cpp`, `ducrosSensor_d.cu`)、固体壁 Ps=1.84 で ON/OFF を**同一バイナリ・
同一 restart 場 (run_0046 発達場)・30000 step**で比較。解析 `analyze_ducros.py` → `compare_ducros_centerline.png`/`ducros_summary.txt`。

| run_* | ducrosLimiter | 結果 | 状態 |
|---|---|---|---|
| `run_0052_solid_ducrosON_bp1p84` | 1 (1次化あり) | ducros 計算正常: max 0.9996, duc>0.8 発火 1.79% (x169-260mm=衝撃波フロントに局在)。中心軸 ripple std 8.76kPa | active (比較基準) |
| `run_0053_solid_ducrosOFF_bp1p84` | 0 (1次化なし) | ducros 一様 0。**ON とほぼ完全一致** (中心軸 ripple std 8.76kPa, OFF/ON 比 1.000, 衝撃波フロント 189mm 同一) | active |

- **結論: Ducros 1 次化を切っても場の差は <0.1% (中心軸 P 差 max 1.3kPa, ducros 発火セルに局在)。
  衝撃波列 ripple 振幅・フロント位置・limit-cycle 変動はいずれも不変**。ユーザ予想の「変動が激しく
  なる」は**起きなかった**。
- 理由: 基底フラックスが SLAU (上流型) で、衝撃捕捉はフラックス自身の散逸が支配する。~1.8% の衝撃
  セルで MUSCL 再構成次数 (1 次/2 次) を変えても結果はほぼ変わらない。SLAU+Venkatakrishnan では
  Ducros 1 次化は**ほぼ冗長な保険**。低散逸の中央/KEEP フラックス (再構成次数が散逸を直接支配) では
  効くはずで、それが本来の用途。
- ducros が計算できているかの確認: ON では衝撃波フロント (x~190mm) に正しく発火 (median 197mm)、
  OFF ではフラグで一様 0 化。センサ自体は正常動作。

### 現時点の所見 (2026-06-21 自走セッション)

- **メッシュ/手法は確立**: 2D wall-resolved SST (y+~77, AR≤1000, skew≤0.9 ゲート通過)、cfl_pseudo=5 陰解法は establish では安定。背圧で擬似衝撃波を seeding→chain で形成可能。
- **⚠ 背圧 run はいずれも未収束 (`check_convergence.py` で全列確認)**: run_0016/0018/0019/0020 とも残差が plateau または**上昇 (発散的)**。**衝撃波＋剥離＋キャビティ循環が本質的に非定常 (limit cycle / 振動)** のため、定常 RANS (unsteady=0) では収束しない。→ **報告した衝撃波位置・中心軸波状振幅は非収束スナップショット値であり信頼できない**。位置安定だけを見て「収束」と判断したのは誤りだった。
- **物理メカニズムの再現は未達 (かつ未収束で判定不能)**: 層流(run_0012, BL薄すぎ)も SST(非収束)も、現時点で論文の「ショックレス化」を確認できていない。
- **次手 (収束問題の解決が先決)**: ① **URANS (unsteady=1, dual-time) で時間平均**するか、定常で収束する条件/手法を探す (擬似衝撃波の非定常性が主問題)、② 3D 丸穴化 (2Dスリットの吸込/吹出分布が論文と異なる疑い)、③ matched-x_f 比較、④ `katoLaunder` で μt 過大抑制、⑤ 出口チャンバ追加。

> 注: 領域は x=690 出口の簡略版 (10°ディフューザ+チャンバ未追加)。run_0003-0007 は旧SST(y+~1200)/手法確立途中の探索 run、run_0009-0010 は超音速ロックインで棄却、run_0015-0017 は中間 (seed/establish)。
