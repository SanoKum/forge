# 37. pintle_nozzle — ピントル可変スロートノズル (3D)

下方から立ち上がる円形給気管が 90°ベンドして収縮拡大ノズルへ入り、ノズル軸上に
ピントルが左端壁から +x に刺さる構成の 3D 解析。ピントル先端形状・位置がノズル効率
(推力・比推力効率) に与える影響を調べる。

参照論文: `papers/pintle_nozzle/Pan_Chen_Ye_2017_pintle_tip_shape_specific_impulse_AtlantisPress.pdf`
(Pan, Chen, Ye 2017)。論文自体は **2D 軸対称・単一理想気体** (γ=1.274, Tc=2889K,
超音速出口外挿) でピントル先端5種・位置3種を比較。本ケースは 90°ベンド給気で軸対称が
破れる実機形状を **3D** で扱う (対称面 1 枚で半割)。

## 形状と座標系

- `x` = ノズル軸 (流れは +x に排出)、`z` = 鉛直 (給気管は下から +z で立ち上がる)。
- 対称面 = `y=0` (x-z 平面)。半割は `y>=0` を残す。
- スロート環状すきま = `Rt - Rp` が律速部。

## ワークフロー

2 経路ある。**prism 境界層が要るなら Salome 経路 (推奨)**、tet のみで良ければ Netgen 経路。

```
[A: Salome 経路 (prism BL 付き, 推奨)]
FreeCAD (cad/build_geom.py) -> pintle_fluid_half.step
   -> Salome headless (cad/mesh_salome.py: NETGEN + SMESH ViscousLayers) -> MED
   -> cad/med_to_msh41.py -> forge.msh -> convertGmshToForge -> forge.h5

[B: Netgen 経路 (pip のみ, tet + 近壁細分)]
FreeCAD (cad/build_geom.py) -> pintle_fluid_half.step
   -> Netgen (cad/mesh_pintle.py: 面命名 + tet + WALL_MAXH 近壁細分 + msh4.1書出)
   -> forge.msh -> convertGmshToForge -> forge.h5
```

- **Salome 9.14 は直リンクで取得可能** (登録不要; `files.salome-platform.org` は **Referer
  ヘッダ必須** — 無いと 403)。`/home/sano/opt/salome/` に導入済み (「環境」参照)。
- **メッシュ設定は `mesh_salome_config.json`** (テンプレ: `cad/`, 単位 mm)。実行 cwd に置くと
  既定値を上書きし、使用値が `mesh_settings_used.json` として出力先に保存される
  (run ごとの設定トレーサビリティ)。ピントル先端の層欠けは 2 要因あった:
  ① **先端キャップ (底面) の symmetry 誤分類** (主因; BC も slip になる実害バグ。分類順の
  修正で解消) ② ノーズ強曲率での層縮退 (`tip_maxh_mm=0.10` の局所細分で解消;
  0.07 は表面自己交差で失敗)。両修正後は wall_pintle 層あり率 100%。
- netgen 生 API (pip) の prism BL はこの形状の凹角で破綻するが、**SMESH の
  StdMeshers_ViscousLayers は凹角の層縮退処理を持ち、同じ形状で成功する** (prism 2.4万生成)。
- 面の境界グループ分けは両経路とも面重心の幾何条件で自動分類 (対話ピック不要・再現可能)。

**単位 (重要)**: CAD/STEP は mm で作るが **forge は SI (m) 前提**。
- Salome 経路: ImportSTEP が mm->m 変換するため MED/forge.h5 は自動的に m (正しい)。
  `mesh_salome.py` は bbox から unit_scale を自動検出してしきい値を合わせる。
- Netgen 経路: `mesh_pintle.py` の `SCALE=0.001` (mm->m) で書き出す (`--scale` で変更可)。
  SCALE 1.0 で書くと 35mm ノズルが 35m になる (初期 run_0001/0002 はこのバグ持ち -> 再生成済)。

- **形状**: `cad/build_geom.py` をパラメトリック編集 (`freecadcmd build_geom.py`)。
  clean な単一閉ソリッドの半割 STEP を出す (OCC boolean の取りこぼし対策に fuzzy fuse +
  `cut(y<0)` 半割 + bbox 検証ループ済)。ノズル輪郭は当面コーン近似、ピントル先端は blunt
  (尖らせると Netgen 表面メッシュが退化頂点で破綻)。給気は垂直円柱でチャンバー下面へ
  (sharp/smooth エルボは sliver/boolean 破綻のため将来改良)。
- **メッシュ**: `cad/mesh_pintle.py` (`netgen-mesher`)。`wall_pintle` を別グループにして
  ピントル軸力を面積分できるようにする。近壁解像度は `WALL_MAXH` (壁 face の局所細分・堅牢)
  を推奨。prism 境界層 (`USE_BL`) は現状この凹角形状では失敗する (下記「現状と次タスク」参照)。
  `MAXH` は throat 環状すきま (Rt-Rp=1mm) を割る値に。
- **取り込み**: Netgen の "Gmsh2 Format" 出力は `$PhysicalNames` を書かず convertGmshToForge が
  読めないため、`mesh_pintle.py` が forge 必須要件 (下記) を満たす **gmsh 4.1 を直接書く**。
- **forge 制約**: float32 の非直交 free-stream 保存が弱い。**1次風上 + 陰解法 (block-DPLUR)
  + `pRef`=チャンバー静圧 + `outlet`=超音速外挿** の安定レシピで投入。投入前に
  `check_mesh_quality.py` で AR<=1000 / skew<=0.9 (PASS) を確認する。

## 環境 (sudo 不要・WSL native)

| ツール | 用途 | 所在 / 起動 |
| --- | --- | --- |
| FreeCAD 1.1.1 | 形状作成 (STEP) | `/home/sano/opt/squashfs-root/usr/bin/freecadcmd` (AppImage展開済, `QT_QPA_PLATFORM=offscreen`) |
| SALOME 9.14.0 | SMESH ViscousLayers (prism BL) | `/home/sano/opt/salome/SALOME-9.14.0-native-UB24.04-SRC/salome` (下記「Salome 起動の注意」) |
| netgen-mesher 6.2 | メッシュ (tet) | venv `/home/sano/work/forge/.venv-mesh` |
| meshio 5.3 + h5py | MED 読取・変換 | 同 venv |
| gmsh | 2D確認・形式変換 | `/usr/bin/gmsh` (既存) |
| convertGmshToForge | msh4.1 -> forge.h5 | `solver_density_cuda/build-native/convertGmshToForge` (tet/prism/pyramid 対応) |

### Salome 起動の注意 (sudo 無し環境)

- DL は Referer 必須: `curl -e "https://www.salome-platform.org/" -A "Mozilla/5.0" <files.salome-platform.org のURL>`。
- native ビルドはシステム python3 とシステム .so を使う。sudo 無しでは:
  - `pip install --user --break-system-packages psutil` (KERNEL 必須)。
  - 不足 .so (boost 等) は `apt-get download <pkg> && dpkg -x` で
    `/home/sano/opt/salome/syslibs/` に展開済み。
- **起動は `--keep-paths` + LD_LIBRARY_PATH が必須** (launcher が環境を上書きするため):

```bash
export LD_LIBRARY_PATH=/home/sano/opt/salome/syslibs/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
/home/sano/opt/salome/SALOME-9.14.0-native-UB24.04-SRC/salome --keep-paths -t cad/mesh_salome.py
```

```bash
# 形状
QT_QPA_PLATFORM=offscreen /home/sano/opt/squashfs-root/usr/bin/freecadcmd cad/build_geom.py
# メッシュ -> forge.msh  (mesh venv で)
source /home/sano/work/forge/.venv-mesh/bin/activate
cd cad && python3 mesh_pintle.py
# 取り込み (要 solverConfig.yaml + bcondConfig.yaml; physID は下表と一致させる)
convertGmshToForge forge.msh forge.h5
python3 ../../../solver_density_cuda/tools/check_mesh_quality.py forge.h5
```

## 寸法 (PARAMS) — 論文スケール仮組み

`cad/build_geom.py` の `P` に設定済 (Pan/Chen/Ye 2017 Fig.2 基準)。**確実な値は L=35mm /
ER=3.0 / 2.5=Rt のみ**。図は径方向が非スケールのため Rc・収縮部は推定、給気管/ベンドは
論文に無く暫定。

| 量 | 値 [mm] | 確度 |
| --- | --- | --- |
| ノズル全長 x_exit | 35 | ◎ 論文 |
| スロート半径 Rt | 2.5 | ○ 論文 2.5 を Rt と解釈 |
| 出口半径 Re | 4.33 (=Rt√3, ER=3) | ◎ 論文 |
| スロート位置 x_throat | 15 (全長の~42%) | ○ 軸方向 |
| チャンバー半径 Rc | 6.0 | △ 推定 (収縮比 2.4) |
| チャンバー直管長 Lc | 5.0 | △ 推定 |
| ピントル胴半径 Rp | 1.5 | △ 2.92 vs 2.5 の閉塞から |
| 給気管半径 Rpipe / Xfeed / z_bot | 3.0 / 4.0 / -20 | ✕ 暫定 (論文に無い; 実機値で要更新) |

γ=1.274 / Tc=2889K / 理論Isp 2212 m/s (solverConfig の物性に反映すること)。

## 現状と次タスク (2026-07-07)

- ✅ 環境構築・パイプライン疎通: FreeCAD→STEP→Netgen→msh4.1→convertGmshToForge→forge.h5
  (`nBconds=5` 認識、全境界 bname 一致) を実証済。
- ✅ 形状: 論文スケール PARAMS で clean な半割 STEP (x[0,35], Rc=6, 給気管 z=-20) を決定的生成。
  MAXH=0.4 で表面自己交差なくメッシュ (tet ~6万)。
- ✅ **メッシュ品質 PASS**: `check_mesh_quality.py` を 3D 対応に拡張 (完全3D は x-y 射影でなく
  面ごと equiangle skew / 全辺 AR で判定; 準2D は従来の x-y 射影を維持)。pintle mesh は
  **3D モードで AR max 2.9 / skew max 0.657 / PASS**。旧 FAIL (skew 36%・AR=inf) は x-y 射影
  アーティファクトで実体でなかった。→ 品質チューニングは不要。
- ✅ **prism 境界層は Salome (SMESH ViscousLayers) で成功** (`run_0003`, prism 2.4万)。
  netgen 生 API (pip) の `boundary_layers` は同じ形状の凹角 (T字給気・スロート・ピントル基部)
  で失敗する (`project_boundaries` 無し=segfault、有りでも体積充填 "too many attempts") ため、
  BL には Salome 経路を使う。Netgen 経路の代替は **近壁 tet 細分** (`WALL_MAXH`, 壁関数前提)。
- ⏳ 実機の給気管/ベンド実寸が分かれば PARAMS を更新 (現状は暫定)。
- ⏳ solverConfig (安定レシピ + γ=1.274/Tc=2889K)・bcondConfig 実値 (全圧/全温) で forge 投入。

## メッシュ取り込みの必須条件 (convertGmshToForge)

リーダ (`solver_density_cuda/mesh/gmshReader.hpp`, `elementType.hpp`) を確認した結果:

- **tetra/prism/pyramid/hex すべて対応**。Salome の tet内部 + prism境界層 + pyramid遷移は
  そのまま読める (連結面マップ・体積計算とも prism/pyramid を処理)。
- エクスポート時の必須条件:
  1. **gmsh 4.1 形式 (msh4)** で書く (リーダは構造的に 4.1 の `$Entities` をパースする)。
     meshio なら `-o gmsh` (4.1)。`gmsh22` は不可。
  2. **1 エンティティ = 物理タグちょうど 1 個** (面/体が 2 個以上だと `exit`)。
     1 つの面を 2 グループに入れない。
  3. **流体ボリュームに 3D 物理グループ `fluid` が必須** (体の物理タグが 1 個でないと `exit`)。
     これが `is3D` 判定のトリガにもなる。
  4. **外表面は全て面グループで覆う** (未タグの外面は BC が付かない)。
     `inlet+outlet+wall+wall_pintle+symmetry` で全外面を覆う。
- **Second Order (中間節点) は不可**。forge は線形要素のみ。NETGEN の "Second Order" は OFF。

## 境界グループ -> forge BC

| 形状の面 | Salome グループ | forge BC |
| --- | --- | --- |
| パイプ下端の円盤 (z=z_bot) | `inlet` | 亜音速圧力入口 (Pt/Tt: 燃焼室全圧・全温) |
| ノズル出口の円盤 (x=x_exit) | `outlet` | 超音速出口 `outflow` (外挿) |
| チャンバー/ノズル/パイプ/ベンドの壁 | `wall` | 断熱 no-slip 壁 |
| ピントル表面 | `wall_pintle` | 断熱 no-slip 壁 (軸力評価用に分離) |
| 半割の平面 (y=0) | `symmetry` | 対称 |

## 方針

ユーザー方針により **3D 半割を直接** 進める (2D 軸対称の段階は踏まない)。
forge 投入で非直交由来の不安定が出て切り分けが要るときは、参考として
`case/23.axi_nozzle` / `27.axi_nozzle_plume` の 2D 軸対称 (γ=1.274・Tc=2889K) で
論文 Table 1 (標準 94.81% / small arc 93.30%) を当て、ソルバ設定の正しさを確認する。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_slau_baseline/` | 論文スケール形状の一様 tet (Netgen maxh=0.4)。solverConfig は仮 (cell-mode placeholder, 物性/BC 未実値) | `forge.msh` (14002節点/60322 tet), `forge.h5` (m 単位・再生成済), 品質 **PASS** (3Dモード AR2.9/skew0.657) | active (**未投入**; ⚠キャップ symmetry 誤分類バグ持ち・使う前に再生成) |
| `run_0002_slau_wallref/` | 近壁 tet 細分版 (`mesh_pintle.py --wall-maxh 0.12 --maxh 0.5`)。等方細分 (prism BL ではない, 壁関数前提) | `forge.msh`/`forge.h5` (~54万 tet, m 単位・再生成済), 品質 **PASS** | active (**未投入**; ⚠キャップ symmetry 誤分類バグ持ち・使う前に再生成) |
| `run_0003_slau_salome_bl/` | Salome SMESH ViscousLayers の prism 境界層付き初版 (総厚0.25mm/4層, maxh0.5)。先端層欠け+**キャップ symmetry 誤分類 (BC slip) バグ持ち** | `pintle_salome.med`, `forge.h5` (tet 40022 + prism 24492), 品質 PASS | ref (比較用; 計算には使わない。本命は run_0004) |
| `run_0004_slau_bl_tipref/` | **本命 BL メッシュ**: 先端局所細分 (`tip_maxh_mm=0.10`; 0.07 は表面自己交差で不可) + **面分類バグ修正** (先端キャップが symmetry 誤分類→BC slip・層除外になっていた; wall_pintle を先に判定)。**先端キャップ含め wall_pintle 層あり率 100%**。設定 `mesh_salome_config.json` (使用値 `mesh_settings_used.json`) | `pintle_salome.med`, `forge.msh`/`forge.h5` (33332節点: tet 64557 + **prism 36384**, m 単位), 品質 **PASS** (AR20.6/skew0.849) | ref (メッシュ正本; 計算は run_0005 以降) |
| `run_0005_slau_trial/` | 試計算 第0段: Euler+slip+1次+陰解法(cfl_pseudo1)+pRef3e6, IC=滞留一様(p=Pt=3e6,T=2889K), outlet=`outflow` | 5000step 完走・残差機械ゼロ = **滞留平衡が完全保存** (outflow+平衡ICは流れを駆動しない)。メッシュ/BC/pRef の保存性 sanity 合格 | ref (平衡保存の確認) |
| `run_0006_slau_trial_diaphragm/` | ダイヤフラム IC (スロートで 3e6/3e5 仕切り) | **step4 NaN** (ピントル先端すきまの薄 prism にジャンプ直撃)。ダイヤフラム起動はこの形状では不可 | ref (失敗診断; res_nan_4.h5) |
| `run_0007_slau_statpress_pr15/` | **段階起動 第1段**: outlet を `outlet_statPress` Ps=2e6 (PR1.5, 非チョーク) に変更, IC=滞留一様 | ✅ 5000step 完走・NaN 0。M max 0.93, P[1.29e6,3.05e6], 場は物理的。VERDICT: NOT CONVERGED (発展途中) | active (第1段の場) |
| `run_0008_slau_statpress_pr30/` | 第2段: Ps=1e6 (PR3) restart | ✅ チョーク+拡大部衝撃 M max 2.32, NaN 0。VERDICT: NOT CONVERGED (衝撃が出口付近で plateau) | active (第2段の場) |
| `run_0009_slau_outflow_full/` | outflow へ早期切替 (出口リップに亜音速ポケット残存のまま) | **step1109 NaN** (出口リップ r=2.7-4.2mm)。亜音速セルの全量外挿は ill-posed → 切替は出口全面超音速後に | ref (失敗診断) |
| `run_0010_slau_statpress_pr100/` | 第3段: Ps=3e5 (PR10) restart → 衝撃を出口外へ | ✅ **出口全面超音速 (M min 2.25)**, M max 2.46, P min 1.41e5≈設計出口圧。全残差下降 (2.3-4.1桁, still converging), NaN 0 | active |
| `run_0011_slau_outflow_super/` | **最終段: outlet=`outflow`** (論文と同じ超音速外挿) restart | ✅ 10000step 完走・NaN 0。場は 0010 と実質不変 (M max 2.46, T max=Tt=2889K 物理的)。VERDICT: NOT CONVERGED (残差は極低絶対値 rms_ro~2e-10 で plateau) | active (cell 試計算の最終場) |
| `run_0012_node_euler/` | **node (median-dual) 化**: node 変換 (dual閉性1.9e-6) + run_0011 場を3D最近傍移植。Euler+slip+1次+陰解法 | ✅ **VERDICT: PASS (converged)** 全列4.2-5.3桁下降。M max 2.38 | active (node ベースライン) |
| `run_0013_node_visc/` | ② 粘性 no-slip (μ=9e-5, k=0.215, `nodeWallDirichlet=1`) restart | ✅ **PASS (converged)** 3.1-4.0桁。**壁ノード4739個 \|U\|=0.000** (no-slip 厳密)。M max 2.34 | active |
| `run_0014_node_sst/` | ③ SST (k/ω 全域非ゼロ init k=4/ω=500, wallTreatmentSST=1) restart | ✅ roK **5.9桁**/roOmega **6.0桁** 下降・NaN 0。μt/μ max 38。VERDICT: NOT CONVERGED は ro/roe が restart 時点で float32 床 (~e-11) のための表示 | active |
| `run_0015_node_sst_2nd/` | ① 2次精度化 (convMethod 1 + Venkat, bndFirstOrder=1) restart | ✅ 全列3-4桁後の低レベル plateau (2次+SSTの微振動)・NaN 0。**M max 2.52** (2次化で膨張シャープ化) | active (ノズル単体の最終場) |
| `run_0016_node_plume/` | **外部プルーム領域追加** (円筒 L60/R25mm, 8境界: +plume_out/plume_far/base)。SMESH 再メッシュ (tet54k+prism25k+pyram176, lip局所0.3mm, SOFT-PASS 0.01%)。IC=ノズル0015移植+プルーム大気 (x>35mm のみ指数ブレンド)。plume_far=`outlet_statPress` (双方向) | ✅ 10000step 完走・NaN 0。**プルームに不足膨張ジェットの衝撃セル構造** (P 0.75-1.75e5 振動)、M max 2.52。VERDICT: NOT CONVERGED (ジェット発達中) | active |
| `run_0017_node_plume_cont/` | 0016 の継続 (20000step) — ジェット発達 | ✅ 完走・NaN 0。M max 2.53, μt/μ max 234 (せん断層)。ただし **rms_roOmega が step~4000 から 1.33e16 で凍結** = 出口リップ (x=35, r=4.3mm; skew0.915+pyramid+BL終端) の少数ノードで ω がリミッタにピン留め。場は局所以外健全 | ref (旧プルーム形状; 後継 run_0018) |
| `run_0019_node_plume_dump/` | 0018 の短継続 (2000step) + **全境界ダンプON** — 力積分の完全収支用 | ✅ **推力 (構造力) Fx=-64.82N (壁+ピントル+base)** vs 出口断面CV 63.70N = **1.8%一致**。対称面 Fy 釣り合い 0.8%。**ピントル軸力 -10.23N (上流向き; 圧力-10.40+摩擦+0.17)**。ṁ=0.031 kg/s | active (力収支のリファレンス) |
| `run_0018_node_plume_big/` | **プルーム改良版**: 領域2倍 (L120×R50mm)・全体1/1.5細分 (maxh2.0/ノズル0.33/リップ0.2; tip0.10は自己交差回避で据え置き)・**base (リップ接続壁) にもレイヤー** (`bl_on_base=True`, BC も wall 化)。IC=0017場移植+拡張域のみ大気ブレンド | ✅ 20000step 完走・NaN 0。**品質 完全PASS (skew max 0.851, 違反0, pyramid 0)** — base レイヤーがリップ角を巻き旧 SOFT-PASS 因子を解消。M max 2.59, μt/μ max 422。tet386k+prism71k (11.3万節点)。ω のリップピン留めは残存 (ω~1.6e8; 微小セル壁 ω として物理オーダー内だが残差 e15 を支配) | **active (最新プルーム場)** |

> メッシュ本体 (`forge.msh`/`forge.h5`) はこの run ディレクトリにある (git には入れない)。
> 再生成は `cad/build_geom.py` → `cad/mesh_pintle.py` → `convertGmshToForge`。
> run を作成・破棄したらこの表を必ず同期する (命名 `run_NNNN_<scheme>_<slug>`)。
