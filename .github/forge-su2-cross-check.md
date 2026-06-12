# forge ↔ SU2 クロスチェック手順

forge の結果(特に軸対称・乱流・近軸挙動)を独立した検証済みソルバ SU2 と突き合わせるための標準手順。
**目的**: forge 固有のスキーム/離散化/数値精度の問題を、同一メッシュ・同一境界条件で SU2 と比較して切り分ける。

> このファイルは [`AGENTS.md`](../AGENTS.md) から参照される。軸対称や乱流の「forge だけ変な値が出る」事象を調べるときは、
> 推測で結論づけず、まず本手順で SU2 と比較すること。

## SU2 の所在と実行

- バイナリ: `.external/su2/bin/SU2_CFD`(v8.5.0 "Harrier", linux64-omp ビルド)。`SU2_SOL` / `SU2_DEF` も同梱。
- 直接実行: `cd <run_dir> && OMP_NUM_THREADS=8 /home/sano/work/forge/.external/su2/bin/SU2_CFD case.cfg`
  - `SU2_RUN` 環境変数は直接 CFD 実行には不要(python ラッパ使用時のみ)。
  - スレッド数は `OMP_NUM_THREADS` で指定。複数 run を並行させるときはコア数(`nproc`)を超えないよう分配する。

## メッシュ生成(forge と同一形状を共有する)

**最重要**: forge と SU2 を**同一の `.geo` から**生成し、「メッシュ違い」を交絡させない。

```bash
# 同一 .geo から両方を出力(gmsh 4.x)
gmsh -2 conical.geo -o conical.msh -format msh2   # forge 用(convertGmshToForge は Gmsh 4.1 形式が必要)
gmsh -2 conical.geo -o conical.su2 -format su2    # SU2 用
```

- **単位スケール**: `.geo` 内の `Mesh.ScalingFactor = 0.001;`(mm→m)は **gmsh が .msh / .su2 両方に適用**するため、
  両者とも SI(メートル)で一致する。SU2 側は `SYSTEM_MEASUREMENTS= SI` のままでよい。
  → 生成後に `.su2` の `NPOIN` 直後の座標が forge メッシュ(`/CELLS/centCoords`)と同レンジか必ず確認する。
- **Physical 名 → SU2 マーカー**: `.geo` の `Physical Curve("inlet"/"outlet"/"wall"/"axis")` がそのまま SU2 のマーカー名になる。
  forge の `bcondConfig.yaml` の physID と対応付けて、下記の cfg マーカーに割り当てる。

## 軸対称 cfg テンプレート(case 29 conical, Pt=4MPa/Tt=1500K/背圧20kPa)

3 本(Euler / 粘性層流 / RANS-SST)を同じメッシュ・BC で回す。差分のみ示す。

共通:
```
AXISYMMETRIC= YES
SYSTEM_MEASUREMENTS= SI
FLUID_MODEL= STANDARD_AIR
GAMMA_VALUE= 1.4
GAS_CONSTANT= 287.058
INIT_OPTION= TD_CONDITIONS
MACH_NUMBER= 0.10
FREESTREAM_PRESSURE= 20000.0
FREESTREAM_TEMPERATURE= 300.0
CONV_NUM_METHOD_FLOW= ROE
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
TIME_DISCRE_FLOW= EULER_IMPLICIT
CFL_NUMBER= 10.0
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.5, 1.5, 1.0, 1.0E6, 0.001, 50 )
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
MARKER_SYM= ( axis )                                   # 軸 = 対称面 (forge: kind: slip/axis)
INLET_TYPE= TOTAL_CONDITIONS
MARKER_INLET= ( inlet, 1500.0, 4000000.0, 1.0, 0.0, 0.0 )   # Tt, Pt, 流れ方向 (x,y,z)
MARKER_OUTLET= ( outlet, 20000.0 )                     # 背圧 Ps (超音速流出時は無視され外挿)
MESH_FORMAT= SU2
MESH_FILENAME= conical.su2
TABULAR_FORMAT= CSV
OUTPUT_FILES= ( RESTART, PARAVIEW, SURFACE_CSV )
VOLUME_OUTPUT= ( COORDINATES, SOLUTION, PRIMITIVE )
```

- **Euler**: `SOLVER= EULER` / `MARKER_EULER= ( wall )`(壁=スリップ)。
- **粘性層流**: `SOLVER= NAVIER_STOKES` / `MARKER_HEATFLUX= ( wall, 0.0 )`(断熱 no-slip)。
  `VISCOSITY_MODEL= SUTHERLAND` / `MU_REF= 1.716e-5` / `MU_T_REF= 273.15` / `SUTHERLAND_CONSTANT= 110.4` / `PRANDTL_LAM= 0.72`。
  (forge `viscMethod: 1` = Sutherland に対応)
- **RANS-SST**: 粘性層流に加えて `SOLVER= RANS` / `KIND_TURB_MODEL= SST` /
  `FREESTREAM_TURBULENCEINTENSITY= 0.01` / `FREESTREAM_TURB2LAMVISCRATIO= 10.0`。
  - 入口乱流: forge の `k`/`omega` 直接指定と、SU2 の(乱流強度 + 粘性比)は厳密一致しない。
    目標 μt/μlam を合わせる(本ケースは ~10)。絶対量は多少ずれても、**軸中心 k の「スパイク有無(形状)」の比較**には十分。

## 収束確認(必須・[AGENTS.md] の収束ルールに従う)

SU2 のノズル流は **リミットサイクル**で残差が下げ止まることが多い。`rms[Rho]` だけでなく:
- `history.csv` の `rms[Rho]`, `rms[Momentum-Y]`, `rms[TKE]`, `rms[Dissipation]` の**トレンド**を見る(単調減 or 振動プラトーか)。
- 積分量(出口 massflux: `MARKER_ANALYZE_AVERAGE= MASSFLUX`)が定常化しているか。
- 早期の `rms[Rho]` だけ低くても、運動量・乱流残差が下げ止まり/上昇していれば**未収束**と判断する。

## VTU の読み取り(注意)

SU2 v8 の `flow.vtu` は `NumberOfComponents= "3"`(= 後の空白)など属性が非標準で、**VTK の `vtkXMLUnstructuredGridReader` が失敗する**ことがある。
その場合は appended-raw バイナリを手動パースする(`header_type="UInt64"`、各 DataArray は 8 byte 長さ接頭辞 + Float32 ペイロード、offset は `<AppendedData>` の `_` 直後からの相対)。
実装例: [`case/29.bell_vs_conical/compare_forge_su2.py`](../case/29.bell_vs_conical/compare_forge_su2.py) の `load_su2()`。

## 既知の比較結果(case 29 conical, 2026-06)

軸対称 SST の「軸中心 k スパイク」を SU2 と比較し、**forge 固有・近軸の数値問題**であることを特定:
- SU2(頂点中心・軸上に節点・倍精度)は軸中心の半径速度 `u_r` が 0 から滑らかに立ち上がり、k は軸で最小(スパイク無し)。
- forge(セル中心)は **float32 の陰解法(block-DPLUR)が近軸第一セルの `u_r` を収束させきれず**(陽解法・倍精度では正しい)、
  偽の `∂u_r/∂r` → 偽ひずみ → SST 生産で k がスパイク。フープ項・Kato–Launder は無関係(下流の対症療法)。
- 詳細: [`.github/plans/architecture-axisym-axis-singularity.md`](plans/architecture-axisym-axis-singularity.md)。
