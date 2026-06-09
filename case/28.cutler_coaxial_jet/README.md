# 28.cutler_coaxial_jet — Cutler 超音速同軸 He/空気ジェット (多成分 TP 検証)

NASA Langley の超音速同軸ジェット (Cutler 他, AIAA-99-3588, NTRS 20040086859) を、多成分
thermally-perfect gas + 組成依存入口 BC + 乱流化学種拡散の検証ケースとして構築する。

## 実験条件 (AIAA-99-3588)

| | 中心ジェット | coflow |
| --- | --- | --- |
| ガス | He 95% + O2 5% (体積) | 空気 |
| Mach | 1.8 (nominal) | 1.8 (nominal) |
| 出口圧 | 1 atm (両者圧力整合) | 1 atm |
| 形状 | 軸対称、中心ジェット出口径 ~10 mm (r≈5.12 mm)、center body lip r≈1.74 mm | |

- 対流マッハ数 $M_c=|u_j-u_\infty|/(c_j+c_\infty)=0.7$。中心ジェットは音速が大きく速度は coflow の 2 倍超。
- 全温 ~300 K (冷流)。検証データ: He モル分率分布、pitot 圧、全温の probe survey (RELIEF/GSAS)。
- **非反応** (He は燃料 simulant)。化学はスコープ外なので好都合。

## 進捗

- **メッシュ** (`mesh/coaxial.geo`): 軸対称 2D 同軸 (中心ジェット/lip/coflow/ambient の 4 バンド構造格子)。
  境界 physID = axis(1)/jet_in(2)/lip(3)/coflow_in(4)/amb_in(5)/top(6)/outflow(7)。38240 セル。
  → **`convertGmshToForge` で一発変換成功、forge も NaN なしで稼働** (`run_0001_meshcheck`, CPG 単成分で
  流れ場が立つことを確認。残差 3.4e-5→1.8e-5、P/T/Ux 有限)。
- **組成依存超音速入口 BC** (plan §10 / M5, 実装・検証完了): `inlet_uniformVelocity` に入口組成 `Y_s`
  を追加 (config `floats: Y0,Y1,...`)。ghost は $\rho=P/(R_{mix}(Y)T)$, $roe=\rho(e_{mix}(Y,T)+e_k)$,
  $\rho Y_s=\rho Y_s^{in}$ (Dirichlet)。あわせて SLAU 対流流束の被移流エンタルピーをセル組成で再構成
  (He/N2 で cp 質量基準が ~5×差のため `sp[0]` 固定だと混合 contact で圧力発散していた)。
  - `run_0003_he_n2_inviscid/` — binary [He,N2], 1次風上・非粘性, 2000 步。NaN/発散なし、$\sum Y=1$、
    中心 He コア (Y_He=1, Ux~1280 m/s = Mach1.8 He)、coflow N2。
  - `run_0008_he_n2_2nd/` — 上を 2次風上でリスタート (計 4000 步)。He ポテンシャルコアが ~0.07m へ伸長。
    成果物: `field_YHe_Mach.png`, `residual_history.png`。
  - `run_0002_cpg_regr/` — CPG 単成分回帰 (修正が単成分/CPG を壊さないことの確認)。
- **startup 注意**: 静止 IC (u=0) + `outlet_statPress` は出口の backflow 分岐で `Pt`/`Tt` 未指定だと NaN。
  IC に下流速度 (u=100) を与え、出口 bcond に `Pt`/`Tt` を指定すること (`gen_binary_ic.py` 参照)。
- **次 (後続)**: 定量照合 (He コア長・濃度減衰 vs RELIEF survey) は粘性 + RANS-SST + 乱流化学種拡散が必要。
  その後 `[He,O2,N2]` 3 成分へ拡張。

## ガス物性 (TP 用)

- 中心ジェット He-O2 (95/5 vol): MW≈5.40 g/mol, γ≈1.65 (He 単原子 + O2 二原子の混合)。
- coflow 空気: N2/O2。簡易には binary He/N2 から始め、後で He-O2/air へ拡張。
