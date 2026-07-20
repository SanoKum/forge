# case/38.channel_wmles — WMLES 周期チャネル (Reτ=550 較正 → 2000)

WMLES 代数壁応力モデル ([turbulence-wmles-wall-stress plan](../../plans/active/turbulence-wmles-wall-stress.md) §6.2)
と IDDES の WMLES モード定量評価 ([turbulence-iddes-sst plan](../../plans/active/turbulence-iddes-sst.md) §4.6)
の共通検証ケース。参照は Hoyas & Jiménez (Reτ=550/2000) の平均速度 log 層。

## 構成 (Reτ=550 最小構成)

- 箱: $2\pi\delta \times 2\delta \times \pi\delta$ (δ=0.005 m)、一様格子 48×24×32 (Δx⁺=72, Δy⁺=46, Δz⁺=54)
- x/z 周期 (node Cartesian periodic)、y 両壁 `wall_isothermal` (Ts=300K) + `wallModelLES: 1`
- 駆動: `bodyForce: [3488.7, 0, 0]` = ρu_τ²/δ → **u_τ=3.85 を直接指定** (Reτ=u_τδ/ν=550, ν=3.5e-5)
- 流体: CPG 空気系 (visc 4.12e-5 定数, Pr=0.72), u_b≈70 m/s (M≈0.2)
- スキーム: 確定版 KEEP-LES スタック (KEEP + matrix ES σ0.05 jump2 precond + advGauge) + WALE, node, unsteady RK3 cfl1
- 1 flow-through ≈ 400 step。メッシュ品質: 一様直方体 (AR≈1.6) につき自明 PASS

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_ret550_smoke` | スモーク (2000 step ≈5 FT): 安定性・WMLES 稼働・周期/体積力/等温ピンの統合確認。IC = Reichardt 平均 + 12% 撹乱 | **完走・NaN 無し・WMLES フル稼働** (utau med 2.76→平衡 3.85 へ発達途上, y+≈25, τ_w=ρu_τ² 整合, q_w<0=散逸熱が壁へ)。撹乱は減衰中 (Uy rms 2.4→0.44) だが減衰は鈍化 | ref (統合スモーク) |
| `run_0002_ret550_dev` | res_2000 から 12000 step 継続 (~35 FT): 乱流自己維持 or 層流化の判定 | **乱流自己維持を確認**: Uy rms 0.50→1.75→2.20・Uz rms 0.95→2.63→3.15 と再成長 (遷移成功)。u_bulk 70→79 (τ_w 平衡へ発達中・オーバーシュート局面)。NaN 無し | ref (遷移確認) |
| `run_0003_ret550_stat` | res_12000 から 30000 step (~75 FT): 統計平衡 (u_τ→3.85) + 時間平均統計収集 (snap 1000 毎) | (実行後に記入) | active |
