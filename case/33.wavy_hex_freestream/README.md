# 33. wavy hex free-stream 検証

正弦波で歪ませた非直交 hex (10^3) で forge の free-stream(一様場)保存を検証する。
非構造/非直交メッシュで forge が発散する原因が「free-stream 保存崩れ(metric誤差×流束)」か
「別機構(エネルギー式等)」かを切り分ける。

- mesh/wavy.h5         : forge 入力メッシュ(閉性 1.28e-9・体積厳密に検証済み, 純幾何は正しい)
- mesh/wavy_fluent.h5  : 合成 Fluent 形式の元データ
- 生成: /tmp/make_synth_fluent.py の gen_wavy (内部ノードを amp*h*sin で歪ませた構造hex)

検証方針(別LLM助言 + 自前): 一様静止場の forge 残差を 1ステップ逆算し、
F(U_c)·r_c (r_c=Σ dir·sv = セル幾何閉じ誤差) と成分ごと・符号付きで突き合わせる。
  - 一致(運動量=p·r, 質量=0, エネルギー=0) → free-stream が主原因
  - エネルギー等に説明できない残差 → 別機構

| `run_0004_keep_pref0` | **KEEP** (pRef なし) 一様静止 | 運動量残差 7.2e-5 → **step11 発散** = float32 free-stream 崩壊を KEEP でも実証 | ref (KEEP 症状記録) |
| `run_0005_keep_pref` | **KEEP + pRef=101325** | 運動量残差 **3.96e-12 (機械精度)**・300step 安定。SLAU と同水準の根治 | ref (KEEP pRef 検証) |
| `run_0006_keep_dissmat_pref` | KEEP + **matrix ES 散逸 (keepDissType=2)** + pRef | 4.07e-12・安定。一様場で Δw=0 → 散逸レイヤは free-stream を厳密に保つ | ref (散逸レイヤ無害確認) |
