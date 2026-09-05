# codex レビュー: 有限速度化学の現状評価 (2026-09-05)

codex-cli 0.153.4 (read-only) による第三者レビュー。依頼元セッションの記録。対象: 引継ぎ文書・Cabra 文献調査・case/47,48 run 表・chemistry plan。

## (1) 到達点の総評

結論は、**反応カーネルの実装検証としては高水準だが、燃焼加熱器を予測的に設計できる検証水準には未到達**です。

探索的な加熱器計算へ進むことはできますが、結果の用途は感度把握・計算成立性確認に限定すべきです。寸法、火炎安定性、壁熱負荷などの設計判断にはまだ使えません。

| 検証 | 評価 |
| --- | --- |
| 0-D | 良好。同一機構の Cantera に対し着火遅れ 32.0/32.2 µs、平衡温度差約1 K。GPU周期境界・反応熱・Strang・ΣY保存も確認され、実装検証として十分 |
| Q1D | 化学非平衡の実装検証として良好。OH/NO は Cantera PFR とよく一致。ただし `check_convergence`: **NOT CONVERGED** で、T/Mの絶対比較も成立していない。[case/46 README](/home/sano/work/forge-chem/case/46.nozzle_h2o2_kinetics/README.md:23) |
| BK | 出口の組成ピーク・位置は良いが、入口 `k×3` 等で合わせた結果。全温ピークは1.08対1.18、着火は約4 cmまで上流ドリフトし実験18–25 cmと大差。`check_convergence`: **NOT CONVERGED**、SST残差プラトー。[run_0028](/home/sano/work/forge-chem/case/47.burrows_kurkov/README.md:43) |
| Cabra | 低Mach反応経路のバグ発見・安定性検証として非常に有益。しかし No-TCI/PaSRとも付着し、全ケース **NOT CONVERGED**。下流半径温度の部分的一致だけで、中心軸種分布や上流温度は不一致。[case/48 README](/home/sano/work/forge-chem/case/48.cabra_h2n2/README.md:47) |

主要最終 HDF5 と残差CSVには NaN/Infがなく、`ro>0`、`T>0`、ΣY≈1でした。一方、公式ツールの再実行では次のとおりです。

- Q1D: `NOT CONVERGED`、`rms_roUy` stalled
- BK: `NOT CONVERGED`、`rms_roK`・`rms_roOmega` stalled
- Cabra混合場: `NOT CONVERGED`
- Cabra PaSR C_mix=1: `NOT CONVERGED`、`rms_roUx`・SST残差が末尾 `RISING`
- Cabra PaSR C_mix=4: `NOT CONVERGED`、`rms_roK` が `RISING`

したがって、計画の「全残差PASS」という完了基準は満たしていません。[plan §6](/home/sano/work/forge-chem/plans/active/chemistry-finite-rate-h2.md:58)

未解決事項の優先順位は少し組み替えるべきです。

1. SSTプラトー／準定常性
2. 反応機構の着火遅れ妥当性
3. CabraのTCIモデル不足
4. BK固有の早着火条件
5. 性能最適化

BK壁・入口条件の追加調整を先行するより、数値定常性と機構依存性を先に閉じる方が合理的です。

## (2) 2択への推奨と理由

**燃焼加熱器の「設計能力」を目標とするなら、TCI実装を推奨します。**  
ただし PaSR/EDC の追加調整ではなく、別計画として **RANS 1st-order CMCを第一候補**にします。Transported PDFはより一般的ですが、粒子統計、混合モデル、ISAT相当、RANSとの整合など実装範囲がかなり大きくなります。

理由は次のとおりです。

- BKとCabraの両方で発熱開始位置が大幅に上流です。これは新しい加熱器でも燃焼長、壁熱流束、圧力損失、局所最高温度を誤る可能性を示します。
- BKは出口組成ピークこそ合っていますが、熱的プロファイルと着火位置は合っておらず、「正しい理由で正しい出口を得た」とは言えません。
- Cabraの下流一致も温度の一部に限られ、未収束です。現状の根拠だけで「下流温度・組成で閉じた」とするのは弱いです。
- 文献調査の主要結論は妥当です。Cabra火炎は化学律速で、10 Kの変化でリフトオフが倍になり得ること、混合モデル定数より機構・温度の影響が大きいことが一次文献でも示されています。[Cao–Pope–Masri 2005](https://tcg.mae.cornell.edu/pubs/Cao_PM_CF_2005.pdf)
- PDF計算では1045 Kの公称値に固執せず、1031–1033 Kで同じリフトオフに合わせて構造を比較しています。ただしこれは自由なチューニングではなく、約31 Kの実験不確かさを条件化した比較です。**単一点合わせではなく H(Tc) の傾きを検証することが本質**です。
- CMCにもCabra H₂/N₂を対象とした実績があります。[1st-order CMC study](https://arxiv.org/abs/2102.04854)

調査メモの「transported PDFかCMCのみ」という表現は少し強すぎます。正確には、**Cabra H₂について、リフトオフと温度応答を説得力ある形で再現した実績が最も強いのがPDF/CMC**です。[調査メモ](/home/sano/work/forge-chem/notes/investigations/cabra-liftoff-model-fidelity-survey.md:14)

納期上TCIを見送る場合は、選択肢2も可能ですが、成果物を「出口状態推定用・探索モデル」と明記し、火炎位置、着火余裕、壁熱負荷、flashback、NOxを適用範囲外にすべきです。

## (3) 優先度付き推奨アクション

### P0 — 設計計算前の必須ゲート

1. **設計QoIと合格基準を確定する**

   出口全温・組成・圧力損失だけか、火炎位置・壁熱流束・安定余裕まで要求するかを明文化する。後者ならTCIは必須です。

2. **SSTプラトーを解消または非定常問題として扱う**

   単なるstep延長ではなく、dual-time/非定常で統計定常化を確認する。`check_quasisteady.py` に少なくとも火炎基部、出口massflux、全エンタルピー流束、主要種流束、Tmax、壁熱流束を追加して判定する。

3. **Jachimowski機構の着火遅れ面を検証する**

   Cabraの最反応混合分率付近、`T_c=1015–1075 K`、1 atm、および実加熱器圧力で Li/Mueller/Keromnes/Burke系と比較する。現行機構は1000 KでGRI系より約2倍速く、BK/Cabra早着火への寄与を排除できていません。[plan](/home/sano/work/forge-chem/plans/active/chemistry-finite-rate-h2.md:102)

4. **TCIを別planで設計する**

   第一候補は1st-order CMC。Cabraについて `T_c±30 K` の応答曲線、統一したOH閾値、平均・RMSの温度／主要種を合格基準にする。「数セッション規模」は楽観的で、独立した機能計画として扱うべきです。

### P1 — 検証の信頼性向上

5. BKのH₂噴射温度、壁温、入口乱流、3D側壁blockageの不確かさを整理する。壁条件をさらに調整する前に、機構感度と着火点の準定常性を確認する。

6. Cabraの入口 `k/ω`、管内長、リップ熱条件、格子感度を確認する。forgeは管内長20 mmですが、文献PDF計算は約50d上流と鋼管熱伝導を含みます。[Gordon et al.](https://web.aeromech.usyd.edu.au/thermofluids/auto_files/Gordon-Masri-Pope-Goldin-Revised.pdf)

7. 実加熱器のend-to-end試行では、質量・元素・全エンタルピー収支、圧力損失、出口一様性、ノズル入口ラジカル、格子／restart依存性を必須出力にする。

### P2

8. `chemistry_source_d` の占有率最適化、`inlet_Pressure_dir` 組成対応。正しさのゲートより後でよいです。

## (4) 指摘事項

- `plans/active` に置き `status: in_progress` としていること自体は正しいです。ただしメタ情報は「BK一次結果」時点で止まり、Cabraが検証仕様・判定基準に追加されていません。[plan冒頭](/home/sano/work/forge-chem/plans/active/chemistry-finite-rate-h2.md:3)
- Phase 2を`done`としながら、計画上の全残差PASS基準をQ1Dが満たしていない点は整理が必要です。
- [methods/index.md](/home/sano/work/forge-chem/methods/index.md:26) が「Phase 0完了、実装未着手」のままです。
- [methods/chemistry.md](/home/sano/work/forge-chem/methods/chemistry.md:152) も、precond版反応熱対角を「未配線」、PaSRを将来扱いとしており現状と不整合です。
- [plans/README.md](/home/sano/work/forge-chem/plans/README.md:26) の概要もPhase 0時点で止まっています。
- run表は存在しますが、複数runを1行にまとめており、AGENTS.mdの「1 run = 1行」に違反しています。case/47には番号なしの`run_dbg_*`が複数、`run_0026`の番号重複もあります。[命名・索引ルール](/home/sano/work/forge-chem/AGENTS.md:50)
- BK出口組成・Cabra火炎基部などの派生量に、対応する `check_quasisteady.py` VERDICTがありません。BKで確認済みなのは`machmax,pmax`だけで、出口組成や着火位置の定常性を保証しません。
- Cabraの数値火炎基部はスクリプトで`Y_OH>1e-4`、調査メモでは標準を約`2e-4`としており、実験OH発光との対応も含め定義が不統一です。
- 計画で要求したcase/28・case/44の`chemistry.enabled: 0`ビット不変検証について、変更ログに明確な完了エビデンスが見当たりません。
- 0-Dテストは手動ビルドで、通常の自動テストには登録されていません。反応機構読込、Jacobian、datum、falloffは回帰テスト化すべきです。
- 外部のClaude artifactを結果の主要参照先にするのは再現性が弱く、repo内PNG・CSV・判定ログを一次情報にすべきです。

ファイル変更は行っていません。
