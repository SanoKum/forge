# ノズル設計ツール Phase 3: ①軸対称風洞 (モード F + 帰還エンジン v1/v2 + 照合)

## メタ

- **area**: `tooling / optimization`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (現在仕様)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) — §4.6①/§4.7 (2026-08-13 改訂 = 帰還 3 段構成) と §5 Phase 3 を実装する
  - 前提: [`../accepted/tooling-nozzle-phase0-foundation.md`](../accepted/tooling-nozzle-phase0-foundation.md) (評価パイプライン) / [`../accepted/tooling-nozzle-moo-loop.md`](../accepted/tooling-nozzle-moo-loop.md) (MOO 基盤)
  - 参考: [`../../notes/investigations/nozzle-top-internal-shock-diagnosis.md`](../../notes/investigations/nozzle-top-internal-shock-diagnosis.md) (forge/SU2 Euler が安価で信頼できる対照であることの確証)
- **created**: `2026-08-13`
- **owner**: `sano` (実装: Claude 自走)

## 1. 目的

①軸対称風洞ノズルを**モード F (中心線マッハ全権)** で設計できる状態にする。核心は
**帰還エンジンの初実装** (親計画 §4.7 の 3 段構成のうち v1・v2):

- **v1 (MOC 逆設計)**: starting line + 軸目標 $M_c(x)$ → 逆マーチで超音速場を充填 →
  **質量流量法**で壁流線を抽出 → **δ\* 経験式**の法線オフセット。CFD なし・ミリ秒
- **v2 (Euler 帰還)**: forge Euler (cell) で検算し、軸 $M(x)$ 残差を v1 特性線網の
  **凍結対応マップ**で壁へ帰還。$\|\Delta M\|_\infty \le 0.5\%\,M_d$ までパス反復

検証は「既知目標分布 → 逆設計 → forge 実測 → 目標回復」の自己一貫 + 一様性
$\varepsilon_M,\varepsilon_\theta$ の実装。MOO (co-kriging 多忠実度)・凝縮確認メニューは
本 plan の後半ステップ (7–8) とし、v1/v2+照合の成立を先に確定する。

## 2. スコープ

- **やる**:
  - `probdef` に `wind_tunnel_axisym` 型 (spec: $P_0,T_0,D_e,M_d$ / ($D_e,r^*,M_d$) 過拘束
    チェック / $r^*$ 閉ループ派生 — 初回は無次元 $r^*{=}1$ 設計 → $D_e$ スケーリング)
  - **逆 MOC マーチ** (`geometry/moc_inverse.py`): 軸データ (M 指定, θ=0) からの横マーチ
    (Goursat) による場充填 + 各 C⁺ フロントの $\int\rho u\,dA$ = スロート流量で壁点抽出。
    軸目標定義域は $M_d$ 到達後「一定」延長 (§4.7(c) — 壁端 C⁻ 足の充足は一様域で自明)
  - モード F 目標 Bézier: starting line の軸点に C2 整合、出口 $(M_d, M'{=}M''{=}0)$ C2
  - δ\* 経験式 (`feedback/deltastar.py`): 圧縮性乱流平板相関の局所適用 (Edenfield 系は option)。
    v1 の役割は初期壁の質だけなので相関精度は帰還収束速度にのみ影響 (B3)
  - **v2 Euler 帰還** (`feedback/euler_loop.py`): forge Euler (cell — node slip 市松の既知
    問題を回避)・R1 同トポロジ再メッシュ・同 index warm restart・v1 特性線網の凍結
    軸↔壁マップ・PM 換算 $\Delta\theta_w=\omega\Delta\nu$・重みランプ・クリップ (§4.7(c))
  - 一様性メトリクス: $\varepsilon_M$ (コア質量流束重み RMS)・$\varepsilon_\theta$ (コア max 流れ角)
  - **照合**: (i) 自己一貫 (既知目標→回復)、(ii) 放射源流・一様流の解析照合 (単体)、
    (iii) 可能なら CONTUR/Sivells 系コンタ (papers/ の AEDC-TR-78-63) との比較
- **後半 (ステップ 7–8, v1/v2 成立後)**:
  - モード F dv ($M_c$ 自由 CP k=3 + $L$ + $R$) の MOO、Euler=低忠実度の co-kriging
  - 凝縮確認メニュー (複数入口条件 → 出口 M・動圧低下表)
- **やらない**:
  - v3 (NS 凍結場フルトレース) — v2 停滞時のみ (§4.7)
  - Kliegel–Levine 高次 — ① は $R_u{=}R_d$ 単一円弧で $R\gtrsim1$ が標準、Sauer で足りる
    見込み (着手時に実値で要否判定、必要なら本 plan に追記)
  - ② 矩形 (側壁 δ\*)・q_peak ランキング検証 (Phase 3 後半〜Phase 4)

## 3. 関連 docs と前提

- 帰還アルゴリズムの正本 = 親計画 §4.7 (a)–(c)。逆マーチの数理は本セッションの議論で
  確定済み: 壁点の依存域は軸の両側 (C⁺ 足上流・C⁻ 足下流)、一様出口では下流延長が自明
- `moc_kernel.py` の実装知見 (楔充填・軸対称 1/r 源項の評価点・前フロント補間) を流用
- forge Euler = `viscMethod: 0, visc: 0.0, thermCond: 0.0` + slip 壁 + cell (夜間バッチ
  run_0085 で検証済みの経路)。SU2 Euler 対照も安価 (rms 10⁻¹¹, ~3 分)

## 4. 設計方針 (差分のみ)

- **スロート近傍の幾何**: ① は $R_u=R_d=R$ 単一円弧 (dv, 初期 1.5–2)。starting line は
  既存 Sauer。壁の所有権: 円弧は「遷音速アンカー」(モード F — 円弧終端 = starting line
  が壁に当たる点、θa という dv は無い)
- **逆マーチの検証規律** (moc_kernel と同格): 一様流 (M 一定 → 円筒壁)・放射源流
  (解析コーン回復) の 2 解析照合 + 格子収束。決定性 (同一入力 → bit 同一壁)
- **v2 の収束判定は実場**: 達成軸 M を固定相対格子に補間し $\|\Delta M\|_\infty$。マスク:
  $M<1.05$ 域・出口近傍ランプ (§4.7(c))
- case は `case/41.wind_tunnel_design/` を新設 (README run 台帳)

## 5. 実装ステップ

1. `geometry/moc_inverse.py`: 逆マーチ + 質量流量法壁抽出 + 単体検証 (一様流/源流/格子収束)
2. モード F 目標 Bézier 構成 (starting line C2 整合) + `probdef` `wind_tunnel_axisym` 型
3. `feedback/deltastar.py` (δ\* 相関) + v1 チェーン結合 (dv → 物理壁 → TFI メッシュ)
4. **v1 E2E**: $M_d=4$ 標準条件で逆設計壁を forge Euler 実測 → 達成軸 M vs 目標 (期待
   ~1–3% 残差) + $\varepsilon_M,\varepsilon_\theta$ 実装
5. `feedback/euler_loop.py` (v2): 凍結マップ帰還のパス反復 → $\|\Delta M\|_\infty \le 0.5\%M_d$
6. 照合: 自己一貫 (別の既知分布から回復) + Sivells 系比較 (可能なら)
7. (後半) モード F MOO (co-kriging: Euler=低 / NS=高)
8. (後半) 凝縮確認メニュー
- 各段で `methods/design/overview.md` を同期

## 6. 検証

- 単体: 逆マーチの解析照合 2 系統 + 格子収束 + 決定性 (§4)
- v1: 逆設計壁の forge Euler 達成度 (定量記録)。メッシュ品質 PASS
- v2: 帰還収束履歴 ($\|\Delta M\|_\infty$ vs パス数) がロバストに単調減、~10 パス内で
  0.5% 到達。R1 restart の過渡長実測 (親計画 §6 Phase 1 の宿題)
- 一様性: 収束解で $\varepsilon_M, \varepsilon_\theta$ を報告 (NS でも 1 本測る)
- 台帳: `case/41` README + 3 ツール VERDICT

## 7. 影響範囲

- 新規: `design/forge_design/geometry/moc_inverse.py`, `design/forge_design/feedback/`,
  `case/41.wind_tunnel_design/`
- 変更: `probdef.py` (問題型追加)・`evaluate/runner.py` (Euler 実行モード・wind tunnel 対応)
- 既存ソルバ無改造

## 8. 完了条件

- [ ] §5 の 1–6 (v1/v2 + 照合)
- [ ] §5 の 7–8 (MOO + 凝縮メニュー) — 先行完了条件を分けて中間コミット可
- [ ] `methods/design/overview.md` 同期
- [ ] status done → `accepted/`、`plans/README.md` 同期

## 9. 変更ログ

- `2026-08-13` (夜間自走) — **§5 ステップ 1–5 完了 (v1+v2 成立)**:
  - **逆 MOC マーチ** (`moc_inverse.py`): 軸 Cauchy 三角充填 + 壁流線抽出。単体 ALL PASS
    (一様流=機械精度 / 源流レイ 0.13% / Md4 ノズル: リップ完成・$r_e/\sqrt\varepsilon$=0.996
    [=√Cd]・質量保存 0.2%・決定性)。軸目標の下流延長則を定量確認 (リップ回復には
    C⁻ 足まで = 素朴 margin の約 2 倍)。
  - **v1 チェーン** (`wall_modef`/`deltastar`/`runner_wt`, `case/41` run_0001): Md=4 設計壁の
    forge Euler 実測 = 出口軸 M 4.024 (+0.6%)・ΔM_max 3.0%Md。誤差構造 = 接合波の軸着地
    スパイク (円弧/設計壁 C1 接合の曲率跳び — マスク領域) + 減速部の滑らかな欠損 (v2 対象)。
    Sauer 線形化の壁 θ 過小は kernel と同じ接線上書きで解消。
  - **v2 Euler 帰還** (`euler_loop.py`, run_0003_v2_tr): 凍結 C⁻ マップ + PM 換算 + 平滑化/
    クリップ + 同一トポロジ再メッシュ + warm restart。固定 ω=0.4 は高 M 末端でリンギング
    (run_0002 — dν/dM 縮小による実効過ゲイン) → **決定的 trust-region (悪化 reject → 最良壁
    巻き戻し + ω 半減) で単調降下 2.89%→**0.45% Md** (13 パス, reject 1回) — §4.7(c) の
    0.5% 基準達成**。図 = `run_0003_v2_tr/v2_convergence.png`。
  - 残り: §5 ステップ 4 の $\varepsilon_M/\varepsilon_\theta$ 実装、ステップ 6 (独立目標での
    自己一貫 + Sivells 照合)、ステップ 7–8 (MOO / 凝縮メニュー)。
- `2026-08-13` — 起票 (親計画 §5 Phase 3。Phase 1 の帰還エンジンを統合)。
