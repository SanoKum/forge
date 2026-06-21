#pragma once

#include "flowFormat.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "yaml-cpp/yaml.h"
#include <stdexcept>


class solverConfig
{
private:
    std::string solConfigFileName;

public:
    std::string meshFormat; 
    std::string meshFileName; 
    std::string valueFileName; 

    int gpu;

    std::string solver;

    //Time
    int endTimeControl; // 0: use dt , 1: cfl
    int nStepOuter;
    int outStepInterval;
    int outStepStart;

    int dtControl; // 0: use dt , 1: cfl
    flow_float totalTime=0.0;
    flow_float dt;
    flow_float dt_pseudo;
    flow_float cfl;
    flow_float cfl_pseudo;
    flow_float implicitRelax = 1.0;
    flow_float implicitRelaxSST = -1.0; // -1: implicitRelax に倒置 (既定動作不変)
    // 軸対称 near-axis 安定化: 擬似時間スペクトル半径に軸項 λ_axis=β·(|u_r|+c)·A_planar を加える。
    // 近軸 (r→0) で Δτ∝CFL·r/(|u_r|+c) を自然に与え半径運動量不安定を抑える。0=不変 (既定)。
    flow_float axisTimestepBeta = 0.0;
    int blockDPLUR = 0;
    // 多成分 face 整合再構成: 0 (既定・ビット不変, 組成は owner セル 1 次=mixed-order)、
    // 1: Y_s を ρ と同じ勾配+limiter で face へ 2 次再構成し thermo/species 流束で同一 face 組成を使う。
    int speciesFaceReconstruction = 0;
    // multispeciesRhoYCommonLimiter: opt-in 診断。0 (既定・ビット不変)、
    // 1: ρ と全 species に共通リミタ ψ_ρY=min(ψ_ρ, min_s ψ_Y_s) を適用し ρ_f=ρ(Y_f) 整合だけを切り分ける
    //    (p・速度は各自のリミタのまま)。nSpecies>1 かつ speciesFaceReconstruction>=1 で有効。
    int multispeciesRhoYCommonLimiter = 0;
    int speciesImplicitCoupling = 0; // 多成分 TP 陰解法 (timeIntegration==11, nSpecies>=2) の化学種更新方式。
                                     // 0: 従来 segregated 点陰的 forward-Euler (既定・ビット不変)。
                                     // 1: 緩和整合 scalar-DPLUR (流れ block と同一 dt_local/implicitRelax/nStepInner
                                     //    sweep で ρY_s を緩和し、要因2 の擬似時間緩和ミスマッチを解消)。
                                     // 詳細: .github/plans/thermophysics-species-implicit-coupling.md。
    int implicitSolvePrecision = 0; // block-DPLUR 線形 solve の内部精度。0: float (既定・高速), 1: double。
                                    // 残差/状態は float のまま、Jacobian 構築+5×5 solve のみ double 化する混合精度
                                    // (iterative refinement)。軸対称 近軸の float 陰解固着 (Uy が −15 でなく
                                    // −0.6 固着) を根治するが double は遅い (RTX で ~×2.6)。詳細:
                                    // .github/plans/precision-mixed-axisym.md。blockDPLUR==1・lowMachPrecond 0/1 でのみ有効。
    int detectNaN = 0;             // 0: off (既定), 1: 保存量+P の非有限値を検査し、
                                   // 見つけたら res_nan_<step>.h5 をダンプして即停止する診断モード。
                                   // off のときは検査を一切行わないため通常実行はビット不変。
                                   // 検査は fused 1 カーネルで device フラグへ集約し、detectNaNInterval
                                   // ステップごとにのみ host へ読み出す (per-step 同期を避ける)。
    int detectNaNInterval = 1;     // detectNaN==1 時にフラグを host 読み出しする間隔 [step] (既定 1=毎ステップ)。
                                   // 大きくすると per-step 同期が減るが NaN 検知が最大 (interval-1) step 遅れる。
    // 残差 RMS を device バッファに常駐させ、この間隔 [step] ごとに 1 回だけ D2H 転送して
    // residual_history.csv へ一括書き出しする。既定 1=毎ステップ flush (従来挙動)。>1 で per-step の
    // 残差 reduction 同期を間引く (毎ステップの値・行構成は不変)。詳細:
    // .github/plans/architecture-residual-monitor-async.md。
    int residualFlushInterval = 1;
    // 定常 (unsteady=0) implicit 経路で `max cfl` の host 読み出し (`thrust::max_element`)・dt 適応・
    // 表示を行う間隔 [step] (既定 1=毎ステップ=従来挙動)。定常では dt_local=cfl_pseudo·dx/λ で cfg.dt が
    // 打ち消され、この host 読み出しは診断のみ (解に不影響) なので >1 で per-step 同期を間引ける。
    // explicit / dual-time 経路には適用しない (cfg.dt が物理的に効くため)。詳細:
    // .github/plans/architecture-perphase-profiling-hotspot.md。
    int cflReportInterval = 1;
    int lowMachPrecond = 0;        // 0: off (従来), 1: Weiss-Smith 低マッハ前処理 (フラックス散逸)
    flow_float precondEps = 0.15;  // 低マッハ前処理の停留点フロア ε (Ur=min(c,max(|u|,ε·c)))。
                                   // ε 小ほど低マッハ振動を強く減衰するが ε≲0.1 は発散 (ε=0.05 で NaN)。
                                   // ε=0.15: M4 ノズルで limit-cycle 振幅 −32% (検証済), ε=0.3: −17%。
    int lowMachThornber = 0;       // 0: off (従来), 1: Thornber 型再構成補正 (SLAU の L/R 速度ジャンプを
                                   // z=min(M,1) で縮約)。lowMachPrecond と直交・併用可。SLAU 経路のみ。
    flow_float dt_max;
    flow_float dt_min;

    int unsteady; // steady , unsteady
    int dualTime; // 0: off , 1: on
    int timeIntegration; // 1: Euler explicit, 3: 3rd Runge explicit, 4: 4th RK, 11: implicit defect-correction

    int isImplicit; // 0: exp , 1:imp;

    // for inner loop
    int nStage;
    int nStepInner;
    int nSubIterDualTime = 20; // dual-time: 物理ステップあたりの擬似時間サブ反復数
    int bdfOrder = 2;          // dual-time: 物理時間 BDF 次数 (1 or 2、初回ステップは BDF1)
    flow_float unsteadyDiagCoef = 0.0; // dual-time: 陰解法対角へ加える物理時間項係数 a/Δt（定常は 0）。driver が毎ステップ設定
    std::vector<flow_float> coef_N;
    std::vector<flow_float> coef_M; 
    std::vector<flow_float> coef_Res;

    std::vector<flow_float> coef_DT_4thRunge;
    std::vector<flow_float> coef_Res_4thRunge;

    flow_float  convMethod;
    int limiter;    // 0: off, 1: Barth-Jespersen, 2: Venkata, -1: legacy

    // free-stream 保存: 対流流束の圧力項を (p_tilde - pRef)*s で組み、非直交メッシュで
    // 大きな p*s を float32 加算する際の桁落ち(metric closure 由来の偽運動量源)を抑える。
    // 既定 0.0 で従来挙動(ビット不変)。一様基準(動作/フリーストリーム)静圧を入れる。
    flow_float  pRef = 0.0;


    int LESorRANS; // 0:no 1:LES 2:RANS
    int LESmodel; // 1:WALE
    int RANSmodel = 0; // 0:none 1:SST
    int scalarDiffusion = 1; // 0:advection-only 1:advection+diffusion
    int dilatationCorrection = 2; // SST生産項の圧縮性補正 0:off 1:deviatoric(A) 2:deviatoric+isotropic(A+B) 既定:2
    int katoLaunder = 0; // SST生産項 Kato-Launder 補正 0:標準 mu_t S^2 1:mu_t S Omega 既定:0 (docs/turbulence §7.5)
    int wallTreatmentSST = 0; // SST壁処理 0:low-Re壁解像(60ν/β₁y²) 1:automatic(y⁺非依存) 既定:0 (docs/turbulence §6.5)
    flow_float turbulentPrandtl = 0.85; // 乱流プラントル数 Pr_t。乱流熱伝導 k_t=cp*mu_t/Pr_t に使用 (既定 0.85)

    int isCompressible;
    int isAxisymmetric = 0;

    // 離散化レイアウト。"cell": cell-centered FVM (既定・従来)、"node": node-centered
    // (中点双対 median-dual) FVM。docs/discretization/ 参照。
    std::string discretization = "cell";

    // 境界隣接 CV の 2 次再構成を 1 次に落とす (0:off 既定, 1:on)。node-centered の壁近傍
    // 高マッハ発散の原因切り分け診断 / ロバス化。docs/discretization/implementation.md §7。
    int bndFirstOrder = 0;

    // node-centered 軸対称: 軸上 CV (R=0) の cell 中心に CV 面積加重重心を使うか (1:既定, ゼロ回転体積回避)。
    // 0 で node 座標 (R=0) を使う。0 は軸ソース OFF と併用前提 (converter で消費)。
    int axisCentroidShift = 1;

    // node-centered 壁 Dirichlet 試作 (0:既定 OFF)。1 で init u=0 + 運動量残差射影。
    // 注: 残差射影は壁圧力寄与も落とすため不正 (検証で壁全域発散)。正解は弱形式半割面 flux (Phase 2,
    // docs/discretization/implementation.md §7.2)。本フラグは切り分け用に残置、既定 OFF。
    int nodeWallDirichlet = 0;

    // 勾配を最小二乗 (LSQ) で計算する (0:既定 Green-Gauss, 1:LSQ)。node-centered の median-dual 近壁で
    // Green-Gauss 面勾配が checkerboard を持ち粘性 BL を振動させるための対策 (over-relaxed 法線項は別途維持)。
    // 境界は bvar 境界値を LSQ 点として含め勾配を閉じる。docs/discretization/implementation.md §7.3。
    int gradLSQ = 0;

    // node-centered の内部双対面で面補間係数を中点 fx=0.5 (φ_f=½(φ_A+φ_B)) に固定する
    // (0:既定 幾何 fx=dualFaceCent 射影比, 1:中点)。標準的な median-dual エッジ補間で近壁 checkerboard を
    // 低減するが、出口リップ近傍 near-wall の解を変えるため SU2 検証まで既定 OFF。cell モードは無関係。
    // docs/discretization/implementation.md §7.4。
    int nodeMidpointFx = 0;

    // node-centered 粘性壁 flux の法線距離 dcc 退化対策 (0:従来ミラーゴースト距離そのまま,
    // 1:既定 退化 dcc をフロアで平滑化)。node モードの壁ノードは壁面上に乗るため、ミラーゴースト
    // cc_ghost=cc+2((pc-cc)·n)n は (pc-cc)·n が雑音 (|dn|∈[1e-8,3e-4], 符号反転) で dcc=2|dn| が
    // 4〜5 桁ばらつき、法線項 (U[ig]-U[ic])/dcc が近壁 twall の特異スパイク・偶奇振動を生む。
    // 1 で dcc_eff=max(dcc, nodeWallDistFloorCoef·vol/sss) と下限クリップする。
    // 【検証所見 2026-06-20】強 no-slip ドラッグ (小 dcc) は安定性に必須だが同時に twall 振動の源でも
    // あり、局所的な距離操作 (純勾配 ∇u·S / 滑らか距離 / floor 大 / 残差射影 Dirichlet) では両者を
    // 分離できない (弱めると no-slip を失い energy が先頭発散; 強いままだと振動が残る)。バルク解
    // (中心線 Mach) は元スキームで Euler/cell と一致し健全。よって既定は 0 (元挙動を保持) とし、
    // 本フラグ+floor は切り分け/将来の正攻法 (内部隣接ノード距離を使うコンパクト法線作用素) 用に
    // 残置する。docs/diffusion/implementation.md, .github/plans/diffusion-node-wall-viscous-distance.md。
    int nodeWallViscGradFlux = 0;

    // 上記フロア係数 c: dcc_eff=max(dcc, c·vol/sss)。小さいほど元挙動に近い (大半の dcc を温存)。
    // 既定 0.05。0 で実質フロア無効 (nodeWallViscGradFlux=0 と同等)。docs/diffusion/implementation.md。
    flow_float nodeWallDistFloorCoef = 0.05;
    int thermalMethod;   // 0: calorically perfect (定数 cp/γ), 2: 多成分 thermally-perfect (NASA-9)
    int viscMethod;      // 0: 定数, 1: Sutherland, 2: kinetic theory (Chapman-Enskog)

    flow_float ro;
    flow_float visc;
    flow_float thermCond;
    flow_float cp;
    flow_float gamma;

    // 多成分 thermally-perfect gas (thermalMethod==2)。calorically-perfect 経路では未使用。
    int nSpecies = 1;                          // 化学種数 (既定 1 = 単成分)
    std::vector<std::string> speciesNames;     // 混合を構成する化学種名。順序が index s を定義
    std::string speciesDBFile = "";            // 任意: NASA-9/LJ 係数の外部 DB (yaml)。空なら内蔵 DB
    int speciesDiffusionMethod = 1;            // 0: 定数 Schmidt, 1: kinetic theory 混合平均拡散
    double thermoHrefTemp = 0.0;               // >0: 各化学種のエンタルピー基準を h_s(thermoHrefTemp)=0 へ
                                               // オフセット (sensible-enthalpy datum)。非反応流では物理不変
                                               // だが、種ごとに桁違いの生成エンタルピー (H2O≈-13.4MJ/kg) を
                                               // 除いて多成分 implicit の roe/roY 緩和ミスマッチ起因の T ジャンプ
                                               // を抑え安定化する。0 (既定) で従来 NASA 絶対基準・ビット不変。
    flow_float Sc = 0.7;                       // 定数 Schmidt 数 (speciesDiffusionMethod==0)
    flow_float Sc_t = 0.7;                     // 乱流 Schmidt 数 (D_t=mu_t/(ro*Sc_t))。
                                               // turbulence.turbulentSchmidt でも設定可 (physProp.Sc_t は後方互換、turbulence 優先)

    // 非平衡凝縮 (4 モーメント方程式 ρg,ρQ2,ρQ1,ρQ0)。docs/condensation/ 参照。
    // Phase 1 はモーメントを受動スカラー (ソース=0) として輸送するのみ。既定 off で従来経路ビット不変。
    int condensation = 0;    // 0: off (既定), 1: on
    int nCondSpecies = 0;    // 凝縮種数。condensation==1 のとき >=1。当面 1 (N2)
    int condModel = 0;       // 凝縮種の物性/核生成/成長モデル。0: N2 (CNT_Iland+Goodheart, CPG),
                             //   1: H2O (Murphy-Koop, CNT+Kantrowitz+Hertz-Knudsen, carrier+TP)
    int condGasSpecies = -1; // carrier+condensible: 凝縮する気相化学種の index (roY{s})。
                             //   -1: pure-condensible (気相=凝縮種, N2 Arthur)。>=0: H2O 等の希薄凝縮 (Wyslouzil)
    int condKantrowitz = 0;  // 核生成の Kantrowitz 非等温補正。0: off (等温 CNT, 既定), 1: on
    int condGrowthModel = 0; // 成長則。0: 既定 (H2O=Hertz-Knudsen, N2=Goodheart), 1: Gyarmathy(熱伝導律速)
    double condGyarmathyC = 3.18; // Gyarmathy の Knudsen 補正係数 1/(1+C·Kn)。標準 3.18 (感度評価用に可変)
    int condTwoTemp = 0;     // 液滴温度 T_d 考慮 (Hertz-Knudsen 経路, 準定常 Hill バランス)。0: 一温度(既定), 1: 二温度
      //          int isCompressible = physProp["isCompressible"].as<int>();
      //          if (isCompressible == 0) flow_float ro = physProp["isCompressible"]["ro"].as<flow_float>();
      //          flow_float visc = physProp["visc"].as<flow_float>();
      //          flow_float thermCond = physProp["thermCond"].as<flow_float>();
      //          flow_float cp = physProp["cp"].as<flow_float>();

    std::string initial;

    flow_float Pref = 1.0/1.4;
    flow_float Tref = 1.0;

    solverConfig();

    void read(std::string);
    void initTimeIntegrationScheme(int timeIntegration);
    int mainLoopCount() const;
    int perStepIterationCount() const;
    const char* perStepIterationLabel() const;
};

