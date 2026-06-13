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
    int blockDPLUR = 0;
    int implicitSolvePrecision = 0; // block-DPLUR 線形 solve の内部精度。0: float (既定・高速), 1: double。
                                    // 残差/状態は float のまま、Jacobian 構築+5×5 solve のみ double 化する混合精度
                                    // (iterative refinement)。軸対称 近軸の float 陰解固着 (Uy が −15 でなく
                                    // −0.6 固着) を根治するが double は遅い (RTX で ~×2.6)。詳細:
                                    // .github/plans/precision-mixed-axisym.md。blockDPLUR==1・lowMachPrecond 0/1 でのみ有効。
    int detectNaN = 0;             // 0: off (既定), 1: 毎ステップ終端で保存量+P の非有限値を検査し、
                                   // 見つけたら res_nan_<step>.h5 をダンプして即停止する診断モード。
                                   // off のときは検査を一切行わないため通常実行はビット不変。
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


    int LESorRANS; // 0:no 1:LES 2:RANS
    int LESmodel; // 1:WALE
    int RANSmodel = 0; // 0:none 1:SST
    int scalarDiffusion = 1; // 0:advection-only 1:advection+diffusion
    int dilatationCorrection = 2; // SST生産項の圧縮性補正 0:off 1:deviatoric(A) 2:deviatoric+isotropic(A+B) 既定:2
    int katoLaunder = 0; // SST生産項 Kato-Launder 補正 0:標準 mu_t S^2 1:mu_t S Omega 既定:0 (docs/turbulence §7.5)

    int isCompressible;
    int isAxisymmetric = 0;
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
    flow_float Sc = 0.7;                       // 定数 Schmidt 数 (speciesDiffusionMethod==0)
    flow_float Sc_t = 0.7;                     // 乱流 Schmidt 数
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

