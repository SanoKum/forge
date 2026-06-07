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

    int isCompressible;
    int isAxisymmetric = 0;
    int thermalMethod;
    int viscMethod;

    flow_float ro;
    flow_float visc;
    flow_float thermCond;
    flow_float cp;
    flow_float gamma;
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

