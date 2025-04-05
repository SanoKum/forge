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
    int nStep;
    int outStepInterval;
    int outStepStart;

    int dtControl; // 0: use dt , 1: cfl
    flow_float totalTime=0.0;
    flow_float dt;
    flow_float dt_pseudo;
    flow_float cfl;
    flow_float cfl_pseudo;
    flow_float dt_max;
    flow_float dt_min;

    int unsteady; // steady , unsteady
    int dualTime; // 0: off , 1: on
    int timeIntegration; // 1: Euler expli , 3: 3rd Runge Explicit

    int isImplicit; // 0: exp , 1:imp;

    // for inner loop
    int nStage;
    int nInnerLoop; 
    std::vector<flow_float> coef_N;
    std::vector<flow_float> coef_M; 
    std::vector<flow_float> coef_Res;

    std::vector<flow_float> coef_DT_4thRunge;
    std::vector<flow_float> coef_Res_4thRunge;

    flow_float  convMethod; 
    int limiter;    // 0: off 1:venkata


    int LESorRANS; // 0:no 1:LES 2:RANS
    int LESmodel; // 1:WALE

    int isCompressible;
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
};

