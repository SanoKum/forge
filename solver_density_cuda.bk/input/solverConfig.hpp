#pragma once

#include "flowFormat.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "yaml-cpp/yaml.h"


class solverConfig
{
private:
    std::string solConfigFileName;

public:
    std::string meshFormat; 
    std::string meshFileName; 

    int gpu;

    std::string solver;

    //Time
    int endTimeControl; // 0: use dt , 1: cfl
    int nStep;
    int outStepInterval;

    int dtControl; // 0: use dt , 1: cfl
    flow_float totalTime=0.0;
    flow_float dt;
    flow_float dt_pseudo;
    flow_float cfl;
    flow_float cfl_pseudo;
    flow_float dt_max;
    flow_float dt_min;

    int timeIntegration; // 1: Euler expli , 3: 3rd Runge Explicit

    int isImplicit; // 0: exp , 1:imp;

    // for inner loop
    int nLoop; 
    std::vector<flow_float> coef_N;
    std::vector<flow_float> coef_M; 
    std::vector<flow_float> coef_Res;

    int convMethod; // 0: 1st Up

    int isCompressible;
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

    flow_float Pref = 101325.0;
    flow_float Tref = 288.15;

    solverConfig();

    void read(std::string);
    void initTimeIntegrationScheme(int timeIntegration);
};

