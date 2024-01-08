#pragma once

#include <string>
#include "flowFormat.hpp"

class solverConfig
{
private:
    std::string solConfigFileName;

public:
    std::string meshFormat; 
    std::string meshFileName; 

    int gpu;

    std::string solver;

    int endTimeControl; // 0: use dt , 1: cfl
    int nStep;
    int outStepInterval;

    int dtControl; // 0: use dt , 1: cfl
    flow_float dt;
    flow_float cfl;

    int convMethod; // 0: 1st Up

    int isCompressible;
    flow_float ro;
    flow_float visc;
    flow_float thermCond;
    flow_float cp;
      //          int isCompressible = physProp["isCompressible"].as<int>();
      //          if (isCompressible == 0) flow_float ro = physProp["isCompressible"]["ro"].as<flow_float>();
      //          flow_float visc = physProp["visc"].as<flow_float>();
      //          flow_float thermCond = physProp["thermCond"].as<flow_float>();
      //          flow_float cp = physProp["cp"].as<flow_float>();

    solverConfig();

    void read(std::string);
};

