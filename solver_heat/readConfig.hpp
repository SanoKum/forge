#ifndef READ_CONFIG_H
#define READ_CONFIG_H

#include <vector>
#include <string>
#include <map>

#include "flowFormat.hpp"
#include "yaml-cpp/yaml.h"

#include <iostream>

class solverConfig
{
private:
    std::string solConfigFileName    = "solverConfig.yaml";

public:
    std::string meshFormat; 
    std::string meshFileName; 

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



    solverConfig()
    {
        try {
            YAML::Node config = YAML::LoadFile(this->solConfigFileName);

            std::string meshFormat   = config["mesh"]["meshFormat"].as<std::string>();
            std::string meshFileName = config["mesh"]["meshFileName"].as<std::string>();
            std::string solver = config["solver"].as<std::string>();

            auto last = config["time"]["last"];
                int endTimeControl = last["control"].as<int>();
                int nStep = last["nStep"].as<int>();
                int endTime = last["time"].as<flow_float>();

            auto deltaT = config["time"]["deltaT"];
                int dtControl = deltaT["control"].as<int>();
                flow_float dt = deltaT["dt"].as<flow_float>();
                flow_float cfl = deltaT["cfl"].as<flow_float>();

            int outStepInterval = config["time"]["outStepInterval"].as<int>();

            auto space = config["space"];
                int convMethod = space["convMethod"].as<int>();


            auto physProp = config["physProp"];
                int isCompressible = physProp["isCompressible"].as<int>();
                flow_float ro = physProp["ro"].as<flow_float>();
                flow_float visc = physProp["visc"].as<flow_float>();
                flow_float thermCond = physProp["thermCond"].as<flow_float>();
                flow_float cp = physProp["cp"].as<flow_float>();

            std::cout << "Mesh Name : " << meshFileName << '\n';
            std::cout << "Mesh Name size : " << meshFileName.size() << '\n';

            std::cout << "Solver Name : " << solver << '\n';
            if        (endTimeControl == 0) { std::cout << "End Step : " << nStep << '\n';
            } else if (endTimeControl == 1) { std::cout << "End Time : " << endTime << '\n';
            } else {
                std::cerr << "something wrong in solver yaml" << std::endl;
                exit(EXIT_FAILURE);
            }

            if        (dtControl == 0) { std::cout << "delta Time : " << dt << '\n'; 
            } else if (dtControl == 1) { std::cout << "CFL : " << cfl << '\n';
            } else {
                std::cerr << "something wrong in solver yaml" << std::endl;
                exit(EXIT_FAILURE);
            }

            std::cout << "convection method : " << convMethod<< '\n';

            this->meshFormat   = meshFormat;
            this->meshFileName = meshFileName;
            this->solver = solver;

            this->endTimeControl = endTimeControl;
            this->nStep = nStep;
            this->outStepInterval = outStepInterval;
    
            this->dtControl = dtControl;
            this->dt = dt;
            this->cfl = cfl;
    
            this->convMethod = convMethod;

            this->isCompressible = isCompressible ;
            this->ro = ro;
            this->visc = visc;
            this->thermCond = thermCond;
            this->cp = cp;

        } catch(const YAML::BadFile& e) {
            std::cerr << e.msg << std::endl;
            exit(EXIT_FAILURE);

        } catch(const YAML::ParserException& e) {
            std::cerr << e.msg << std::endl;
            exit(EXIT_FAILURE);
        }
    }
};


#endif