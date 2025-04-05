#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "yaml-cpp/yaml.h"

#include "cuda_forge/fluct_variables_d.cuh"


using namespace std;

struct bcondConfFormat{
    int physID;
    std::string physName;
    std::string kind;
    int outputHDFflg;
    std::map<std::string, int> inputInts ;
    std::map<std::string, flow_float> inputFloats;
    std::map<std::string, int> neuDirFlag ;

    map<string, map<string,int>> valueTypesOfBC
    {
        { // 0: float variables, 1: uniform float variables (read), 3: int variables , 
          //10:float information, 11:int information
          {"wall" , {
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,1},
              {"Uy"  ,1},
              {"Uz"  ,1},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,0},
              {"Ps"  ,0},
              {"ypls",0},
              {"twall_x",0},
              {"twall_y",0},
              {"twall_z",0},
          }},

          {"wall_isothermal", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,1},
              {"Uy"  ,1},
              {"Uz"  ,1},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,1},
              {"Ps"  ,0},
              {"ypls",0},
              {"twall_x",0},
              {"twall_y",0},
              {"twall_z",0},

          }},

          {"inlet_uniformVelocity", { 
              {"ro"  ,1},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,1},
              {"Uy"  ,1},
              {"Uz"  ,1},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,0},
              {"Ps"  ,1},

          }},

          {"inlet_fluctVelocity", { 
              {"ro"  ,1},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Ux0" ,1},
              {"Uy0" ,1},
              {"Uz0" ,1},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,0},
              {"Ps"  ,1},

          }},


          {"inlet_Pressure", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Tt"  ,1},
              {"Pt"  ,1},
              {"Ts"  ,0},
              {"Ps"  ,0},
          }},

          {"inlet_Pressure_dir", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,1},
              {"Uy"  ,1},
              {"Uz"  ,1},
              {"Tt"  ,1},
              {"Pt"  ,1},
              {"Ts"  ,0},
              {"Ps"  ,0},
          }},

          {"outlet_statPress", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Tt"  ,1},
              {"Pt"  ,1},
              {"Ts"  ,0},
              {"Ps"  ,1},
          }},

          {"outflow", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Tt"  ,1},
              {"Pt"  ,1},
              {"Ts"  ,0},
              {"Ps"  ,0},

          }},
 
          {"slip", { 
              {"ro"  ,-1},
              {"roUx",-1},
              {"roUy",-1},
              {"roUz",-1},
              {"roe" ,-1},
              {"Ux"  ,-1},
              {"Uy"  ,-1},
              {"Uz"  ,-1},
              {"Tt"  ,-1},
              {"Pt"  ,-1},
              {"Ts"  ,-1},
              {"Ps"  ,-1},
          }},

          {"periodic", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,0},
              {"Ps"  ,0},
              {"dtheta", 10},
              {"dx", 10},
              {"dy", 10},
              {"dz", 10},
              {"type" , 11},
              {"partnerBCID"  , 11},
              {"partnerPlnID" , 3},
              {"partnerCellID", 3},
          }},
        }
    };

    bcondConfFormat();
};

//void setBcondsValue(solverConfig& cfg , mesh& msh , variables& var , matrix& mat_p);
void applyBconds(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p , fluct_variables& fluct);

void copyBcondsGradient(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p);

void readBcondConfig(solverConfig& , vector<bcond>& );

void wall_isothermal(solverConfig& , bcond& , mesh& , variables& , matrix& );

void inlet_uniformVelocity(solverConfig& , bcond& , mesh& , variables& , matrix& );

void inlet_Pressure(solverConfig& , bcond& , mesh& , variables& , matrix& );

void outlet_statPress(solverConfig& , bcond& , mesh& , variables& , matrix& );

void outflow(solverConfig& , bcond& , mesh& , variables& , matrix& );

void slip(solverConfig& , bcond& , mesh& , variables& , matrix& );

void periodic(solverConfig& , bcond& , mesh& , variables& , matrix& );
