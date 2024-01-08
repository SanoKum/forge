#ifndef BOUNDARY_COND_H
#define BOUNDARY_COND_H

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

using namespace std;

struct bcondConfFormat{
    int physID;
    std::string kind;
    std::map<std::string, int> inputInts ;
    std::map<std::string, flow_float> inputFloats;
    std::map<std::string, int> neuDirFlag ;

    map<string, map<string,int>> valueTypesOfBC
    {
        { // 0: no gradient, 1: dirichlet, -1:slip, -2:special
          {"wall" , {
              {"Ux", 1},
              {"Uy", 1},
              {"Uz", 1},
              {"P" , 0},
              {"T" , 0},
              {"ro", 0}, }},

          {"wall_isothermal", { 
              {"Ux", 1},
              {"Uy", 1},
              {"Uz", 1},
              {"P" , 0},
              {"T" , 1},
              {"ro", 0}, }},

          {"inlet_uniformVelocity", { 
              {"Ux", 1},
              {"Uy", 1},
              {"Uz", 1},
              {"P" , 0},
              {"T" , 1},
              {"ro", 1}, }},

          {"outlet_statPress", { 
              {"Ux", 0},
              {"Uy", 0},
              {"Uz", 0},
              {"P" , 1},
              {"T" , 0},
              {"ro", 0}, }},

          {"slip", { 
              {"Ux", -1},
              {"Uy", -1},
              {"Uz", -1},
              {"P" , -1},
              {"T" , -1},
              {"ro", -1}, }},
        }
    };

    bcondConfFormat();
};

//void setBcondsValue(solverConfig& cfg , mesh& msh , variables& var , matrix& mat_p);
void setBcondsValue(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p);

void readAndSetBcondConfig(vector<bcond>& );

void wall_isothermal(solverConfig& , bcond& , mesh& , variables& , matrix& );

void inlet_uniformVelocity(solverConfig& , bcond& , mesh& , variables& , matrix& );

void outlet_statPress(solverConfig& , bcond& , mesh& , variables& , matrix& );

void slip(solverConfig& , bcond& , mesh& , variables& , matrix& );

#endif

