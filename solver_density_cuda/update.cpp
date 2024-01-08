#pragma once

#include "update.hpp"
#include <vector>

using namespace std;

void updateVariablesOuter(solverConfig& cfg , mesh& msh , variables& v , matrix& mat_p)
//void updateVariablesForNextLoop(solverConfig& cfg , mesh& msh , variables& v)
{
    vector<flow_float>& ro   = v.c["ro"];
    vector<flow_float>& roUx = v.c["roUx"];
    vector<flow_float>& roUy = v.c["roUy"];
    vector<flow_float>& roUz = v.c["roUz"];
    vector<flow_float>& roe  = v.c["roe"];

    vector<flow_float>& roN  = v.c["roN"];
    vector<flow_float>& roUxN= v.c["roUxN"];
    vector<flow_float>& roUyN= v.c["roUyN"];
    vector<flow_float>& roUzN= v.c["roUzN"];
    vector<flow_float>& roeN = v.c["roeN"];

    vector<flow_float>& roM  = v.c["roM"];
    vector<flow_float>& roUxM= v.c["roUxM"];
    vector<flow_float>& roUyM= v.c["roUyM"];
    vector<flow_float>& roUzM= v.c["roUzM"];
    vector<flow_float>& roeM = v.c["roeM"];

    //vector<flow_float>& convx = v.c["convx"];
    //vector<flow_float>& convy = v.c["convy"];
    //vector<flow_float>& convz = v.c["convz"];

    //vector<flow_float>& diffx = v.c["diffx"];
    //vector<flow_float>& diffy = v.c["diffy"];
    //vector<flow_float>& diffz = v.c["diffz"];

    vector<flow_float>& divU = v.c["divU"];

    for (geom_int ic=0 ; ic<msh.nCells_all; ic++)
    {

        roN[ic] = ro[ic];
        roUxN[ic] = roUx[ic];
        roUyN[ic] = roUy[ic];
        roUzN[ic] = roUz[ic];
        roeN[ic] = roe[ic];

        roM[ic] = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic] = roe[ic];

        divU[ic] = 0.0;
    }

}

void updateVariablesInner(solverConfig& cfg , mesh& msh , variables& v , matrix& mat_p)
//void updateVariablesForNextLoop(solverConfig& cfg , mesh& msh , variables& v)
{
    vector<flow_float>& ro   = v.c["ro"];
    vector<flow_float>& roUx = v.c["roUx"];
    vector<flow_float>& roUy = v.c["roUy"];
    vector<flow_float>& roUz = v.c["roUz"];
    vector<flow_float>& roe  = v.c["roe"];

    vector<flow_float>& roM  = v.c["roM"];
    vector<flow_float>& roUxM= v.c["roUxM"];
    vector<flow_float>& roUyM= v.c["roUyM"];
    vector<flow_float>& roUzM= v.c["roUzM"];
    vector<flow_float>& roeM = v.c["roeM"];

    vector<flow_float>& divU = v.c["divU"];

    for (geom_int ic=0 ; ic<msh.nCells_all; ic++)
    {

        roM[ic] = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic] = roe[ic];

        divU[ic] = 0.0;
    }

}
