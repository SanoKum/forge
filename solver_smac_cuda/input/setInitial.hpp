#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "solverConfig.hpp"

#include <cmath>

void setInitial(solverConfig& cfg , mesh& msh , variables& v)
{
    // *** set initial value ***
    for (geom_int i = 0 ; i<msh.nCells ; i++)
    {
        v.c["T"][i] = 300.0;
        v.c["Ux"][i] = 10.0;
        v.c["Uy"][i] = 0.0;
        v.c["Uz"][i] = 0.0;
        v.c["P"][i] = 0.0;
        v.c["ro"][i] = 1.0;
    }
};
