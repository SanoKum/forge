#ifndef UPDATE_H
#define UPDATE_H

#include <mesh.hpp>
#include <variables.hpp>

void updateVariablesForNextLoop(solverConfig &cfg , mesh &msh , variables &v , matrix& mat_p)
{
    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        v.c["TN"][ic] = v.c["T"][ic];

        v.c["diffT"][ic] = 0.0;
    }

    geom_int ic0, ic1;
    geom_int ip_loc0, ip_loc1;

    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        ic0 = msh.planes[ip].iCells[0];
        ic1 = msh.planes[ip].iCells[1];

        ip_loc0  = mat_p.localPlnOfCell[ip][0];
        ip_loc1  = mat_p.localPlnOfCell[ip][1];

        mat_p.lhs[ic0][0] = 0.0;
        mat_p.lhs[ic0][ip_loc0] = 0.0;
        mat_p.rhs[ic0] = 0.0;

        mat_p.lhs[ic1][0] = 0.0;
        mat_p.lhs[ic1][ip_loc1] = 0.0;
        mat_p.rhs[ic1] = 0.0;
    }
}

#endif