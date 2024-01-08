#ifndef GRADIENT_H
#define GRADIENT_H

#include <mesh.hpp>
#include <variables.hpp>

#include <cmath>

void calcGradient(solverConfig &cfg , mesh &msh , variables &v)
{
    //geom_float ss;
    vector<geom_float> sv(3);

    vector<geom_float> dccv(3);
    geom_float dcc;

    vector<geom_float> dc1pv(3);
    geom_float dc1p;

    vector<geom_float> dc2pv(3);
    geom_float dc2p;

    vector<geom_float> pcent(3);
    vector<geom_float> c1cent(3);
    vector<geom_float> c2cent(3);
    geom_int ic1;
    geom_int ic2;

    geom_float f;

    geom_float Tf;

    for (geom_int ie=0 ; ie<msh.nCells; ie++)
    {
        v.c["dTdx"][ie] =  0.0;
        v.c["dTdy"][ie] =  0.0;
        v.c["dTdz"][ie] =  0.0;
    }

    // normal plane
    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        ic1     = msh.planes[ip].iCells[0];
        ic2     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;
        //ss      = msh.planes[ip].surfArea;
        pcent   = msh.planes[ip].centCoords;

        c1cent  = msh.cells[ic1].centCoords;
        c2cent  = msh.cells[ic2].centCoords;

        dccv[0] = c2cent[0] - c1cent[0];
        dccv[1] = c2cent[1] - c1cent[1];
        dccv[2] = c2cent[2] - c1cent[2];
        dcc     = sqrt( pow(dccv[0], 2.0) + pow(dccv[1], 2.0) + pow(dccv[2], 2.0));

        dc1pv[0] = pcent[0] - c1cent[0];
        dc1pv[1] = pcent[1] - c1cent[1];
        dc1pv[2] = pcent[2] - c1cent[2];
        dc1p    = sqrt( pow(dc1pv[0], 2.0) + pow(dc1pv[1], 2.0) + pow(dc1pv[2], 2.0));

        dc2pv[0] = pcent[0] - c2cent[0];
        dc2pv[1] = pcent[1] - c2cent[1];
        dc2pv[2] = pcent[2] - c2cent[2];
        dc2p    = sqrt( pow(dc2pv[0], 2.0) + pow(dc2pv[1], 2.0) + pow(dc2pv[2], 2.0));

        f = dc2p/dcc;

        Tf = f*v.c["T"][ic1] + (1.0-f)*v.c["T"][ic2] ;

        v.c["dTdx"][ic1] +=  sv[0]*Tf;
        v.c["dTdy"][ic1] +=  sv[1]*Tf;
        v.c["dTdz"][ic1] +=  sv[2]*Tf;

        v.c["dTdx"][ic2] += -sv[0]*Tf;
        v.c["dTdy"][ic2] += -sv[1]*Tf;
        v.c["dTdz"][ic2] += -sv[2]*Tf;
    }

    // boundary plane
    for (geom_int ip=msh.nNormalPlanes ; ip<msh.nPlanes ; ip++)
    {
        ic1     = msh.planes[ip].iCells[0];
        sv      = msh.planes[ip].surfVect;
    
        v.c["dTdx"][ic1] +=  sv[0]*v.p["T"][ip];
        v.c["dTdy"][ic1] +=  sv[1]*v.p["T"][ip];
        v.c["dTdz"][ic1] +=  sv[2]*v.p["T"][ip];
    }

    geom_float volume;
    for (geom_int ie=0 ; ie<msh.nCells; ie++)
    {
        volume = msh.cells[ie].volume;
        v.c["dTdx"][ie] = v.c["dTdx"][ie]/volume;
        v.c["dTdy"][ie] = v.c["dTdy"][ie]/volume;
        v.c["dTdz"][ie] = v.c["dTdz"][ie]/volume;
    }

}

#endif