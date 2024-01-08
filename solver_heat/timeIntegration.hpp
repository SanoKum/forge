#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include <mesh.hpp>
#include <variables.hpp>

void timeIntegration(solverConfig &cfg , mesh &msh , variables &v)
{
    geom_float dt = 0.00001;

    // normal plane
    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        geom_float vol = msh.cells[ic].volume;

        v.c["T"][ic] = v.c["T"][ic] + v.c["diffT"][ic]*cfg.dt/vol;

    }
}

#endif