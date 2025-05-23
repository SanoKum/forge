#include "setStructualVariables.hpp"
#include "flowFormat.hpp"

#include <cmath>

void setStructualVariables(solverConfig& cfg , cudaConfig& cuda , mesh& msh , variables& v)
{
    //geom_float ss;
    std::vector<geom_float> sv(3);

    std::vector<geom_float> dccv(3);
    geom_float dcc;

    std::vector<geom_float> dc1pv(3);
    geom_float dc1p;

    std::vector<geom_float> dc2pv(3);
    geom_float dc2p;

    std::vector<geom_float> pcent(3);
    std::vector<geom_float> c1cent(3);
    std::vector<geom_float> c2cent(3);

    geom_float volume;

    geom_int ic1;
    geom_int ic2;

    geom_float f;
    std::vector<flow_float>& fxp = v.p["fx"]; std::vector<flow_float>& dccp = v.p["dcc"];
    std::vector<flow_float>& vvol = v.c["volume"];
    std::vector<flow_float>& vccx = v.c["ccx"];
    std::vector<flow_float>& vccy = v.c["ccy"];
    std::vector<flow_float>& vccz = v.c["ccz"];

    for (geom_int ip=0 ; ip<msh.nPlanes ; ip++)
    {
        ic1     = msh.planes[ip].iCells[0];
        ic2     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;

        pcent   = msh.planes[ip].centCoords;

        c1cent  = msh.cells[ic1].centCoords;
        c2cent  = msh.cells[ic2].centCoords;

        dccv[0] = c2cent[0] - c1cent[0];
        dccv[1] = c2cent[1] - c1cent[1];
        dccv[2] = c2cent[2] - c1cent[2];
        dcc     = sqrt( pow(dccv[0], 2.0) + pow(dccv[1], 2.0) + pow(dccv[2], 2.0));

        dc2pv[0] = pcent[0] - c2cent[0];
        dc2pv[1] = pcent[1] - c2cent[1];
        dc2pv[2] = pcent[2] - c2cent[2];
        dc2p     = sqrt( pow(dc2pv[0], 2.0) + pow(dc2pv[1], 2.0) + pow(dc2pv[2], 2.0));

        fxp[ip]  = dc2p/dcc;
        dccp[ip] = dcc;

    }

    // cell
    for (geom_int ic=0 ; ic<msh.nCells ; ic++)
    {
        //if (cfg.initial == "Taylor") {
        //    msh.cells[ic].volume = (2*M_PI)*(2*M_PI)*(2*M_PI)/(pow(32.0,3));
        //}
        c1cent   = msh.cells[ic].centCoords;
        vvol[ic] = msh.cells[ic].volume;
        vccx[ic] = c1cent[0];
        vccy[ic] = c1cent[1];
        vccz[ic] = c1cent[2];
    }

}
