#include "calcCFL.hpp"

#include <cmath>

using namespace std;

void calcCFL(solverConfig& cfg, mesh& msh, variables& v)
{
    vector<flow_float>& cfl = v.c["cfl"];
    vector<flow_float>& US  = v.p["US"];

    // 界面の流速でCFLを見る
    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        geom_int ic0 = msh.planes[ip].iCells[0];
        geom_int ic1 = msh.planes[ip].iCells[1];

        geom_float vol0 = msh.cells[ic0].volume;
        geom_float vol1 = msh.cells[ic1].volume;

        geom_float surfArea = msh.planes[ip].surfArea;

        geom_float dx0 = vol0/surfArea;
        geom_float dx1 = vol0/surfArea;

        cfl[ic0] = max(cfl[ic0] , cfg.dt*abs(US[ip])/dx0/surfArea);
        cfl[ic1] = max(cfl[ic1] , cfg.dt*abs(US[ip])/dx1/surfArea);
    }

    for (geom_int ip=msh.nNormalPlanes ; ip<msh.nPlanes ; ip++)
    {
        geom_int ic0 = msh.planes[ip].iCells[0];
        geom_float vol0 = msh.cells[ic0].volume;
        geom_float surfArea = msh.planes[ip].surfArea;
        geom_float dx0 = vol0/surfArea;
        cfl[ic0] = max(cfl[ic0] , cfg.dt*US[ip]/dx0/surfArea);
    }


};