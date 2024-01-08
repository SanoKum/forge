#ifndef UPDATE_H
#define UPDATE_H

#include "update.hpp"
#include <vector>

using namespace std;

void updateVariablesForNextLoop(solverConfig& cfg , mesh& msh , variables& v , matrix& mat_p)
{
    vector<flow_float>& US = v.p["US"];
    vector<flow_float>& USN = v.p["USN"];

    vector<flow_float>& Ux = v.c["Ux"];
    vector<flow_float>& Uy = v.c["Uy"];
    vector<flow_float>& Uz = v.c["Uz"];

    vector<flow_float>& UxN = v.c["UxN"];
    vector<flow_float>& UyN = v.c["UyN"];
    vector<flow_float>& UzN = v.c["UzN"];

    vector<flow_float>& T = v.c["T"];
    vector<flow_float>& TN = v.c["TN"];

    vector<flow_float>& dUxdx = v.c["dUxdx"];
    vector<flow_float>& dUydx = v.c["dUydx"];
    vector<flow_float>& dUzdx = v.c["dUzdx"];

    vector<flow_float>& dUxdy = v.c["dUxdy"];
    vector<flow_float>& dUydy = v.c["dUydy"];
    vector<flow_float>& dUzdy = v.c["dUzdy"];

    vector<flow_float>& dUxdz = v.c["dUxdz"];
    vector<flow_float>& dUydz = v.c["dUydz"];
    vector<flow_float>& dUzdz = v.c["dUzdz"];

    vector<flow_float>& dPdx = v.c["dPdx"];
    vector<flow_float>& dPdy = v.c["dPdy"];
    vector<flow_float>& dPdz = v.c["dPdz"];

    vector<flow_float>& convx = v.c["convx"];
    vector<flow_float>& convy = v.c["convy"];
    vector<flow_float>& convz = v.c["convz"];

    vector<flow_float>& diffx = v.c["diffx"];
    vector<flow_float>& diffy = v.c["diffy"];
    vector<flow_float>& diffz = v.c["diffz"];

    vector<flow_float>& src = v.c["src"];
    vector<flow_float>& res = v.c["res"];

    vector<flow_float>& diffT = v.c["diffT"];
    vector<flow_float>& dP = v.c["dP"];
    vector<flow_float>& dPPdx = v.c["dPPdx"];
    vector<flow_float>& dPPdy = v.c["dPPdy"];
    vector<flow_float>& dPPdz = v.c["dPPdz"];
    vector<flow_float>& divU = v.c["divU"];
    vector<flow_float>& divU_vol = v.c["divU*vol"];
    vector<flow_float>& divU_star = v.c["divU_star"];

    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {

        UxN[ic] = Ux[ic];
        UyN[ic] = Uy[ic];
        UzN[ic] = Uz[ic];

        TN[ic] = T[ic];

        convx[ic] = 0.0;
        convy[ic] = 0.0;
        convz[ic] = 0.0;

        diffx[ic] = 0.0;
        diffy[ic] = 0.0;
        diffz[ic] = 0.0;

        dP[ic] = 0.0;
        dPPdx[ic] = 0.0;
        dPPdy[ic] = 0.0;
        dPPdz[ic] = 0.0;

        divU[ic] = 0.0;
        divU_vol[ic] = 0.0;
        divU_star[ic] = 0.0;
    }

    geom_int ic0, ic1;
    geom_int ip_loc0, ip_loc1;

    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        USN = US;

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