#include "gradient.hpp"
#include "flowFormat.hpp"

#include <cmath>

void calcGradient(solverConfig& cfg , mesh& msh , variables& v)
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
    geom_int ic1;
    geom_int ic2;

    geom_float f;

    flow_float Pf;
    flow_float Tf;
    flow_float Uxf , Uyf , Uzf;

    std::vector<flow_float>& dUxdx = v.c["dUxdx"];
    std::vector<flow_float>& dUxdy = v.c["dUxdy"];
    std::vector<flow_float>& dUxdz = v.c["dUxdz"];

    std::vector<flow_float>& dUydx = v.c["dUydx"];
    std::vector<flow_float>& dUydy = v.c["dUydy"];
    std::vector<flow_float>& dUydz = v.c["dUydz"];

    std::vector<flow_float>& dUzdx = v.c["dUzdx"];
    std::vector<flow_float>& dUzdy = v.c["dUzdy"];
    std::vector<flow_float>& dUzdz = v.c["dUzdz"];

    std::vector<flow_float>& Ux = v.c["Ux"];
    std::vector<flow_float>& Uy = v.c["Uy"];
    std::vector<flow_float>& Uz = v.c["Uz"];

    std::vector<flow_float>& P = v.c["P"];
    std::vector<flow_float>& T = v.c["T"];

    std::vector<flow_float>& Uxp = v.p["Ux"];
    std::vector<flow_float>& Uyp = v.p["Uy"];
    std::vector<flow_float>& Uzp = v.p["Uz"];


    std::vector<flow_float>& Pp = v.p["P"];
    std::vector<flow_float>& Tp = v.p["T"];

    std::vector<flow_float>& dPdx = v.c["dPdx"];
    std::vector<flow_float>& dPdy = v.c["dPdy"];
    std::vector<flow_float>& dPdz = v.c["dPdz"];

    std::vector<flow_float>& dTdx = v.c["dTdx"];
    std::vector<flow_float>& dTdy = v.c["dTdy"];
    std::vector<flow_float>& dTdz = v.c["dTdz"];

    for (geom_int ie=0 ; ie<msh.nCells; ie++)
    {
        dUxdx[ie] =  0.0;
        dUxdy[ie] =  0.0;
        dUxdz[ie] =  0.0;

        dUydx[ie] =  0.0;
        dUydy[ie] =  0.0;
        dUydz[ie] =  0.0;

        dUzdx[ie] =  0.0;
        dUzdy[ie] =  0.0;
        dUzdz[ie] =  0.0;

        dPdx[ie] =  0.0;
        dPdy[ie] =  0.0;
        dPdz[ie] =  0.0;

        dTdx[ie] =  0.0;
        dTdy[ie] =  0.0;
        dTdz[ie] =  0.0;
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

        Uxf= f*Ux[ic1]+ (1.0-f)*Ux[ic2] ;
        Uyf= f*Uy[ic1]+ (1.0-f)*Uy[ic2] ;
        Uzf= f*Uz[ic1]+ (1.0-f)*Uz[ic2] ;

        Pf = f*P[ic1] + (1.0-f)*P[ic2] ;
        Tf = f*T[ic1] + (1.0-f)*T[ic2] ;

        dUxdx[ic1] +=  sv[0]*Uxf;
        dUxdy[ic1] +=  sv[1]*Uxf;
        dUxdz[ic1] +=  sv[2]*Uxf;

        dUxdx[ic2] += -sv[0]*Uxf;
        dUxdy[ic2] += -sv[1]*Uxf;
        dUxdz[ic2] += -sv[2]*Uxf;


        dUydx[ic1] +=  sv[0]*Uyf;
        dUydy[ic1] +=  sv[1]*Uyf;
        dUydz[ic1] +=  sv[2]*Uyf;

        dUydx[ic2] += -sv[0]*Uyf;
        dUydy[ic2] += -sv[1]*Uyf;
        dUydz[ic2] += -sv[2]*Uyf;


        dUzdx[ic1] +=  sv[0]*Uzf;
        dUzdy[ic1] +=  sv[1]*Uzf;
        dUzdz[ic1] +=  sv[2]*Uzf;

        dUzdx[ic2] += -sv[0]*Uzf;
        dUzdy[ic2] += -sv[1]*Uzf;
        dUzdz[ic2] += -sv[2]*Uzf;


        dTdx[ic1] +=  sv[0]*Tf;
        dTdy[ic1] +=  sv[1]*Tf;
        dTdz[ic1] +=  sv[2]*Tf;

        dTdx[ic2] += -sv[0]*Tf;
        dTdy[ic2] += -sv[1]*Tf;
        dTdz[ic2] += -sv[2]*Tf;

        dPdx[ic1] +=  sv[0]*Pf;
        dPdy[ic1] +=  sv[1]*Pf;
        dPdz[ic1] +=  sv[2]*Pf;

        dPdx[ic2] += -sv[0]*Pf;
        dPdy[ic2] += -sv[1]*Pf;
        dPdz[ic2] += -sv[2]*Pf;

    }

    // boundary plane
    for (geom_int ip=msh.nNormalPlanes ; ip<msh.nPlanes ; ip++)
    {
        ic1     = msh.planes[ip].iCells[0];
        sv      = msh.planes[ip].surfVect;

        dUxdx[ic1] +=  sv[0]*Uxp[ip];
        dUxdy[ic1] +=  sv[1]*Uxp[ip];
        dUxdz[ic1] +=  sv[2]*Uxp[ip];

        dUydx[ic1] +=  sv[0]*Uyp[ip];
        dUydy[ic1] +=  sv[1]*Uyp[ip];
        dUydz[ic1] +=  sv[2]*Uyp[ip];

        dUzdx[ic1] +=  sv[0]*Uzp[ip];
        dUzdy[ic1] +=  sv[1]*Uzp[ip];
        dUzdz[ic1] +=  sv[2]*Uzp[ip];
    
        dTdx[ic1] +=  sv[0]*Tp[ip];
        dTdy[ic1] +=  sv[1]*Tp[ip];
        dTdz[ic1] +=  sv[2]*Tp[ip];

        dPdx[ic1] +=  sv[0]*Pp[ip];
        dPdy[ic1] +=  sv[1]*Pp[ip];
        dPdz[ic1] +=  sv[2]*Pp[ip];

    }

    geom_float volume;
    for (geom_int ie=0 ; ie<msh.nCells; ie++)
    {
        volume = msh.cells[ie].volume;
        dUxdx[ie] = dUxdx[ie]/volume;
        dUxdy[ie] = dUxdy[ie]/volume;
        dUxdz[ie] = dUxdz[ie]/volume;

        dUydx[ie] = dUydx[ie]/volume;
        dUydy[ie] = dUydy[ie]/volume;
        dUydz[ie] = dUydz[ie]/volume;

        dUzdx[ie] = dUzdx[ie]/volume;
        dUzdy[ie] = dUzdy[ie]/volume;
        dUzdz[ie] = dUzdz[ie]/volume;

        dPdx[ie] = dPdx[ie]/volume;
        dPdy[ie] = dPdy[ie]/volume;
        dPdz[ie] = dPdz[ie]/volume;

        dTdx[ie] = dTdx[ie]/volume;
        dTdy[ie] = dTdy[ie]/volume;
        dTdz[ie] = dTdz[ie]/volume;
    }

}
