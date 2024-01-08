#include "gradient.hpp"
#include "flowFormat.hpp"

#include <cmath>

void gradientGauss(solverConfig& cfg , mesh& msh , variables& v)
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

    std::vector<flow_float>& dPdx = v.c["dPdx"];
    std::vector<flow_float>& dPdy = v.c["dPdy"];
    std::vector<flow_float>& dPdz = v.c["dPdz"];

    std::vector<flow_float>& dTdx = v.c["dTdx"];
    std::vector<flow_float>& dTdy = v.c["dTdy"];
    std::vector<flow_float>& dTdz = v.c["dTdz"];

    std::vector<flow_float>& Ux = v.c["Ux"];
    std::vector<flow_float>& Uy = v.c["Uy"];
    std::vector<flow_float>& Uz = v.c["Uz"];
    std::vector<flow_float>& P = v.c["P"];
    std::vector<flow_float>& T = v.c["T"];

    std::vector<flow_float>& fxp = v.p["fx"];

    for (geom_int ic=0 ; ic<msh.nCells_all; ic++)
    {
        dUxdx[ic] =  0.0;
        dUxdy[ic] =  0.0;
        dUxdz[ic] =  0.0;

        dUydx[ic] =  0.0;
        dUydy[ic] =  0.0;
        dUydz[ic] =  0.0;

        dUzdx[ic] =  0.0;
        dUzdy[ic] =  0.0;
        dUzdz[ic] =  0.0;

        dPdx[ic] =  0.0;
        dPdy[ic] =  0.0;
        dPdz[ic] =  0.0;

        dTdx[ic] =  0.0;
        dTdy[ic] =  0.0;
        dTdz[ic] =  0.0;
    }

    // normal plane
    for (geom_int ip=0 ; ip<msh.nPlanes ; ip++)
    {
        ic1     = msh.planes[ip].iCells[0];
        ic2     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;

        f = fxp[ip];

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

    geom_float volume;
    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        volume = msh.cells[ic].volume;
        dUxdx[ic] = dUxdx[ic]/volume;
        dUxdy[ic] = dUxdy[ic]/volume;
        dUxdz[ic] = dUxdz[ic]/volume;

        dUydx[ic] = dUydx[ic]/volume;
        dUydy[ic] = dUydy[ic]/volume;
        dUydz[ic] = dUydz[ic]/volume;

        dUzdx[ic] = dUzdx[ic]/volume;
        dUzdy[ic] = dUzdy[ic]/volume;
        dUzdz[ic] = dUzdz[ic]/volume;

        dPdx[ic] = dPdx[ic]/volume;
        dPdy[ic] = dPdy[ic]/volume;
        dPdz[ic] = dPdz[ic]/volume;

        dTdx[ic] = dTdx[ic]/volume;
        dTdy[ic] = dTdy[ic]/volume;
        dTdz[ic] = dTdz[ic]/volume;
    }

}
