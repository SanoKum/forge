#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <mesh.hpp>
#include <variables.hpp>

#include <cmath>

void calcDiffusion(solverConfig &cfg , mesh &msh , variables &v , matrix& mat_p)
{
    geom_float ss;
    vector<geom_float> sv(3);
    vector<geom_float> nv(3); // normal

    vector<geom_float> nv_dia(3); // diagonal
    vector<geom_float> nv_nodia(3); // non-diagonal

    vector<geom_float> pcent(3);
    vector<geom_float> c1cent(3);
    vector<geom_float> c2cent(3);
    geom_int ic0;
    geom_int ic1;
    geom_float dn;
    flow_float temp;


    geom_int ip_loc0;
    geom_int ip_loc1;
    
    vector<geom_float> dccv(3);
    geom_float dcc;

    vector<geom_float> dc1pv(3);
    geom_float dc1p;

    vector<geom_float> dc2pv(3);
    geom_float dc2p;

    geom_float dccv_dot_nv;

    geom_float f;

    geom_float temp_ndia;

    //v.cel["diffT"][ic1] =  0.0;
    //v.cel["diffT"][ic2] =  0.0;

    // normal plane
    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        ic0     = msh.planes[ip].iCells[0];
        ic1     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;
        ss      = msh.planes[ip].surfArea;

        nv[0]   = sv[0]/ss;
        nv[1]   = sv[1]/ss;
        nv[2]   = sv[2]/ss;

        pcent   = msh.planes[ip].centCoords;

        c1cent  = msh.cells[ic0].centCoords;
        c2cent  = msh.cells[ic1].centCoords;

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

        dccv_dot_nv = dccv[0]*nv[0] + dccv[1]*nv[1] + dccv[2]*nv[2];

        nv_nodia[0] = nv[0] - dccv[0]/dccv_dot_nv;
        nv_nodia[1] = nv[1] - dccv[1]/dccv_dot_nv;
        nv_nodia[2] = nv[2] - dccv[2]/dccv_dot_nv;

        f = dc2p/dcc;

        // over-relaxed approach
        temp_ndia = +(   (f*v.c["dTdx"][ic0]+(1.0-f)*v.c["dTdx"][ic1]) *nv_nodia[0]   // non-diagonal term 
                       + (f*v.c["dTdy"][ic0]+(1.0-f)*v.c["dTdy"][ic1]) *nv_nodia[1] 
                       + (f*v.c["dTdz"][ic0]+(1.0-f)*v.c["dTdz"][ic1]) *nv_nodia[2] )*ss ;

        temp = (v.c["T"][ic1] - v.c["T"][ic0])/dcc*ss // diagonal term
                + temp_ndia;

        v.c["diffT"][ic0] +=  temp;
        v.c["diffT"][ic1] += -temp;

        ip_loc0  = mat_p.localPlnOfCell[ip][0];
        ip_loc1  = mat_p.localPlnOfCell[ip][1];

        mat_p.lhs[ic0][0] += ss/dcc;
        mat_p.lhs[ic0][ip_loc0] = -ss/dcc;
        mat_p.rhs[ic0] += +temp_ndia;

        mat_p.lhs[ic1][0] += ss/dcc;
        mat_p.lhs[ic1][ip_loc1] = -ss/dcc;
        mat_p.rhs[ic1] += -temp_ndia;
    }

//    for (bcond& bc : msh.bconds)
//    {
//        for (geom_int& ip : bc.iPlanes)
//        {
//            ic0     = msh.planes[ip].iCells[0];
//            sv      = msh.planes[ip].surfVect;
//            ss      = msh.planes[ip].surfArea;
//            pcent   = msh.planes[ip].centCoords;
//    
//            c1cent  = msh.cells[ic0].centCoords;
//    
//            dn = ( (pcent[0] - c1cent[0])*sv[0]
//                  +(pcent[1] - c1cent[1])*sv[1]
//                  +(pcent[2] - c1cent[2])*sv[2] )/ss;
//        
//            v.c["diffT"][ic0] += (v.p["T"][ip] - v.c["T"][ic0])/dn*ss;
//
//            mat_p.lhs[ic0][0] += ss/dn;
//            mat_p.rhs[ic0] += v.p["T"][ip]/dn*ss; 
//        }
//    }

}

#endif