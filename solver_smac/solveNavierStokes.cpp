#include "solveNavierStokes.hpp"
#include <cmath>

using namespace std;

void solveNavierStokes(solverConfig &cfg , mesh &msh , variables &v , matrix& mat_ns)
{
    geom_float ss;
    vector<geom_float> sv(3);
    vector<geom_float> nv(3); // normal

    vector<geom_float> sv_dia(3); // diagonal
    vector<geom_float> sv_nodia(3); // non-diagonal

    vector<geom_float> pcent(3);
    vector<geom_float> c0cent(3);
    vector<geom_float> c1cent(3);
    geom_int ic0;
    geom_int ic1;
    geom_float dn;
    flow_float temp;


    geom_int ip_loc0;
    geom_int ip_loc1;
    
    vector<geom_float> dccv(3);
    geom_float dcc;

    vector<geom_float> dc0pv(3);
    geom_float dc0p;

    vector<geom_float> dc1pv(3);
    geom_float dc1p;

    geom_float dccv_dot_sv;

    vector<geom_float> deltav(3);
    geom_float delta;

    geom_float temp_ndia;
    geom_float temp_ndia_x;
    geom_float temp_ndia_y;
    geom_float temp_ndia_z;

    flow_float US;
    flow_float rof;
    flow_float volf;
    flow_float mass;

    geom_float cosT;

    vector<flow_float>& Ux = v.c["Ux"];
    vector<flow_float>& Uy = v.c["Uy"];
    vector<flow_float>& Uz = v.c["Uz"];

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

    vector<flow_float>& UxN = v.c["UxN"];
    vector<flow_float>& UyN = v.c["UyN"];
    vector<flow_float>& UzN = v.c["UzN"];

    vector<flow_float>& divU_vol = v.c["divU*vol"];
    vector<flow_float>& divU_star = v.c["divU_star"];

    vector<flow_float>& P = v.c["P"];
    vector<flow_float>& T = v.c["T"];
    vector<flow_float>& ro = v.c["ro"];

    vector<flow_float>& convx = v.c["convx"];
    vector<flow_float>& convy = v.c["convy"];
    vector<flow_float>& convz = v.c["convz"];

    vector<flow_float>& diffx = v.c["diffx"];
    vector<flow_float>& diffy = v.c["diffy"];
    vector<flow_float>& diffz = v.c["diffz"];


    vector<flow_float>& Uxp = v.p["Ux"];
    vector<flow_float>& Uyp = v.p["Uy"];
    vector<flow_float>& Uzp = v.p["Uz"];

    vector<flow_float>& rop = v.p["ro"];
    vector<flow_float>& USp = v.p["US"];
    vector<flow_float>& Pp  = v.p["P"];

    vector<flow_float>& fxp  = v.p["fx"];


    // normal plane
    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        ic0     = msh.planes[ip].iCells[0];
        ic1     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;
        ss      = msh.planes[ip].surfArea;

        //nv[0]   = sv[0]/ss;
        //nv[1]   = sv[1]/ss;
        //nv[2]   = sv[2]/ss;

        pcent   = msh.planes[ip].centCoords;

        c0cent  = msh.cells[ic0].centCoords;
        c1cent  = msh.cells[ic1].centCoords;

        dccv[0] = c1cent[0] - c0cent[0];
        dccv[1] = c1cent[1] - c0cent[1];
        dccv[2] = c1cent[2] - c0cent[2];
        dcc     = sqrt( pow(dccv[0], 2.0) + pow(dccv[1], 2.0) + pow(dccv[2], 2.0));

        dccv_dot_sv = dccv[0]*sv[0] + dccv[1]*sv[1] + dccv[2]*sv[2];

        deltav[0] = dccv[0]*pow(ss, 2.0)/dccv_dot_sv;
        deltav[1] = dccv[1]*pow(ss, 2.0)/dccv_dot_sv;
        deltav[2] = dccv[2]*pow(ss, 2.0)/dccv_dot_sv;
        delta = dcc*pow(ss, 2.0)/dccv_dot_sv;

        sv_nodia[0] = sv[0] - deltav[0];
        sv_nodia[1] = sv[1] - deltav[1];
        sv_nodia[2] = sv[2] - deltav[2];


        //dc0pv[0] = pcent[0] - c0cent[0];
        //dc0pv[1] = pcent[1] - c0cent[1];
        //dc0pv[2] = pcent[2] - c0cent[2];
        //dc0p    = sqrt( pow(dc0pv[0], 2.0) + pow(dc0pv[1], 2.0) + pow(dc0pv[2], 2.0));

        //dc1pv[0] = pcent[0] - c1cent[0];
        //dc1pv[1] = pcent[1] - c1cent[1];
        //dc1pv[2] = pcent[2] - c1cent[2];
        //dc1p    = sqrt( pow(dc1pv[0], 2.0) + pow(dc1pv[1], 2.0) + pow(dc1pv[2], 2.0));

        //f = dc1p/dcc;

        // conve u
        rof = fxp[ip]*ro[ic0] +(1.0-fxp[ip])*ro[ic1];
        US = (fxp[ip]*Ux[ic0] +(1.0-fxp[ip])*Ux[ic1])*sv[0]
            +(fxp[ip]*Uy[ic0] +(1.0-fxp[ip])*Uy[ic1])*sv[1]
            +(fxp[ip]*Uz[ic0] +(1.0-fxp[ip])*Uz[ic1])*sv[2];
        //mass = rof*US;

        convx[ic0] += 0.5*rof*( (US+abs(US))*Ux[ic0] + (US-abs(US))*Ux[ic1] );
        convy[ic0] += 0.5*rof*( (US+abs(US))*Uy[ic0] + (US-abs(US))*Uy[ic1] );
        convz[ic0] += 0.5*rof*( (US+abs(US))*Uz[ic0] + (US-abs(US))*Uz[ic1] );

        convx[ic1] -= 0.5*rof*( (US+abs(US))*Ux[ic0] + (US-abs(US))*Ux[ic1] );
        convy[ic1] -= 0.5*rof*( (US+abs(US))*Uy[ic0] + (US-abs(US))*Uy[ic1] );
        convz[ic1] -= 0.5*rof*( (US+abs(US))*Uz[ic0] + (US-abs(US))*Uz[ic1] );

        temp_ndia_x = (fxp[ip]*dUxdx[ic0]+(1.0-fxp[ip])*dUxdx[ic1]) *sv_nodia[0]   // non-diagonal term 
                    + (fxp[ip]*dUxdy[ic0]+(1.0-fxp[ip])*dUxdy[ic1]) *sv_nodia[1] 
                    + (fxp[ip]*dUxdz[ic0]+(1.0-fxp[ip])*dUxdz[ic1]) *sv_nodia[2] ;

        temp_ndia_y = (fxp[ip]*dUydx[ic0]+(1.0-fxp[ip])*dUydx[ic1]) *sv_nodia[0]   // non-diagonal term 
                    + (fxp[ip]*dUydy[ic0]+(1.0-fxp[ip])*dUydy[ic1]) *sv_nodia[1] 
                    + (fxp[ip]*dUydz[ic0]+(1.0-fxp[ip])*dUydz[ic1]) *sv_nodia[2] ;

        temp_ndia_z =  (fxp[ip]*dUzdx[ic0]+(1.0-fxp[ip])*dUzdx[ic1]) *sv_nodia[0]   // non-diagonal term 
                     + (fxp[ip]*dUzdy[ic0]+(1.0-fxp[ip])*dUzdy[ic1]) *sv_nodia[1] 
                     + (fxp[ip]*dUzdz[ic0]+(1.0-fxp[ip])*dUzdz[ic1]) *sv_nodia[2] ;


        diffx[ic0] +=  cfg.visc*( (Ux[ic1] - Ux[ic0])/dcc*delta + temp_ndia_x);
        diffy[ic0] +=  cfg.visc*( (Uy[ic1] - Uy[ic0])/dcc*delta + temp_ndia_y);
        diffz[ic0] +=  cfg.visc*( (Uz[ic1] - Uz[ic0])/dcc*delta + temp_ndia_z);

        diffx[ic1] -=  cfg.visc*( (Ux[ic1] - Ux[ic0])/dcc*delta + temp_ndia_x);
        diffy[ic1] -=  cfg.visc*( (Uy[ic1] - Uy[ic0])/dcc*delta + temp_ndia_y);
        diffz[ic1] -=  cfg.visc*( (Uz[ic1] - Uz[ic0])/dcc*delta + temp_ndia_z);
    }

    // boundary conditions
    for (auto& bc : msh.bconds)
    {
        if (bc.valueTypes["Ux"] == -1) continue; // ignore slip boundaries.

        // calculate convect & diffusion term
        for (geom_int& ip : bc.iPlanes)
        {
            ic0     = msh.planes[ip].iCells[0];
            sv      = msh.planes[ip].surfVect;
            ss      = msh.planes[ip].surfArea;
            pcent   = msh.planes[ip].centCoords;

            c0cent  = msh.cells[ic0].centCoords;

            US = (Uxp[ip]*sv[0] + Uyp[ip]*sv[1] + Uzp[ip]*sv[2]);

            convx[ic0] += US*Uxp[ip]*rop[ip];
            convy[ic0] += US*Uyp[ip]*rop[ip];
            convz[ic0] += US*Uzp[ip]*rop[ip];

            dn = ( (pcent[0] - c0cent[0])*sv[0]
                  +(pcent[1] - c0cent[1])*sv[1]
                  +(pcent[2] - c0cent[2])*sv[2] )/ss;

            diffx[ic0] += cfg.visc*(Uxp[ip] - Ux[ic0])/dn*ss;
            diffy[ic0] += cfg.visc*(Uyp[ip] - Uy[ic0])/dn*ss;
            diffz[ic0] += cfg.visc*(Uzp[ip] - Uz[ic0])/dn*ss;
        }
    }

    // time integration
    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        geom_float vol = msh.cells[ic].volume;

        Ux[ic] = UxN[ic] -dPdx[ic]/ro[ic]*cfg.dt + (-convx[ic] + diffx[ic] )*cfg.dt/vol/ro[ic];
        Uy[ic] = UyN[ic] -dPdy[ic]/ro[ic]*cfg.dt + (-convy[ic] + diffy[ic] )*cfg.dt/vol/ro[ic];
        Uz[ic] = UzN[ic] -dPdz[ic]/ro[ic]*cfg.dt + (-convz[ic] + diffz[ic] )*cfg.dt/vol/ro[ic];

    }

    geom_float vol0;
    geom_float vol1;


    flow_float Uxf, Uyf, Uzf;


    // rhi-chow interporation
    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        ic0     = msh.planes[ip].iCells[0];
        ic1     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;
        ss      = msh.planes[ip].surfArea;

        vol0 = msh.cells[ic0].volume;
        vol1 = msh.cells[ic1].volume;

        pcent   = msh.planes[ip].centCoords;

        c0cent  = msh.cells[ic0].centCoords;
        c1cent  = msh.cells[ic1].centCoords;

        dccv[0] = c1cent[0] - c0cent[0];
        dccv[1] = c1cent[1] - c0cent[1];
        dccv[2] = c1cent[2] - c0cent[2];
        dcc     = sqrt( pow(dccv[0], 2.0) + pow(dccv[1], 2.0) + pow(dccv[2], 2.0));

        //dc0pv[0] = pcent[0] - c0cent[0];
        //dc0pv[1] = pcent[1] - c0cent[1];
        //dc0pv[2] = pcent[2] - c0cent[2];
        //dc0p    = sqrt( pow(dc0pv[0], 2.0) + pow(dc0pv[1], 2.0) + pow(dc0pv[2], 2.0));

        //dc1pv[0] = pcent[0] - c1cent[0];
        //dc1pv[1] = pcent[1] - c1cent[1];
        //dc1pv[2] = pcent[2] - c1cent[2];
        //dc1p    = sqrt( pow(dc1pv[0], 2.0) + pow(dc1pv[1], 2.0) + pow(dc1pv[2], 2.0));

        dccv_dot_sv = dccv[0]*sv[0] + dccv[1]*sv[1] + dccv[2]*sv[2];

        deltav[0] = dccv[0]*pow(ss, 2.0)/dccv_dot_sv;
        deltav[1] = dccv[1]*pow(ss, 2.0)/dccv_dot_sv;
        deltav[2] = dccv[2]*pow(ss, 2.0)/dccv_dot_sv;
        delta = dcc*pow(ss, 2.0)/dccv_dot_sv;

        sv_nodia[0] = sv[0] - deltav[0];
        sv_nodia[1] = sv[1] - deltav[1];
        sv_nodia[2] = sv[2] - deltav[2];

        // conve u
        rof  = fxp[ip]*ro[ic0] +(1.0-fxp[ip])*ro[ic1];
        rop[ip] = rof;

        Uxf = fxp[ip]*(Ux[ic0] + dPdx[ic0]*cfg.dt/ro[ic0]) +(1.0-fxp[ip])*(Ux[ic1] + dPdx[ic1]*cfg.dt/ro[ic1]);
        Uyf = fxp[ip]*(Uy[ic0] + dPdy[ic0]*cfg.dt/ro[ic0]) +(1.0-fxp[ip])*(Uy[ic1] + dPdy[ic1]*cfg.dt/ro[ic1]);
        Uzf = fxp[ip]*(Uz[ic0] + dPdz[ic0]*cfg.dt/ro[ic0]) +(1.0-fxp[ip])*(Uz[ic1] + dPdz[ic1]*cfg.dt/ro[ic1]);

        US = Uxf*sv[0] + Uyf*sv[1] + Uzf*sv[2];

        temp_ndia = (fxp[ip]*dPdx[ic0]+(1.0-fxp[ip])*dPdx[ic1]) *sv_nodia[0]   // non-diagonal term 
                  + (fxp[ip]*dPdy[ic0]+(1.0-fxp[ip])*dPdy[ic1]) *sv_nodia[1] 
                  + (fxp[ip]*dPdz[ic0]+(1.0-fxp[ip])*dPdz[ic1]) *sv_nodia[2]  ;

        US += (-(P[ic1] - P[ic0])/dcc*delta - temp_ndia)*cfg.dt/rof;

        USp[ip] = US;

        divU_vol[ic0] += US;
        divU_vol[ic1] -= US;
    }

    // boundary conditions
    for (auto& bc : msh.bconds)
    {
        if (bc.bcondKind == "outlet_statPress" ) 
        {
            for (geom_int& ip : bc.iPlanes)
            {
                ic0     = msh.planes[ip].iCells[0];
                sv      = msh.planes[ip].surfVect;
                ss      = msh.planes[ip].surfArea;
                pcent   = msh.planes[ip].centCoords;

                c0cent  = msh.cells[ic0].centCoords;

                flow_float dn = ( (pcent[0] - c0cent[0])*sv[0]
                                 +(pcent[1] - c0cent[1])*sv[1]
                                 +(pcent[2] - c0cent[2])*sv[2] )/ss;

                flow_float US;

                US = (Ux[ic0] + dPdx[ic0]*cfg.dt/ro[ic0])*sv[0]
                    +(Uy[ic0] + dPdy[ic0]*cfg.dt/ro[ic0])*sv[1]
                    +(Uz[ic0] + dPdz[ic0]*cfg.dt/ro[ic0])*sv[2];

                US += (-(Pp[ip] - P[ic0])/dn*ss)*cfg.dt/ro[ic0];
                USp[ip] = US;
                divU_vol[ic0] += US;
            }
        } else if (bc.bcondKind == "inlet_uniformVelocity" ) 
        {
            // calculate convect & diffusion term
            for (geom_int& ip : bc.iPlanes)
            {
                ic0     = msh.planes[ip].iCells[0];
                sv      = msh.planes[ip].surfVect;
                ss      = msh.planes[ip].surfArea;
                pcent   = msh.planes[ip].centCoords;

                //c0cent  = msh.cells[ic0].centCoords;

                US = Uxp[ip]*sv[0] +Uyp[ip]*sv[1] +Uzp[ip]*sv[2];
                USp[ip] = US;
                divU_vol[ic0] += US;
            }
        } else {
            for (geom_int& ip : bc.iPlanes)
            {
                USp[ip] = 0.0;
            }
        }
    }


    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        geom_float vol = msh.cells[ic].volume;

        divU_star[ic]= divU_vol[ic] / vol;
    }

}
