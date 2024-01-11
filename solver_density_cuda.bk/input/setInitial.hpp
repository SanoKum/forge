#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "solverConfig.hpp"

#include <cmath>

void setInitial(solverConfig& cfg , mesh& msh , variables& v)
{
    // *** set initial value ***
    for (geom_int i = 0 ; i<msh.nCells ; i++)
    {
//        flow_float x = msh.cells[i].centCoords[0];
//
//        flow_float ro = 1.2;
//        flow_float T  = 270;
//        flow_float Ux = 80.0;
//        flow_float Uy = 0.0;
//        flow_float Uz = 0.0;
//        flow_float P  = ro*T*(cfg.gamma-1.0)/cfg.gamma*cfg.cp;
//        flow_float roe  = ro*(0.5*(Ux*Ux+Uy*Uy+Uz*Uz) +cfg.cp/cfg.gamma*T);
//
//        v.c["T"][i] = T;
//        v.c["Ux"][i] = Ux;
//        v.c["Uy"][i] = Uy;
//        v.c["Uz"][i] = Uz;
//        v.c["P"][i] = P;
//        v.c["ro"][i] = ro;
//        v.c["roUx"][i] = ro*Ux;
//        v.c["roUy"][i] = ro*Uy;
//        v.c["roUz"][i] = ro*Uz;
//        v.c["roe"][i]  = roe;
        
        // 1D_shock_tube

        if (cfg.initial == "sod") 
        {
            flow_float gam  = 1.4;
            flow_float roL  = 1.0;
            flow_float prsL = 1.0;
            flow_float velL = 0.0;
            flow_float roR  = 0.125;
            flow_float prsR = 0.1;
            flow_float velR = 0.0;

            flow_float x = msh.cells[i].centCoords[0];

            if (x<0.5) {
                v.c["ro"][i] = roL;
                v.c["roUx"][i] = roL*velL;
                v.c["roUy"][i] = 0.0;
                v.c["roUz"][i] = 0.0;
                v.c["roe"][i] = prsL/(gam-1.0) + 0.5*roL*velL*velL;
            } else {
                v.c["ro"][i] = roR;
                v.c["roUx"][i] = roR*velR;
                v.c["roUy"][i] = 0.0;
                v.c["roUz"][i] = 0.0;
                v.c["roe"][i] = prsR/(gam-1.0) + 0.5*roR*velR*velR;
            }
        }

        if (cfg.initial == "mach3")
        {
            // mach3 
            flow_float gam  = 1.4;
            flow_float roL  = 1.4;
            flow_float prsL = 1.0;
            flow_float velL = 3.0;

            flow_float x = msh.cells[i].centCoords[0];

            v.c["ro"][i] = roL;
            v.c["roUx"][i] = roL*velL;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i] = prsL/(gam-1.0) + 0.5*roL*velL*velL;
        }

        if (cfg.initial == "bump")
        {
            // bump
            flow_float gam  = 1.4;
            flow_float Pt = 120193.0;
            flow_float Tt = 302.557;
            flow_float M  = 0.5;

            flow_float Ps = Pt*pow((1.0+0.5*(gam-1.0)*M*M), -gam/(gam-1.0));
            flow_float Ts = Tt*pow((1.0+0.5*(gam-1.0)*M*M), -1.0);
            flow_float ro = gam*Ps/(cfg.cp*(gam-1.0)*Ts);
            flow_float a = sqrt(gam*Ps/ro);

            flow_float x = msh.cells[i].centCoords[0];

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = M*a*ro;
            v.c["roUy"][i] = 0.0*ro;
            v.c["roUz"][i] = 0.0*ro;
            v.c["roe"][i] =  pow(1.0+0.5*(gam-1.0)*M*M,-gam/(gam-1.0))*Pt/(gam-1.0) + 0.5*ro*pow((M*a),2.0);
        }

        if (cfg.initial == "Taylor-Green")
        {
            flow_float gam = 1.4;
            flow_float cp  = gam-1.0;
            flow_float R   = cp*(gam-1.0)/gam;

            flow_float L   = 1.0; // length

            flow_float ro0 = 1.0;
            flow_float P0  = 1.0/gam;
            //flow_float T0  = 1.0; // P0/(ro0*R);
            flow_float M0  = 0.4;
            flow_float c0  = 1.0;//sqrt(gam*P0/ro0);

            flow_float V0  = M0*c0;

            flow_float x = msh.cells[i].centCoords[0];
            flow_float y = msh.cells[i].centCoords[1];
            flow_float z = msh.cells[i].centCoords[2];

            flow_float u0 = V0*sin(x/L)*cos(y/L)*cos(z/L);
            flow_float v0 =-V0*cos(x/L)*sin(y/L)*cos(z/L);
            flow_float w0 = 0.0;

           flow_float P1 = P0 +ro0*V0*V0/16.0*(cos(2*x/L)+cos(2*y/L))*(cos(2*z/L)+2.0);

            //flow_float ro1 = P1/(R*T0);

            v.c["ro"][i]   = ro0;
            v.c["roUx"][i] = ro0*u0;
            v.c["roUy"][i] = -ro0*v0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P1/(gam-1.0) + 0.5*ro0*(u0*u0+v0*v0+w0*w0);
        }

        if (cfg.initial == "half_sphere")
        {
            flow_float cp  = cfg.cp;
            flow_float gam = cfg.gamma;
            flow_float R   = cp*(gam-1.0)/gam;

            flow_float ro0 = 0.02026;
            flow_float T0  = 63.73;
            flow_float P0  = ro0*R*T0;

            flow_float u0  = 1296.22;


            v.c["ro"][i]   = ro0;
            v.c["roUx"][i] = ro0*u0;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P0/(gam-1.0) + 0.5*ro0*(u0*u0);
        }

    }
};
