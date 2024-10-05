#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "solverConfig.hpp"
#include "calcWallDistance_kdtree.hpp"

#include <cmath>

void setInitial(solverConfig& cfg , mesh& msh , variables& v)
{
    // *** set initial value ***
    for (geom_int i = 0 ; i<msh.nCells ; i++)
    {
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
        } else if (cfg.initial == "mach3") {
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

        } else if (cfg.initial == "bump") {
            // bump
            flow_float gam  = 1.4;
            flow_float Pt = 100000.0;
            flow_float Tt = 293.15 ;
            flow_float M  = 0.675;
            //flow_float M  = 1.65;

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


        } else if (cfg.initial == "bump_M1.65") {
            // bump
            flow_float gam  = 1.4;
            flow_float Pt = 100000.0;
            flow_float Tt = 293.15 ;
            //flow_float M  = 0.675;
            flow_float M  = 1.65;

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

        } else if (cfg.initial == "Taylor-Green" or cfg.initial == "Taylor-Green_M0.1") {
            flow_float gam = 1.4;
            flow_float cp  = gam-1.0;
            flow_float R   = cp*(gam-1.0)/gam;

            flow_float L   = 1.0; // length

            flow_float ro0 = 1.0;
            flow_float P0  = 1.0/gam;
            flow_float M0  = 0.4;

            if (cfg.initial == "Taylor-Green_M0.1") M0 = 0.1;

            flow_float x = msh.cells[i].centCoords[0];
            flow_float y = msh.cells[i].centCoords[1];
            flow_float z = msh.cells[i].centCoords[2];

            flow_float u0 = M0*sin(x/L)*cos(y/L)*cos(z/L);
            flow_float v0 =-M0*cos(x/L)*sin(y/L)*cos(z/L);
            flow_float w0 = 0.0;

            flow_float P1 = P0 +ro0*M0*M0/16.0*(cos(2*x/L)+cos(2*y/L))*(cos(2*z/L)+2.0);

            //flow_float ro1 = P1/(R*T0);

            v.c["ro"][i]   = ro0;
            v.c["roUx"][i] = ro0*u0;
            v.c["roUy"][i] = ro0*v0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P1/(gam-1.0) + 0.5*ro0*(u0*u0+v0*v0+w0*w0);
        } else if (cfg.initial == "half_sphere") {
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

        } else if (cfg.initial == "supercritical") {
            flow_float cp  = cfg.cp;
            flow_float gam = cfg.gamma;
            flow_float R   = cp*(gam-1.0)/gam;

            flow_float M = 0.3;
            flow_float P = 101325.0;
            flow_float T = 288.15;

            flow_float ro= P/(R*T);
            flow_float c = sqrt(gam*R*T);
            flow_float u = c*M;

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);

        } else if (cfg.initial == "laval") {
            flow_float cp  = cfg.cp;
            flow_float gam = cfg.gamma;
            flow_float R   = cp*(gam-1.0)/gam;

            flow_float M = 0.80;
            flow_float P0 = 59100.0;
            flow_float T0 = 286.65;

            flow_float Ps = P0/(pow(1.0+0.5*(gam-1.0)*M*M, gam/(gam-1.0)));
            flow_float Ts = T0/(1.0+0.5*(gam-1.0)*M*M);

            flow_float ro= Ps/(R*Ts);
            flow_float c = sqrt(gam*R*Ts);
            flow_float u = c*M;

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = Ps/(gam-1.0) + 0.5*ro*(u*u);

        } else if (cfg.initial == "hifire") {
            flow_float cp  = cfg.cp;
            flow_float gam = cfg.gamma;
            flow_float R   = cp*(gam-1.0)/gam;

            flow_float M = 7.16;
            flow_float P = 4620.0;
            flow_float T = 231.7;

            flow_float ro= P/(R*T);
            flow_float c = sqrt(gam*R*T);
            flow_float u = c*M;

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);

        } else if (cfg.initial == "flare") {
            flow_float cp  = cfg.cp;
            flow_float gam = cfg.gamma;
            flow_float R   = cp*(gam-1.0)/gam;

            flow_float M = 6.0;
            //flow_float P = 930.9*1000;
            //flow_float T = 433.0;
            flow_float T = 52.805;

            flow_float ro= 0.03888;
            flow_float P= ro*R*T;
            flow_float c = sqrt(gam*R*T);
            flow_float u = c*M;

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);

        } else if (cfg.initial == "poiseuille") {
            flow_float cp  = cfg.cp;
            flow_float gam = cfg.gamma;
            flow_float R   = cp*(gam-1.0)/gam;

            flow_float P = 100000.0;
            flow_float T = 288.15;

            flow_float ro= P/(R*T);
            flow_float u = 10;

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);



        } else {
            std::cout << "Error: Unknown initial" << std::endl;
            exit;
        }
    }

    std::list<std::string> names = {"ro", "roUx", "roUy", "roUz", "roe"};
    v.copyVariables_cell_H2D(names);

    calcWallDistance_kdtree(cfg , msh , v);
};


