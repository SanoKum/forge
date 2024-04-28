#include "setDT.hpp"

void setDT(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& v)
{
    if (cfg.gpu==1) {
        setDT_d_wrapper(cfg , cuda_cfg , msh , v);
        return;
    }
    std::vector<flow_float>& cfl        = v.c["cfl"];
    std::vector<flow_float>& cfl_pseudo = v.c["cfl_pseudo"];
    std::vector<flow_float>& sonic = v.c["sonic"];
    std::vector<flow_float>& Ux    = v.c["Ux"];
    std::vector<flow_float>& Uy    = v.c["Uy"];
    std::vector<flow_float>& Uz    = v.c["Uz"];
    std::vector<flow_float>& fxp   = v.p["fx"];

    for (geom_int ic=0 ; ic<msh.nCells_all ; ic++)
    {
        cfl[ic]        = -1e+30;
        cfl_pseudo[ic] = -1e+30;
    }


    // 界面の流速でCFLを見る
    for (geom_int ip=0 ; ip<msh.nPlanes ; ip++)
    {
        geom_int ic0 = msh.planes[ip].iCells[0];
        geom_int ic1 = msh.planes[ip].iCells[1];


        geom_float vol0 = msh.cells[ic0].volume;
        geom_float vol1 = msh.cells[ic1].volume;

        flow_float sx = msh.planes[ip].surfVect[0];
        flow_float sy = msh.planes[ip].surfVect[1];
        flow_float sz = msh.planes[ip].surfVect[2];

        geom_float surfArea = msh.planes[ip].surfArea;

        geom_float dx0 = vol0/surfArea;
        geom_float dx1 = vol0/surfArea;

        flow_float Ux0 = Ux[ic0];
        flow_float Uy0 = Uy[ic0];
        flow_float Uz0 = Uz[ic0];

        flow_float Ux1 = Ux[ic1];
        flow_float Uy1 = Uy[ic1];
        flow_float Uz1 = Uz[ic1];

        flow_float US  = (fxp[ip]*Ux0 + (1.0-fxp[ip])*Ux1)*sx
                        +(fxp[ip]*Uy0 + (1.0-fxp[ip])*Uy1)*sy
                        +(fxp[ip]*Uz0 + (1.0-fxp[ip])*Uz1)*sz;

        cfl[ic0] = std::max(cfl[ic0] , cfg.dt*(abs(US)/surfArea+sonic[ic0])/dx0);
        cfl[ic1] = std::max(cfl[ic1] , cfg.dt*(abs(US)/surfArea+sonic[ic1])/dx1);

        cfl_pseudo[ic0] = std::max(cfl_pseudo[ic0] , cfg.dt_pseudo*(abs(US)/surfArea+sonic[ic0])/dx0);
        cfl_pseudo[ic1] = std::max(cfl_pseudo[ic1] , cfg.dt_pseudo*(abs(US)/surfArea+sonic[ic1])/dx1);
    }

    // get max cfl and change DT
    flow_float cfl_max        = *max_element(begin(cfl), end(cfl));
    flow_float cfl_pseudo_max = *max_element(begin(cfl_pseudo), end(cfl_pseudo));
    cfg.totalTime = cfg.totalTime + cfg.dt;

    std::cout << "  DT              : " << cfg.dt  << std::endl;
    std::cout << "  Max CFL         : " << cfl_max << std::endl;
    if (cfg.isImplicit == 1) {
        std::cout << "  Max CFL(pseudo) : " << cfl_pseudo_max << std::endl;
    }
    std::cout << "  Total Time      : " << cfg.totalTime << std::endl;

    if (cfg.dtControl == 1)
    {
        flow_float dt_new        = cfg.dt*cfg.cfl/cfl_max;
        flow_float dt_pseudo_new = cfg.dt*cfg.cfl_pseudo/cfl_max;
        cfg.dt        = std::max(std::min(dt_new       , cfg.dt_max), cfg.dt_min);
        cfg.dt_pseudo = std::max(std::min(dt_pseudo_new, cfg.dt_max), cfg.dt_min);
    }
};