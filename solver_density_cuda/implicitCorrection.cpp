#include "implicitCorrection.hpp"

using namespace std;

void implicitCorrection(int loop , solverConfig &cfg , mesh &msh , variables &v , matrix& mat_ns)
{
    if (cfg.isImplicit != 1) return; // function only for Implicit scheme

    vector<flow_float>& ro = v.c["ro"];
    vector<flow_float>& roUx = v.c["roUx"];
    vector<flow_float>& roUy = v.c["roUy"];
    vector<flow_float>& roUz = v.c["roUz"];
    vector<flow_float>& roe  = v.c["roe"];

    vector<flow_float>& roN = v.c["roN"];
    vector<flow_float>& roUxN = v.c["roUxN"];
    vector<flow_float>& roUyN = v.c["roUyN"];
    vector<flow_float>& roUzN = v.c["roUzN"];
    vector<flow_float>& roeN = v.c["roeN"];

    vector<flow_float>& roM = v.c["roM"];
    vector<flow_float>& roUxM = v.c["roUxM"];
    vector<flow_float>& roUyM = v.c["roUyM"];
    vector<flow_float>& roUzM = v.c["roUzM"];
    vector<flow_float>& roeM = v.c["roeM"];

    vector<flow_float>& res_ro   = v.c["res_ro"];
    vector<flow_float>& res_roUx = v.c["res_roUx"];
    vector<flow_float>& res_roUy = v.c["res_roUy"];
    vector<flow_float>& res_roUz = v.c["res_roUz"];
    vector<flow_float>& res_roe  = v.c["res_roe"];

    vector<flow_float>& res_ro_dual   = v.c["res_ro_dual"];
    vector<flow_float>& res_roUx_dual = v.c["res_roUx_dual"];
    vector<flow_float>& res_roUy_dual = v.c["res_roUy_dual"];
    vector<flow_float>& res_roUz_dual = v.c["res_roUz_dual"];
    vector<flow_float>& res_roe_dual  = v.c["res_roe_dual"];

    // time integration
    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        geom_float vol = msh.cells[ic].volume;

        if (cfg.timeIntegration == 10) {// dual-time stepping & explicit scheme
            // see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
            // N: previous outer step , M: previous inner loop
            res_ro_dual[ic]   = -(ro[ic]-roN[ic])*vol/cfg.dt     + res_ro[ic];
            res_roUx_dual[ic] = -(roUx[ic]-roUxN[ic])*vol/cfg.dt + res_roUx[ic];
            res_roUy_dual[ic] = -(roUy[ic]-roUyN[ic])*vol/cfg.dt + res_roUy[ic];
            res_roUz_dual[ic] = -(roUz[ic]-roUzN[ic])*vol/cfg.dt + res_roUz[ic];
            res_roe_dual[ic]  = -(roe[ic]-roeN[ic])*vol/cfg.dt   + res_roe[ic];
        }
    }

    // get max cfl and change DT
    flow_float res_ro_dual_max   = *max_element(begin(res_ro_dual)  , end(res_ro_dual));
    flow_float res_roUx_dual_max = *max_element(begin(res_roUx_dual), end(res_roUx_dual));
    flow_float res_roUy_dual_max = *max_element(begin(res_roUy_dual), end(res_roUy_dual));
    flow_float res_roUz_dual_max = *max_element(begin(res_roUz_dual), end(res_roUz_dual));
    flow_float res_roe_dual_max  = *max_element(begin(res_roe_dual) , end(res_roe_dual));

    std::cout << "    Residual mass max       : " << res_ro_dual_max << endl;
    std::cout << "    Residual momentum x max : " << res_roUx_dual_max << endl;
    std::cout << "    Residual momentum y max : " << res_roUy_dual_max << endl;
    std::cout << "    Residual momentum z max : " << res_roUz_dual_max << endl;
    std::cout << "    Residual Total Energy   : " << res_roe_dual_max << endl;
}