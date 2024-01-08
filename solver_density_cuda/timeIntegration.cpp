#include "timeIntegration.hpp"

using namespace std;

void timeIntegration(int loop , solverConfig &cfg , mesh &msh , variables &v , matrix& mat_ns)
{
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

        if (cfg.isImplicit == 0) { // explicit
            // N: previous outer step , M: previous inner loop
            ro[ic]   = cfg.coef_N[loop]*roN[ic]   + cfg.coef_M[loop]*roM[ic]   + cfg.coef_Res[loop]*res_ro[ic]*cfg.dt/vol;
            roUx[ic] = cfg.coef_N[loop]*roUxN[ic] + cfg.coef_M[loop]*roUxM[ic] + cfg.coef_Res[loop]*res_roUx[ic]*cfg.dt/vol;
            roUy[ic] = cfg.coef_N[loop]*roUyN[ic] + cfg.coef_M[loop]*roUyM[ic] + cfg.coef_Res[loop]*res_roUy[ic]*cfg.dt/vol;
            roUz[ic] = cfg.coef_N[loop]*roUzN[ic] + cfg.coef_M[loop]*roUzM[ic] + cfg.coef_Res[loop]*res_roUz[ic]*cfg.dt/vol;
            roe[ic]  = cfg.coef_N[loop]*roeN[ic]  + cfg.coef_M[loop]*roeM[ic]  + cfg.coef_Res[loop]*res_roe[ic]*cfg.dt/vol;

        } else if (cfg.timeIntegration == 10) { // implicit (dual-time stepping & explicit scheme)
            ro[ic]   = cfg.coef_N[loop]*roN[ic]   + cfg.coef_Res[loop]*res_ro_dual[ic]*cfg.dt/vol;        
            roUx[ic] = cfg.coef_N[loop]*roUxN[ic] + cfg.coef_Res[loop]*res_roUx_dual[ic]*cfg.dt/vol;        
            roUy[ic] = cfg.coef_N[loop]*roUyN[ic] + cfg.coef_Res[loop]*res_roUy_dual[ic]*cfg.dt/vol;    
            roUz[ic] = cfg.coef_N[loop]*roUzN[ic] + cfg.coef_Res[loop]*res_roUz_dual[ic]*cfg.dt/vol;
            roe[ic]  = cfg.coef_N[loop]*roeN[ic]  + cfg.coef_Res[loop]*res_roe_dual[ic]*cfg.dt/vol;
        }
    }
}
