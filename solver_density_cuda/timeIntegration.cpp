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

    vector<flow_float>& res_ro_m   = v.c["res_ro_m"];
    vector<flow_float>& res_roUx_m = v.c["res_roUx_m"];
    vector<flow_float>& res_roUy_m = v.c["res_roUy_m"];
    vector<flow_float>& res_roUz_m = v.c["res_roUz_m"];
    vector<flow_float>& res_roe_m  = v.c["res_roe_m"];

    // time integration
    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        geom_float vol = msh.cells[ic].volume;

        if (cfg.timeIntegration == 4) { // 4th order
            if (loop==0) { // first
                res_ro_m[ic]   = 0.0;
                res_roUx_m[ic] = 0.0;
                res_roUy_m[ic] = 0.0;
                res_roUz_m[ic] = 0.0;
                res_roe_m[ic]  = 0.0;
            }

            res_ro_m[ic]   += cfg.coef_Res_4thRunge[loop]*res_ro[ic]*cfg.dt/vol;
            res_roUx_m[ic] += cfg.coef_Res_4thRunge[loop]*res_roUx[ic]*cfg.dt/vol;
            res_roUy_m[ic] += cfg.coef_Res_4thRunge[loop]*res_roUy[ic]*cfg.dt/vol;
            res_roUz_m[ic] += cfg.coef_Res_4thRunge[loop]*res_roUz[ic]*cfg.dt/vol;
            res_roe_m[ic]  += cfg.coef_Res_4thRunge[loop]*res_roe[ic]*cfg.dt/vol;


            if (loop<3) { 
                ro[ic]   = roN[ic]   + cfg.coef_DT_4thRunge[loop]*res_ro[ic]*cfg.dt/vol;
                roUx[ic] = roUxN[ic] + cfg.coef_DT_4thRunge[loop]*res_roUx[ic]*cfg.dt/vol;
                roUy[ic] = roUyN[ic] + cfg.coef_DT_4thRunge[loop]*res_roUy[ic]*cfg.dt/vol;
                roUz[ic] = roUzN[ic] + cfg.coef_DT_4thRunge[loop]*res_roUz[ic]*cfg.dt/vol;
                roe[ic]  = roeN[ic]  + cfg.coef_DT_4thRunge[loop]*res_roe[ic]*cfg.dt/vol;
            } else {
                res_ro[ic]   = res_ro_m[ic]*vol/cfg.dt ;
                res_roUx[ic] = res_roUx_m[ic]*vol/cfg.dt ;
                res_roUy[ic] = res_roUy_m[ic]*vol/cfg.dt ;
                res_roUz[ic] = res_roUz_m[ic]*vol/cfg.dt ;
                res_roe[ic]  = res_roe_m[ic]*vol/cfg.dt ;

                ro[ic]   = roN[ic]  + res_ro_m[ic] ;
                roUx[ic] = roUxN[ic]+ res_roUx_m[ic] ;
                roUy[ic] = roUyN[ic]+ res_roUy_m[ic] ;
                roUz[ic] = roUzN[ic]+ res_roUz_m[ic] ;
                roe[ic]  = roeN[ic] + res_roe_m[ic] ;

            }

        } else if (cfg.isImplicit == 0) { // explicit
            // N: previous outer step , M: previous inner loop
            ro[ic]   = cfg.coef_N[loop]*roN[ic]   + cfg.coef_M[loop]*roM[ic]   + cfg.coef_Res[loop]*res_ro[ic]*cfg.dt/vol;
            roUx[ic] = cfg.coef_N[loop]*roUxN[ic] + cfg.coef_M[loop]*roUxM[ic] + cfg.coef_Res[loop]*res_roUx[ic]*cfg.dt/vol;
            roUy[ic] = cfg.coef_N[loop]*roUyN[ic] + cfg.coef_M[loop]*roUyM[ic] + cfg.coef_Res[loop]*res_roUy[ic]*cfg.dt/vol;
            roUz[ic] = cfg.coef_N[loop]*roUzN[ic] + cfg.coef_M[loop]*roUzM[ic] + cfg.coef_Res[loop]*res_roUz[ic]*cfg.dt/vol;
            roe[ic]  = cfg.coef_N[loop]*roeN[ic]  + cfg.coef_M[loop]*roeM[ic]  + cfg.coef_Res[loop]*res_roe[ic]*cfg.dt/vol;

        } else if (cfg.timeIntegration == 10) { // implicit (m-time stepping & explicit scheme)
            ro[ic]   = cfg.coef_N[loop]*roN[ic]   + cfg.coef_Res[loop]*res_ro_m[ic]*cfg.dt/vol;        
            roUx[ic] = cfg.coef_N[loop]*roUxN[ic] + cfg.coef_Res[loop]*res_roUx_m[ic]*cfg.dt/vol;        
            roUy[ic] = cfg.coef_N[loop]*roUyN[ic] + cfg.coef_Res[loop]*res_roUy_m[ic]*cfg.dt/vol;    
            roUz[ic] = cfg.coef_N[loop]*roUzN[ic] + cfg.coef_Res[loop]*res_roUz_m[ic]*cfg.dt/vol;
            roe[ic]  = cfg.coef_N[loop]*roeN[ic]  + cfg.coef_Res[loop]*res_roe_m[ic]*cfg.dt/vol;
        }
    }
}
