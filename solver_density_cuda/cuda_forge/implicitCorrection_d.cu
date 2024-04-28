#include "implicitCorrection_d.cuh"

__global__ void dualtime_explicit_d
// see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
// N: previous outer step , M: previous inner loop

( 
 geom_int dt ,

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol ,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe,

 flow_float* res_ro_dual ,
 flow_float* res_roUx_dual ,
 flow_float* res_roUy_dual ,
 flow_float* res_roUz_dual ,
 flow_float* res_roe_dual 

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    geom_float v = vol[ic];

    if (ic < nCells) {
        // N: previous outer step , M: previous inner loop
        res_ro_dual[ic]   = -(ro[ic]-roN[ic])*v/dt     + res_ro[ic];
        res_roUx_dual[ic] = -(roUx[ic]-roUxN[ic])*v/dt + res_roUx[ic];
        res_roUy_dual[ic] = -(roUy[ic]-roUyN[ic])*v/dt + res_roUy[ic];
        res_roUz_dual[ic] = -(roUz[ic]-roUzN[ic])*v/dt + res_roUz[ic];
        res_roe_dual[ic]  = -(roe[ic]-roeN[ic])*v/dt   + res_roe[ic];
    }
    __syncthreads();
}


void implicitCorrection_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isImplicit != 1) return;

    if (cfg.timeIntegration == 10) {
        dualtime_explicit_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            cfg.dt , 

            // mesh structure
            msh.nCells_all , msh.nCells ,
            var.c_d["volume"],

            // basic variables
            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] ,
            var.c_d["res_ro_m"], var.c_d["res_roUx_m"], var.c_d["res_roUy_m"], var.c_d["res_roUz_m"] , var.c_d["res_roe_m"] 
        ) ;


    }
}