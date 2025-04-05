#include "gasProperties_d.cuh"

__global__ void gasProperties_d
( 
 // gas properties
 int thermalMethod , int viscMethod , 

 // gas properties
 flow_float gamma , flow_float cp , flow_float visc_lam,

 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* T   ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 flow_float* vis_lam_array,
 flow_float* gam_array,
 flow_float* cp_array

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    flow_float ek;
    flow_float intE;

    if (ic < nCells_all) {

        if (thermalMethod == 0) {
            gam_array[ic] = gamma;
            cp_array[ic] = cp;
        }

        if (viscMethod == 0) {
            vis_lam_array[ic] = visc_lam;

        } else if (viscMethod == 1) { // sutherland
            flow_float T0  = 273.0;
            flow_float mu0 = 1.716e-5;
            flow_float Smu = 111.0;
            vis_lam_array[ic] = mu0*pow(T[ic]/T0,3.0/2.0)*(T0+Smu)/(T[ic]+Smu);
        }
    }
}


void gasProperties_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    gasProperties_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        cfg.thermalMethod,  
        cfg.viscMethod ,  

        // gas properties
        cfg.gamma , cfg.cp , cfg.visc,

        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"], var.c_d["roUy"] , var.c_d["roUz"], var.c_d["roe"] ,
        var.c_d["P"]   , var.c_d["Ht"]  , var.c_d["sonic"], var.c_d["T"], 
        var.c_d["Ux"]  , var.c_d["Uy"]  , var.c_d["Uz"] ,

        var.c_d["vis_lam"] , var.c_d["gamma"] , var.c_d["cp"] 
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}