#include "calcGradient_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"


__global__ void ducrosSensor_d
( 
 // mesh structure
 geom_int nCells,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* Ht  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* ducros
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    //https://doi.org/10.2514/6.2018-3710
    if (ic < nCells ) {
        flow_float Ux_0 = Ux[ic];
        flow_float Uy_0 = Uy[ic];
        flow_float Uz_0 = Uz[ic];

        flow_float dUxdx_0 = dUxdx[ic];
        flow_float dUxdy_0 = dUxdy[ic];
        flow_float dUxdz_0 = dUxdz[ic];
        flow_float dUydx_0 = dUydx[ic];
        flow_float dUydy_0 = dUydy[ic];
        flow_float dUydz_0 = dUydz[ic];
        flow_float dUzdx_0 = dUzdx[ic];
        flow_float dUzdy_0 = dUzdy[ic];
        flow_float dUzdz_0 = dUzdz[ic];

        geom_float volume = vol[ic];

        flow_float divu = dUxdx_0 + dUydy_0 + dUzdz_0;
        flow_float divu2 = divu*divu;
        flow_float Umag = sqrt(pow(Ux_0,2.0) + pow(Uy_0,2.0) + pow(Uz_0,2.0));
        flow_float eps = 1e-12;
        flow_float nu = 0.1;
        flow_float omega2 = pow(nu*Umag/pow(volume, 1.0/3.0),2.0);

        flow_float theta = divu2/(divu2 + omega2 + eps);

        flow_float vort = pow(dUzdy[ic]-dUydz[ic],2.0)
                        + pow(dUxdz[ic]-dUzdx[ic],2.0)
                        + pow(dUydx[ic]-dUxdy[ic],2.0);

        flow_float D = min((4.0/3.0)*divu2/(divu2+vort+eps), 1.0);

        flow_float fai = 0.01;

        ducros[ic] = (max((D-fai),0.0))*theta+fai;
    }
}


void ducrosSensor_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    ducrosSensor_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        var.c_d["ro"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"]  , 
        var.c_d["T"]  ,
        var.c_d["Ht"] ,

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["ducros"]  
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}