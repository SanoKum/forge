#include "cuda_forge/calcStructualVariables_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "iostream"

__global__ 
void calcStructualVariables_d 
( 
 geom_int nPlanes, geom_int nNormalPlanes,
 geom_int* plane_cells ,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz ,  geom_float* ss,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz ,  
 geom_float* ccx ,  geom_float* ccy ,  geom_float* ccz ,  
 geom_float* fx  ,  geom_float* dcc 
)
{
    geom_int  ip  = blockDim.x*blockIdx.x + threadIdx.x;

    //if (ip < nNormalPlanes) {
    if (ip < nPlanes) { // include ghost cell
        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float ccx0 = ccx[ic0];
        geom_float ccy0 = ccy[ic0];
        geom_float ccz0 = ccz[ic0];

        geom_float ccx1 = ccx[ic1];
        geom_float ccy1 = ccy[ic1];
        geom_float ccz1 = ccz[ic1];

        geom_float dcx = ccx1 - ccx0;
        geom_float dcy = ccy1 - ccy0;
        geom_float dcz = ccz1 - ccz0;
        geom_float dc  = sqrtf( powf(dcx, 2.0) + pow(dcy, 2.0) + pow(dcz, 2.0));

        geom_float dc0px = sx[ip]*(pcx[ip] - ccx0)/ss[ip];
        geom_float dc0py = sy[ip]*(pcy[ip] - ccy0)/ss[ip];
        geom_float dc0pz = sz[ip]*(pcz[ip] - ccz0)/ss[ip];
        geom_float dc0p  = sqrtf( powf(dc0px, 2.0) + pow(dc0py, 2.0) + pow(dc0pz, 2.0));

        geom_float dc1px = sx[ip]*(pcx[ip] - ccx1)/ss[ip];
        geom_float dc1py = sy[ip]*(pcy[ip] - ccy1)/ss[ip];
        geom_float dc1pz = sz[ip]*(pcz[ip] - ccz1)/ss[ip];
        geom_float dc1p  = sqrtf( powf(dc1px, 2.0) + pow(dc1py, 2.0) + pow(dc1pz, 2.0));

        fx [ip] = dc1p/(dc0p + dc1p);
        dcc[ip] = dc;
    }
};

void calcStructualVariables_d_wrapper(cudaConfig& cuda_cfg , mesh& msh,  variables& v)
{
    calcStructualVariables_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>>(
        msh.nPlanes, msh.nNormalPlanes,
        msh.map_plane_cells_d,
        v.p_d["sx"] , v.p_d["sy"] , v.p_d["sz"], v.p_d["ss"],
        v.p_d["pcx"], v.p_d["pcy"], v.p_d["pcz"],
        v.c_d["ccx"], v.c_d["ccy"], v.c_d["ccz"],
        v.p_d["fx"] , v.p_d["dcc"]
    );

    //for (auto& bc : msh.bconds)
    //{
    //    calcStructualVariables_bp_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
    //        bc.iPlanes.size() ,
    //        bc.map_bplane_plane_d,  bc.map_bplane_cell_d,
    //        v.p_d["sx"] , v.p_d["sy"] , v.p_d["sz"], v.p_d["ss"],
    //        v.p_d["pcx"], v.p_d["pcy"], v.p_d["pcz"],
    //        v.c_d["ccx"], v.c_d["ccy"], v.c_d["ccz"],
    //        v.p_d["fx"] , v.p_d["dcc"]
    //    );
    //}

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

};

