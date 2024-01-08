#include "cuda_nagare/calcStructualVariables_d.cuh"

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

    if (ip < nNormalPlanes) {
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

        geom_float dc2px = pcx[ip] - ccx1;
        geom_float dc2py = pcy[ip] - ccy1;
        geom_float dc2pz = pcz[ip] - ccz1;
        geom_float dc2p  = sqrtf( powf(dc2px, 2.0) + pow(dc2py, 2.0) + pow(dc2pz, 2.0));

        fx [ip] = dc2p/dc;
        dcc[ip] = dc;
        //printf("gpu ip=%d, fx=%f, dcc=%f\n", ip, fx[ip], dcc[ip]);
        __syncthreads();

    }
};

__global__ 
void calcStructualVariables_bp_d 
( 
 // mesh structure
 geom_int nb, 
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  

 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz ,  geom_float* ss,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz ,  
 geom_float* ccx ,  geom_float* ccy ,  geom_float* ccz ,  
 geom_float* fx  ,  geom_float* dcc 
)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip  = bplane_plane[ib];
        geom_int  ic0 = bplane_cell[ib];

        geom_float ccx0 = ccx[ic0];
        geom_float ccy0 = ccy[ic0];
        geom_float ccz0 = ccz[ic0];

        geom_float pcx0 = pcx[ip];
        geom_float pcy0 = pcy[ip];
        geom_float pcz0 = pcz[ip];

        geom_float sx0  = sx[ip];
        geom_float sy0  = sy[ip];
        geom_float sz0  = sz[ip];
        geom_float ss0  = sz[ip];

        geom_float dn   = (  (pcx0 - ccx0)*sx0 
                           + (pcy0 - ccy0)*sy0 
                           + (pcz0 - ccz0)*sz0 )/ss0;

        fx [ip] = 0.0;
        dcc[ip] = dn;
        __syncthreads();
    }
};

void calcStructualVariables_d_wrapper(cudaConfig& cuda_cfg , mesh& msh,  variables& v)
{
    calcStructualVariables_d<<<cuda_cfg.dimGrid_nplane , cuda_cfg.dimBlock>>>(
        msh.nPlanes, msh.nNormalPlanes,
        msh.map_nplane_cells_d,
        v.p_d["sx"] , v.p_d["sy"] , v.p_d["sz"], v.p_d["ss"],
        v.p_d["pcx"], v.p_d["pcy"], v.p_d["pcz"],
        v.c_d["ccx"], v.c_d["ccy"], v.c_d["ccz"],
        v.p_d["fx"] , v.p_d["dcc"]
    );

    for (auto& bc : msh.bconds)
    {
        calcStructualVariables_bp_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
            bc.iPlanes.size() ,
            bc.map_bplane_plane_d,  bc.map_bplane_cell_d,
            v.p_d["sx"] , v.p_d["sy"] , v.p_d["sz"], v.p_d["ss"],
            v.p_d["pcx"], v.p_d["pcy"], v.p_d["pcz"],
            v.c_d["ccx"], v.c_d["ccy"], v.c_d["ccz"],
            v.p_d["fx"] , v.p_d["dcc"]
        );
    }

};

