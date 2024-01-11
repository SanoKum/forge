#include "calcGradient_d.cuh"
#include "cuda_nagare/cudaWrapper.cuh"

__global__ void calcGradient_np_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* nplane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;

    __syncthreads();
    if (ip < nNormalPlanes) {
        geom_int  ic0 = nplane_cells[2*ip+0];
        geom_int  ic1 = nplane_cells[2*ip+1];

        geom_float f   = fx[ip];
        flow_float Uxf = f*Ux[ic0] + (1.0-f)*Ux[ic1];
        flow_float Uyf = f*Uy[ic0] + (1.0-f)*Uy[ic1];
        flow_float Uzf = f*Uz[ic0] + (1.0-f)*Uz[ic1];

        flow_float Pf  = f*P[ic0] + (1.0-f)*P[ic1];
        flow_float Tf  = f*T[ic0] + (1.0-f)*T[ic1];
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];

        atomicAdd(&dUxdx[ic0], sxx*Uxf);
        atomicAdd(&dUxdy[ic0], syy*Uxf);
        atomicAdd(&dUxdz[ic0], szz*Uxf);

        atomicAdd(&dUydx[ic0], sxx*Uyf);
        atomicAdd(&dUydy[ic0], syy*Uyf);
        atomicAdd(&dUydz[ic0], szz*Uyf);

        atomicAdd(&dUzdx[ic0], sxx*Uzf);
        atomicAdd(&dUzdy[ic0], syy*Uzf);
        atomicAdd(&dUzdz[ic0], szz*Uzf);

        atomicAdd(&dTdx[ic0], sxx*Tf);
        atomicAdd(&dTdy[ic0], syy*Tf);
        atomicAdd(&dTdz[ic0], szz*Tf);

        atomicAdd(&dPdx[ic0], sxx*Pf);
        atomicAdd(&dPdy[ic0], syy*Pf);
        atomicAdd(&dPdz[ic0], szz*Pf);


        atomicAdd(&dUxdx[ic1], -sxx*Uxf);
        atomicAdd(&dUxdy[ic1], -syy*Uxf);
        atomicAdd(&dUxdz[ic1], -szz*Uxf);

        atomicAdd(&dUydx[ic1], -sxx*Uyf);
        atomicAdd(&dUydy[ic1], -syy*Uyf);
        atomicAdd(&dUydz[ic1], -szz*Uyf);

        atomicAdd(&dUzdx[ic1], -sxx*Uzf);
        atomicAdd(&dUzdy[ic1], -syy*Uzf);
        atomicAdd(&dUzdz[ic1], -szz*Uzf);

        atomicAdd(&dTdx[ic1], -sxx*Tf);
        atomicAdd(&dTdy[ic1], -syy*Tf);
        atomicAdd(&dTdz[ic1], -szz*Tf);

        atomicAdd(&dPdx[ic1], -sxx*Pf);
        atomicAdd(&dPdy[ic1], -syy*Pf);
        atomicAdd(&dPdz[ic1], -szz*Pf);
    }
    __syncthreads();
}


__global__ void calcGradient_bp_d
( 
 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz
)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip  = bplane_plane[ib];
        geom_int  ic0 = bplane_cell[ib];

        flow_float Uxf = Uxp[ip];
        flow_float Uyf = Uyp[ip];
        flow_float Uzf = Uzp[ip];

        //TEMP flow_float Pf  = Pp[ip];
        flow_float Pf  = P[ic0];
        flow_float Tf  = Tp[ip];

        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];

        __syncthreads();
        atomicAdd(&dUxdx[ic0], sxx*Uxf);
        atomicAdd(&dUxdy[ic0], syy*Uxf);
        atomicAdd(&dUxdz[ic0], szz*Uxf);

        atomicAdd(&dUydx[ic0], sxx*Uyf);
        atomicAdd(&dUydy[ic0], syy*Uyf);
        atomicAdd(&dUydz[ic0], szz*Uyf);

        atomicAdd(&dUzdx[ic0], sxx*Uzf);
        atomicAdd(&dUzdy[ic0], syy*Uzf);
        atomicAdd(&dUzdz[ic0], szz*Uzf);

        atomicAdd(&dTdx[ic0], sxx*Tf);
        atomicAdd(&dTdy[ic0], syy*Tf);
        atomicAdd(&dTdz[ic0], szz*Tf);

        atomicAdd(&dPdx[ic0], sxx*Pf);
        atomicAdd(&dPdy[ic0], syy*Pf);
        atomicAdd(&dPdz[ic0], szz*Pf);
        __syncthreads();
    }
}

__global__ void calcGradient_vol_d
( 
 // mesh structure
 geom_int nCells,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells ) {

        geom_float volume = vol[ic];

        dUxdx[ic] = dUxdx[ic]/volume;
        dUxdy[ic] = dUxdy[ic]/volume;
        dUxdz[ic] = dUxdz[ic]/volume;
                                    
        dUydx[ic] = dUydx[ic]/volume;
        dUydy[ic] = dUydy[ic]/volume;
        dUydz[ic] = dUydz[ic]/volume;

        dUzdx[ic] = dUzdx[ic]/volume;
        dUzdy[ic] = dUzdy[ic]/volume;
        dUzdz[ic] = dUzdz[ic]/volume;

        dPdx[ic] = dPdx[ic]/volume;
        dPdy[ic] = dPdy[ic]/volume;
        dPdz[ic] = dPdz[ic]/volume;

        dTdx[ic] = dTdx[ic]/volume;
        dTdy[ic] = dTdy[ic]/volume;
        dTdz[ic] = dTdz[ic]/volume;
    }
}

void calcGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    // initialize
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUxdx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUxdy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUxdz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUydx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUydy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUydz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUzdx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUzdy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUzdz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdz"] , 0.0, msh.nCells*sizeof(flow_float)));

    // sum over normal planes
    calcGradient_np_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_nplane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        var.p_d["T"]  , var.c_d["T"] ,
        var.p_d["Ux"] , var.c_d["Ux"] ,
        var.p_d["Uy"] , var.c_d["Uy"] ,
        var.p_d["Uz"] , var.c_d["Uz"] ,
        var.p_d["ro"] , var.c_d["ro"] ,
        var.p_d["P"]  , var.c_d["P"]  , 

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"]
    ) ;

    // sum over boundary plane
     
    // --------------------------------
    // *** sum over boudnary planes ***
    // --------------------------------
    // boundary conditions
    for (auto& bc : msh.bconds)
    {
        calcGradient_bp_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
            // mesh structure
            bc.iPlanes.size() , 
            bc.map_bplane_plane_d,  bc.map_bplane_cell_d,
            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
            var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

            // basic variables
            var.p_d["T"]  , var.c_d["T"] ,
            var.p_d["Ux"] , var.c_d["Ux"] ,
            var.p_d["Uy"] , var.c_d["Uy"] ,
            var.p_d["Uz"] , var.c_d["Uz"] ,
            var.p_d["ro"] , var.c_d["ro"] ,
            var.p_d["P"]  , var.c_d["P"]  , 

            // gradient
            var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
            var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
            var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
            var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
            var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"]
        ) ;
    }

    // 
    calcGradient_vol_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        var.p_d["T"]  , var.c_d["T"] ,
        var.p_d["Ux"] , var.c_d["Ux"] ,
        var.p_d["Uy"] , var.c_d["Uy"] ,
        var.p_d["Uz"] , var.c_d["Uz"] ,
        var.p_d["ro"] , var.c_d["ro"] ,
        var.p_d["P"]  , var.c_d["P"]  , 

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"]
    ) ;


}