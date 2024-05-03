#include "calcGradient_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"

__global__ void calcGradient_1_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
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
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,
 flow_float* divU
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;

    __syncthreads();

    if (ip < nPlanes) {
        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float f   = fx[ip];
        flow_float Uxf = f*Ux[ic0] + (1.0-f)*Ux[ic1];
        flow_float Uyf = f*Uy[ic0] + (1.0-f)*Uy[ic1];
        flow_float Uzf = f*Uz[ic0] + (1.0-f)*Uz[ic1];

        flow_float Pf  = f*P[ic0] + (1.0-f)*P[ic1];
        flow_float Tf  = f*T[ic0] + (1.0-f)*T[ic1];
        flow_float rof = f*ro[ic0]+ (1.0-f)*ro[ic1];
        flow_float Htf = f*Ht[ic0]+ (1.0-f)*Ht[ic1];
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];

        flow_float US  = Uxf*sxx +Uyf*syy +Uzf*szz;

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

        atomicAdd(&drodx[ic0], sxx*rof);
        atomicAdd(&drody[ic0], syy*rof);
        atomicAdd(&drodz[ic0], szz*rof);

        atomicAdd(&dHtdx[ic0], sxx*Htf);
        atomicAdd(&dHtdy[ic0], syy*Htf);
        atomicAdd(&dHtdz[ic0], szz*Htf);

        atomicAdd(&divU[ic0], US);


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

        atomicAdd(&drodx[ic1], -sxx*rof);
        atomicAdd(&drody[ic1], -syy*rof);
        atomicAdd(&drodz[ic1], -szz*rof);

        atomicAdd(&dHtdx[ic1], -sxx*Htf);
        atomicAdd(&dHtdy[ic1], -syy*Htf);
        atomicAdd(&dHtdz[ic1], -szz*Htf);

        atomicAdd(&divU[ic1], -US);
    }
    __syncthreads();
}

__global__ void calcGradient_2_d
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
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,
 flow_float* divU
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

        drodx[ic] = drodx[ic]/volume;
        drody[ic] = drody[ic]/volume;
        drodz[ic] = drodz[ic]/volume;

        dTdx[ic] = dTdx[ic]/volume;
        dTdy[ic] = dTdy[ic]/volume;
        dTdz[ic] = dTdz[ic]/volume;

        dHtdx[ic] = dHtdx[ic]/volume;
        dHtdy[ic] = dHtdy[ic]/volume;
        dHtdz[ic] = dHtdz[ic]/volume;


        divU[ic] = divU[ic]/volume;
    }
}

__global__ void limiter_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells, 
 geom_int* cell_planes_index, geom_int* cell_planes,  

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
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,
 flow_float* divU   ) 
 {
    geom_int ic0 = blockDim.x*blockIdx.x + threadIdx.x;

    __syncthreads();

    if (ic0 < nCells) {
        geom_int ip; 
        geom_int ic1; 
        geom_int ic0_temp; 
        geom_int ic1_temp; 
        geom_float dcp_x; 
        geom_float dcp_y; 
        geom_float dcp_z; 

        geom_int index_st = cell_planes_index[ic0];
        geom_int index_en = cell_planes_index[ic0+1];
        geom_int np = index_en - index_st;

        flow_float dUx_max=-1e+30;
        flow_float dUx_min=1e+30;
        flow_float dUy_max=-1e+30;
        flow_float dUy_min=1e+30;
        flow_float dUz_max=-1e+30;
        flow_float dUz_min=1e+30;
        flow_float dro_max=-1e+30;
        flow_float dro_min=1e+30;
        flow_float dP_max =-1e+30;
        flow_float dP_min =1e+30;
        flow_float dT_max =-1e+30;
        flow_float dT_min =1e+30;
        flow_float dHt_max=-1e+30;
        flow_float dHt_min=1e+30;

        flow_float denomi;

        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            ip  = cell_planes[ilp];

            if (ip >= nNormalPlanes) continue;

            ic0_temp = plane_cells[2*ip+0];
            ic1_temp = plane_cells[2*ip+1];
            ic1 = ic0_temp + ic1_temp -ic0;

            dUx_max = max(dUx_max, Ux[ic1]-Ux[ic0]);
            dUy_max = max(dUy_max, Uy[ic1]-Uy[ic0]);
            dUz_max = max(dUz_max, Uz[ic1]-Uz[ic0]);
            dro_max = max(dro_max, ro[ic1]-ro[ic0]);
            dP_max  = max(dP_max , P[ic1] -P[ic0]);
            dT_max  = max(dT_max , T[ic1] -T[ic0]);
            dHt_max = max(dHt_max, Ht[ic1]-Ht[ic0]);

            dUx_min = min(dUx_min, Ux[ic1]-Ux[ic0]);
            dUy_min = min(dUy_min, Uy[ic1]-Uy[ic0]);
            dUz_min = min(dUz_min, Uz[ic1]-Uz[ic0]);
            dro_min = min(dro_min, ro[ic1]-ro[ic0]);
            dP_min  = min(dP_min , P[ic1] -P[ic0]);
            dT_min  = min(dT_min , T[ic1] -T[ic0]);
            dHt_min = min(dHt_min, Ht[ic1]-Ht[ic0]);
        }

        dUx_max = max(dUx_max, 0.0 );
        dUy_max = max(dUy_max, 0.0 );
        dUz_max = max(dUz_max, 0.0 );
        dro_max = max(dro_max, 0.0 );
        dP_max  = max(dP_max , 0.0 );
        dT_max  = max(dT_max , 0.0 );
        dHt_max = max(dHt_max, 0.0 );

        dUx_min = min(dUx_min, 0.0 );
        dUy_min = min(dUy_min, 0.0 );
        dUz_min = min(dUz_min, 0.0 );
        dro_min = min(dro_min, 0.0 );
        dP_min  = min(dP_min , 0.0 );
        dT_min  = min(dT_min , 0.0 );
        dHt_min = min(dHt_min, 0.0 );


        flow_float limiter_Ux = 1e+30;
        flow_float limiter_Uy = 1e+30;
        flow_float limiter_Uz = 1e+30;
        flow_float limiter_ro = 1e+30;
        flow_float limiter_P  = 1e+30;
        flow_float limiter_T  = 1e+30;
        flow_float limiter_Ht = 1e+30;
        flow_float delta;
        flow_float limiter_temp;
        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            ip  = cell_planes[ilp];

            if (ip >= nNormalPlanes) continue;

            dcp_x = pcx[ip] - ccx[ic0];
            dcp_y = pcy[ip] - ccy[ic0];
            dcp_z = pcz[ip] - ccz[ic0];

            delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
                                dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUx_max , dUx_min);
            limiter_temp = venkata_limiter(delta, vol[ic0]);
            limiter_Ux = min(limiter_Ux, limiter_temp);
            
            delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
                                dUydx[ic0], dUydy[ic0], dUydz[ic0], dUy_max , dUy_min);
            limiter_temp = venkata_limiter(delta, vol[ic0]);
            limiter_Uy = min(limiter_Uy, limiter_temp);
 
            delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
                                dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUz_max , dUz_min);
            limiter_temp = venkata_limiter(delta, vol[ic0]);
            limiter_Uz = min(limiter_Uz, limiter_temp);
 
            delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
                                drodx[ic0], drody[ic0], drodz[ic0], dro_max , dro_min);
            limiter_temp = venkata_limiter(delta, vol[ic0]);
            limiter_ro = min(limiter_ro, limiter_temp);
 
            delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
                                dPdx[ic0], dPdy[ic0], dPdz[ic0], dP_max , dP_min);
            limiter_temp = venkata_limiter(delta, vol[ic0]);
            limiter_P = min(limiter_P, limiter_temp);
 
            delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
                                dTdx[ic0], dTdy[ic0], dTdz[ic0], dT_max , dT_min);
            limiter_temp = venkata_limiter(delta, vol[ic0]);
            limiter_T = min(limiter_T, limiter_temp);
  
            delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
                                dHtdx[ic0], dHtdy[ic0], dHtdz[ic0], dHt_max , dHt_min);
            limiter_temp = venkata_limiter(delta, vol[ic0]);
            limiter_Ht = min(limiter_Ht, limiter_temp);
        }

        dUxdx[ic0] = dUxdx[ic0]*limiter_Ux;
        dUxdy[ic0] = dUxdy[ic0]*limiter_Ux;
        dUxdz[ic0] = dUxdz[ic0]*limiter_Ux;

        dUydx[ic0] = dUydx[ic0]*limiter_Uy;
        dUydy[ic0] = dUydy[ic0]*limiter_Uy;
        dUydz[ic0] = dUydz[ic0]*limiter_Uy;

        dUzdx[ic0] = dUzdx[ic0]*limiter_Uz;
        dUzdy[ic0] = dUzdy[ic0]*limiter_Uz;
        dUzdz[ic0] = dUzdz[ic0]*limiter_Uz;

        drodx[ic0] = drodx[ic0]*limiter_ro;
        drody[ic0] = drody[ic0]*limiter_ro;
        drodz[ic0] = drodz[ic0]*limiter_ro;

        dPdx[ic0] = dPdx[ic0]*limiter_P;
        dPdy[ic0] = dPdy[ic0]*limiter_P;
        dPdz[ic0] = dPdz[ic0]*limiter_P;

        dTdx[ic0] = dTdx[ic0]*limiter_T;
        dTdy[ic0] = dTdy[ic0]*limiter_T;
        dTdz[ic0] = dTdz[ic0]*limiter_T;

        dHtdx[ic0] = dHtdx[ic0]*limiter_Ht;
        dHtdy[ic0] = dHtdy[ic0]*limiter_Ht;
        dHtdz[ic0] = dHtdz[ic0]*limiter_Ht;
    }
    __syncthreads();
}

__device__ flow_float calcDeltaIJ(geom_float pcx , geom_float pcy, geom_float pcz, 
                                  flow_float dudx, flow_float dudy,flow_float dudz,
                                  flow_float delu_max, flow_float delu_min ) {
    flow_float denomi = dudx*pcx + dudy*pcx + dudz*pcz;
    flow_float delta;

    if (denomi > 1e-20) {
        delta = delu_max/denomi;
    } else if (denomi < -1e-20) {
        delta = delu_min/denomi;
    } else {
        delta = 0.0;
    }

    return delta;
}

__device__ flow_float venkata_limiter(flow_float x , flow_float volume) {
    flow_float K = 10.0;
    flow_float eps = sqrt(K*K*K*volume);
    return (x*x + 2.0*x + eps*eps)/(x*x + x + 2.0 + eps*eps);
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
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drodx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drody"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drodz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["divU"] , 0.0, msh.nCells*sizeof(flow_float)));

    // sum over planes
    calcGradient_1_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        var.c_d["ro"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"]  , 
        var.c_d["T"] ,
        var.c_d["Ht"] ,

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],
        var.c_d["dHtdx"] , var.c_d["dHtdy"] , var.c_d["dHtdz"],
        var.c_d["divU"]  
    ) ;

    calcGradient_2_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
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
        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],
        var.c_d["dHtdx"] , var.c_d["dHtdy"] , var.c_d["dHtdz"],
        var.c_d["divU"]  
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    // sum over planes
    if (cfg.limiter == 1) {
        limiter_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            // mesh structure
            msh.nCells,
            msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
            msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
            var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

            // basic variables
            var.c_d["ro"] ,
            var.c_d["Ux"] ,
            var.c_d["Uy"] ,
            var.c_d["Uz"] ,
            var.c_d["P"]  , 
            var.c_d["T"] ,
            var.c_d["Ht"] ,

            // gradient
            var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
            var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
            var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
            var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
            var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
            var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],
            var.c_d["dHtdx"] , var.c_d["dHtdy"] , var.c_d["dHtdz"],
            var.c_d["divU"]  
        ) ;
    }


    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}