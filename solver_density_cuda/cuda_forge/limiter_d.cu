#include "calcGradient_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"

__global__ void limiter_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells, 
 geom_int* cell_planes_index, geom_int* cell_planes,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,

 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* Ht  ,
 flow_float* limiter  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz
 ) 
 {
    geom_int ic0 = blockDim.x*blockIdx.x + threadIdx.x;

    __syncthreads();

    if (ic0 < nCells) {
        geom_int ip; 
        geom_int ic1; 
        geom_float dcp_x; 
        geom_float dcp_y; 
        geom_float dcp_z; 

        geom_int index_st = cell_planes_index[ic0];
        geom_int index_en = cell_planes_index[ic0+1];
        geom_int np = index_en - index_st;

        //flow_float dUx_max=-1e+30;
        //flow_float dUx_min=1e+30;
        //flow_float dUy_max=-1e+30;
        //flow_float dUy_min=1e+30;
        //flow_float dUz_max=-1e+30;
        //flow_float dUz_min=1e+30;
        flow_float dro_max=-1e+30;
        flow_float dro_min=1e+30;
        //flow_float dP_max =-1e+30;
        //flow_float dP_min =1e+30;
        //flow_float dT_max =-1e+30;
        //flow_float dT_min =1e+30;
        //flow_float dHt_max=-1e+30;
        //flow_float dHt_min=1e+30;

        flow_float denomi;

        //flow_float limiter_Ux ;
        //flow_float limiter_Uy ;
        //flow_float limiter_Uz ;
        flow_float limiter_ro ;
        //flow_float limiter_P  ;
        //flow_float limiter_T  ;
        //flow_float limiter_Ht ;
        flow_float delta;
        flow_float deltap;
        flow_float deltam;
        flow_float volume = vol[ic0];


        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            ip  = cell_planes[ilp];

            if (ip >= nNormalPlanes) continue;

            ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] -ic0;

            //dUx_max = max(dUx_max, Ux[ic1]-Ux[ic0]);
            //dUy_max = max(dUy_max, Uy[ic1]-Uy[ic0]);
            //dUz_max = max(dUz_max, Uz[ic1]-Uz[ic0]);
            dro_max = max(dro_max, ro[ic1]-ro[ic0]);
            //dP_max  = max(dP_max , P[ic1] -P[ic0]);
            //dT_max  = max(dT_max , T[ic1] -T[ic0]);
            //dHt_max = max(dHt_max, Ht[ic1]-Ht[ic0]);

            //dUx_min = min(dUx_min, Ux[ic1]-Ux[ic0]);
            //dUy_min = min(dUy_min, Uy[ic1]-Uy[ic0]);
            //dUz_min = min(dUz_min, Uz[ic1]-Uz[ic0]);
            dro_min = min(dro_min, ro[ic1]-ro[ic0]);
            //dP_min  = min(dP_min , P[ic1] -P[ic0]);
            //dT_min  = min(dT_min , T[ic1] -T[ic0]);
            //dHt_min = min(dHt_min, Ht[ic1]-Ht[ic0]);
        }

        //dUx_max = max(dUx_max, 0.0 );
        //dUy_max = max(dUy_max, 0.0 );
        //dUz_max = max(dUz_max, 0.0 );
        dro_max = max(dro_max, 0.0 );
        //dP_max  = max(dP_max , 0.0 );
        //dT_max  = max(dT_max , 0.0 );
        //dHt_max = max(dHt_max, 0.0 );

        //dUx_min = min(dUx_min, 0.0 );
        //dUy_min = min(dUy_min, 0.0 );
        //dUz_min = min(dUz_min, 0.0 );
        dro_min = min(dro_min, 0.0 );
        //dP_min  = min(dP_min , 0.0 );
        //dT_min  = min(dT_min , 0.0 );
        //dHt_min = min(dHt_min, 0.0 );


        //limiter_Ux = 1e+30;
        //limiter_Uy = 1e+30;
        //limiter_Uz = 1e+30;
        limiter_ro = 1e+30;
        //limiter_P  = 1e+30;
        //limiter_T  = 1e+30;
        //limiter_Ht = 1e+30;

        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            ip  = cell_planes[ilp];

            if (ip >= nNormalPlanes) continue;

            dcp_x = pcx[ip] - ccx[ic0];
            dcp_y = pcy[ip] - ccy[ic0];
            dcp_z = pcz[ip] - ccz[ic0];

            //delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
            //                    dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUx_max , dUx_min);
            //limiter_Ux = min(limiter_Ux, venkata_limiter(delta, volume  ));
            
            //delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
            //                    dUydx[ic0], dUydy[ic0], dUydz[ic0], dUy_max , dUy_min);
            //limiter_Uy = min(limiter_Uy, venkata_limiter(delta, volume  ));
 
            //delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
            //                    dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUz_max , dUz_min);
            //limiter_Uz = min(limiter_Uz, venkata_limiter(delta, volume  ));
 
            //calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
            //            drodx[ic0], drody[ic0], drodz[ic0], dro_max , dro_min,
            //            &delta, &deltap, &deltam);

            flow_float denomi = drodx[ic0]*dcp_x + drody[ic0]*dcp_y + drodz[ic0]*dcp_z;

            if (denomi > 1e-12) {
                delta = dro_max/denomi;
            } else if (denomi < -1e-12) {
                delta = dro_min/denomi;
            } else {
                delta = 1.0;
            }

            //barth&Jespersen
            limiter_ro = min(limiter_ro, min(1.0,delta));

            //venkata
            //deltap = dro_max;
            //deltam = dro_min;

            //limiter_ro = min(limiter_ro, venkata_limiter(deltap, deltam, volume));

            //printf("ip=%d, ic0=%d\n", ip, ic0);
            //printf("dx=%e, dy=%e, ez=%e\n", dcp_x, dcp_y, dcp_z);
            //printf("drodx=%e, drody=%e, drodz=%e\n", drodx[ic0], drodx[ic0], drodx[ic0]);
            //printf("delta=%e, deltap=%e, deltam=%e\n", delta, deltap, deltam);
            //printf("limiter_ro=%e\n", limiter_ro);

 
            //delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
            //                    dPdx[ic0], dPdy[ic0], dPdz[ic0], dP_max , dP_min);
            //limiter_P = min(limiter_P, venkata_limiter(delta, volume  ));
 
            //delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
            //                    dTdx[ic0], dTdy[ic0], dTdz[ic0], dT_max , dT_min);
            //limiter_T = min(limiter_T, venkata_limiter(delta, volume  ));
  
            //delta = calcDeltaIJ(dcp_x, dcp_y, dcp_z, 
            //                    dHtdx[ic0], dHtdy[ic0], dHtdz[ic0], dHt_max , dHt_min);
            //limiter_Ht = min(limiter_Ht, venkata_limiter(delta, volume  ));
        }

        limiter[ic0] = limiter_ro;
        
        //printf("limiter=%e\n", limiter_ro);

        //dUxdx[ic0] *= limiter_Ux;
        //dUxdy[ic0] *= limiter_Ux;
        //dUxdz[ic0] *= limiter_Ux;

        //dUydx[ic0] *= limiter_Uy;
        //dUydy[ic0] *= limiter_Uy;
        //dUydz[ic0] *= limiter_Uy;

        //dUzdx[ic0] *= limiter_Uz;
        //dUzdy[ic0] *= limiter_Uz;
        //dUzdz[ic0] *= limiter_Uz;

        //drodx[ic0] *= limiter_ro;
        //drody[ic0] *= limiter_ro;
        //drodz[ic0] *= limiter_ro;

        //dPdx[ic0] *= limiter_P;
        //dPdy[ic0] *= limiter_P;
        //dPdz[ic0] *= limiter_P;

        //dTdx[ic0] *= limiter_T;
        //dTdy[ic0] *= limiter_T;
        //dTdz[ic0] *= limiter_T;

        //dHtdx[ic0] *= limiter_Ht;
        //dHtdy[ic0] *= limiter_Ht;
        //dHtdz[ic0] *= limiter_Ht;
    }
    __syncthreads();
}

//__device__ void calcDeltaIJ(geom_float pcx , geom_float pcy, geom_float pcz, 
//                            flow_float dudx, flow_float dudy,flow_float dudz,
//                            flow_float delu_max, flow_float delu_min ,
//                            flow_float* delta, flow_float* deltap, flow_float* deltam
//                            ) {
//    flow_float denomi = dudx*pcx + dudy*pcy + dudz*pcz;
//    flow_float zero=0.0;
//    flow_float res1;
//    flow_float res2;
//    flow_float res3;
//
//    if (denomi > 1e-12) {
//        res1 = delu_max/denomi;
//    } else if (denomi < -1e-12) {
//        res1 = delu_min/denomi;
//    } else {
//        res1 = zero;
//    }
//
//    res2 = delu_max/denomi;
//    res3 = delu_min/denomi;
//
//    *delta = res1;
//    *deltap= res2;
//    *deltam= res3;
//
//}
//
// Limiters for Unstructured Higher-Order Accurate Solutions of the Euler Equations
// Krzysztof Michalak
__device__ flow_float venkata_limiter(flow_float delp, flow_float delm , flow_float volume) {
    flow_float K = 10.0;
    flow_float eps2 = K*K*K*volume;
    //return (x*x + 2.0*x + eps*eps)/(x*x + x + 2.0 + eps*eps);
    flow_float res;
    if (abs(delm)>1e-12) {
        res = ((delp*delp+eps2)*delm + 2.0*delm*delm*delp)/(delp*delp+2.0*delm*delm+delm*delp+eps2)/delm;
    } else {
        res = 0.0;
    }
    return res;
}


void limiter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["limiter"] , 1.0, msh.nCells_all*sizeof(flow_float)));

    // sum over planes
    if (cfg.limiter == 1) {
        limiter_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> ( 
            // mesh structure
            msh.nCells,
            msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
            msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],

            // basic variables
            var.c_d["ro"] ,
            var.c_d["Ux"] ,
            var.c_d["Uy"] ,
            var.c_d["Uz"] ,
            var.c_d["P"]  , 
            var.c_d["T"] ,
            var.c_d["Ht"] ,
            var.c_d["limiter"] ,

            // gradient
            var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
            var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
            var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
            var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
            var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
            var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],
            var.c_d["dHtdx"] , var.c_d["dHtdy"] , var.c_d["dHtdz"]
        ) ;
    }


    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}