#include "limiter_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"

__global__ void fill_limiter_d(flow_float* values, geom_int nValues, flow_float value)
{
    geom_int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < nValues) {
        values[index] = value;
    }
}

// Limiters for Unstructured Higher-Order Accurate Solutions of the Euler Equations
// Krzysztof Michalak

__device__ flow_float venkata_limiter(flow_float delta_p_max, flow_float delta_p_min, 
                                      flow_float delta_m, flow_float volume) {

    flow_float K = 1.;
    flow_float eps2 = K*K*K*volume;
    //return (x*x + 2.0*x + eps*eps)/(x*x + x + 2.0 + eps*eps);
    flow_float res;

    // K11: 元は /(...)/delta_m と除算2回。/(denom*delta_m) に統合して除算1回に（compute律速の limiter 向け）。
    if (delta_m > 1e-20) {
        flow_float delta_p = delta_p_max;
        res = ((delta_p*delta_p+eps2)*delta_m +2*delta_m*delta_m*delta_p)
              /((delta_p*delta_p +2.0*delta_m*delta_m +delta_p*delta_m +eps2)*delta_m);
    } else if (delta_m < -1e-20) {
        flow_float delta_p = delta_p_min;
        res = ((delta_p*delta_p+eps2)*delta_m +2*delta_m*delta_m*delta_p)
              /((delta_p*delta_p +2.0*delta_m*delta_m +delta_p*delta_m +eps2)*delta_m);
    } else {
        res = 1.0;
    }

    return res;
}

__device__ flow_float barth_Jespersen_limiter(flow_float delta_p_max, flow_float delta_p_min, 
                                              flow_float delta_m, flow_float volume) {

    flow_float res;

    if (delta_m > 1e-20) {
        res = min(1.0, delta_p_max/delta_m);
    } else if (delta_m < -1e-20) {
        res = min(1.0, delta_p_min/delta_m);
    } else {
        res = 1.0;
    }

    return min(res, 1.0);
}


//__device__ flow_float nishikawa_r1_limiter(deltas delta_dash) {
//    flow_float res;
//    flow_float del = delta_dash.del;
//    flow_float rik = delta_dash.rik;
//
//    if (del < 0.0) {
//        res = 0.0;
//    } else if (del >= 0.0 && del < 1.0/rik) {
//        res = rik*del*(1.0+0.5*(1.0-rik)*del);
//    } else if (del < 1.0 && del >= 1.0/rik) {
//        res = 1.0 + rik*0.5/(1.0-rik)*(1.0-del)*(1.0-del);
//    } else if (del >= 1.0) {
//        res = 1.0;
//    }
//
//    return res;
//}

// Modified multi-dimensional limiting process with enhanced shock stability on unstructured grids
__global__ void limiter_r1_d
( 
 int limiter_scheme,
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells, 
 geom_int* cell_planes_index, geom_int* cell_planes,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,

 // variables
 flow_float* Q  ,

 flow_float* limiter_Q  ,

 flow_float* dQdx  , flow_float* dQdy , flow_float* dQdz
) 
{
    geom_int ic0 = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic0 < nCells) {
        if (limiter_scheme == 0) {
            limiter_Q[ic0] = 1.0;
            return;
        }

        geom_int ip; 
        geom_int ic1; 
        geom_float dcp_x; 
        geom_float dcp_y; 
        geom_float dcp_z; 

        geom_float dcp2_x; 
        geom_float dcp2_y; 
        geom_float dcp2_z; 

        geom_int index_st = cell_planes_index[ic0];
        geom_int index_en = cell_planes_index[ic0+1];
        geom_int np = index_en - index_st;

        flow_float Q_max =Q[ic0];
        flow_float Q_min =Q[ic0];

        flow_float denomi;

        flow_float limiter_Q_temp  = 1.0;
        flow_float limiter_Q_temp2 = 1.0;

        deltas delta;
        deltas delta_dash;
        flow_float volume = vol[ic0];

        int int_one=1;
        int int_two=2;
        int int_three=3;

        flow_float (*limiter_function)(flow_float, flow_float, flow_float, flow_float);

        if (limiter_scheme == int_one) { //barth
            limiter_function = barth_Jespersen_limiter;
        } else if (limiter_scheme == int_two or limiter_scheme == -1) { //venkata
            limiter_function = venkata_limiter;
        //} else if (limiter_scheme == int_three or limiter_scheme == -1) { //
        //    limiter_function = nishikawa_r1_limiter;
        } else {
            printf("Error: something wrong");
        }

        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            ip  = cell_planes[ilp];

            if (ip >= nNormalPlanes) continue;

            ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] -ic0;

            Q_max = max(Q_max, Q[ic1]);
            Q_min = min(Q_min, Q[ic1]);
        }

        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            ip  = cell_planes[ilp];

            if (ip >= nNormalPlanes) continue;

            ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] -ic0;

            dcp_x = pcx[ip] - ccx[ic0];
            dcp_y = pcy[ip] - ccy[ic0];
            dcp_z = pcz[ip] - ccz[ic0];

            // K9: dcp2_*, ri, rk, rik は未使用（dead）。compute律速の limiter から sqrt 2回/面/変数を除去。

            flow_float delta_p_max;
            flow_float delta_p_min;
            flow_float delta_m;

            flow_float Qt = Q[ic0] + dQdx[ic0]*dcp_x + dQdy[ic0]*dcp_y + dQdz[ic0]*dcp_z;

            delta_p_max = Q_max - Q[ic0];
            delta_p_min = Q_min - Q[ic0];
            delta_m     = Qt    - Q[ic0];

            limiter_Q_temp = limiter_function(delta_p_max, delta_p_min, delta_m, volume);

            limiter_Q_temp2 = min(limiter_Q_temp2, limiter_Q_temp);
        }

        limiter_Q[ic0] = min(max(limiter_Q_temp2, 0.0),1.0);

            //ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] -ic0;
        //limiter_Q[ic1] = dQdx[ic0];
    }
}


void limiter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    fill_limiter_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(var.c_d["limiter_ro"], msh.nCells_all, 1.0);
    fill_limiter_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(var.c_d["limiter_Ux"], msh.nCells_all, 1.0);
    fill_limiter_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(var.c_d["limiter_Uy"], msh.nCells_all, 1.0);
    fill_limiter_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(var.c_d["limiter_Uz"], msh.nCells_all, 1.0);
    fill_limiter_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(var.c_d["limiter_P"], msh.nCells_all, 1.0);

    if (cfg.limiter == 0) {
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        return;
    }

    // ro
    limiter_r1_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> ( 
        cfg.limiter,
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],

        // basic variables
        var.c_d["ro"] ,
        var.c_d["limiter_ro"] ,

        // gradient
        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"]
    ) ;

    // Ux
    limiter_r1_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> ( 
        cfg.limiter,
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],

        // basic variables
        var.c_d["Ux"] ,
        var.c_d["limiter_Ux"] ,

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"]
    ) ;

    // Uy
    limiter_r1_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> ( 
        cfg.limiter,
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],

        // basic variables
        var.c_d["Uy"] ,
        var.c_d["limiter_Uy"] ,

        // gradient
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"]
    ) ;

    // Uz
    limiter_r1_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> ( 
        cfg.limiter,
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],

        // basic variables
        var.c_d["Uz"] ,
        var.c_d["limiter_Uz"] ,

        // gradient
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"]
    ) ;


    // Ht
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["limiter_Ht"] , 1.0, msh.nCells_all*sizeof(flow_float)));
    //limiter_r1_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> ( 
    //    cfg.limiter,
    //    // mesh structure
    //    msh.nCells,
    //    msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
    //    msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
    //    var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
    //    var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],

    //    // basic variables
    //    var.c_d["Ht"] ,
    //    var.c_d["limiter_Ht"] ,

    //    // gradient
    //    var.c_d["dHtdx"] , var.c_d["dHtdy"] , var.c_d["dHtdz"]
    //) ;

    // P
    limiter_r1_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> ( 
        cfg.limiter,
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],

        // basic variables
        var.c_d["P"] ,
        var.c_d["limiter_P"] ,

        // gradient
        var.c_d["dPdx"] , var.c_d["dPdy"] , var.c_d["dPdz"]
    ) ;


    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}