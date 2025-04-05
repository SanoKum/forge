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
 flow_float* roe  ,
 flow_float* Ht  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx , flow_float* droUxdy, flow_float* droUxdz,
 //flow_float* droUydx , flow_float* droUydy, flow_float* droUydz,
 //flow_float* droUzdx , flow_float* droUzdy, flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,
 //flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,

 flow_float* divU
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;

    __syncthreads();

    //if (ip < nNormalPlanes) {
    if (ip < nPlanes) {
        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float f   = fx[ip];
        flow_float Uxf = f*Ux[ic0] + (1.0-f)*Ux[ic1];
        flow_float Uyf = f*Uy[ic0] + (1.0-f)*Uy[ic1];
        flow_float Uzf = f*Uz[ic0] + (1.0-f)*Uz[ic1];

        flow_float Pf  = f*P[ic0]  + (1.0-f)*P[ic1];
        flow_float Tf  = f*T[ic0]  + (1.0-f)*T[ic1];
        flow_float rof = f*ro[ic0] + (1.0-f)*ro[ic1];
        flow_float roef= f*roe[ic0]+ (1.0-f)*roe[ic1];
        flow_float Htf= f*Ht[ic0]+ (1.0-f)*Ht[ic1];
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];

        flow_float US  = Uxf*sxx +Uyf*syy +Uzf*szz;

//        if (ic0 == 0 or ic1 == 0) {
//            printf("ic0 = %d, ic1 = %d, f = %f, Uxf = %f, Uyf = %f, Uzf = %f\n", ic0, ic1, f, Uxf, Uyf, Uzf);
//            printf("sxx = %f, syy = %f, szz = %f\n", sxx, syy, szz);
//            printf("Pf  = %f\n", Pf);    
//            printf("US  = %f\n", US);    
//        }
//
//        if (ic0 == 1 or ic1 == 1) {
//            printf("ic0 = %d, ic1 = %d, f = %f, Uxf = %f, Uyf = %f, Uzf = %f\n", ic0, ic1, f, Uxf, Uyf, Uzf);
//        }



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

        //atomicAdd(&dHtdx[ic0], sxx*Htf);
        //atomicAdd(&dHtdy[ic0], syy*Htf);
        //atomicAdd(&dHtdz[ic0], szz*Htf);



        atomicAdd(&drodx[ic0], sxx*rof);
        atomicAdd(&drody[ic0], syy*rof);
        atomicAdd(&drodz[ic0], szz*rof);


        //atomicAdd(&droUxdx[ic0], sxx*rof*Uxf);
        //atomicAdd(&droUxdy[ic0], syy*rof*Uxf);
        //atomicAdd(&droUxdz[ic0], szz*rof*Uxf);

        //atomicAdd(&droUydx[ic0], sxx*rof*Uyf);
        //atomicAdd(&droUydy[ic0], syy*rof*Uyf);
        //atomicAdd(&droUydz[ic0], szz*rof*Uyf);

        //atomicAdd(&droUzdx[ic0], sxx*rof*Uzf);
        //atomicAdd(&droUzdy[ic0], syy*rof*Uzf);
        //atomicAdd(&droUzdz[ic0], szz*rof*Uzf);

        //atomicAdd(&droedx[ic0], sxx*roef);
        //atomicAdd(&droedy[ic0], syy*roef);
        //atomicAdd(&droedz[ic0], szz*roef);


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

        //atomicAdd(&dHtdx[ic1], -sxx*Pf);
        //atomicAdd(&dHtdy[ic1], -syy*Pf);
        //atomicAdd(&dHtdz[ic1], -szz*Pf);

        atomicAdd(&drodx[ic1], -sxx*rof);
        atomicAdd(&drody[ic1], -syy*rof);
        atomicAdd(&drodz[ic1], -szz*rof);

        //atomicAdd(&droUxdx[ic1], -sxx*rof*Uxf);
        //atomicAdd(&droUxdy[ic1], -syy*rof*Uxf);
        //atomicAdd(&droUxdz[ic1], -szz*rof*Uxf);

        //atomicAdd(&droUydx[ic1], -sxx*rof*Uyf);
        //atomicAdd(&droUydy[ic1], -syy*rof*Uyf);
        //atomicAdd(&droUydz[ic1], -szz*rof*Uyf);

        //atomicAdd(&droUzdx[ic1], -sxx*rof*Uzf);
        //atomicAdd(&droUzdy[ic1], -syy*rof*Uzf);
        //atomicAdd(&droUzdz[ic1], -szz*rof*Uzf);

        //atomicAdd(&droedx[ic1], -sxx*roef);
        //atomicAdd(&droedy[ic1], -syy*roef);
        //atomicAdd(&droedz[ic1], -szz*roef);

        atomicAdd(&divU[ic1], -US);
    }
    __syncthreads();
}

__global__ void calcGradient_b_d
( 
  // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

// variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* Tt   ,
 flow_float* Pt   ,
 flow_float* Ts   ,
 flow_float* Ps   ,
 flow_float* Ht   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx  , flow_float* droUxdy , flow_float* droUxdz,
 //flow_float* droUydx  , flow_float* droUydy , flow_float* droUydz,
 //flow_float* droUzdx  , flow_float* droUzdy , flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,
 //flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,

 flow_float* divU
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];
       
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        flow_float rof   = ro[ib];
        flow_float roUxf = roUx[ib];
        flow_float roUyf = roUy[ib];
        flow_float roUzf = roUz[ib];
        flow_float roef  = roe[ib];
        flow_float Uxf = roUxf/rof;
        flow_float Uyf = roUyf/rof;
        flow_float Uzf = roUzf/rof;
        flow_float ek = 0.5*(Uxf*Uxf + Uyf*Uyf + Uzf*Uzf);
        flow_float Pf = Ps[ib];
        flow_float Tf = Ts[ib];
        flow_float Htf= roef/rof + Pf/rof;

        flow_float US  = Uxf*sxx +Uyf*syy +Uzf*szz;

        atomicAdd(&dUxdx[ic], sxx*Uxf);
        atomicAdd(&dUxdy[ic], syy*Uxf);
        atomicAdd(&dUxdz[ic], szz*Uxf);

        atomicAdd(&dUydx[ic], sxx*Uyf);
        atomicAdd(&dUydy[ic], syy*Uyf);
        atomicAdd(&dUydz[ic], szz*Uyf);

        atomicAdd(&dUzdx[ic], sxx*Uzf);
        atomicAdd(&dUzdy[ic], syy*Uzf);
        atomicAdd(&dUzdz[ic], szz*Uzf);

        atomicAdd(&dTdx[ic], sxx*Tf);
        atomicAdd(&dTdy[ic], syy*Tf);
        atomicAdd(&dTdz[ic], szz*Tf);

        atomicAdd(&dPdx[ic], sxx*Pf);
        atomicAdd(&dPdy[ic], syy*Pf);
        atomicAdd(&dPdz[ic], szz*Pf);

        //atomicAdd(&dHtdx[ic], sxx*Htf);
        //atomicAdd(&dHtdy[ic], syy*Htf);
        //atomicAdd(&dHtdz[ic], szz*Htf);

        atomicAdd(&drodx[ic], sxx*rof);
        atomicAdd(&drody[ic], syy*rof);
        atomicAdd(&drodz[ic], szz*rof);

        //atomicAdd(&droUxdx[ic], sxx*rof*Uxf);
        //atomicAdd(&droUxdy[ic], syy*rof*Uxf);
        //atomicAdd(&droUxdz[ic], szz*rof*Uxf);

        //atomicAdd(&droUydx[ic], sxx*rof*Uyf);
        //atomicAdd(&droUydy[ic], syy*rof*Uyf);
        //atomicAdd(&droUydz[ic], szz*rof*Uyf);

        //atomicAdd(&droUzdx[ic], sxx*rof*Uzf);
        //atomicAdd(&droUzdy[ic], syy*rof*Uzf);
        //atomicAdd(&droUzdz[ic], szz*rof*Uzf);

        //atomicAdd(&droedx[ic], sxx*roef);
        //atomicAdd(&droedy[ic], syy*roef);
        //atomicAdd(&droedz[ic], szz*roef);

        atomicAdd(&divU[ic], US);
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
 flow_float* roe  ,
 flow_float* Ht  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx  , flow_float* droUxdy , flow_float* droUxdz,
 //flow_float* droUydx  , flow_float* droUydy , flow_float* droUydz,
 //flow_float* droUzdx  , flow_float* droUzdy , flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,
 //flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,

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

        //droUxdx[ic] = droUxdx[ic]/volume;
        //droUxdy[ic] = droUxdy[ic]/volume;
        //droUxdz[ic] = droUxdz[ic]/volume;

        //droUydx[ic] = droUydx[ic]/volume;
        //droUydy[ic] = droUydy[ic]/volume;
        //droUydz[ic] = droUydz[ic]/volume;

        //droUzdx[ic] = droUzdx[ic]/volume;
        //droUzdy[ic] = droUzdy[ic]/volume;
        //droUzdz[ic] = droUzdz[ic]/volume;

        //droedx[ic] = droedx[ic]/volume;
        //droedy[ic] = droedy[ic]/volume;
        //droedz[ic] = droedz[ic]/volume;

        //dHtdx[ic] = dHtdx[ic]/volume;
        //dHtdy[ic] = dHtdy[ic]/volume;
        //dHtdz[ic] = dHtdz[ic]/volume;

        divU[ic] = divU[ic]/volume;
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
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drodx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drody"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drodz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["divU"] , 0.0, msh.nCells*sizeof(flow_float)));

    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUxdx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUxdy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUxdz"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUydx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUydy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUydz"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUzdx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUzdy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUzdz"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droedx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droedy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droedz"], 0.0, msh.nCells*sizeof(flow_float)));


    // sum over planes
    //calcGradient_1_d<<<cuda_cfg.dimGrid_nplane , cuda_cfg.dimBlock>>> ( 
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
        var.c_d["roe"] ,
        var.c_d["Ht"] ,

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],

        //var.c_d["droUxdx"] , var.c_d["droUxdy"] , var.c_d["droUxdz"],
        //var.c_d["droUydx"] , var.c_d["droUydy"] , var.c_d["droUydz"],
        //var.c_d["droUzdx"] , var.c_d["droUzdy"] , var.c_d["droUzdz"],
        //var.c_d["droedx"]  , var.c_d["droedy"]  , var.c_d["droedz"],
        //var.c_d["dHtdx"]  , var.c_d["dHtdy"]  , var.c_d["dHtdz"],
 
        var.c_d["divU"]  
    ) ;


    //for (auto& bc : msh.bconds)
    //{
    //    calcGradient_b_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
    //        // mesh structure
    //        bc.iPlanes.size(),
    //        bc.map_bplane_plane_d,  
    //        bc.map_bplane_cell_d,  
    //        bc.map_bplane_cell_ghst_d,

    //        // mesh structure
    //        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
    //        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
    //        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

    //        // boundary variables
    //        bc.bvar_d["ro"],
    //        bc.bvar_d["roUx"],
    //        bc.bvar_d["roUy"],
    //        bc.bvar_d["roUz"],
    //        bc.bvar_d["roe"],
    //        bc.bvar_d["Ux"],
    //        bc.bvar_d["Uy"],
    //        bc.bvar_d["Uz"],
    //        bc.bvar_d["Tt"],
    //        bc.bvar_d["Pt"],
    //        bc.bvar_d["Ts"],
    //        bc.bvar_d["Ps"],
    //        bc.bvar_d["Ht"],

    //        // gradient
    //        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
    //        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
    //        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
    //        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
    //        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
    //        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],

    //        //var.c_d["droUxdx"] , var.c_d["droUxdy"] , var.c_d["droUxdz"],
    //        //var.c_d["droUydx"] , var.c_d["droUydy"] , var.c_d["droUydz"],
    //        //var.c_d["droUzdx"] , var.c_d["droUzdy"] , var.c_d["droUzdz"],
    //        //var.c_d["droedx"]  , var.c_d["droedy"]  , var.c_d["droedz"],
    //        //var.c_d["dHtdx"]  , var.c_d["dHtdy"]  , var.c_d["dHtdz"],
 
    //        var.c_d["divU"]  
    //    ) ;
    //}
 

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
        var.c_d["roe"] ,
        var.c_d["Ht"] ,

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],

        //var.c_d["droUxdx"] , var.c_d["droUxdy"] , var.c_d["droUxdz"],
        //var.c_d["droUydx"] , var.c_d["droUydy"] , var.c_d["droUydz"],
        //var.c_d["droUzdx"] , var.c_d["droUzdy"] , var.c_d["droUzdz"],
        //var.c_d["droedx"]  , var.c_d["droedy"]  , var.c_d["droedz"],
        //var.c_d["dHtdx"]  , var.c_d["dHtdy"]  , var.c_d["dHtdz"],
 
        var.c_d["divU"]  
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}