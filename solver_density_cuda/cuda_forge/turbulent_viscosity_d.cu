#include "calcGradient_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"


__global__ void WALE_d
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
 flow_float* vis_turb , flow_float* wall_dist
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells ) {

        flow_float karman = 0.41;
        //flow_float Cw = 0.5;
        flow_float Cw = 0.325;

        geom_float volume = vol[ic];


        flow_float g11 = dUxdx[ic];
        flow_float g12 = dUxdy[ic];
        flow_float g13 = dUxdz[ic];
        flow_float g21 = dUydx[ic];
        flow_float g22 = dUydy[ic];
        flow_float g23 = dUydz[ic];
        flow_float g31 = dUzdx[ic];
        flow_float g32 = dUzdy[ic];
        flow_float g33 = dUzdz[ic];

        flow_float Sd11 = 0.5*(g11*g11+g11*g11) - 1.0/3.0*g11*g11;
        flow_float Sd12 = 0.5*(g12*g12+g21*g21);
        flow_float Sd13 = 0.5*(g13*g13+g31*g31);
        flow_float Sd21 = 0.5*(g21*g21+g12*g12);
        flow_float Sd22 = 0.5*(g22*g22+g22*g22) - 1.0/3.0*g22*g22;
        flow_float Sd23 = 0.5*(g23*g23+g32*g32);
        flow_float Sd31 = 0.5*(g31*g31+g13*g13);
        flow_float Sd32 = 0.5*(g32*g32+g23*g23);
        flow_float Sd33 = 0.5*(g33*g33+g33*g33) - 1.0/3.0*g33*g33;


        flow_float S11 = 0.5*(g11+g11);
        flow_float S12 = 0.5*(g12+g21);
        flow_float S13 = 0.5*(g13+g31);
        flow_float S21 = 0.5*(g21+g12);
        flow_float S22 = 0.5*(g22+g22);
        flow_float S23 = 0.5*(g23+g32);
        flow_float S31 = 0.5*(g31+g13);
        flow_float S32 = 0.5*(g32+g23);
        flow_float S33 = 0.5*(g33+g33);

        flow_float SdijSdij = Sd11*Sd11 + Sd12*Sd12 + Sd13*Sd13
                            + Sd21*Sd21 + Sd22*Sd22 + Sd23*Sd23
                            + Sd31*Sd31 + Sd32*Sd32 + Sd33*Sd33;

        flow_float SijSij = S11*S11 + S12*S12 + S13*S13
                          + S21*S21 + S22*S22 + S23*S23
                          + S31*S31 + S32*S32 + S33*S33;

        geom_float d = wall_dist[ic];
        flow_float Ls = min(karman*d , Cw*pow(vol[ic],1.0/3.0));

        vis_turb[ic] =  ro[ic]*Ls*Ls*pow(SdijSdij, 3.0/2.0)/(pow(SijSij, 5.0/2.0) + pow(SdijSdij, 5.0/4.0));
    }
}


void turbulent_viscosity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.LESorRANS == 1) { // LES

        if (cfg.LESmodel = 1) {
            WALE_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
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
                var.c_d["vis_turb"] , var.c_d["wall_dist"]
            ) ;
        }

    } else {
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["vis_turb"], 0.0, msh.nCells_all*sizeof(flow_float)));
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}