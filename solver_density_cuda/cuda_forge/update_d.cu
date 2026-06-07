#include "dependentVariables_d.cuh"

__global__ void updateVariablesOuter_d
( 
 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* roK  ,
 flow_float* roOmega  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,
 flow_float* roKN ,
 flow_float* roOmegaN ,
 
 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM 

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells_all) {

        roN[ic]   = ro[ic];
        roUxN[ic] = roUx[ic];
        roUyN[ic] = roUy[ic];
        roUzN[ic] = roUz[ic];
        roeN[ic]  = roe[ic];
        roKN[ic] = roK[ic];
        roOmegaN[ic] = roOmega[ic];

        roM[ic]   = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic]  = roe[ic];
        roKM[ic] = roK[ic];
        roOmegaM[ic] = roOmega[ic];

    }
}


void updateVariablesOuter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    updateVariablesOuter_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"]  , var.c_d["roe"]  , var.c_d["roK"] , var.c_d["roOmega"] ,
        var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"] , var.c_d["roeN"] , var.c_d["roKN"] , var.c_d["roOmegaN"] ,
        var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"] , var.c_d["roeM"] , var.c_d["roKM"] , var.c_d["roOmegaM"] 
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

__global__ void updateVariablesInner_d
( 
 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* roK  ,
 flow_float* roOmega  ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM 

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells_all) {

        roM[ic]   = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic]  = roe[ic];
        roKM[ic] = roK[ic];
        roOmegaM[ic] = roOmega[ic];
    }
}


void updateVariablesInner_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    updateVariablesInner_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"]  , var.c_d["roe"]  , var.c_d["roK"] , var.c_d["roOmega"] ,
        var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"] , var.c_d["roeM"] , var.c_d["roKM"] , var.c_d["roOmegaM"] 
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

__global__ void applyScalarImplicitCorrection_d
(
 geom_int nCells_all , geom_int nCells,

 flow_float* ro,
 flow_float* roUx,
 flow_float* roUy,
 flow_float* roUz,
 flow_float* roe,

 flow_float* roN,
 flow_float* roUxN,
 flow_float* roUyN,
 flow_float* roUzN,
 flow_float* roeN,

 flow_float* dq_ro,
 flow_float* dq_roUx,
 flow_float* dq_roUy,
 flow_float* dq_roUz,
 flow_float* dq_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        ro[ic] = roN[ic] + dq_ro[ic];
        roUx[ic] = roUxN[ic] + dq_roUx[ic];
        roUy[ic] = roUyN[ic] + dq_roUy[ic];
        roUz[ic] = roUzN[ic] + dq_roUz[ic];
        roe[ic] = roeN[ic] + dq_roe[ic];
    }
}

void applyScalarImplicitCorrection_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    applyScalarImplicitCorrection_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells_all , msh.nCells,
        var.c_d["ro"],
        var.c_d["roUx"],
        var.c_d["roUy"],
        var.c_d["roUz"],
        var.c_d["roe"],
        var.c_d["roN"],
        var.c_d["roUxN"],
        var.c_d["roUyN"],
        var.c_d["roUzN"],
        var.c_d["roeN"],
        var.c_d["dq_ro_old"],
        var.c_d["dq_roUx_old"],
        var.c_d["dq_roUy_old"],
        var.c_d["dq_roUz_old"],
        var.c_d["dq_roe_old"]
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}