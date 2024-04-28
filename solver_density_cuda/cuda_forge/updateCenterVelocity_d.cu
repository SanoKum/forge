#include "updateCenterVelocity_d.cuh"

__global__ void updateCenterVelocity_d
( 
 flow_float dt,
 // mesh structure
 geom_int nCells,
 geom_float* vol , 

 // variables
 flow_float* convx , flow_float* convy , flow_float* convz,
 flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* Ux  , flow_float* UxN ,
 flow_float* Uy  , flow_float* UyN ,
 flow_float* Uz  , flow_float* UzN ,
 flow_float* ro  ,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        __syncthreads();
        geom_float volume = vol[ic];

        Ux[ic] = UxN[ic] -dPdx[ic]/ro[ic]*dt + (-convx[ic]+diffx[ic])*dt/volume/ro[ic];
        Uy[ic] = UyN[ic] -dPdy[ic]/ro[ic]*dt + (-convy[ic]+diffy[ic])*dt/volume/ro[ic];
        Uz[ic] = UzN[ic] -dPdz[ic]/ro[ic]*dt + (-convz[ic]+diffz[ic])*dt/volume/ro[ic];

        //printf("%d UxN=%f, dt=%f, vol=%f, ro=%f, convx=%f, diffx=%f\n", ic, UxN[ic], dt, volume, ro[ic], convx[ic], diffx[ic]);
        __syncthreads();
    }
}

void updateCenterVelocity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    // ------------------------------
    // *** sum over normal planes ***
    // ------------------------------
    updateCenterVelocity_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        cfg.dt, 
        
        msh.nCells,
        var.c_d["volume"], 

        // basic variables
        var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
        var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
        var.c_d["Ux"] , var.c_d["UxN"] ,
        var.c_d["Uy"] , var.c_d["UyN"] ,
        var.c_d["Uz"] , var.c_d["UzN"] ,
        var.c_d["ro"] ,

        // gradient
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"]
    ) ;

}