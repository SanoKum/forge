#include "cuda_nagare/boundaryCond_d.cuh"

__global__ 
void slip_d 
( 
 geom_int nb,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];

        geom_float Un =  (sx[ip]*Ux[ic] + sy[ip]*Uy[ic] + sz[ip]*Uz[ic])/ss[ip];

        Tp[ip]  = T[ic];
        Uxp[ip] = Ux[ic] - Un*sx[ip]/ss[ip];
        Uyp[ip] = Uy[ic] - Un*sy[ip]/ss[ip];
        Uzp[ip] = Uz[ic] - Un*sz[ip]/ss[ip];
        rop[ip] = ro[ic];
        Pp[ip]  = P[ic];
    }
};

void slip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    slip_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>( 
        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  bc.map_bplane_cell_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.p_d["T"]  , var.c_d["T"] ,
        var.p_d["Ux"] , var.c_d["Ux"] ,
        var.p_d["Uy"] , var.c_d["Uy"] ,
        var.p_d["Uz"] , var.c_d["Uz"] ,
        var.p_d["ro"] , var.c_d["ro"] ,
        var.p_d["P"]  , var.c_d["P"] 
    );
}

__global__ 
void wall_d 
( 
 geom_int nb,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,
 // bvar
 flow_float* Uxb ,  
 flow_float* Uyb ,  
 flow_float* Uzb   
)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];

        Uxp[ip] = Uxb[ib];
        Uyp[ip] = Uyb[ib];
        Uzp[ip] = Uzb[ib];
        
        Tp[ip]  = T[ic];
        rop[ip] = ro[ic];
        Pp[ip]  = P[ic];


        //printf("*** WALL ***\n");
        //printf("%d , %d , Uxb  = %f\n", ibl , ib , Uxb[ibl]);
        //printf("%d , %d , Uyb  = %f\n", ibl , ib , Uyb[ibl]);
        //printf("%d , %d , Uzb  = %f\n", ibl , ib , Uzb[ibl]);
        //printf("%d , %d , P[ic]= %f\n", ibl , ib , P[ic]);
        __syncthreads();
    }
};

void wall_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    wall_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  bc.map_bplane_cell_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.p_d["T"]  , var.c_d["T"] ,
        var.p_d["Ux"] , var.c_d["Ux"] ,
        var.p_d["Uy"] , var.c_d["Uy"] ,
        var.p_d["Uz"] , var.c_d["Uz"] ,
        var.p_d["ro"] , var.c_d["ro"] ,
        var.p_d["P"]  , var.c_d["P"]  , 

        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"]
    ) ;
}


__global__ 
void outlet_statPress_d 
( 
 geom_int nb,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,
 // bvar
 flow_float* Pb 
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];

        //printf("ip=%d, ib=%d, Pb=%f\n", ip, ib, Pb[ib]);

        Pp[ip]  = Pb[ib];

        rop[ip] = ro[ic];
        Tp [ip] = T [ic];
        Uxp[ip] = Ux[ic];
        Uyp[ip] = Uy[ic];
        Uzp[ip] = Uz[ic];
        
        __syncthreads();
    }
};

void outlet_statPress_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    outlet_statPress_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  bc.map_bplane_cell_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.p_d["T"]  , var.c_d["T"] ,
        var.p_d["Ux"] , var.c_d["Ux"] ,
        var.p_d["Uy"] , var.c_d["Uy"] ,
        var.p_d["Uz"] , var.c_d["Uz"] ,
        var.p_d["ro"] , var.c_d["ro"] ,
        var.p_d["P"]  , var.c_d["P"]  , 

        bc.bvar_d["P"]
        //bc.bvar_d["Uy"],
        //bc.bvar_d["Uz"]
    ) ;
}




__global__ 
void inlet_uniformVelocity_d
( 
 geom_int nb,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,
 // bvar
 flow_float* Uxb ,  
 flow_float* Uyb ,  
 flow_float* Uzb , 
 flow_float* Tb  , 
 flow_float* rob  
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];

        Uxp[ip] = Uxb[ib];
        Uyp[ip] = Uyb[ib];
        Uzp[ip] = Uzb[ib];
        Tp [ip] = Tb [ib];
        //rop[ip] = rob[ib];
        rop[ip] = ro[ic];
 
        Pp[ip]  = P[ic];

        //printf("*** inlet ***\n");
        //printf("%d , %d , Uxb  = %f\n", ib , Uxb[ib]);
        //printf("%d , %d , Uyb  = %f\n", ib , Uyb[ib]);
        //printf("%d , %d , Uzb  = %f\n", ib , Uzb[ib]);
        //printf("%d , %d , P[ic]= %f\n", ib , P[ic]);

        __syncthreads();
    }
};

void inlet_uniformVelocity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    inlet_uniformVelocity_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  bc.map_bplane_cell_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.p_d["T"]  , var.c_d["T"] ,
        var.p_d["Ux"] , var.c_d["Ux"] ,
        var.p_d["Uy"] , var.c_d["Uy"] ,
        var.p_d["Uz"] , var.c_d["Uz"] ,
        var.p_d["ro"] , var.c_d["ro"] ,
        var.p_d["P"]  , var.c_d["P"]  , 

        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["T"],
        bc.bvar_d["ro"]
    ) ;
}






