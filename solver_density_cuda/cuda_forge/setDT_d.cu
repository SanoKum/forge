//#if defined(__CUDA_ARCH__)
//#pragma message "content of __CUDA_ARCH__: " __CUDA_ARCH__
//#endif

#include "setDT_d.cuh"


__global__ void setDT_d
( 
 flow_float dt,
 flow_float dt_pseudo,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* cfl ,
 flow_float* cfl_pseudo ,
 flow_float* sonic,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    if (ip < nPlanes) {

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        //__syncthreads();

        geom_float f = fx[ip];

        geom_float vol0 = vol[ic0];
        geom_float vol1 = vol[ic1];
        
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];


        geom_float dx0 = vol0/sss;
        geom_float dx1 = vol0/sss;

        flow_float Ux0 = Ux[ic0];
        flow_float Uy0 = Uy[ic0];
        flow_float Uz0 = Uz[ic0];

        flow_float Ux1 = Ux[ic1];
        flow_float Uy1 = Uy[ic1];
        flow_float Uz1 = Uz[ic1];

        flow_float US  = (fx[ip]*Ux0 + (1.0-fx[ip])*Ux1)*sxx
                        +(fx[ip]*Uy0 + (1.0-fx[ip])*Uy1)*syy
                        +(fx[ip]*Uz0 + (1.0-fx[ip])*Uz1)*szz;
        
        //flow_float cfl_temp0 = dt*(abs(US)/sss+sonic[ic0])/dx0;
        //flow_float cfl_temp1 = dt*(abs(US)/sss+sonic[ic1])/dx1;

        cfl[ic0] = max(cfl[ic0] , dt*(abs(US)/sss+sonic[ic0])/dx0);
        cfl[ic1] = max(cfl[ic1] , dt*(abs(US)/sss+sonic[ic1])/dx1);

        cfl_pseudo[ic0] = max(cfl_pseudo[ic0] , dt_pseudo*(abs(US)/sss+sonic[ic0])/dx0);
        cfl_pseudo[ic1] = max(cfl_pseudo[ic1] , dt_pseudo*(abs(US)/sss+sonic[ic1])/dx1);        

        //atomicMax(&cfl[ic0]  , dt*(abs(US)/sss+sonic[ic0])/dx0);
        //atomicMax(&cfl[ic1]  , dt*(abs(US)/sss+sonic[ic1])/dx1);

        //atomicMax(&cfl_pseudo[ic0]  , dt_pseudo*(abs(US)/sss+sonic[ic0])/dx0);
        //atomicMax(&cfl_pseudo[ic1]  , dt_pseudo*(abs(US)/sss+sonic[ic1])/dx1);
    }

}


void setDT_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    setDT_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        cfg.dt,
        cfg.dt_pseudo,

        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        //var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
        //var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["cfl"]  , 
        var.c_d["cfl_pseudo"]  , 
        var.c_d["sonic"]  , 
        var.c_d["Ux"]  , 
        var.c_d["Uy"]  , 
        var.c_d["Uz"]  
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}