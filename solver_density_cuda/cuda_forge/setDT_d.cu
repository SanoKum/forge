//#if defined(__CUDA_ARCH__)
//#pragma message "content of __CUDA_ARCH__: " __CUDA_ARCH__
//#endif

#include "setDT_d.cuh"


__global__ void setCFL_pln_d
( 
 flow_float dt,
 flow_float dt_pseudo,
 flow_float visc,

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
 flow_float* Uz  ,

 //plane variables
 flow_float* cfl_pln ,
 flow_float* cfl_pseudo_pln

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
        geom_float dx1 = vol1/sss;
        geom_float dx_min = min(dx0,dx1);

        flow_float Ux0 = Ux[ic0];
        flow_float Uy0 = Uy[ic0];
        flow_float Uz0 = Uz[ic0];

        flow_float Ux1 = Ux[ic1];
        flow_float Uy1 = Uy[ic1];
        flow_float Uz1 = Uz[ic1];

        flow_float US  = (f*Ux0 + (1.0-f)*Ux1)*sxx
                        +(f*Uy0 + (1.0-f)*Uy1)*syy
                        +(f*Uz0 + (1.0-f)*Uz1)*szz;

        flow_float rof = f*ro[ic0] + (1.0-f)*ro[ic1];
        flow_float lambda = abs(US)/sss + sonic[ic0] + 2.0*visc/(rof*dx_min);

        cfl_pln[ip] = dt*lambda/dx_min;

        cfl_pseudo_pln[ip] = dt_pseudo*lambda/dx_min;
    }
}

__global__ void setCFL_cell_d
( 
 int dtControl, flow_float cfg_target,
 flow_float dt,
 flow_float dt_pseudo,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* cell_planes_index, geom_int* cell_planes,  
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
 flow_float* Uz  ,

 //plane variables
 flow_float* cfl_pln ,
 flow_float* cfl_pseudo_pln

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {

        geom_int index_st = cell_planes_index[ic];
        geom_int index_en = cell_planes_index[ic+1];
        geom_int np = index_en - index_st;

        cfl[ic] = 0.0;
        cfl_pseudo[ic] = 0.0;

        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            geom_int ip = cell_planes[ilp];
            
            cfl[ic]        = max(cfl[ic], cfl_pln[ip]);
            cfl_pseudo[ic] = max(cfl_pseudo[ic], cfl_pseudo_pln[ip]);
        }

        //if (dtControl == 1){
        //    dt = 
        //}
    }
}


void setDT_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    setCFL_pln_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        cfg.dt,
        cfg.dt_pseudo,
        cfg.visc,

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
        var.c_d["cfl_pseudo"] , 
        var.c_d["sonic"]  , 
        var.c_d["Ux"]  , 
        var.c_d["Uy"]  , 
        var.c_d["Uz"]  ,

        var.p_d["cfl_pln"]  , 
        var.p_d["cfl_pseudo_pln"]  
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    setCFL_cell_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        cfg.dtControl, cfg.cfl,
        cfg.dt,
        cfg.dt_pseudo,

        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes, msh.map_cell_planes_index_d , msh.map_cell_planes_d,
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
        var.c_d["cfl_pseudo"] , 
        var.c_d["sonic"]  , 
        var.c_d["Ux"]  , 
        var.c_d["Uy"]  , 
        var.c_d["Uz"]  ,

        var.p_d["cfl_pln"]  , 
        var.p_d["cfl_pseudo_pln"]  
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    //cublasHandle_t handle;
    //cublasStatus_t stat;
    //cublasCreate(&handle);
    //int max_ic;
    //stat = cublasIsamax(handle, msh.nCells, var.c_d["cfl"], 1, &max_ic);
    //if (stat != CUBLAS_STATUS_SUCCESS) {
    //    printf("Max failed\n");
    //}

    flow_float cfl_max;

    thrust::device_ptr<flow_float> d_ptr = thrust::device_pointer_cast(var.c_d["cfl"]);
    cfl_max = *(thrust::max_element(d_ptr, d_ptr + msh.nCells));

    if (cfg.dtControl == 1) { // cfl based time control
        flow_float cfl_target = cfg.cfl;
        cfg.dt = cfg.dt*cfl_target/cfl_max;

        cfg.dt = max(cfg.dt, cfg.dt_min);
        cfg.dt = min(cfg.dt, cfg.dt_max);
    }

    printf("  max cfl : %f    \n", cfl_max);
    printf("  dt      : %e [s]\n", cfg.dt);

}