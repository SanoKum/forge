#include "convectiveFlux_d.cuh"

__global__ void viscousFlux_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 flow_float mu ,  flow_float thermCond,

 // variables
//flow_float* convx , flow_float* convy , flow_float* convz,
// flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* ro   ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* Ps  ,
 flow_float* Ht  ,
 flow_float* sonic,
 
 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    if (ip < nPlanes) {

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float f = fx[ip];
        
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        flow_float ccx_0 = ccx[ic0];
        flow_float ccy_0 = ccy[ic0];
        flow_float ccz_0 = ccz[ic0];

        flow_float ccx_1 = ccx[ic1];
        flow_float ccy_1 = ccy[ic1];
        flow_float ccz_1 = ccz[ic1];

        flow_float dcc_x = ccx_1 - ccx_0;
        flow_float dcc_y = ccy_1 - ccy_0;
        flow_float dcc_z = ccz_1 - ccz_0;
        flow_float dcc   = sqrt(dcc_x*dcc_x +dcc_y*dcc_y +dcc_z*dcc_z) ;

        flow_float Uxf = f*Ux[ic0] + (1.0-f)*Ux[ic1];
        flow_float Uyf = f*Uy[ic0] + (1.0-f)*Uy[ic1];
        flow_float Uzf = f*Uz[ic0] + (1.0-f)*Uz[ic1];

        flow_float dUxdxf = f*dUxdx[ic0] + (1.0-f)*dUxdx[ic1];
        flow_float dUxdyf = f*dUxdy[ic0] + (1.0-f)*dUxdy[ic1];
        flow_float dUxdzf = f*dUxdz[ic0] + (1.0-f)*dUxdz[ic1];

        flow_float dUydxf = f*dUydx[ic0] + (1.0-f)*dUydx[ic1];
        flow_float dUydyf = f*dUydy[ic0] + (1.0-f)*dUydy[ic1];
        flow_float dUydzf = f*dUydz[ic0] + (1.0-f)*dUydz[ic1];

        flow_float dUzdxf = f*dUzdx[ic0] + (1.0-f)*dUzdx[ic1];
        flow_float dUzdyf = f*dUzdy[ic0] + (1.0-f)*dUzdy[ic1];
        flow_float dUzdzf = f*dUzdz[ic0] + (1.0-f)*dUzdz[ic1];

        flow_float dTdxf = f*dTdx[ic0] + (1.0-f)*dTdx[ic1];
        flow_float dTdyf = f*dTdy[ic0] + (1.0-f)*dTdy[ic1];
        flow_float dTdzf = f*dTdz[ic0] + (1.0-f)*dTdz[ic1];


        flow_float alpha = sss*sss/(dcc_x*sxx +dcc_y*syy +dcc_z*szz); // over relaxed

        flow_float tau_x = mu*((Ux[ic1] -Ux[ic0])/dcc)*alpha*dcc;
        tau_x += mu*(dUxdxf*(sxx-alpha*dcc_x) +dUydxf*(syy-alpha*dcc_y) +dUzdxf*(szz-alpha*dcc_z));
        tau_x += mu*(dUxdxf*sxx + dUydxf*syy + dUzdxf*szz);
        tau_x += -mu*2.0/3.0*(dUxdxf+dUydyf+dUzdzf)*sxx;

        flow_float tau_y = mu*((Uy[ic1] -Uy[ic0])/dcc)*alpha*dcc;
        tau_y += mu*(dUxdyf*(sxx-alpha*dcc_x) +dUydyf*(syy-alpha*dcc_y) +dUzdyf*(szz-alpha*dcc_z));
        tau_y += mu*(dUxdyf*sxx + dUydyf*syy + dUzdyf*szz);
        tau_y += -mu*2.0/3.0*(dUxdxf+dUydyf+dUzdzf)*syy;

        flow_float tau_z = mu*((Uz[ic1] -Uz[ic0])/dcc)*alpha*dcc;
        tau_z += mu*(dUxdzf*(sxx-alpha*dcc_x) +dUydzf*(syy-alpha*dcc_y) +dUzdzf*(szz-alpha*dcc_z));
        tau_z += mu*(dUxdzf*sxx + dUydzf*syy + dUzdzf*szz);
        tau_z += -mu*2.0/3.0*(dUxdxf+dUydyf+dUzdzf)*szz;

        flow_float res_ro_temp   = 0.0;
        flow_float res_roUx_temp = tau_x;
        flow_float res_roUy_temp = tau_y;
        flow_float res_roUz_temp = tau_z;
        flow_float res_roe_temp  = tau_x*Uxf +tau_y*Uyf +tau_z*Uzf; 
        res_roe_temp += thermCond*(dTdxf*sxx +dTdyf*syy +dTdzf*szz);

        atomicAdd(&res_ro[ic0]  , res_ro_temp);
        atomicAdd(&res_roUx[ic0], res_roUx_temp);
        atomicAdd(&res_roUy[ic0], res_roUy_temp);
        atomicAdd(&res_roUz[ic0], res_roUz_temp);
        atomicAdd(&res_roe[ic0] , res_roe_temp);

        atomicAdd(&res_ro[ic1]  , -res_ro_temp);
        atomicAdd(&res_roUx[ic1], -res_roUx_temp);
        atomicAdd(&res_roUy[ic1], -res_roUy_temp);
        atomicAdd(&res_roUz[ic1], -res_roUz_temp);
        atomicAdd(&res_roe[ic1] , -res_roe_temp);
    }

    __syncthreads();
}


void viscousFlux_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    // ------------------------------
    // *** sum over normal planes ***
    // ------------------------------
    viscousFlux_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        cfg.visc , cfg.thermCond,

        // basic variables
        //var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
        //var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"]  , 
        var.c_d["Uy"]  , 
        var.c_d["Uz"]  , 
        var.c_d["P"]  , 
        var.c_d["Ht"]  , 
        var.c_d["sonic"]  , 

        var.c_d["res_ro"] ,
        var.c_d["res_roUx"] ,
        var.c_d["res_roUy"] ,
        var.c_d["res_roUz"] ,
        var.c_d["res_roe"]  ,
       
        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],

        var.c_d["dTdx"] , var.c_d["dTdy"] , var.c_d["dTdz"]
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}