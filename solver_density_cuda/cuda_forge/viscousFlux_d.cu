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
 flow_float* vis_lam   , flow_float* vis_turb  ,

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
 flow_float* Ts  ,
 
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


    //if (ip < nPlanes) { 
    if (ip < nNormalPlanes) {

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

        flow_float delta   = dcc*sss*sss/(dcc_x*sxx +dcc_y*syy +dcc_z*szz); // over relaxed
        flow_float delta_x = dcc_x*sss*sss/(dcc_x*sxx +dcc_y*syy +dcc_z*szz); 
        flow_float delta_y = dcc_y*sss*sss/(dcc_x*sxx +dcc_y*syy +dcc_z*szz); 
        flow_float delta_z = dcc_z*sss*sss/(dcc_x*sxx +dcc_y*syy +dcc_z*szz); 
        flow_float k_x = sxx - delta_x; 
        flow_float k_y = syy - delta_y; 
        flow_float k_z = szz - delta_z; 
        flow_float divu = dUxdxf+dUydyf+dUzdzf;

        flow_float v_lam  = f*vis_lam [ic0] + (1.0-f)*vis_lam [ic1] ;
        flow_float v_turb = f*vis_turb[ic0] + (1.0-f)*vis_turb[ic1] ;
        flow_float mu_total = v_lam + v_turb;

        flow_float tau_x = mu_total*((Ux[ic1] -Ux[ic0])/dcc)*delta_x;
        //flow_float tau_x = mu_total*((Ux[ic1] -Ux[ic0])/dcc)*sxx;
        tau_x += mu_total*(dUxdxf*k_x +dUydxf*k_y +dUzdxf*k_z);
        //tau_x += mu_total*(dUxdxf*sxx + dUydxf*syy + dUzdxf*szz);
        tau_x += -mu_total*2.0/3.0*(divu)*sxx;

        flow_float tau_y = mu_total*((Uy[ic1] -Uy[ic0])/dcc)*delta_y;
        //flow_float tau_y = mu_total*((Uy[ic1] -Uy[ic0])/dcc)*syy;
        tau_y += mu_total*(dUxdyf*k_x +dUydyf*k_y +dUzdyf*k_z);
        //tau_y += mu_total*(dUxdyf*sxx + dUydyf*syy + dUzdyf*szz);
        tau_y += -mu_total*2.0/3.0*(divu)*syy;

        flow_float tau_z = mu_total*((Uz[ic1] -Uz[ic0])/dcc)*delta_z;
        //flow_float tau_z = mu_total*((Uz[ic1] -Uz[ic0])/dcc)*szz;
        tau_z += mu_total*(dUxdzf*k_x +dUydzf*k_y +dUzdzf*k_z);
        //tau_z += mu_total*(dUxdzf*sxx + dUydzf*syy + dUzdzf*szz);
        tau_z += -mu_total*2.0/3.0*(divu)*szz;

        flow_float heatflux = thermCond*((Ts[ic1] -Ts[ic0])/dcc)*delta;
        heatflux += thermCond*(dTdxf*k_x +dTdyf*k_y +dTdzf*k_z);

        flow_float res_ro_temp   = 0.0;
        flow_float res_roUx_temp = tau_x;
        flow_float res_roUy_temp = tau_y;
        flow_float res_roUz_temp = tau_z;
        flow_float res_roe_temp  = tau_x*Uxf +tau_y*Uyf +tau_z*Uzf; 
        res_roe_temp += heatflux;

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


__global__ void viscousFlux_wall_d
( 
  // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 flow_float mu ,  flow_float thermCond,
 flow_float* vis_lam   , flow_float* vis_turb  ,

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
 flow_float* Ts  ,
 
 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 // bvar
 flow_float* ypls_b ,
 flow_float* twall_x_b ,
 flow_float* twall_y_b ,
 flow_float* twall_z_b ,
 flow_float* Ux_b ,
 flow_float* Uy_b ,
 flow_float* Uz_b ,
 flow_float* Ts_b 
//flow_float* sx_b ,
 //flow_float* sy_b ,
 //flow_float* sz_b 

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

        flow_float ccx_0 = ccx[ic];
        flow_float ccy_0 = ccy[ic];
        flow_float ccz_0 = ccz[ic];

        flow_float pcx_1 = pcx[ip];
        flow_float pcy_1 = pcy[ip];
        flow_float pcz_1 = pcz[ip];

        flow_float ccx_1 = ccx[ig];
        flow_float ccy_1 = ccy[ig];
        flow_float ccz_1 = ccz[ig];

        //flow_float dcc_x = pcx_1 - ccx_0;
        //flow_float dcc_y = pcy_1 - ccy_0;
        //flow_float dcc_z = pcz_1 - ccz_0;

        flow_float dcc_x = ccx_1 - ccx_0;
        flow_float dcc_y = ccy_1 - ccy_0;
        flow_float dcc_z = ccz_1 - ccz_0;
 
        flow_float dcc   = sqrt(dcc_x*dcc_x +dcc_y*dcc_y +dcc_z*dcc_z) ;

        flow_float dUxdxf = dUxdx[ic] ;
        flow_float dUxdyf = dUxdy[ic] ;
        flow_float dUxdzf = dUxdz[ic] ;

        flow_float dUydxf = dUydx[ic] ;
        flow_float dUydyf = dUydy[ic] ;
        flow_float dUydzf = dUydz[ic] ;

        flow_float dUzdxf = dUzdx[ic] ;
        flow_float dUzdyf = dUzdy[ic] ;
        flow_float dUzdzf = dUzdz[ic] ;

        flow_float v_turb = vis_turb[ic] ;
        flow_float mu_total = vis_lam[ic] + v_turb;

        flow_float divu = dUxdxf+dUydyf+dUzdzf;

        //flow_float tau_x = mu_total*((Ux_b[ib] - Ux[ic])/dcc)*sxx;
        flow_float tau_x = mu_total*((Ux[ig] - Ux[ic])/dcc)*sxx;
        //tau_x += mu_total*(dUxdxf*sxx + dUydxf*syy + dUzdxf*szz);
        tau_x += -mu_total*2.0/3.0*(divu)*sxx;

        flow_float tau_y = mu_total*((Uy[ig] - Uy[ic])/dcc)*syy;
        //tau_y += mu_total*(dUxdyf*sxx + dUydyf*syy + dUzdyf*szz);
        tau_y += -mu_total*2.0/3.0*(divu)*syy;

        flow_float tau_z = mu_total*((Uz[ig] - Uz[ic])/dcc)*szz;
        //tau_z += mu_total*(dUxdzf*sxx + dUydzf*syy + dUzdzf*szz);
        tau_z += -mu_total*2.0/3.0*(divu)*szz;

        flow_float heatflux = thermCond*((Ts[ig]- Ts[ic])/dcc)*sss;

        flow_float res_ro_temp   = 0.0;
        flow_float res_roUx_temp = tau_x;
        flow_float res_roUy_temp = tau_y;
        flow_float res_roUz_temp = tau_z;
        flow_float res_roe_temp  = tau_x*Ux_b[ib] +tau_y*Uy_b[ib] +tau_z*Uz_b[ib]; 

        //printf("res_roe_temp=%e\n", res_roe_temp);
        //printf("Uyb[ib]=%e\n", Uy_b[ib]);
        //printf("Uzb[ib]=%e\n", Uz_b[ib]);

        res_roe_temp += heatflux;

        //res_roe_temp = 0.0;

        atomicAdd(&res_ro[ic]  , res_ro_temp);
        atomicAdd(&res_roUx[ic], res_roUx_temp);
        atomicAdd(&res_roUy[ic], res_roUy_temp);
        atomicAdd(&res_roUz[ic], res_roUz_temp);
        atomicAdd(&res_roe[ic] , res_roe_temp);

        twall_x_b[ib] = tau_x/sss;
        twall_y_b[ib] = tau_y/sss;
        twall_z_b[ib] = tau_z/sss;

        flow_float twall = sqrt(tau_x*tau_x + tau_y*tau_y + tau_z*tau_z)/sss;
        flow_float utau = sqrt(twall/ro[ic]);

        ypls_b[ib] = ro[ic]*utau*dcc/mu_total;
        //if (ib == 1) {
        //    printf("ib = %d\n", ib);
        //    printf("ip = %d\n", ip);
        //    printf("tau_x = %e, tau_y = %e, tau_z = %e\n", tau_x, tau_y, tau_z);
        //    printf("c Ux = %e, Uy = %e, Uz = %e\n", Ux[ic], Uy[ic], Uz[ic]);
        //    printf("g Ux = %e, Uy = %e, Uz = %e\n", Ux[ig], Uy[ig], Uz[ig]);
        //    printf("b Ux = %e, Uy = %e, Uz = %e\n", Ux_b[ib], Uy_b[ib], Uz_b[ib]);
        //}
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
        var.c_d["vis_lam"], var.c_d["vis_turb"],

        // basic variables
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
        var.c_d["T"]  , 

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

    for (auto& bc : msh.bconds)
    {
        if (bc.bcondKind == "wall" or bc.bcondKind == "wall_isothermal") {
            viscousFlux_wall_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
                // mesh structure
                bc.iPlanes.size(),
                bc.map_bplane_plane_d,  
                bc.map_bplane_cell_d,  
                bc.map_bplane_cell_ghst_d,

                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

                cfg.visc , cfg.thermCond,
                var.c_d["vis_lam"], var.c_d["vis_turb"],

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
                var.c_d["T"]  , 

                var.c_d["res_ro"] ,
                var.c_d["res_roUx"] ,
                var.c_d["res_roUy"] ,
                var.c_d["res_roUz"] ,
                var.c_d["res_roe"]  ,
       
                // gradient
                var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
                var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
                var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
                var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],

                bc.bvar_d["ypls"],
                bc.bvar_d["twall_x"],
                bc.bvar_d["twall_y"],
                bc.bvar_d["twall_z"],
                bc.bvar_d["Ux"],
                bc.bvar_d["Uy"],
                bc.bvar_d["Uz"],
                bc.bvar_d["Ts"]
 
                //bc.bvar_d["sx"],
                //bc.bvar_d["sy"],
                //bc.bvar_d["sz"]
            ) ;
        }
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}