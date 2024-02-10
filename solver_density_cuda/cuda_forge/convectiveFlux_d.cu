#include "convectiveFlux_d.cuh"

__global__ void SLAU_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

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
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    if (ip < nPlanes) {

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        //__syncthreads();

        geom_float f = fx[ip];
        
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        flow_float ro_p= ro[ic0];
        flow_float ro_m= ro[ic1];

        flow_float ro_L= ro_p;
        flow_float ro_R= ro_m;

        flow_float u_p = Ux[ic0];
        flow_float v_p = Uy[ic0];
        flow_float w_p = Uz[ic0];
        flow_float h_p = Ht[ic0];
        flow_float u_m = Ux[ic1];
        flow_float v_m = Uy[ic1];
        flow_float w_m = Uz[ic1];
        flow_float h_m = Ht[ic1];


        flow_float P_p = Ps[ic0];
        flow_float P_m = Ps[ic1];
        flow_float Vn_p = ((Ux[ic0])*sxx +(Uy[ic0])*syy +(Uz[ic0])*szz)/sss;
        flow_float Vn_m = ((Ux[ic1])*sxx +(Uy[ic1])*syy +(Uz[ic1])*szz)/sss;

        flow_float roVn_p = (roUx[ic0]*sxx + roUy[ic0]*syy + roUz[ic0]*szz)/sss;
        flow_float roVn_m = (roUx[ic1]*sxx + roUy[ic1]*syy + roUz[ic1]*szz)/sss;

        flow_float VnL = ((Ux[ic0])*sxx +(Uy[ic0])*syy +(Uz[ic0])*szz)/sss;
        flow_float VnR = ((Ux[ic1])*sxx +(Uy[ic1])*syy +(Uz[ic1])*szz)/sss;

        flow_float VnL_abs = abs(VnL);
        flow_float VnR_abs = abs(VnR);

        flow_float c_hat = 0.5*(sonic[ic0] + sonic[ic1]);

        flow_float M_p = Vn_p/c_hat;
        flow_float M_m = Vn_m/c_hat;

        flow_float ro_del= ro_m - ro_p;
        flow_float Vn_del= Vn_m - Vn_p;
        flow_float P_del = P_m  - P_p;

        flow_float beta_p, beta_m;

        if (abs(M_p)>=1.0){
            beta_p = 0.5*(M_p + abs(M_p))/M_p;
        } else {
            beta_p = 0.25*pow((M_p+1.0),2.0)*(2.0-M_p);
        }

        if (abs(M_m)>=1.0){
            beta_m = 0.5*(M_m - abs(M_m))/M_m;
        } else {
            beta_m = 0.25*pow((M_m-1.0),2.0)*(2.0+M_m);
        }

        flow_float zero = 0.0;
        flow_float one  = 1.0;
        flow_float half = 0.5;

        flow_float g = -max(min(M_p,zero),-one)*min(max(M_m,zero),one);
        flow_float Vn_hat_abs   = (ro_p*abs(Vn_p) + ro_m*abs(Vn_m))/(ro_p + ro_m);
        flow_float Vn_hat_p_abs = (1.0-g)*Vn_hat_abs + g*VnL_abs;
        flow_float Vn_hat_m_abs = (1.0-g)*Vn_hat_abs + g*VnR_abs;

        flow_float M_hat = min(one, sqrt(half*(u_m*u_m+v_m*v_m+w_m*w_m +u_p*u_p+v_p*v_p+w_p*w_p))/c_hat);
        flow_float chi = (1.0-M_hat)*(1.0-M_hat);

        flow_float p_tilde = half*(P_p+P_m) +half*(beta_p-beta_m)*(P_p-P_m)
                             +(one-chi)*(beta_p+beta_m-one)*half*(P_m+P_p); 

        flow_float mdot = sss*0.5*((ro_L*(VnL+Vn_hat_p_abs)+ro_R*(VnR-Vn_hat_m_abs)) -chi/(c_hat)*P_del);

        flow_float res_ro_temp   = mdot;
        flow_float res_roUx_temp = 0.5*(mdot+abs(mdot))*u_p +0.5*(mdot-abs(mdot))*u_m +p_tilde*sxx;
        flow_float res_roUy_temp = 0.5*(mdot+abs(mdot))*v_p +0.5*(mdot-abs(mdot))*v_m +p_tilde*syy;
        flow_float res_roUz_temp = 0.5*(mdot+abs(mdot))*w_p +0.5*(mdot-abs(mdot))*w_m +p_tilde*szz;
        flow_float res_roe_temp  = 0.5*(mdot+abs(mdot))*h_p +0.5*(mdot-abs(mdot))*h_m ;

        atomicAdd(&res_ro[ic0]  , -res_ro_temp);
        atomicAdd(&res_roUx[ic0], -res_roUx_temp);
        atomicAdd(&res_roUy[ic0], -res_roUy_temp);
        atomicAdd(&res_roUz[ic0], -res_roUz_temp);
        atomicAdd(&res_roe[ic0] , -res_roe_temp);

        atomicAdd(&res_ro[ic1]  , res_ro_temp);
        atomicAdd(&res_roUx[ic1], res_roUx_temp);
        atomicAdd(&res_roUy[ic1], res_roUy_temp);
        atomicAdd(&res_roUz[ic1], res_roUz_temp);
        atomicAdd(&res_roe[ic1] , res_roe_temp);
    }

    __syncthreads();
}


void convectiveFlux_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    // initialize
    cudaMemset(var.c_d["res_ro"]  , 0.0, msh.nCells*sizeof(flow_float));
    cudaMemset(var.c_d["res_roUx"], 0.0, msh.nCells*sizeof(flow_float));
    cudaMemset(var.c_d["res_roUy"], 0.0, msh.nCells*sizeof(flow_float));
    cudaMemset(var.c_d["res_roUz"], 0.0, msh.nCells*sizeof(flow_float));
    cudaMemset(var.c_d["res_roe"] , 0.0, msh.nCells*sizeof(flow_float));

    // ------------------------------
    // *** sum over normal planes ***
    // ------------------------------
    SLAU_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
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
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"]
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}