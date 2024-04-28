#include "convectiveFlux_d.cuh"
#include "AUSM_d.cuh"

__device__ flow_float betaPls(flow_float M, flow_float alpha) // ok
{
    if (abs(M) >= 1.0) {
        return 0.5*(1.0+sign_sano(M));
    } else {
        return 0.25*pow(M+1.0,2.0)*(2.0-M) +alpha*M*pow(M*M-1, 2.0);
    }
}

__device__ flow_float betaMns(flow_float M, flow_float alpha) // ok
{
    if (abs(M) >= 1.0) {
        return 0.5*(1.0-sign_sano(M));
    } else {
        return 0.25*pow(M-1.0,2.0)*(2.0+M) -alpha*M*pow(M*M-1, 2.0);
    }
}

__device__ flow_float MPls(flow_float M) // ok
{
    flow_float beta = 1.0/8.0;
    if (abs(M) >= 1.0) {
        return 0.5*(M+abs(M));
    } else {
        return +0.25*pow(M+1.0,2.0) +beta*pow(M*M-1.0, 2.0);
    }
}

__device__ flow_float MMns(flow_float M) // ok
{
    flow_float beta = 1.0/8.0;
    if (abs(M) >= 1.0) {
        return 0.5*(M-abs(M));
    } else {
        return -0.25*pow(M-1.0,2.0) -beta*pow(M*M-1.0, 2.0);
    }
}

__global__ void AUSMp_d
( 
 flow_float ga,

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
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz

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

        flow_float ccx_0 = ccx[ic0];
        flow_float ccy_0 = ccy[ic0];
        flow_float ccz_0 = ccz[ic0];

        flow_float ccx_1 = ccx[ic1];
        flow_float ccy_1 = ccy[ic1];
        flow_float ccz_1 = ccz[ic1];

        flow_float dcc_x = ccx_1 - ccx_0;
        flow_float dcc_y = ccy_1 - ccy_0;
        flow_float dcc_z = ccz_1 - ccz_0;

        flow_float ro_L= MUSCL(ro[ic0], ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], dcc_x, dcc_y, dcc_z);
        flow_float Ux_L= MUSCL(Ux[ic0], Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dcc_x, dcc_y, dcc_z);
        flow_float Uy_L= MUSCL(Uy[ic0], Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dcc_x, dcc_y, dcc_z);
        flow_float Uz_L= MUSCL(Uz[ic0], Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dcc_x, dcc_y, dcc_z);
        flow_float P_L = MUSCL(Ps[ic0], Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dcc_x, dcc_y, dcc_z);
        flow_float Ht_L= MUSCL(Ht[ic0], Ht[ic1], dHtdx[ic0], dHtdy[ic0], dHtdz[ic0], dcc_x, dcc_y, dcc_z);

        flow_float ro_R= MUSCL(ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], -dcc_x, -dcc_y, -dcc_z); 
        flow_float Ux_R= MUSCL(Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], -dcc_x, -dcc_y, -dcc_z);
        flow_float Uy_R= MUSCL(Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], -dcc_x, -dcc_y, -dcc_z);
        flow_float Uz_R= MUSCL(Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], -dcc_x, -dcc_y, -dcc_z);
        flow_float P_R = MUSCL(Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , -dcc_x, -dcc_y, -dcc_z);
        flow_float Ht_R= MUSCL(Ht[ic1], Ht[ic0], dHtdx[ic1], dHtdy[ic1], dHtdz[ic1], -dcc_x, -dcc_y, -dcc_z);

        flow_float U_L = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float U_R = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        //flow_float c_L = sqrt(ga*P_L/ro_L);
        //flow_float c_R = sqrt(ga*P_R/ro_R);
        //flow_float c_half = 0.5*(c_L+c_R);

        Ht_L = ga/(ga-1.0)*P_L/ro_L + 0.5*(Ux_L*Ux_L + Uy_L*Uy_L + Uz_L*Uz_L);
        Ht_R = ga/(ga-1.0)*P_R/ro_R + 0.5*(Ux_R*Ux_R + Uy_R*Uy_R + Uz_R*Uz_R);

        flow_float c_star_L = sqrt(2.0*(ga-1.0)/(ga+1.0)*Ht_L); //ok
        flow_float c_star_R = sqrt(2.0*(ga-1.0)/(ga+1.0)*Ht_R); //ok
        flow_float c_tilde_L = pow(c_star_L,2.0)/max(c_star_L, abs(U_L)); //ok
        flow_float c_tilde_R = pow(c_star_R,2.0)/max(c_star_R, abs(U_R)); //ok
        flow_float c_half = min(c_tilde_L, c_tilde_R); //ok

        flow_float M_L = U_L/c_half;
        flow_float M_R = U_R/c_half;

        flow_float M_half = MPls(M_L) + MMns(M_R);
        flow_float alpha = 3.0/16.0;
        flow_float p_tilde = betaPls(M_L,alpha)*P_L + betaMns(M_R,alpha)*P_R;

        flow_float m;
        if (M_half > 0) {
            m = M_half*c_half*ro_L;
        } else {
            m = M_half*c_half*ro_R;
        }

        flow_float res_ro_temp   = m*sss;
        flow_float res_roUx_temp = (0.5*(m+abs(m))*Ux_L +0.5*(m-abs(m))*Ux_R)*sss +p_tilde*sxx;
        flow_float res_roUy_temp = (0.5*(m+abs(m))*Uy_L +0.5*(m-abs(m))*Uy_R)*sss +p_tilde*syy;
        flow_float res_roUz_temp = (0.5*(m+abs(m))*Uz_L +0.5*(m-abs(m))*Uz_R)*sss +p_tilde*szz;
        flow_float res_roe_temp  = (0.5*(m+abs(m))*Ht_L +0.5*(m-abs(m))*Ht_R)*sss;

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


__global__ void AUSMp_UP_d
( 
 flow_float ga,

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
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz

)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    if (ip < nPlanes) {

        flow_float Kp = 0.25;
        flow_float Ku = 0.75;
        flow_float sig= 1.0;
        flow_float Mco= 0.1;

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        //__syncthreads();

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

        flow_float ro_L= MUSCL(ro[ic0], ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], dcc_x, dcc_y, dcc_z);
        flow_float Ux_L= MUSCL(Ux[ic0], Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dcc_x, dcc_y, dcc_z);
        flow_float Uy_L= MUSCL(Uy[ic0], Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dcc_x, dcc_y, dcc_z);
        flow_float Uz_L= MUSCL(Uz[ic0], Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dcc_x, dcc_y, dcc_z);
        flow_float P_L = MUSCL(Ps[ic0], Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dcc_x, dcc_y, dcc_z);
        flow_float Ht_L= MUSCL(Ht[ic0], Ht[ic1], dHtdx[ic0], dHtdy[ic0], dHtdz[ic0], dcc_x, dcc_y, dcc_z);

        flow_float ro_R= MUSCL(ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], -dcc_x, -dcc_y, -dcc_z); 
        flow_float Ux_R= MUSCL(Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], -dcc_x, -dcc_y, -dcc_z);
        flow_float Uy_R= MUSCL(Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], -dcc_x, -dcc_y, -dcc_z);
        flow_float Uz_R= MUSCL(Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], -dcc_x, -dcc_y, -dcc_z);
        flow_float P_R = MUSCL(Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , -dcc_x, -dcc_y, -dcc_z);
        flow_float Ht_R= MUSCL(Ht[ic1], Ht[ic0], dHtdx[ic1], dHtdy[ic1], dHtdz[ic1], -dcc_x, -dcc_y, -dcc_z);

        flow_float U_L = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float U_R = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        flow_float c_L = sqrt(ga*P_L/ro_L);
        flow_float c_R = sqrt(ga*P_R/ro_R);
        flow_float c_half = 0.5*(c_L+c_R);

        flow_float M_bar = sqrt(U_L*U_L + U_R*U_R)/(2.0*c_half*c_half);
        flow_float M_o = sqrt(min(1.0, max(M_bar*M_bar, Mco*Mco)));
        flow_float fa = M_o*(2.0-M_o);
        flow_float alpha = 3.0/16.0*(-4.0+5.0*fa*fa);

        flow_float M_L = U_L/c_half;
        flow_float M_R = U_R/c_half;

        flow_float p_u = -Ku*betaPls(M_L,alpha)*betaMns(M_R, alpha)*(ro_L+ro_R)*(fa*c_half)*(U_L-U_R);
        flow_float p_tilde = betaPls(M_L,alpha)*P_L + betaMns(M_R,alpha)*P_R + p_u;

        flow_float ro_half = 0.5*(ro_L + ro_R);

        flow_float M_p = -Kp/fa*max(1.0-sig*M_bar*M_bar, 0.0)*(P_R-P_L)/(ro_half*c_half*c_half);

        flow_float M_half = MPls(M_L) + MMns(M_R) + M_p;


        flow_float m;
        if (M_half > 0) {
            m = M_half*c_half*ro_L;
        } else {
            m = M_half*c_half*ro_R;
        }

//        if (ip==2) {
//            printf("m=%f\n", m);
//            printf("M_p=%f\n", M_p);
//            printf("M_half=%f\n", M_half);
//            printf("M_L=%f\n", M_L);
//            printf("M_R=%f\n", M_R);
//            printf("MPls(M_L)=%f\n", MPls(M_L));
//            printf("MMns(M_L)=%f\n", MMns(M_L));
//            printf("MPls(M_R)=%f\n", MPls(M_R));
//            printf("MMns(M_R)=%f\n", MMns(M_R));
//            printf("Ux_L=%f\n", Ux_L);
//            printf("Uy_L=%f\n", Uy_L);
//            printf("Uz_L=%f\n", Uz_L);
//            printf("Ht_L=%f\n", Ht_L);
//            printf("Ht_R=%f\n", Ht_R);
//            printf("P_L=%f\n", P_L);
//            printf("P_R=%f\n", P_R);
//            printf("c_L=%f\n", c_L);
//            printf("c_R=%f\n", c_R);
//            printf("ro_L=%f\n", ro_L);
//            printf("ro_R=%f\n", ro_R);
//
//        }

        flow_float res_ro_temp   = m*sss;
        flow_float res_roUx_temp = (0.5*(m+abs(m))*Ux_L +0.5*(m-abs(m))*Ux_R)*sss +p_tilde*sxx;
        flow_float res_roUy_temp = (0.5*(m+abs(m))*Uy_L +0.5*(m-abs(m))*Uy_R)*sss +p_tilde*syy;
        flow_float res_roUz_temp = (0.5*(m+abs(m))*Uz_L +0.5*(m-abs(m))*Uz_R)*sss +p_tilde*szz;
        flow_float res_roe_temp  = (0.5*(m+abs(m))*Ht_L +0.5*(m-abs(m))*Ht_R)*sss;

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

