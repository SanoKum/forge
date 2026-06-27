// =============================================================================
// convectiveFlux_d.cu の dispatch から到達不能なスキーム実装の退避先。
//   - wrapper の cfg.solver 分岐に存在せず launch されない (AUSM+/AUSM+-UP/KEEP)。
//   - 単一 TU を維持するため convectiveFlux_d.cu に textual include される断片であり、
//     **スタンドアロン TU ではない**。CMake の source に .cu/.cuh として足さないこと
//     (足すと __device__ グローバルが多重定義/別モジュール化し cudaMemcpyToSymbol が壊れる)。
//   - 復活させる場合は wrapper に分岐を戻し、参照ヘルパの可視性を確認する。
// =============================================================================

// interp_general: dead な AUSM/KEEP のみが関数ポインタ経由で参照する汎用補間。
// live の interp_dispatch は 1stUp/MUSCL_2nd/3rd/MINMOD のみ使うため共通ヘッダから退避。
__device__ flow_float interp_general(int scheme, int limit_scheme,
                            flow_float phiC, flow_float phiD,
                            flow_float dphidxC, flow_float dphidyC, flow_float dphidzC,
                            flow_float dphidxD, flow_float dphidyD, flow_float dphidzD,
                            flow_float dx , flow_float dy , flow_float dz,
                            flow_float cpdx, flow_float cpdy, flow_float cpdz,
                            flow_float f, flow_float limiter
                           )
{
    flow_float phif;
    flow_float k;
    flow_float r;
    flow_float psi_r;

    flow_float DD2dx;
    flow_float DD2dy;
    flow_float DD2dz;
    flow_float phiDD;
    flow_float phiU;
    flow_float limit;

    if (limit_scheme == -1) { // minmod  not completed
        r = (phiD - phiC)/(phiC - phiD);

        DD2dx = (2.0*f-1.0)*dx;
        DD2dy = (2.0*f-1.0)*dy;
        DD2dz = (2.0*f-1.0)*dz;

        phiDD = phiD - limiter*(DD2dx*dphidxD + DD2dy*dphidyD + DD2dz*dphidzD );
        phiU  = phiDD -4.0*(1.0-f)*limiter*(dx*dphidxC + dy*dphidyC + dz*dphidzC );

        r = (phiC - phiU)/(phiD - phiC);

        limiter = max(0.0, min(1.0,r));
    }

    if (scheme == 0) {
        phif = phiC;
    } else if (scheme == 1) { // 2nd order
        phif = phiC + limiter*(dphidxC*cpdx +dphidyC*cpdy +dphidzC*cpdz);
    } else if (scheme == 2) {// 3rd order
        k = 1.0/3.0;
        phif = phiC + limiter*(0.5*k*(phiD-phiC) +(1.0-k)*(dphidxC*cpdx +dphidyC*cpdy +dphidzC*cpdz));
    } else if (scheme == -1) {// ghost
        phif = phiC;
    }

    return phif;
};

inline __device__ flow_float betaPls(flow_float M, flow_float alpha) // ok
{
    if (abs(M) >= 1.0) {
        return 0.5*(1.0+sign_sano(M));
    } else {
        return 0.25*sq(M+1.0)*(2.0-M) +alpha*M*sq(M*M-1);
    }
}

inline __device__ flow_float betaMns(flow_float M, flow_float alpha) // ok
{
    if (abs(M) >= 1.0) {
        return 0.5*(1.0-sign_sano(M));
    } else {
        return 0.25*sq(M-1.0)*(2.0+M) -alpha*M*sq(M*M-1);
    }
}

inline __device__ flow_float MPls(flow_float M) // ok
{
    flow_float beta = 1.0/8.0;
    if (abs(M) >= 1.0) {
        return 0.5*(M+abs(M));
    } else {
        return +0.25*sq(M+1.0) +beta*sq(M*M-1.0);
    }
}

inline __device__ flow_float MMns(flow_float M) // ok
{
    flow_float beta = 1.0/8.0;
    if (abs(M) >= 1.0) {
        return 0.5*(M-abs(M));
    } else {
        return -0.25*sq(M-1.0) -beta*sq(M*M-1.0);
    }
}

__global__ void AUSMp_d
( 
 int conv_scheme, int limit_scheme,

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

 flow_float* limiter_ro,
 flow_float* limiter_Ux,
 flow_float* limiter_Uy,
 flow_float* limiter_Uz,
 flow_float* limiter_P,

 flow_float* ducros,

 flow_float* drodx    , flow_float* drody   , flow_float* drodz,
 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz
 

)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    //if (ip < nPlanes) {
    if (ip < nNormalPlanes) {

        flow_float (*interp)(int , int , flow_float  , flow_float , 
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float );

        if (ip >= nNormalPlanes) conv_scheme = -1; // ghost

        if (conv_scheme == 0 || conv_scheme == -1) {
            interp = &interp_1stUp;
        } else if (conv_scheme == 1 && limit_scheme >= 0) {
            interp = &interp_MUSCL_2nd;
        } else if (conv_scheme == 2 && limit_scheme >= 0) {
            interp = &interp_MUSCL_3rd;
        } else if (limit_scheme == -1) { 
            interp = &interp_MINMOD;
        }

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

        flow_float dc0p_x = pcx[ip] - ccx_0;
        flow_float dc0p_y = pcy[ip] - ccy_0;
        flow_float dc0p_z = pcz[ip] - ccz_0;

        flow_float dc1p_x = pcx[ip] - ccx_1;
        flow_float dc1p_y = pcy[ip] - ccy_1;
        flow_float dc1p_z = pcz[ip] - ccz_1;

        flow_float ro_L = interp(conv_scheme, limit_scheme, ro[ic0] , ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], drodx[ic1], drody[ic1], drodz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_ro[ic0]);
        flow_float Ux_L = interp(conv_scheme, limit_scheme, Ux[ic0] , Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Ux[ic0]);
        flow_float Uy_L = interp(conv_scheme, limit_scheme, Uy[ic0] , Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uy[ic0]);
        flow_float Uz_L = interp(conv_scheme, limit_scheme, Uz[ic0] , Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uz[ic0]);
        flow_float P_L  = interp(conv_scheme, limit_scheme, Ps[ic0] , Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_P[ic0]);
        //flow_float Ux_L = roUx_L/ro_L;
        //flow_float Uy_L = roUy_L/ro_L;
        //flow_float Uz_L = roUz_L/ro_L;
        flow_float roe_L = P_L/(ga-1.0) + 0.5*ro_L*(Ux_L*Ux_L +Uy_L*Uy_L +Uz_L*Uz_L);
        //flow_float P_L = (ga-1.0)*roe_L - 0.5*(ga-1.0)*ro_L*(Ux_L*Ux_L +Uy_L*Uy_L +Uz_L*Uz_L);
        flow_float Ht_L= roe_L/ro_L + P_L/ro_L;

        flow_float ro_R  = interp(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro[ic1]); 
        flow_float Ux_R  = interp(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Ux[ic1]);
        flow_float Uy_R  = interp(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uy[ic1]);
        flow_float Uz_R  = interp(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uz[ic1]);
        flow_float P_R   = interp(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dPdx[ic0] , dPdy[ic0] , dPdz[ic0] ,-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_P[ic1]);
        flow_float roe_R = P_R/(ga-1.0) + 0.5*ro_R*(Ux_R*Ux_R +Uy_R*Uy_R +Uz_R*Uz_R);
        //flow_float P_R = (ga-1.0)*roe_R - 0.5*(ga-1.0)*ro_R*(Ux_R*Ux_R +Uy_R*Uy_R +Uz_R*Uz_R);
        flow_float Ht_R= roe_R/ro_R + P_R/ro_R;




        flow_float U_L = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float U_R = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        //flow_float c_L = sqrt(ga*P_L/ro_L);
        //flow_float c_R = sqrt(ga*P_R/ro_R);
        //flow_float c_half = 0.5*(c_L+c_R);

        Ht_L = ga/(ga-1.0)*P_L/ro_L + 0.5*(Ux_L*Ux_L + Uy_L*Uy_L + Uz_L*Uz_L);
        Ht_R = ga/(ga-1.0)*P_R/ro_R + 0.5*(Ux_R*Ux_R + Uy_R*Uy_R + Uz_R*Uz_R);

        flow_float c_star_L = sqrt(2.0*(ga-1.0)/(ga+1.0)*Ht_L); //ok
        flow_float c_star_R = sqrt(2.0*(ga-1.0)/(ga+1.0)*Ht_R); //ok
        flow_float c_tilde_L = sq(c_star_L)/max(c_star_L, abs(U_L)); //ok
        flow_float c_tilde_R = sq(c_star_R)/max(c_star_R, abs(U_R)); //ok
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
 int conv_scheme, int limit_scheme,

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

 flow_float* limiter_ro,
 flow_float* limiter_Ux,
 flow_float* limiter_Uy,
 flow_float* limiter_Uz,
 flow_float* limiter_P,

 flow_float* ducros,

 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz

)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    //if (ip < nPlanes) {
    if (ip < nNormalPlanes) {


        flow_float (*interp)(int , int , flow_float  , flow_float , 
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float );


        if (ip >= nNormalPlanes) conv_scheme = -1; // ghost

        if (conv_scheme == 0 || conv_scheme == -1) {
            interp = &interp_1stUp;
        } else if (conv_scheme == 1 && limit_scheme >= 0) {
            interp = &interp_MUSCL_2nd;
        } else if (conv_scheme == 2 && limit_scheme >= 0) {
            interp = &interp_MUSCL_3rd;
        } else if (limit_scheme == -1) { 
            interp = &interp_MINMOD;
        }




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

        flow_float dc0p_x = pcx[ip] - ccx_0;
        flow_float dc0p_y = pcy[ip] - ccy_0;
        flow_float dc0p_z = pcz[ip] - ccz_0;

        flow_float dc1p_x = pcx[ip] - ccx_1;
        flow_float dc1p_y = pcy[ip] - ccy_1;
        flow_float dc1p_z = pcz[ip] - ccz_1;


        flow_float ro_L = interp(conv_scheme, limit_scheme, ro[ic0] , ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], drodx[ic1], drody[ic1], drodz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_ro[ic0]);
        flow_float Ux_L = interp(conv_scheme, limit_scheme, Ux[ic0] , Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Ux[ic0]);
        flow_float Uy_L = interp(conv_scheme, limit_scheme, Uy[ic0] , Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uy[ic0]);
        flow_float Uz_L = interp(conv_scheme, limit_scheme, Uz[ic0] , Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uz[ic0]);
        flow_float P_L  = interp(conv_scheme, limit_scheme, Ps[ic0] , Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_P[ic0]);
        //flow_float Ux_L = roUx_L/ro_L;
        //flow_float Uy_L = roUy_L/ro_L;
        //flow_float Uz_L = roUz_L/ro_L;
        flow_float roe_L = P_L/(ga-1.0) + 0.5*ro_L*(Ux_L*Ux_L +Uy_L*Uy_L +Uz_L*Uz_L);
        //flow_float P_L = (ga-1.0)*roe_L - 0.5*(ga-1.0)*ro_L*(Ux_L*Ux_L +Uy_L*Uy_L +Uz_L*Uz_L);
        flow_float Ht_L= roe_L/ro_L + P_L/ro_L;

        flow_float ro_R  = interp(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro[ic1]); 
        flow_float Ux_R  = interp(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Ux[ic1]);
        flow_float Uy_R  = interp(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uy[ic1]);
        flow_float Uz_R  = interp(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uz[ic1]);
        flow_float P_R   = interp(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dPdx[ic0] , dPdy[ic0] , dPdz[ic0] ,-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_P[ic1]);
        flow_float roe_R = P_R/(ga-1.0) + 0.5*ro_R*(Ux_R*Ux_R +Uy_R*Uy_R +Uz_R*Uz_R);
        //flow_float P_R = (ga-1.0)*roe_R - 0.5*(ga-1.0)*ro_R*(Ux_R*Ux_R +Uy_R*Uy_R +Uz_R*Uz_R);
        flow_float Ht_R= roe_R/ro_R + P_R/ro_R;

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

__global__ void KEEP_FVS_d
( 
 int conv_scheme, int limit_scheme,

 // gas property
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

 flow_float* limiter_ro,
 flow_float* limiter_Ux,
 flow_float* limiter_Uy,
 flow_float* limiter_Uz,
 flow_float* limiter_P,


 flow_float* ducros,

 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz
 
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;

    //if (ip < nPlanes) {
    if (ip < nNormalPlanes) {

        flow_float dUx;
        flow_float dUy;
        flow_float dUz;
        flow_float dU;

        flow_float A[5][5];
        //flow_float lam[5][5];
        flow_float lam[5];
        flow_float R[5][5];
        flow_float L[5][5];
        flow_float dW[5];
        flow_float difQ1[5];
        flow_float difQ2[5];

        flow_float Q[5];

        int zero = 0;
        int one  = 1;
        int two  = 2;
        int mone = -1;

        int i , j;

        flow_float (*interp)(int , int , flow_float  , flow_float , 
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float , flow_float ,
                             flow_float , flow_float );

        //if (limit_scheme == mone) { 
        //    interp = &interp_MINMOD;
        //}

        interp = &interp_general;
        //interp = &interp_MUSCL_2nd;
        //interp = &interp_MINMOD;

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        //__syncthreads();

        geom_float f = fx[ip];
        
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        geom_float nx = sxx/sss;
        geom_float ny = syy/sss;
        geom_float nz = szz/sss;

        flow_float ccx_0 = ccx[ic0];
        flow_float ccy_0 = ccy[ic0];
        flow_float ccz_0 = ccz[ic0];

        flow_float ccx_1 = ccx[ic1];
        flow_float ccy_1 = ccy[ic1];
        flow_float ccz_1 = ccz[ic1];

        flow_float dcc_x = ccx_1 - ccx_0;
        flow_float dcc_y = ccy_1 - ccy_0;
        flow_float dcc_z = ccz_1 - ccz_0;

        flow_float dc0p_x = pcx[ip] - ccx_0;
        flow_float dc0p_y = pcy[ip] - ccy_0;
        flow_float dc0p_z = pcz[ip] - ccz_0;

        flow_float dc1p_x = pcx[ip] - ccx_1;
        flow_float dc1p_y = pcy[ip] - ccy_1;
        flow_float dc1p_z = pcz[ip] - ccz_1;

        flow_float ro_L = interp(conv_scheme, limit_scheme, ro[ic0] , ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], drodx[ic1], drody[ic1], drodz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_ro[ic0]);
        flow_float Ux_L = interp(conv_scheme, limit_scheme, Ux[ic0] , Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Ux[ic0]);
        flow_float Uy_L = interp(conv_scheme, limit_scheme, Uy[ic0] , Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uy[ic0]);
        flow_float Uz_L = interp(conv_scheme, limit_scheme, Uz[ic0] , Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uz[ic0]);
        flow_float P_L  = interp(conv_scheme, limit_scheme, Ps[ic0] , Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_P[ic0]);
        flow_float roe_L = P_L/(ga-1.0) + 0.5*ro_L*(Ux_L*Ux_L +Uy_L*Uy_L +Uz_L*Uz_L);
        flow_float Ht_L= roe_L/ro_L + P_L/ro_L;
        flow_float ca_L  = sqrt(ga*P_L/ro_L);

        flow_float ro_R  = interp(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro[ic1]); 
        flow_float Ux_R  = interp(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Ux[ic1]);
        flow_float Uy_R  = interp(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uy[ic1]);
        flow_float Uz_R  = interp(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uz[ic1]);
        flow_float P_R   = interp(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dPdx[ic0] , dPdy[ic0] , dPdz[ic0] ,-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_P[ic1]);
        flow_float roe_R = P_R/(ga-1.0) + 0.5*ro_R*(Ux_R*Ux_R +Uy_R*Uy_R +Uz_R*Uz_R);

        flow_float U_L = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float U_R = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;
        flow_float Ht_R= roe_R/ro_R + P_R/ro_R;

        // calc roe average //ok
        flow_float roa = ro_L;
        flow_float ua  = Ux_L;
        flow_float va  = Uy_L;
        flow_float wa  = Uz_L;
        flow_float Ha  = Ht_L;
        flow_float ca  = sqrt((ga-1.0)*(Ha-0.5*(ua*ua +va*va +wa*wa)));
        flow_float Ua  = U_L;

        flow_float P_ro = 0.5*(ga-1.0)*(ua*ua + va*va + wa*wa);
        flow_float P_roe = ga-1.0;
        flow_float P_rou = -ua*(ga-1.0);
        flow_float P_rov = -va*(ga-1.0);
        flow_float P_row = -wa*(ga-1.0);
        flow_float z = ua*ua+va*va+wa*wa - P_ro/P_roe;
        flow_float fai = P_ro-ca*ca;

        // see  Numerical simulation of dense gas ows on unstructured grids, Colonna
        // right eigen matrix //ok
        R[0][0] = 1.0;      R[0][1] = nx;               R[0][2] = ny;              R[0][3] = nz;              R[0][4] = 1.0;
        R[1][0] = Ha-ca*Ua; R[1][1] = z*nx+va*nz-wa*ny; R[1][2] = z*ny+wa*nx-ua*nz;R[1][3] = z*nz+ua*ny-va*nx;R[1][4] = Ha+ca*Ua;
        R[2][0] = ua-ca*nx; R[2][1] = ua*nx;            R[2][2] = ua*ny-nz        ;R[2][3] = ua*nz+ny        ;R[2][4] = ua+ca*nx;
        R[3][0] = va-ca*ny; R[3][1] = va*nx+nz;         R[3][2] = va*ny           ;R[3][3] = va*nz-nx        ;R[3][4] = va+ca*ny;
        R[4][0] = wa-ca*nz; R[4][1] = wa*nx-ny;         R[4][2] = wa*ny+nx        ;R[4][3] = wa*nz;           R[4][4] = wa+ca*nz;

        // left eigen matrix
        L[0][0] = 0.5*(P_ro+ca*Ua);            L[0][1] = P_roe*0.5; L[0][2] = (P_rou-ca*nx)*0.5 ; L[0][3] = (P_rov-ca*ny)*0.5 ; L[0][4] = (P_row-ca*nz)*0.5;
        L[1][0] = -fai*nx-ca*ca*(va*nz-wa*ny); L[1][1] = -P_roe*nx; L[1][2] = -P_rou*nx         ; L[1][3] = -P_rov*nx+ca*ca*nz; L[1][4] = -P_row*nx-ca*ca*ny;
        L[2][0] = -fai*ny-ca*ca*(wa*nx-ua*nz); L[2][1] = -P_roe*ny; L[2][2] = -P_rou*ny-ca*ca*nz; L[2][3] = -P_rov*ny         ; L[2][4] = -P_row*ny+ca*ca*nx;
        L[3][0] = -fai*nz-ca*ca*(ua*ny-va*nx); L[3][1] = -P_roe*nz; L[3][2] = -P_rou*nz+ca*ca*ny; L[3][3] = -P_rov*nz-ca*ca*nx; L[3][4] = -P_row*nz;
        L[4][0] = 0.5*(P_ro-ca*Ua); 
        L[4][1] = 0.5*P_roe;
        L[4][2] = 0.5*(P_rou+ca*nx);
        L[4][3] = 0.5*(P_rov+ca*ny);
        L[4][4] = 0.5*(P_row+ca*nz);

        Q[0] = ro_L;
        Q[1] = ro_L*Ht_L-P_L ;
        Q[2] = ro_L*Ux_L;
        Q[3] = ro_L*Uy_L;
        Q[4] = ro_L*Uz_L;

        for (i=0; i<5; i++) {
            for (j=0; j<5; j++) {
                L[i][j] = L[i][j]/(ca*ca);
            }
        }

        lam[0] = abs(Ua - ca);
        lam[1] = abs(Ua     );
        lam[2] = abs(Ua     );
        lam[3] = abs(Ua     );
        lam[4] = abs(Ua + ca);


        dW[0] = 0.0; 
        dW[1] = 0.0; 
        dW[2] = 0.0; 
        dW[3] = 0.0; 
        dW[4] = 0.0; 

        for (i=0; i<5; i++) {
            for (j=0; j<5; j++) {
                dW[i] += L[i][j]*Q[j];
            }
        }

        difQ1[0] = 0.0;
        difQ1[1] = 0.0;
        difQ1[2] = 0.0;
        difQ1[3] = 0.0;
        difQ1[4] = 0.0;

        for (i=0; i<5; i++) {
            difQ1[i] = lam[i]*dW[i];
        }
 
        difQ2[0] = 0.0;
        difQ2[1] = 0.0;
        difQ2[2] = 0.0;
        difQ2[3] = 0.0;
        difQ2[4] = 0.0;

        for (i=0; i<5; i++) {
            for (j=0; j<5; j++) {
                difQ2[i] += R[i][j]*difQ1[j];
            }
        }

        flow_float res_ro_temp   = 0.5*difQ2[0]*sss;
        flow_float res_roe_temp  = 0.5*difQ2[1]*sss;
        flow_float res_roUx_temp = 0.5*difQ2[2]*sss;
        flow_float res_roUy_temp = 0.5*difQ2[3]*sss;
        flow_float res_roUz_temp = 0.5*difQ2[4]*sss;

        roa= ro_R;
        ua = Ux_R;
        va = Uy_R;
        wa = Uz_R;
        Ha = Ht_R;
        ca  = sqrt((ga-1.0)*(Ha-0.5*(ua*ua +va*va +wa*wa)));
        Ua = U_R;

        P_ro = 0.5*(ga-1.0)*(ua*ua + va*va + wa*wa);
        P_roe = ga-1.0;
        P_rou = -ua*(ga-1.0);
        P_rov = -va*(ga-1.0);
        P_row = -wa*(ga-1.0);
        z = ua*ua+va*va+wa*wa - P_ro/P_roe;
        fai = P_ro-ca*ca;

        // see  Numerical simulation of dense gas ows on unstructured grids, Colonna
        // right eigen matrix //ok
        R[0][0] = 1.0;      R[0][1] = nx;               R[0][2] = ny;              R[0][3] = nz;              R[0][4] = 1.0;
        R[1][0] = Ha-ca*Ua; R[1][1] = z*nx+va*nz-wa*ny; R[1][2] = z*ny+wa*nx-ua*nz;R[1][3] = z*nz+ua*ny-va*nx;R[1][4] = Ha+ca*Ua;
        R[2][0] = ua-ca*nx; R[2][1] = ua*nx;            R[2][2] = ua*ny-nz        ;R[2][3] = ua*nz+ny        ;R[2][4] = ua+ca*nx;
        R[3][0] = va-ca*ny; R[3][1] = va*nx+nz;         R[3][2] = va*ny           ;R[3][3] = va*nz-nx        ;R[3][4] = va+ca*ny;
        R[4][0] = wa-ca*nz; R[4][1] = wa*nx-ny;         R[4][2] = wa*ny+nx        ;R[4][3] = wa*nz;           R[4][4] = wa+ca*nz;

        // left eigen matrix
        L[0][0] = 0.5*(P_ro+ca*Ua);            L[0][1] = P_roe*0.5; L[0][2] = (P_rou-ca*nx)*0.5 ; L[0][3] = (P_rov-ca*ny)*0.5 ; L[0][4] =  (P_row-ca*nz)*0.5;
        L[1][0] = -fai*nx-ca*ca*(va*nz-wa*ny); L[1][1] = -P_roe*nx; L[1][2] = -P_rou*nx         ; L[1][3] = -P_rov*nx+ca*ca*nz; L[1][4] = -P_row*nx-ca*ca*ny;
        L[2][0] = -fai*ny-ca*ca*(wa*nx-ua*nz); L[2][1] = -P_roe*ny; L[2][2] = -P_rou*ny-ca*ca*nz; L[2][3] = -P_rov*ny         ; L[2][4] = -P_row*ny+ca*ca*nx;
        L[3][0] = -fai*nz-ca*ca*(ua*ny-va*nx); L[3][1] = -P_roe*nz; L[3][2] = -P_rou*nz+ca*ca*ny; L[3][3] = -P_rov*nz-ca*ca*nx; L[3][4] = -P_row*nz;
        L[4][0] = 0.5*(P_ro-ca*Ua); 
        L[4][1] = 0.5*P_roe;
        L[4][2] = 0.5*(P_rou+ca*nx);
        L[4][3] = 0.5*(P_rov+ca*ny);
        L[4][4] = 0.5*(P_row+ca*nz);

        Q[0] = ro_R;
        Q[1] = ro_R*Ht_R-P_R ;
        Q[2] = ro_R*Ux_R;
        Q[3] = ro_R*Uy_R;
        Q[4] = ro_R*Uz_R;

        for (i=0; i<5; i++) {
            for (j=0; j<5; j++) {
                L[i][j] = L[i][j]/(ca*ca);
            }
        }

        lam[0] = abs(Ua - ca);
        lam[1] = abs(Ua     );
        lam[2] = abs(Ua     );
        lam[3] = abs(Ua     );
        lam[4] = abs(Ua + ca);


        dW[0] = 0.0; 
        dW[1] = 0.0; 
        dW[2] = 0.0; 
        dW[3] = 0.0; 
        dW[4] = 0.0; 

        for (i=0; i<5; i++) {
            for (j=0; j<5; j++) {
                dW[i] += L[i][j]*Q[j];
            }
        }

        difQ1[0] = 0.0;
        difQ1[1] = 0.0;
        difQ1[2] = 0.0;
        difQ1[3] = 0.0;
        difQ1[4] = 0.0;

        //for (int i=0; i<5; i++) {
        //    for (int j=0; j<5; j++) {
        //        difQ1[i] += lam[i][j]*dW[j];
        //    }
        //}

        for (i=0; i<5; i++) {
            difQ1[i] = lam[i]*dW[i];
        }
 
        difQ2[0] = 0.0;
        difQ2[1] = 0.0;
        difQ2[2] = 0.0;
        difQ2[3] = 0.0;
        difQ2[4] = 0.0;

        for (i=0; i<5; i++) {
            for (j=0; j<5; j++) {
                difQ2[i] += R[i][j]*difQ1[j];
            }
        }

        res_ro_temp   = res_ro_temp   -0.5*difQ2[0]*sss;
        res_roe_temp  = res_roe_temp  -0.5*difQ2[1]*sss;
        res_roUx_temp = res_roUx_temp -0.5*difQ2[2]*sss;
        res_roUy_temp = res_roUy_temp -0.5*difQ2[3]*sss;
        res_roUz_temp = res_roUz_temp -0.5*difQ2[4]*sss;


        if (false) {
            // KEEP
            flow_float duc = max(min(max(ducros[ic0], ducros[ic1]),1.0),0.0);

            flow_float Ctilde  = 0.5*(ro[ic0]+ro[ic1])*0.5*( (Ux[ic0]+Ux[ic1])*nx
                                                            +(Uy[ic0]+Uy[ic1])*ny
                                                            +(Uz[ic0]+Uz[ic1])*nz );
            flow_float Mtildex = Ctilde*(Ux[ic0]+Ux[ic1])*0.5;
            flow_float Mtildey = Ctilde*(Uy[ic0]+Uy[ic1])*0.5;
            flow_float Mtildez = Ctilde*(Uz[ic0]+Uz[ic1])*0.5;

            flow_float Ktilde = Ctilde*0.5*(Ux[ic0]*Ux[ic1] +Uy[ic0]*Uy[ic1] +Uz[ic0]*Uz[ic1]);
            flow_float Itilde = Ctilde*0.5*(Ps[ic0]/ro[ic0] +Ps[ic1]/ro[ic1])/(ga-1.0);

            flow_float Gtildex = 0.5*(Ps[ic0]+Ps[ic1])*nx;
            flow_float Gtildey = 0.5*(Ps[ic0]+Ps[ic1])*ny;
            flow_float Gtildez = 0.5*(Ps[ic0]+Ps[ic1])*nz;

            flow_float Ptilde = 0.5*((Ux[ic0]*Ps[ic1] + Ux[ic1]*Ps[ic0])*nx
                                    +(Uy[ic0]*Ps[ic1] + Uy[ic1]*Ps[ic0])*ny
                                    +(Uz[ic0]*Ps[ic1] + Uz[ic1]*Ps[ic0])*nz);

            res_ro_temp   = Ctilde*sss                      + res_ro_temp;
            res_roe_temp  = (Ktilde + Itilde + Ptilde)*sss  + res_roe_temp;
            res_roUx_temp = ((Mtildex + Gtildex)*sss)       + res_roUx_temp;
            res_roUy_temp = ((Mtildey + Gtildey)*sss)       + res_roUy_temp;
            res_roUz_temp = ((Mtildez + Gtildez)*sss)       + res_roUz_temp;

        } else {

            flow_float mdot = 0.5*(ro_L*(Ux_L*nx+Uy_L*ny+Uz_L*nz) + ro_R*(Ux_R*nx+Uy_R*ny+Uz_R*nz))*sss ;

            res_ro_temp   = mdot                                        + res_ro_temp;
            res_roe_temp  = mdot*0.5*(Ht_L + Ht_R)                      + res_roe_temp;
            res_roUx_temp = mdot*0.5*(Ux_L + Ux_R) +0.5*(P_L + P_R)*sxx + res_roUx_temp;
            res_roUy_temp = mdot*0.5*(Uy_L + Uy_R) +0.5*(P_L + P_R)*syy + res_roUy_temp;
            res_roUz_temp = mdot*0.5*(Uz_L + Uz_R) +0.5*(P_L + P_R)*szz + res_roUz_temp;
        }

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
