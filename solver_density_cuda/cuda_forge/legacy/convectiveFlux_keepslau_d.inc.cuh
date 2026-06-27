// =============================================================================
// convectiveFlux_d.cu の dispatch から到達不能なスキーム実装の退避先。
//   - wrapper の cfg.solver 分岐に存在せず launch されない (AUSM+/AUSM+-UP/KEEP)。
//   - 単一 TU を維持するため convectiveFlux_d.cu に textual include される断片であり、
//     **スタンドアロン TU ではない**。CMake の source に .cu/.cuh として足さないこと
//     (足すと __device__ グローバルが多重定義/別モジュール化し cudaMemcpyToSymbol が壊れる)。
//   - 復活させる場合は wrapper に分岐を戻し、参照ヘルパの可視性を確認する。
// =============================================================================

__global__ void KEEP_SLAU_d
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
 //flow_float* limiter   ,

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

    flow_float (*interp)(int , int , flow_float  , flow_float , 
                         flow_float , flow_float , flow_float ,
                         flow_float , flow_float , flow_float ,
                         flow_float , flow_float , flow_float ,
                         flow_float , flow_float , flow_float ,
                         flow_float , flow_float );

    //if (ip < nPlanes) {
    if (ip < nNormalPlanes) {

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
        //flow_float Ux_L = roUx_L/ro_L;
        //flow_float Uy_L = roUy_L/ro_L;
        //flow_float Uz_L = roUz_L/ro_L;
        flow_float roe_L = P_L/(ga-1.0) + 0.5*ro_L*(Ux_L*Ux_L +Uy_L*Uy_L +Uz_L*Uz_L);
        //flow_float P_L = (ga-1.0)*roe_L - 0.5*(ga-1.0)*ro_L*(Ux_L*Ux_L +Uy_L*Uy_L +Uz_L*Uz_L);
        flow_float Ht_L= roe_L/ro_L + P_L/ro_L;
        flow_float ca_L  = sqrt(ga*P_L/ro_L);

        flow_float ro_R  = interp(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro[ic1]); 
        flow_float Ux_R  = interp(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Ux[ic1]);
        flow_float Uy_R  = interp(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uy[ic1]);
        flow_float Uz_R  = interp(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uz[ic1]);
        flow_float P_R   = interp(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dPdx[ic0] , dPdy[ic0] , dPdz[ic0] ,-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_P[ic1]);
        flow_float roe_R = P_R/(ga-1.0) + 0.5*ro_R*(Ux_R*Ux_R +Uy_R*Uy_R +Uz_R*Uz_R);
        //flow_float P_R = (ga-1.0)*roe_R - 0.5*(ga-1.0)*ro_R*(Ux_R*Ux_R +Uy_R*Uy_R +Uz_R*Uz_R);
        flow_float Ht_R= roe_R/ro_R + P_R/ro_R;

        flow_float ro_p= ro_L;
        flow_float ro_m= ro_R;
        flow_float u_p = Ux_L;
        flow_float v_p = Uy_L;
        flow_float w_p = Uz_L;
        flow_float h_p = Ht_L;
        flow_float u_m = Ux_R;
        flow_float v_m = Uy_R;
        flow_float w_m = Uz_R;
        flow_float h_m = Ht_R;
        flow_float P_p = P_L;
        flow_float P_m = P_R;

        //flow_float Vn_p = ((Ux[ic0])*sxx +(Uy[ic0])*syy +(Uz[ic0])*szz)/sss;
        //flow_float Vn_m = ((Ux[ic1])*sxx +(Uy[ic1])*syy +(Uz[ic1])*szz)/sss;

        //flow_float roVn_p = (roUx[ic0]*sxx + roUy[ic0]*syy + roUz[ic0]*szz)/sss;
        //flow_float roVn_m = (roUx[ic1]*sxx + roUy[ic1]*syy + roUz[ic1]*szz)/sss;

        //flow_float VnL = ((Ux[ic0])*sxx +(Uy[ic0])*syy +(Uz[ic0])*szz)/sss;
        //flow_float VnR = ((Ux[ic1])*sxx +(Uy[ic1])*syy +(Uz[ic1])*szz)/sss;
        flow_float Vn_p = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float Vn_m = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        flow_float roVn_p = (ro_L*Ux_L*sxx + ro_L*Uy_L*syy + ro_L*Uz_L*szz)/sss;
        flow_float roVn_m = (ro_R*Ux_R*sxx + ro_R*Uy_R*syy + ro_R*Uz_R*szz)/sss;

        flow_float VnL = Vn_p;
        flow_float VnR = Vn_m;

        flow_float VnL_abs = abs(VnL);
        flow_float VnR_abs = abs(VnR);

//TODO: change c
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
            beta_p = 0.25*sq((M_p+1.0))*(2.0-M_p);
        }

        if (abs(M_m)>=1.0){
            beta_m = 0.5*(M_m - abs(M_m))/M_m;
        } else {
            beta_m = 0.25*sq((M_m-1.0))*(2.0+M_m);
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
        //flow_float lim_ro = min(limiter_ro[ic0], limiter_ro[ic1]);

        flow_float lim = min(limiter_Ux[ic0],limiter_Ux[ic1]);
        lim =  min(lim,  min(limiter_Uy[ic0],limiter_Uy[ic1]));
        lim =  min(lim,  min(limiter_Uz[ic0],limiter_Uz[ic1]));

        duc = max(duc, 1.0-lim);

        res_ro_temp   = (1.0-duc)*Ctilde*sss + duc*res_ro_temp;
        res_roUx_temp = (1.0-duc)*((Mtildex + Gtildex)*sss) + duc*res_roUx_temp;
        res_roUy_temp = (1.0-duc)*((Mtildey + Gtildey)*sss) + duc*res_roUy_temp;
        res_roUz_temp = (1.0-duc)*((Mtildez + Gtildez)*sss) + duc*res_roUz_temp;
        res_roe_temp  = (1.0-duc)*(Ktilde + Itilde + Ptilde)*sss + duc*res_roe_temp;

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


//__global__ void convectiveFlux_boundary_d // roe
//( 
// flow_float ga,
//
//  // mesh structure
// geom_int nb,
// geom_int* bplane_plane,  
// geom_int* bplane_cell,  
// geom_int* bplane_cell_ghst,  
//
// geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
// geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
// geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
//
// // variables
// flow_float* ro   ,
// flow_float* roUx  ,
// flow_float* roUy  ,
// flow_float* roUz  ,
// flow_float* roe ,
// flow_float* Ux  ,
// flow_float* Uy  ,
// flow_float* Uz  ,
// flow_float* Ps  ,
// flow_float* Ht  ,
// flow_float* sonic,
// flow_float* Ts  ,
//
// // bvar
// flow_float* rob ,
// flow_float* roUxb ,
// flow_float* roUyb ,
// flow_float* roUzb ,
// flow_float* roeb ,
// flow_float* Uxb ,
// flow_float* Uyb ,
// flow_float* Uzb ,
// flow_float* Ttb ,
// flow_float* Ptb ,
// flow_float* Tsb ,
// flow_float* Psb ,
// 
// flow_float* res_ro   ,
// flow_float* res_roUx  ,
// flow_float* res_roUy  ,
// flow_float* res_roUz  ,
// flow_float* res_roe   
//)
//{
//    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;
//
//    if (ib < nb) {
//        geom_int  ip = bplane_plane[ib];
//        geom_int  ic = bplane_cell[ib];
//        geom_int  ig = bplane_cell_ghst[ib];
//
//        flow_float dro;
//        flow_float dP;
//        flow_float dUx;
//        flow_float dUy;
//        flow_float dUz;
//        flow_float dU;
//
//        flow_float A[5][5];
//        flow_float lam[5][5];
//        flow_float R[5][5];
//        flow_float L[5][5];
//        flow_float dW[5];
//        flow_float difQ1[5];
//        flow_float difQ2[5];
//
//        flow_float delQ[5];
//       
//        geom_float sxx = sx[ip];
//        geom_float syy = sy[ip];
//        geom_float szz = sz[ip];
//        geom_float sss = ss[ip];
//
//        geom_float nx = sxx/sss;
//        geom_float ny = syy/sss;
//        geom_float nz = szz/sss;
//
//        flow_float ccx_0 = ccx[ic];
//        flow_float ccy_0 = ccy[ic];
//        flow_float ccz_0 = ccz[ic];
//
//        flow_float pcx_1 = pcx[ip];
//        flow_float pcy_1 = pcy[ip];
//        flow_float pcz_1 = pcz[ip];
//
//        flow_float dcc_x = pcx_1 - ccx_0;
//        flow_float dcc_y = pcy_1 - ccy_0;
//        flow_float dcc_z = pcz_1 - ccz_0;
//        flow_float dcc   = sqrt(dcc_x*dcc_x +dcc_y*dcc_y +dcc_z*dcc_z) ;
//
//        flow_float ro_L = ro[ic];
//        flow_float Ux_L = Ux[ic];
//        flow_float Uy_L = Uy[ic];
//        flow_float Uz_L = Uz[ic];
//        flow_float P_L  = Ps[ic];
//        flow_float Ht_L = Ht[ic];
//
//        flow_float ro_R= ro[ig]; 
//        flow_float Ux_R= Ux[ig];
//        flow_float Uy_R= Uy[ig];
//        flow_float Uz_R= Uz[ig];
//        flow_float P_R = Ps[ig];
//        flow_float Ht_R= Ht[ig];
//
//        // calc delta value
//        dro = ro_R - ro_L;
//        dP  = P_R  - P_L;
//        dUx = Ux_R - Ux_L;
//        dUy = Uy_R - Uy_L;
//        dUz = Uz_R - Uz_L;
//        dU  = dUx*sxx +dUy*syy +dUz*szz;
//
//        // calc roe average //ok
//        flow_float roa = sqrt(ro_R*ro_L);
//        flow_float ua = (sqrt(ro_L)*Ux_L + sqrt(ro_R)*Ux_R)/(sqrt(ro_R)+sqrt(ro_L));
//        flow_float va = (sqrt(ro_L)*Uy_L + sqrt(ro_R)*Uy_R)/(sqrt(ro_R)+sqrt(ro_L));
//        flow_float wa = (sqrt(ro_L)*Uz_L + sqrt(ro_R)*Uz_R)/(sqrt(ro_R)+sqrt(ro_L));
//        flow_float Ha = (sqrt(ro_L)*Ht_L + sqrt(ro_R)*Ht_R)/(sqrt(ro_R)+sqrt(ro_L));
//        flow_float ca  = sqrt((ga-1.0)*(Ha-0.5*(ua*ua +va*va +wa*wa)));
//        flow_float Ua  = ua*nx +va*ny +wa*nz;
//
//        flow_float P_ro = 0.0;
//        flow_float P_roe = ga-1.0;
//        flow_float P_rou = ua*(ga-1.0);
//        flow_float P_rov = va*(ga-1.0);
//        flow_float P_row = wa*(ga-1.0);
//        flow_float z = ua*ua+va*va+wa*wa - P_ro/P_roe;
//        flow_float fai = P_ro-ca*ca;
//
//        // see  Numerical simulation of dense gas ows on unstructured grids, Colonna
//        // right eigen matrix //ok
//        R[0][0] = 1.0;      R[0][1] = nx;               R[0][2] = ny;              R[0][3] = nz;              R[0][4] = 1.0;
//        R[1][0] = Ha-ca*Ua; R[1][1] = z*nx+va*nz-wa*ny; R[1][2] = z*ny+wa*nx-ua*nz;R[1][3] = z*nz+ua*ny-va*nx;R[1][4] = Ha+ca*Ua;
//        R[2][0] = ua-ca*nx; R[2][1] = ua*nx;            R[2][2] = ua*ny-nz        ;R[2][3] = ua*nz+ny        ;R[2][4] = ua+ca*nx;
//        R[3][0] = va-ca*ny; R[3][1] = va*nx+nz;         R[3][2] = va*ny           ;R[3][3] = va*nz-nx        ;R[3][4] = va+ca*ny;
//        R[4][0] = wa-ca*nz; R[4][1] = wa*nx-ny;         R[4][2] = wa*ny+nx        ;R[4][3] = wa*nz;           R[4][4] = wa+ca*nz;
//
//        // left eigen matrix
//        L[0][0] = 0.5*(P_ro+ca*Ua);            L[0][1] = P_roe*0.5; L[0][2] = (P_rou-ca*nx)*0.5 ; L[0][3] = (P_rov-ca*ny)*0.5 ; L[0][4] =  (P_row-ca*nz)*0.5;
//        L[1][0] = -fai*nx-ca*ca*(va*nz-wa*ny); L[1][1] = -P_roe*nx; L[1][2] = -P_rou*nx         ; L[1][3] = -P_rov*nx+ca*ca*nz; L[1][4] = -P_row*nx-ca*ca*ny;
//        L[2][0] = -fai*ny-ca*ca*(wa*nx-ua*nz); L[2][1] = -P_roe*ny; L[2][2] = -P_rou*ny-ca*ca*nz; L[2][3] = -P_rov*ny         ; L[2][4] = -P_row*ny+ca*ca*nx;
//        L[3][0] = -fai*nz-ca*ca*(ua*ny-va*nx); L[3][1] = -P_roe*nz; L[3][2] = -P_rou*nz+ca*ca*ny; L[3][3] = -P_rov*nz-ca*ca*nx; L[3][4] = -P_row*nz;
//        L[4][0] = 0.5*(P_ro-ca*Ua); 
//        L[4][1] = 0.5*P_roe;
//        L[4][2] = 0.5*(P_rou+ca*nx);
//        L[4][3] = 0.5*(P_rov+ca*ny);
//        L[4][4] = 0.5*(P_row+ca*nz);
//
//        delQ[0] = ro_R      - ro_L;
//        delQ[1] = (ro_R*Ht_R-P_R) - (ro_L*Ht_L-P_L) ;
//        delQ[2] = ro_R*Ux_R - ro_L*Ux_L;
//        delQ[3] = ro_R*Uy_R - ro_L*Uy_L;
//        delQ[4] = ro_R*Uz_R - ro_L*Uz_L;
//
//        //printf("delro=%e\n", ro_R);
//
//        for (int i=0; i<5; i++) {
//            for (int j=0; j<5; j++) {
//                L[i][j] = L[i][j]/(ca*ca);
//            }
//        }
//
//        // diagonal matrix;
//        for (int i=0; i<5; i++) {
//            for (int j=0; j<5; j++) {
//                lam[i][j] =0.0;
//            }
//        }
//
//        lam[0][0] = abs(Ua - ca);
//        lam[1][1] = abs(Ua     );
//        lam[2][2] = abs(Ua     );
//        lam[3][3] = abs(Ua     );
//        lam[4][4] = abs(Ua + ca);
//
//        dW[0] = 0.0; 
//        dW[1] = 0.0; 
//        dW[2] = 0.0; 
//        dW[3] = 0.0; 
//        dW[4] = 0.0; 
//
//        for (int i=0; i<5; i++) {
//            for (int j=0; j<5; j++) {
//                dW[i] += L[i][j]*delQ[j];
//            }
//        }
//
//        //dW[0] = 0.5*(dP/ca -roa*dU)/ca;
//        //dW[1] = (dro-dP/(ca*ca))*nx + roa*(nz*dUy-ny*dUz);
//        //dW[2] = (dro-dP/(ca*ca))*ny + roa*(nx*dUz-nz*dUx);
//        //dW[3] = (dro-dP/(ca*ca))*ny + roa*(ny*dUx-nx*dUy);
//        //dW[4] = 0.5*(dP/ca +roa*dU)/ca;
//
// 
//        difQ1[0] = 0.0;
//        difQ1[1] = 0.0;
//        difQ1[2] = 0.0;
//        difQ1[3] = 0.0;
//        difQ1[4] = 0.0;
//
//        for (int i=0; i<5; i++) {
//            for (int j=0; j<5; j++) {
//                difQ1[i] += lam[i][j]*dW[j];
//            }
//        }
//
//        difQ2[0] = 0.0;
//        difQ2[1] = 0.0;
//        difQ2[2] = 0.0;
//        difQ2[3] = 0.0;
//        difQ2[4] = 0.0;
//
//        for (int i=0; i<5; i++) {
//            for (int j=0; j<5; j++) {
//                difQ2[i] += R[i][j]*difQ1[j];
//            }
//        }
//
//        flow_float mdot = 0.5*(ro_L*(Ux_L*nx+Uy_L*ny+Uz_L*nz) + ro_R*(Ux_R*nx+Uy_R*ny+Uz_R*nz))*sss ;
//
//        flow_float res_ro_temp   = mdot                                        -0.5*difQ2[0]*sss;
//        flow_float res_roe_temp  = mdot*0.5*(Ht_L + Ht_R)                      -0.5*difQ2[1]*sss;
//        flow_float res_roUx_temp = mdot*0.5*(Ux_L + Ux_R) +0.5*(P_L + P_R)*sxx -0.5*difQ2[2]*sss;
//        flow_float res_roUy_temp = mdot*0.5*(Uy_L + Uy_R) +0.5*(P_L + P_R)*syy -0.5*difQ2[3]*sss;
//        flow_float res_roUz_temp = mdot*0.5*(Uz_L + Uz_R) +0.5*(P_L + P_R)*szz -0.5*difQ2[4]*sss;
//
//        atomicAdd(&res_ro[ic]  , -res_ro_temp);
//        atomicAdd(&res_roUx[ic], -res_roUx_temp);
//        atomicAdd(&res_roUy[ic], -res_roUy_temp);
//        atomicAdd(&res_roUz[ic], -res_roUz_temp);
//        atomicAdd(&res_roe[ic] , -res_roe_temp);
//    }
//
//    __syncthreads();
//}
