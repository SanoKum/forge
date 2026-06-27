// =============================================================================
// 対流フラックス HLLE_d の実装断片。
//   - convectiveFlux_d.cu に textual include される断片であり **スタンドアロン TU ではない**。
//     単一 TU を維持するため (rdc 不要 / __device__ グローバル・cudaMemcpyToSymbol の共有を保つ)、
//     CMake の source には .cu/.cuh として足さないこと。
//   - 参照する共通ヘルパ (interp_dispatch 等) と診断グローバルは include 元の上方で定義済み。
// =============================================================================

__global__ void HLLE_d
(
 int conv_scheme, int limit_scheme,

 // gas property
 flow_float ga,

 // 非平衡凝縮 (二相): エネルギー流束を二相全エンタルピーに補正。g_total==nullptr で従来 (ビット不変)。
 flow_float cp_cpg, flow_float* g_total, flow_float* T_cell, int condModel,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,
 geom_int nNormal_halo_Planes, geom_int* normal_halo_planes_d,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 flow_float* massflux,

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
    geom_int ip_orig = blockDim.x*blockIdx.x + threadIdx.x;

    if (ip_orig < nNormal_halo_Planes) {
        geom_int ip = normal_halo_planes_d[ip_orig];

        geom_int ic0 = plane_cells[2*ip+0];
        geom_int ic1 = plane_cells[2*ip+1];

        // Force 1st-order upwind at non-periodic boundary planes (ghost on R side).
        // Gradients/limiters are not maintained for ghost cells, so reconstructing
        // from them would corrupt the face state. The ghost value itself encodes
        // the boundary condition and must be used as-is.
        if (ic1 >= nCells) {
            conv_scheme = -1;
        }

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

        flow_float ro_L = interp_dispatch(conv_scheme, limit_scheme, ro[ic0] , ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], drodx[ic1], drody[ic1], drodz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_ro[ic0]);
        flow_float Ux_L = interp_dispatch(conv_scheme, limit_scheme, Ux[ic0] , Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Ux[ic0]);
        flow_float Uy_L = interp_dispatch(conv_scheme, limit_scheme, Uy[ic0] , Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uy[ic0]);
        flow_float Uz_L = interp_dispatch(conv_scheme, limit_scheme, Uz[ic0] , Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uz[ic0]);
        flow_float P_L  = interp_dispatch(conv_scheme, limit_scheme, Ps[ic0] , Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_P[ic0]);

        flow_float ro_R = interp_dispatch(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0], -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro[ic1]);
        flow_float Ux_R = interp_dispatch(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Ux[ic1]);
        flow_float Uy_R = interp_dispatch(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uy[ic1]);
        flow_float Uz_R = interp_dispatch(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uz[ic1]);
        flow_float P_R  = interp_dispatch(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_P[ic1]);

        flow_float small_rho = 1.0e-4f;
        flow_float small_p   = 1.0f;
        flow_float small_a2  = 1.0f;

        ro_L = max(ro_L, small_rho);
        ro_R = max(ro_R, small_rho);
        P_L = max(P_L, small_p);
        P_R = max(P_R, small_p);

        flow_float velocity2_L = Ux_L*Ux_L + Uy_L*Uy_L + Uz_L*Uz_L;
        flow_float velocity2_R = Ux_R*Ux_R + Uy_R*Uy_R + Uz_R*Uz_R;

        flow_float roe_L = P_L/(ga-1.0) + 0.5*ro_L*velocity2_L;
        flow_float roe_R = P_R/(ga-1.0) + 0.5*ro_R*velocity2_R;
        // 非平衡凝縮 (二相): 単相 roe/Ht を二相 (内部エネルギー e=(cv+gR)T-gL) に補正。差 g(cpT-L)/質量。
        // HLLE は roe (S_L S_R(roe_R-roe_L)) と Ht (ro Ht U) を両方使うため両方補正。g_total==nullptr で従来。
        if (g_total != nullptr) {
            const CondSpeciesProps cprC = (condModel == 1) ? condProps_H2O() : condProps_N2();
            const flow_float dL0 = g_total[ic0]*(cp_cpg*T_cell[ic0] - (flow_float)cond_latent(cprC, (double)T_cell[ic0]));
            const flow_float dR0 = g_total[ic1]*(cp_cpg*T_cell[ic1] - (flow_float)cond_latent(cprC, (double)T_cell[ic1]));
            roe_L += ro_L*dL0;
            roe_R += ro_R*dR0;
        }
        flow_float Ht_L = (roe_L + P_L)/ro_L;
        flow_float Ht_R = (roe_R + P_R)/ro_R;

        flow_float U_L = Ux_L*nx + Uy_L*ny + Uz_L*nz;
        flow_float U_R = Ux_R*nx + Uy_R*ny + Uz_R*nz;

        flow_float ca_L = sqrt(max(ga*P_L/ro_L, small_a2));
        flow_float ca_R = sqrt(max(ga*P_R/ro_R, small_a2));

        flow_float sqrt_ro_L = sqrt(ro_L);
        flow_float sqrt_ro_R = sqrt(ro_R);
        flow_float sqrt_ro_sum = max(sqrt_ro_L + sqrt_ro_R, small_rho);
        flow_float ua = (sqrt_ro_L*Ux_L + sqrt_ro_R*Ux_R)/sqrt_ro_sum;
        flow_float va = (sqrt_ro_L*Uy_L + sqrt_ro_R*Uy_R)/sqrt_ro_sum;
        flow_float wa = (sqrt_ro_L*Uz_L + sqrt_ro_R*Uz_R)/sqrt_ro_sum;
        flow_float Ha = (sqrt_ro_L*Ht_L + sqrt_ro_R*Ht_R)/sqrt_ro_sum;
        flow_float q2a = ua*ua + va*va + wa*wa;
        flow_float ca = sqrt(max((ga-1.0)*(Ha-0.5*q2a), small_a2));
        flow_float Ua = ua*nx + va*ny + wa*nz;

        flow_float S_L = min(U_L - ca_L, min(U_R - ca_R, Ua - ca));
        flow_float S_R = max(U_L + ca_L, max(U_R + ca_R, Ua + ca));

        flow_float F_L_0 = ro_L*U_L;
        flow_float F_L_1 = ro_L*Ux_L*U_L + P_L*nx;
        flow_float F_L_2 = ro_L*Uy_L*U_L + P_L*ny;
        flow_float F_L_3 = ro_L*Uz_L*U_L + P_L*nz;
        flow_float F_L_4 = ro_L*Ht_L*U_L;

        flow_float F_R_0 = ro_R*U_R;
        flow_float F_R_1 = ro_R*Ux_R*U_R + P_R*nx;
        flow_float F_R_2 = ro_R*Uy_R*U_R + P_R*ny;
        flow_float F_R_3 = ro_R*Uz_R*U_R + P_R*nz;
        flow_float F_R_4 = ro_R*Ht_R*U_R;

        flow_float res_ro_temp;
        flow_float res_roUx_temp;
        flow_float res_roUy_temp;
        flow_float res_roUz_temp;
        flow_float res_roe_temp;

        if (0.0 <= S_L) {
            res_ro_temp = F_L_0*sss;
            res_roUx_temp = F_L_1*sss;
            res_roUy_temp = F_L_2*sss;
            res_roUz_temp = F_L_3*sss;
            res_roe_temp = F_L_4*sss;
        } else if (S_R <= 0.0) {
            res_ro_temp = F_R_0*sss;
            res_roUx_temp = F_R_1*sss;
            res_roUy_temp = F_R_2*sss;
            res_roUz_temp = F_R_3*sss;
            res_roe_temp = F_R_4*sss;
        } else {
            flow_float denom = max(S_R - S_L, 1.0e-8);
            res_ro_temp = ((S_R*F_L_0 - S_L*F_R_0 + S_L*S_R*(ro_R - ro_L))/denom)*sss;
            res_roUx_temp = ((S_R*F_L_1 - S_L*F_R_1 + S_L*S_R*(ro_R*Ux_R - ro_L*Ux_L))/denom)*sss;
            res_roUy_temp = ((S_R*F_L_2 - S_L*F_R_2 + S_L*S_R*(ro_R*Uy_R - ro_L*Uy_L))/denom)*sss;
            res_roUz_temp = ((S_R*F_L_3 - S_L*F_R_3 + S_L*S_R*(ro_R*Uz_R - ro_L*Uz_L))/denom)*sss;
            res_roe_temp = ((S_R*F_L_4 - S_L*F_R_4 + S_L*S_R*(roe_R - roe_L))/denom)*sss;
        }

        atomicAdd(&res_ro[ic0], -res_ro_temp);
        atomicAdd(&res_roUx[ic0], -res_roUx_temp);
        atomicAdd(&res_roUy[ic0], -res_roUy_temp);
        atomicAdd(&res_roUz[ic0], -res_roUz_temp);
        atomicAdd(&res_roe[ic0], -res_roe_temp);

        atomicAdd(&res_ro[ic1], res_ro_temp);
        atomicAdd(&res_roUx[ic1], res_roUx_temp);
        atomicAdd(&res_roUy[ic1], res_roUy_temp);
        atomicAdd(&res_roUz[ic1], res_roUz_temp);
        atomicAdd(&res_roe[ic1], res_roe_temp);
    }

    __syncthreads();
}
