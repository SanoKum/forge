#include <cstdlib>
#include "convectiveFlux_d.cuh"
#include "lowMachPrecond_d.cuh"
#include "speciesTransport_d.cuh"  // species_roY_device_ptr()
#include "condensationProperties_d.cuh"  // n2_latent (二相エネルギー流束の潜熱補正)

// free-stream 保存: 対流流束の圧力項を (p_tilde - d_pRef)*s で組むための基準静圧。
// wrapper で cfg.pRef を cudaMemcpyToSymbol。既定 0.0 で従来挙動 (ビット不変)。
// 非直交メッシュで大きな p*s を float32 加算する際の桁落ち(metric closure 由来の
// 偽運動量源)を抑える。詳細: .github/plans/freestream-preserving-flux.md
__constant__ flow_float d_pRef;

// D2a/D4 診断 (env ゲート, 既定 off = ビット不変): 多成分 TP の contact 混合層 limit-cycle 切り分け。
//   FORGE_CONTACT_1ST=1 : 面組成センサ s_Y=max_s|Y_sR-Y_sL|/(Y_sR+Y_sL+ε) が g_contactThresh を超える
//                         内部面で flow(ρ,p,u) の MUSCL 再構成を無効化し1次 (セル値) にする (D2a)。species は不変。
//   FORGE_CONTACT_LOG=1 : s_Y>g_contactLogThresh の面で L/R 状態を printf (mixed-order/limiter chatter 観察, D4)。
__device__ int        g_contact1st = 0;
__device__ flow_float g_contactThresh = 0.05f;
__device__ int        g_contactLog = 0;
__device__ flow_float g_contactLogThresh = 0.3f;
// 連続ブレンド版 (chatter-free): ω_Y=min(1, s_Y/g_contactBlend) で flow 再構成を 2次→1次 へ滑らかに寄せる。
// 0=off。hard 切替 (FORGE_CONTACT_1ST) と違いセンサ閾値の bang-bang chatter を起こさない (D2a 清浄版/製品候補)。
__device__ flow_float g_contactBlend = 0.0f;
// S2 診断 (FORGE_FACE_THERMOY=1, 既定0): TP face thermo の組成 mixed-order を解消。
// R_mix/T/h を **face 補間組成** Y_face=f·Y_ic0+(1-f)·Y_ic1 (正規化) で評価し、ρ_L,P_L (2次再構成) と
// 同じ face 位置の組成に整合させる (現行は owner セル組成=1次 → ΔT_f^MO~30K)。species 流束は不変。
__device__ int g_faceThermoY = 0;
// speciesFaceReconstruction==1: Y_s を ρ と同じ勾配+limiter で face へ再構成し thermo/species 流束で
// 同一 face 組成を使う (proper S2/S3)。wrapper で cfg.speciesFaceReconstruction を設定。
__device__ int g_speciesFaceRecon = 0;

// K7: pow(x,2.0) は exp(2*log(x)) に展開され重い。2乗は乗算 1 命令で済むため sq() に置換。
__device__ __forceinline__ flow_float sq(flow_float x) { return x * x; }

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

__device__ flow_float interp_MUSCL_2nd(int scheme, int limit_scheme,
                                       flow_float phiC, flow_float phiD, 
                                       flow_float dphidx, flow_float dphidy, flow_float dphidz,
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

    phif = phiC + limiter*(dphidx*cpdx +dphidy*cpdy +dphidz*cpdz);

    return phif;
};

__device__ flow_float interp_MUSCL_3rd(int scheme, int limit_scheme,
                                       flow_float phiC, flow_float phiD, 
                                       flow_float dphidx, flow_float dphidy, flow_float dphidz,
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

    k = 1.0/3.0;
    phif = phiC + limiter*(0.5*k*(phiD-phiC) +(1.0-k)*(dphidx*cpdx +dphidy*cpdy +dphidz*cpdz));

    return phif;
};


__device__ flow_float interp_1stUp(int scheme, int limit_scheme,
                                 flow_float phiC, flow_float phiD, 
                                 flow_float dphidx, flow_float dphidy, flow_float dphidz,
                                 flow_float dphidxD, flow_float dphidyD, flow_float dphidzD,
                                 flow_float dx , flow_float dy , flow_float dz,
                                 flow_float cpdx, flow_float cpdy, flow_float cpdz,
                                 flow_float f, flow_float limiter
                                )
{
    return phiC;
};

__device__ flow_float interp_MINMOD(int scheme, int limit_scheme,
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

    DD2dx = (2.0*f-1.0)*dx;
    DD2dy = (2.0*f-1.0)*dy;
    DD2dz = (2.0*f-1.0)*dz;

    phiDD = phiD - (DD2dx*dphidxD + DD2dy*dphidyD + DD2dz*dphidzD );
    phiU  = phiDD -4.0*(1.0-f)*(dx*dphidxC + dy*dphidyC + dz*dphidzC );

    r = (phiC - phiU)/(phiD - phiC);

    //limit = sign_sano(r)*max(0.0, (min(abs(r), sign_sano(r))));
    limit = max(0.0, min(1.0,r));

    phif = phiC + limit*(dphidxC*cpdx +dphidyC*cpdy +dphidzC*cpdz);

    return phif;
};

__device__ __forceinline__ flow_float interp_dispatch(int scheme, int limit_scheme,
                                                      flow_float phiC, flow_float phiD,
                                                      flow_float dphidxC, flow_float dphidyC, flow_float dphidzC,
                                                      flow_float dphidxD, flow_float dphidyD, flow_float dphidzD,
                                                      flow_float dx, flow_float dy, flow_float dz,
                                                      flow_float cpdx, flow_float cpdy, flow_float cpdz,
                                                      flow_float f, flow_float limiter)
{
    if (scheme == 0 || scheme == -1) {
        return interp_1stUp(scheme, limit_scheme,
                            phiC, phiD,
                            dphidxC, dphidyC, dphidzC,
                            dphidxD, dphidyD, dphidzD,
                            dx, dy, dz,
                            cpdx, cpdy, cpdz,
                            f, limiter);
    }

    if (scheme == 1 && limit_scheme >= 0) {
        return interp_MUSCL_2nd(scheme, limit_scheme,
                                phiC, phiD,
                                dphidxC, dphidyC, dphidzC,
                                dphidxD, dphidyD, dphidzD,
                                dx, dy, dz,
                                cpdx, cpdy, cpdz,
                                f, limiter);
    }

    if (scheme == 2 && limit_scheme >= 0) {
        return interp_MUSCL_3rd(scheme, limit_scheme,
                                phiC, phiD,
                                dphidxC, dphidyC, dphidzC,
                                dphidxD, dphidyD, dphidzD,
                                dx, dy, dz,
                                cpdx, cpdy, cpdz,
                                f, limiter);
    }

    return interp_MINMOD(scheme, limit_scheme,
                         phiC, phiD,
                         dphidxC, dphidyC, dphidzC,
                         dphidxD, dphidyD, dphidzD,
                         dx, dy, dz,
                         cpdx, cpdy, cpdz,
                         f, limiter);
}

__device__ __forceinline__ flow_float apply_ducros_limiter(flow_float limiter, flow_float duc)
{
    if (duc <= 0.8) {
        return limiter;
    }

    return max(0.0, (1.0 - duc) * limiter);
}



__device__ flow_float betaPls_slau(flow_float M)
{
    if (abs(M) >= 1.0) {
        return 0.25*(2.0-M)*sq(M+1.0);
    } else {
        return 0.5*(1.0+sign_sano(+M));
    }
}

__device__ flow_float betaMns_slau(flow_float M)
{
    if (abs(M) >= 1.0) {
        return 0.25*(2.0+M)*sq(M-1.0);
    } else {
        return 0.5*(1.0+sign_sano(-M));
    }
}

__global__ void SLAU_d
(
 int conv_scheme, int limit_scheme,
 int slauVariant,   // 1: SLAU, 2: SLAU2 (圧力束第3項のみ低マッハ改良)
 int lowMachPrecond, flow_float precondEps,   // 低マッハ前処理 (1: 散逸スケールを c'、0: 従来 c_hat)
 int lowMachThornber,                         // Thornber 再構成補正 (1: L/R 速度ジャンプを z=min(M,1) で縮約)

 flow_float ga,
 int thermalMethod,                           // 0: calorically perfect, 2: thermally-perfect (NASA-9)
 const SpeciesThermo* sp, int nSpecies,       // thermally-perfect 用化学種データ
 flow_float** roY,                            // 多成分: セル毎 ρY_s (nullptr で単成分 sp[0])
 flow_float** Yd_recon,                       // face 整合再構成用 Y{s} とセル勾配 ∇Y{s} (g_speciesFaceRecon 時)
 flow_float** dYdx_recon, flow_float** dYdy_recon, flow_float** dYdz_recon,
 flow_float* Rmix_cell,                       // M6: per-cell 混合比気体定数 R[ic] (面エンタルピー用キャッシュ)

 // 非平衡凝縮 (二相): エネルギー流束を二相全エンタルピーに補正する。g_total==nullptr で従来 (ビット不変)。
 // cp_cpg=CPG 比熱 (pure)。condModel=凝縮種 (0:N2, 1:H2O carrier の潜熱選択)。
 flow_float cp_cpg, flow_float* g_total, flow_float* T_cell, int condModel,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,
 geom_int nNormal_ghst_Planes , geom_int* normal_ghst_planes_d,
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


    //if (ip < nPlanes) {
    if (ip_orig < nNormal_ghst_Planes) {

        geom_int ip = normal_ghst_planes_d[ip_orig];

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        // At non-periodic boundary planes, ic1 is a ghost cell (index >= nCells)
        // whose state has been set by applyBconds(). Gradients/limiters are not
        // computed for ghost cells, so we must NOT reconstruct from ghost; force
        // 1st-order upwind here so that the face state equals the ghost state
        // exactly (which encodes the boundary condition).
        if (ic1 >= nCells) {
            conv_scheme = -1;
        }

        // D2a/D4 診断: 多成分 TP 内部面の組成センサ s_Y = max_s|Y_sR-Y_sL|/(Y_sR+Y_sL+ε)。
        flow_float sY_dbg = 0.0;
        if ((g_contact1st || g_contactLog || g_contactBlend > 0.0) && nSpecies > 1 && roY != nullptr && ic1 < nCells) {
            const flow_float ro0 = max(ro[ic0], (flow_float)1.0e-30);
            const flow_float ro1 = max(ro[ic1], (flow_float)1.0e-30);
            for (int s = 0; s < nSpecies; ++s) {
                const flow_float yl = roY[s][ic0] / ro0;
                const flow_float yr = roY[s][ic1] / ro1;
                sY_dbg = max(sY_dbg, fabsf(yr - yl) / (yr + yl + (flow_float)1.0e-12));
            }
            // D2a: 組成勾配が強い面で flow(ρ,p,u) を1次化 (既存 ghost 用 conv_scheme=-1 を流用)。species 流束は不変。
            if (g_contact1st && sY_dbg > g_contactThresh) conv_scheme = -1;
        }

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

        //flow_float lim_ro = min(limiter_ro[ic0], limiter_ro[ic1]);
        //flow_float lim_Ux = min(limiter_Ux[ic0], limiter_Ux[ic1]);
        //flow_float lim_Uy = min(limiter_Uy[ic0], limiter_Uy[ic1]);
        //flow_float lim_Uz = min(limiter_Uz[ic0], limiter_Uz[ic1]);
        //flow_float lim_P  = min(limiter_P [ic0], limiter_P[ic1]);

        flow_float ro_L = interp_dispatch(conv_scheme, limit_scheme, ro[ic0] , ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], drodx[ic1], drody[ic1], drodz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_ro[ic0]);
        flow_float Ux_L = interp_dispatch(conv_scheme, limit_scheme, Ux[ic0] , Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Ux[ic0]);
        flow_float Uy_L = interp_dispatch(conv_scheme, limit_scheme, Uy[ic0] , Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uy[ic0]);
        flow_float Uz_L = interp_dispatch(conv_scheme, limit_scheme, Uz[ic0] , Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uz[ic0]);
        flow_float P_L  = interp_dispatch(conv_scheme, limit_scheme, Ps[ic0] , Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_P[ic0]);
        //flow_float Ux_L = roUx_L/ro_L;
        //flow_float Uy_L = roUy_L/ro_L;
        //flow_float Uz_L = roUz_L/ro_L;
        // velocity2_L / h_p はブレンド後に算出するため後段へ移動 (lowMachThornber 対応)。

        flow_float ro_R  = interp_dispatch(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro[ic1]); 
        flow_float Ux_R  = interp_dispatch(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Ux[ic1]);
        flow_float Uy_R  = interp_dispatch(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uy[ic1]);
        flow_float Uz_R  = interp_dispatch(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uz[ic1]);
        flow_float P_R   = interp_dispatch(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dPdx[ic0] , dPdy[ic0] , dPdz[ic0] ,-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_P[ic1]);

        // D2a 連続ブレンド (chatter-free): 強組成勾配ほど flow 再構成をセル値(1次)へ滑らかに寄せる。
        if (g_contactBlend > 0.0 && sY_dbg > 0.0) {
            const flow_float w = min((flow_float)1.0, sY_dbg / g_contactBlend);
            if (w > (flow_float)0.0) {
                const flow_float w1 = (flow_float)1.0 - w;
                ro_L = w1*ro_L + w*ro[ic0]; Ux_L = w1*Ux_L + w*Ux[ic0]; Uy_L = w1*Uy_L + w*Uy[ic0]; Uz_L = w1*Uz_L + w*Uz[ic0]; P_L = w1*P_L + w*Ps[ic0];
                ro_R = w1*ro_R + w*ro[ic1]; Ux_R = w1*Ux_R + w*Ux[ic1]; Uy_R = w1*Uy_R + w*Uy[ic1]; Uz_R = w1*Uz_R + w*Uz[ic1]; P_R = w1*P_R + w*Ps[ic1];
            }
        }
        // 低マッハ Thornber 再構成補正: L/R 速度ジャンプを z=min(M,1) で縮約し、低マッハで
        // O(1/M) に増大する速度ジャンプ由来の散逸を抑える。lowMachThornber==0 で恒等 (ビット不変)、
        // M>=1 で z=1 (超音速域不変)。圧力 P_L/R・密度 ro_L/R は不変。理論は docs/convection/theory.md。
        if (lowMachThornber == 1) {
            flow_float c_hat_th = 0.5*(sonic[ic0] + sonic[ic1]);
            flow_float v2L_th = Ux_L*Ux_L + Uy_L*Uy_L + Uz_L*Uz_L;
            flow_float v2R_th = Ux_R*Ux_R + Uy_R*Uy_R + Uz_R*Uz_R;
            flow_float z_th = min(static_cast<flow_float>(1.0), sqrt(0.5*(v2L_th + v2R_th))/c_hat_th);
            flow_float um, du;
            um = 0.5*(Ux_L+Ux_R); du = 0.5*(Ux_L-Ux_R); Ux_L = um + z_th*du; Ux_R = um - z_th*du;
            um = 0.5*(Uy_L+Uy_R); du = 0.5*(Uy_L-Uy_R); Uy_L = um + z_th*du; Uy_R = um - z_th*du;
            um = 0.5*(Uz_L+Uz_R); du = 0.5*(Uz_L-Uz_R); Uz_L = um + z_th*du; Uz_R = um - z_th*du;
        }
        // velocity2_L / h_p はブレンド後に算出 (L 再構成直後から移動)。
        flow_float velocity2_L = Ux_L*Ux_L + Uy_L*Uy_L + Uz_L*Uz_L;
        flow_float velocity2_R = Ux_R*Ux_R + Uy_R*Uy_R + Uz_R*Uz_R;
        flow_float h_p, h_m;
        if (thermalMethod == 2) {
            // TP gas: 被移流全エンタルピーを NASA 絶対エンタルピーで再構成する。
            //   ga/(ga-1)·P/ρ = cp·T であって h(T)=∫cp dT' ではないため誤り。
            //   面温度 T_face = P_face/(ρ_face·R_mix) → h_mix(T_face) を NASA で評価し ek を加える。
            //   多成分 (roY!=nullptr): L 側は ic0, R 側は ic1 のセル組成 Y を使う (組成は 1 次)。
            //   He/N2 の質量基準 cp は ~5 倍違うため、sp[0] 固定だと混合 contact で
            //   エネルギー束が不整合になり圧力発散する。
            if (nSpecies > 1 && roY != nullptr) {
                double YL[THERMO_MAX_SPECIES], YR[THERMO_MAX_SPECIES];
                const double roL_c = (double)max(ro[ic0], (flow_float)1.0e-30);
                const double roR_c = (double)max(ro[ic1], (flow_float)1.0e-30);
                double sL=0.0, sR=0.0;
                for (int s=0;s<nSpecies;s++){
                    double yl=(double)roY[s][ic0]/roL_c; if(yl<0.0)yl=0.0; YL[s]=yl; sL+=yl;
                    double yr=(double)roY[s][ic1]/roR_c; if(yr<0.0)yr=0.0; YR[s]=yr; sR+=yr;
                }
                const double iL=1.0/(sL>1.0e-30?sL:1.0e-30), iR=1.0/(sR>1.0e-30?sR:1.0e-30);
                for (int s=0;s<nSpecies;s++){ YL[s]*=iL; YR[s]*=iR; }
                double RgL, RgR;
                if (g_speciesFaceRecon && Yd_recon != nullptr) {
                    // proper S2/S3: Y_s を ρ と同じ勾配+limiter (limiter_ro) で face へ再構成。
                    // positivity は clamp(≥0)+正規化 (ΣY=1) で閉じる (簡易版; 共通 θ_Y は今後)。
                    double sLr=0.0, sRr=0.0;
                    for (int s=0;s<nSpecies;s++){
                        flow_float ylr = interp_dispatch(conv_scheme, limit_scheme, Yd_recon[s][ic0], Yd_recon[s][ic1],
                            dYdx_recon[s][ic0], dYdy_recon[s][ic0], dYdz_recon[s][ic0],
                            dYdx_recon[s][ic1], dYdy_recon[s][ic1], dYdz_recon[s][ic1],
                            dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_ro[ic0]);
                        flow_float yrr = interp_dispatch(conv_scheme, limit_scheme, Yd_recon[s][ic1], Yd_recon[s][ic0],
                            dYdx_recon[s][ic1], dYdy_recon[s][ic1], dYdz_recon[s][ic1],
                            dYdx_recon[s][ic0], dYdy_recon[s][ic0], dYdz_recon[s][ic0],
                            -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro[ic1]);
                        if (ylr < 0.0) ylr = 0.0; if (yrr < 0.0) yrr = 0.0;
                        YL[s]=(double)ylr; YR[s]=(double)yrr; sLr+=(double)ylr; sRr+=(double)yrr;
                    }
                    const double iLr=1.0/(sLr>1.0e-30?sLr:1.0e-30), iRr=1.0/(sRr>1.0e-30?sRr:1.0e-30);
                    for (int s=0;s<nSpecies;s++){ YL[s]*=iLr; YR[s]*=iRr; }
                    RgL = thermo_R_mix(sp, nSpecies, YL);
                    RgR = thermo_R_mix(sp, nSpecies, YR);
                } else if (g_faceThermoY) {
                    // S2: 整合 face 組成 Y_face=f·Y_ic0+(1-f)·Y_ic1 (正規化) を L/R thermo に使い、
                    // ρ_L,P_L (2次) と同じ face 位置の組成に揃える。R_mix も Y_face から作る。
                    double Yf[THERMO_MAX_SPECIES], sf=0.0;
                    for (int s=0;s<nSpecies;s++){ Yf[s] = (double)f*YL[s] + (1.0-(double)f)*YR[s]; sf+=Yf[s]; }
                    const double inv = 1.0/(sf>1.0e-30?sf:1.0e-30);
                    for (int s=0;s<nSpecies;s++){ Yf[s]*=inv; YL[s]=Yf[s]; YR[s]=Yf[s]; }
                    RgL = RgR = thermo_R_mix(sp, nSpecies, Yf);
                } else {
                    // R_mix はセル組成のみの定数 (現行 mixed-order)。dependentVariables で計算済の per-cell 値。
                    RgL = (double)Rmix_cell[ic0];
                    RgR = (double)Rmix_cell[ic1];
                }
                const double Tl = (double)P_L/((double)ro_L*RgL);
                const double Tr = (double)P_R/((double)ro_R*RgR);
                h_p = (flow_float)(thermo_h_mix(sp, nSpecies, YL, Tl) + 0.5*(double)velocity2_L);
                h_m = (flow_float)(thermo_h_mix(sp, nSpecies, YR, Tr) + 0.5*(double)velocity2_R);
            } else {
                const double Rg = thermo_R_species(sp[0]);
                const double Tl = (double)P_L/((double)ro_L*Rg);
                const double Tr = (double)P_R/((double)ro_R*Rg);
                h_p = (flow_float)(thermo_h_mass(sp[0], Tl) + 0.5*(double)velocity2_L);
                h_m = (flow_float)(thermo_h_mass(sp[0], Tr) + 0.5*(double)velocity2_R);
            }
            // carrier+condensible (H2O in N2) 二相補正: 凝縮した水の潜熱を引く
            //   h_2phase = h_gas^全蒸気 - g L_w + ek。気相混合 h は全蒸気なので -g L_w を足すだけ。
            if (g_total != nullptr) {
                const CondSpeciesProps cprL = (condModel == 1) ? condProps_H2O() : condProps_N2();
                h_p -= (flow_float)((double)g_total[ic0]*cond_latent(cprL, (double)T_cell[ic0]));
                h_m -= (flow_float)((double)g_total[ic1]*cond_latent(cprL, (double)T_cell[ic1]));
            }
        } else {
            h_p = ga*P_L/((ga-1.0)*ro_L) + 0.5*velocity2_L;
            h_m = ga*P_R/((ga-1.0)*ro_R) + 0.5*velocity2_R;
            // 非平衡凝縮 (二相): 単相 h=cp(1-g)T+ek を二相全エンタルピー h=cpT-gL+ek に補正。
            //   差 = g(cpT - L)。これを落とすとエネルギー流束が潜熱分を運ばず全エンタルピー非保存になる。
            //   c_hat (音速) は気相近似で単相のまま。g・T はセル値(1次, g は元来1次風上移流で整合)。
            //   g_total==nullptr (凝縮 off) で従来 (ビット不変)。
            if (g_total != nullptr) {
                const CondSpeciesProps cprC = (condModel == 1) ? condProps_H2O() : condProps_N2();
                const flow_float gL0 = g_total[ic0], gR0 = g_total[ic1];
                const flow_float TL0 = T_cell[ic0],  TR0 = T_cell[ic1];
                h_p += gL0*(cp_cpg*TL0 - (flow_float)cond_latent(cprC, (double)TL0));
                h_m += gR0*(cp_cpg*TR0 - (flow_float)cond_latent(cprC, (double)TR0));
            }
        }

        //flow_float Vn_p = ((Ux[ic0])*sxx +(Uy[ic0])*syy +(Uz[ic0])*szz)/sss;
        //flow_float Vn_m = ((Ux[ic1])*sxx +(Uy[ic1])*syy +(Uz[ic1])*szz)/sss;

        //flow_float roVn_p = (roUx[ic0]*sxx + roUy[ic0]*syy + roUz[ic0]*szz)/sss;
        //flow_float roVn_m = (roUx[ic1]*sxx + roUy[ic1]*syy + roUz[ic1]*szz)/sss;

        //flow_float VnL = ((Ux[ic0])*sxx +(Uy[ic0])*syy +(Uz[ic0])*szz)/sss;
        //flow_float VnR = ((Ux[ic1])*sxx +(Uy[ic1])*syy +(Uz[ic1])*szz)/sss;
        flow_float Vn_p = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float Vn_m = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        flow_float Vn_p_abs = abs(Vn_p);
        flow_float Vn_m_abs = abs(Vn_m);

//TODO: change c
        flow_float c_hat = 0.5*(sonic[ic0] + sonic[ic1]);

        // D4 診断: 強組成勾配面の L/R 状態を出力 (mixed-order face-state / limiter chatter 観察)。
        // 一面=一行/step。後処理で ip でフィルタすれば pseudo-step 時系列になる。
        if (g_contactLog && sY_dbg > g_contactLogThresh) {
            const flow_float RgL = (nSpecies > 1 && Rmix_cell != nullptr) ? Rmix_cell[ic0] : (flow_float)0.0;
            const flow_float RgR = (nSpecies > 1 && Rmix_cell != nullptr) ? Rmix_cell[ic1] : (flow_float)0.0;
            const flow_float TL  = (RgL > 0.0) ? P_L / (ro_L * RgL) : (flow_float)0.0;
            const flow_float TR  = (RgR > 0.0) ? P_R / (ro_R * RgR) : (flow_float)0.0;
            // mixed-order 誤差: 整合 face 組成 (面平均 Y) で R_mix を作った場合の T との差 ΔT_f^MO。
            flow_float Rmix_face = 0.0;
            if (nSpecies > 1 && roY != nullptr) {
                const flow_float ro0c = max(ro[ic0],(flow_float)1e-30), ro1c = max(ro[ic1],(flow_float)1e-30);
                for (int s = 0; s < nSpecies; ++s) {
                    const flow_float yf = (flow_float)0.5*(roY[s][ic0]/ro0c + roY[s][ic1]/ro1c);
                    Rmix_face += yf * (flow_float)thermo_R_species(sp[s]);
                }
            }
            const flow_float TL_cons = (Rmix_face > 0.0) ? P_L/(ro_L*Rmix_face) : (flow_float)0.0;
            const flow_float dT_MO   = TL - TL_cons;
            // limiter (ψ_ρ,ψ_p,ψ_u) と cell 値も出し、再構成が実際に効いていたか (sharp で ψ≈0 か) を判定可能に。
            printf("CLOG ip=%d ic0=%d ic1=%d sY=%.4f roL=%.5f roR=%.5f PL=%.1f PR=%.1f UxL=%.4f UxR=%.4f limP=%.4f limRo=%.4f limUx=%.4f roC=%.5f PC=%.1f UxC=%.4f RgL=%.2f RgR=%.2f TL=%.2f TR=%.2f hL=%.6e cL=%.2f\n",
                   (int)ip, (int)ic0, (int)ic1, (double)sY_dbg, (double)ro_L, (double)ro_R,
                   (double)P_L, (double)P_R, (double)Ux_L, (double)Ux_R, (double)limiter_P[ic0],
                   (double)limiter_ro[ic0], (double)limiter_Ux[ic0],
                   (double)ro[ic0], (double)Ps[ic0], (double)Ux[ic0],
                   (double)RgL, (double)RgR, (double)TL, (double)TR, (double)h_p, (double)sonic[ic0]);
            printf("CMO ip=%d sY=%.4f TL=%.3f TLcons=%.3f dT_MO=%.4f RgL=%.2f Rmix_face=%.2f\n",
                   (int)ip, (double)sY_dbg, (double)TL, (double)TL_cons, (double)dT_MO, (double)RgL, (double)Rmix_face);
        }

        flow_float M_p = Vn_p/c_hat;
        flow_float M_m = Vn_m/c_hat;

        flow_float P_del = P_R - P_L;

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
        flow_float Vn_hat_abs   = (ro_L*Vn_p_abs + ro_R*Vn_m_abs)/(ro_L + ro_R);
        flow_float Vn_hat_p_abs = (1.0-g)*Vn_hat_abs + g*Vn_p_abs;
        flow_float Vn_hat_m_abs = (1.0-g)*Vn_hat_abs + g*Vn_m_abs;

        flow_float M_hat = min(one, sqrt(half*(velocity2_R + velocity2_L))/c_hat);
        flow_float chi = (1.0-M_hat)*(1.0-M_hat);

        flow_float pressure_sum = P_L + P_R;
        // 圧力束の第3項のみ slauVariant で分岐 (mdot は SLAU/SLAU2 共通)。
        // SLAU : (1-chi)(beta_p+beta_m-1) * (P_L+P_R)/2   ... 低マッハで消失し圧力散逸が乏しい
        // SLAU2: (beta_p+beta_m-1) * sqrt((|u_L|^2+|u_R|^2)/2) * roBar * c_hat  ... M に比例した低マッハ散逸
        flow_float p_third;
        if (slauVariant == 2) {
            flow_float roBar = half*(ro_L + ro_R);
            p_third = (beta_p+beta_m-one)*sqrt(half*(velocity2_L+velocity2_R))*roBar*c_hat;
        } else {
            p_third = (one-chi)*(beta_p+beta_m-one)*half*pressure_sum;
        }
        flow_float p_tilde = half*pressure_sum + half*(beta_p-beta_m)*(P_L-P_R) + p_third;

        // 低マッハ前処理: 圧力散逸項のスケール c_hat を前処理音速 c_diss に置き換える。
        // lowMachPrecond==0 では c_diss==c_hat でビット不変。M>=1 でも c'=c_hat に復帰。
        flow_float c_diss = c_hat;
        if (lowMachPrecond >= 1) {   // 1: フラックス散逸前処理のみ / 2: + LHS 完全前処理 (どちらも c' 散逸を使う)
            flow_float velMag_face = sqrt(half*(velocity2_L + velocity2_R));
            flow_float Un_face     = half*(Vn_p + Vn_m);
            c_diss = lowMachCprime(c_hat, velMag_face, Un_face, precondEps);
        }

        flow_float mdot = sss*0.5*((ro_L*(Vn_p+Vn_hat_p_abs)+ro_R*(Vn_m-Vn_hat_m_abs)) -chi/(c_diss)*P_del);
        massflux[ip] = mdot;

        flow_float res_ro_temp   = mdot;
        flow_float p_tilde_r = p_tilde - d_pRef;   // free-stream 保存: 基準静圧を差し引いて float32 桁落ちを抑制
        flow_float res_roUx_temp = 0.5*(mdot+abs(mdot))*Ux_L +0.5*(mdot-abs(mdot))*Ux_R +p_tilde_r*sxx;
        flow_float res_roUy_temp = 0.5*(mdot+abs(mdot))*Uy_L +0.5*(mdot-abs(mdot))*Uy_R +p_tilde_r*syy;
        flow_float res_roUz_temp = 0.5*(mdot+abs(mdot))*Uz_L +0.5*(mdot-abs(mdot))*Uz_R +p_tilde_r*szz;
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

inline __device__ flow_float sign_sano(flow_float x)
{
    return x > 0 ? 1 : (x<0 ? -1 : 0);
}

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



__global__ void ROE_d
(
 int conv_scheme, int limit_scheme,

 // gas property
 flow_float ga,
 int thermalMethod,                           // 0: calorically perfect, 2: thermally-perfect (NASA-9)
 const SpeciesThermo* sp, int nSpecies,       // thermally-perfect 用化学種データ

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

    //if (ip < nPlanes) {
    if (ip_orig < nNormal_halo_Planes) {

        geom_int ip = normal_halo_planes_d[ip_orig];

        flow_float dro;
        flow_float dP;
        flow_float dUx;
        flow_float dUy;
        flow_float dUz;
        flow_float dU;

        flow_float lam[5];
        flow_float R[5][5];
        flow_float L[5][5];
        flow_float dW[5];
        flow_float difQ1[5];
        flow_float difQ2[5];
        flow_float delQ[5];

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        // Force 1st-order upwind at non-periodic boundary planes (ghost on R side).
        // Gradients/limiters are not maintained for ghost cells, so reconstructing
        // from them would corrupt the face state. The ghost value itself encodes
        // the boundary condition and must be used as-is.
        if (ic1 >= nCells) {
            conv_scheme = -1;
        }

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

        flow_float small_rho = 1.0e-4f;
        flow_float small_p   = 1.0f;
        flow_float small_a2  = 1.0f;

        flow_float duc = max(min(max(ducros[ic0], ducros[ic1]), 1.0), 0.0);
        flow_float limiter_ro_L = apply_ducros_limiter(limiter_ro[ic0], duc);
        flow_float limiter_Ux_L = apply_ducros_limiter(limiter_Ux[ic0], duc);
        flow_float limiter_Uy_L = apply_ducros_limiter(limiter_Uy[ic0], duc);
        flow_float limiter_Uz_L = apply_ducros_limiter(limiter_Uz[ic0], duc);
        flow_float limiter_P_L = apply_ducros_limiter(limiter_P[ic0], duc);
        flow_float limiter_ro_R = apply_ducros_limiter(limiter_ro[ic1], duc);
        flow_float limiter_Ux_R = apply_ducros_limiter(limiter_Ux[ic1], duc);
        flow_float limiter_Uy_R = apply_ducros_limiter(limiter_Uy[ic1], duc);
        flow_float limiter_Uz_R = apply_ducros_limiter(limiter_Uz[ic1], duc);
        flow_float limiter_P_R = apply_ducros_limiter(limiter_P[ic1], duc);

        flow_float ro_L = interp_dispatch(conv_scheme, limit_scheme, ro[ic0] , ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], drodx[ic1], drody[ic1], drodz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_ro_L);
        flow_float Ux_L = interp_dispatch(conv_scheme, limit_scheme, Ux[ic0] , Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Ux_L);
        flow_float Uy_L = interp_dispatch(conv_scheme, limit_scheme, Uy[ic0] , Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uy_L);
        flow_float Uz_L = interp_dispatch(conv_scheme, limit_scheme, Uz[ic0] , Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uz_L);
        flow_float P_L  = interp_dispatch(conv_scheme, limit_scheme, Ps[ic0] , Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_P_L);

        ro_L = max(ro_L, small_rho);
        P_L = max(P_L, small_p);

        //flow_float Ux_L = roUx_L/ro_L;
        //flow_float Uy_L = roUy_L/ro_L;
        //flow_float Uz_L = roUz_L/ro_L;
        flow_float v2_L = Ux_L*Ux_L +Uy_L*Uy_L +Uz_L*Uz_L;
        flow_float roe_L, Ht_L, ca_L;
        if (thermalMethod == 2) {
            // TP gas: 内部エネルギー・全エンタルピー・音速を NASA で再構成 (sp[0], 単成分)
            const double Rg = thermo_R_species(sp[0]);
            const double Tl = (double)P_L/((double)ro_L*Rg);
            const double hl = thermo_h_mass(sp[0], Tl);
            const double cpl= thermo_cp_mass(sp[0], Tl);
            const double gl = cpl/((cpl-Rg) > 1.0e-6 ? (cpl-Rg) : 1.0e-6);
            roe_L = (flow_float)((double)ro_L*((hl - Rg*Tl) + 0.5*(double)v2_L));
            Ht_L  = (flow_float)(hl + 0.5*(double)v2_L);
            ca_L  = (flow_float)sqrt(max(gl*Rg*Tl, (double)small_a2));
        } else {
            roe_L = P_L/(ga-1.0) + 0.5*ro_L*v2_L;
            Ht_L  = roe_L/ro_L + P_L/ro_L;
            ca_L  = sqrt(max(ga*P_L/ro_L, small_a2));
            // 非平衡凝縮 (二相): roe/Ht を二相に補正 (差 g(cpT-L)/質量)。g_total==nullptr で従来。
            if (g_total != nullptr) {
                const CondSpeciesProps cprC = (condModel == 1) ? condProps_H2O() : condProps_N2();
                const flow_float dL0 = g_total[ic0]*(cp_cpg*T_cell[ic0] - (flow_float)cond_latent(cprC, (double)T_cell[ic0]));
                roe_L += ro_L*dL0;  Ht_L += dL0;
            }
        }

        flow_float ro_R  = interp_dispatch(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro_R); 
        flow_float Ux_R  = interp_dispatch(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Ux_R);
        flow_float Uy_R  = interp_dispatch(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uy_R);
        flow_float Uz_R  = interp_dispatch(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uz_R);
        flow_float P_R   = interp_dispatch(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dPdx[ic0] , dPdy[ic0] , dPdz[ic0] ,-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_P_R);

        ro_R = max(ro_R, small_rho);
        P_R = max(P_R, small_p);

        flow_float v2_R = Ux_R*Ux_R +Uy_R*Uy_R +Uz_R*Uz_R;
        flow_float roe_R, Ht_R, ca_R;
        if (thermalMethod == 2) {
            const double Rg = thermo_R_species(sp[0]);
            const double Tr = (double)P_R/((double)ro_R*Rg);
            const double hr = thermo_h_mass(sp[0], Tr);
            const double cpr= thermo_cp_mass(sp[0], Tr);
            const double gr = cpr/((cpr-Rg) > 1.0e-6 ? (cpr-Rg) : 1.0e-6);
            roe_R = (flow_float)((double)ro_R*((hr - Rg*Tr) + 0.5*(double)v2_R));
            Ht_R  = (flow_float)(hr + 0.5*(double)v2_R);
            ca_R  = (flow_float)sqrt(max(gr*Rg*Tr, (double)small_a2));
        } else {
            roe_R = P_R/(ga-1.0) + 0.5*ro_R*v2_R;
            Ht_R  = roe_R/ro_R + P_R/ro_R;
            ca_R  = sqrt(max(ga*P_R/ro_R, small_a2));
            if (g_total != nullptr) {
                const CondSpeciesProps cprC = (condModel == 1) ? condProps_H2O() : condProps_N2();
                const flow_float dR0 = g_total[ic1]*(cp_cpg*T_cell[ic1] - (flow_float)cond_latent(cprC, (double)T_cell[ic1]));
                roe_R += ro_R*dR0;  Ht_R += dR0;
            }
        }

        flow_float U_L = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float U_R = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        // calc delta value
        dro = ro_R - ro_L;
        dP  = P_R  - P_L;
        dUx = Ux_R - Ux_L;
        dUy = Uy_R - Uy_L;
        dUz = Uz_R - Uz_L;
        dU  = (dUx*sxx +dUy*syy +dUz*szz)/sss;

        // calc roe average //ok
        flow_float sqrt_ro_L = sqrt(ro_L);
        flow_float sqrt_ro_R = sqrt(ro_R);
        flow_float sqrt_ro_sum = max(sqrt_ro_R + sqrt_ro_L, small_rho);
        flow_float roa = sqrt_ro_R*sqrt_ro_L;
        flow_float ua = (sqrt_ro_L*Ux_L + sqrt_ro_R*Ux_R)/sqrt_ro_sum;
        flow_float va = (sqrt_ro_L*Uy_L + sqrt_ro_R*Uy_R)/sqrt_ro_sum;
        flow_float wa = (sqrt_ro_L*Uz_L + sqrt_ro_R*Uz_R)/sqrt_ro_sum;
        flow_float Ha = (sqrt_ro_L*Ht_L + sqrt_ro_R*Ht_R)/sqrt_ro_sum;
        flow_float Ua  = ua*nx +va*ny +wa*nz;

        // Roe 平均状態の有効比熱比 ga_eff と音速 ca。
        //   CPG: ga_eff=ga, ca=sqrt((ga-1)(Ha-ek))。
        //   TP : Roe 平均静エンタルピー h~=Ha-ek から T~ を反転し ga_eff=cp~/(cp~-R),
        //        ca=sqrt(ga_eff·R·T~) (= sqrt(ga_eff·P/ρ))。圧力微分ブロックは κ=ga_eff-1 を
        //        採用 (χ=∂P/∂ρ|_{ρe} 項は stage-A で無視; defect-correction で収束解は不変)。
        flow_float ga_eff = ga;
        flow_float ca;
        const flow_float ek_a = 0.5*(ua*ua + va*va + wa*wa);
        if (thermalMethod == 2) {
            const double Yone = 1.0;   // M1 単成分の質量分率
            const double Rg = thermo_R_species(sp[0]);
            const double ha = (double)Ha - (double)ek_a;
            const double Ta = thermo_T_from_h(sp, 1, &Yone, ha, 300.0, 50.0, 6000.0);
            const double cpa = thermo_cp_mass(sp[0], Ta);
            ga_eff = (flow_float)(cpa/((cpa-Rg) > 1.0e-6 ? (cpa-Rg) : 1.0e-6));
            ca = (flow_float)sqrt(max((double)ga_eff*Rg*Ta, (double)small_a2));
        } else {
            ca = sqrt(max((ga-1.0)*(Ha-ek_a), small_a2));
        }

        flow_float P_ro = 0.5*(ga_eff-1.0)*(ua*ua + va*va + wa*wa);
        flow_float P_roe = ga_eff-1.0;
        flow_float P_rou = -ua*(ga_eff-1.0);
        flow_float P_rov = -va*(ga_eff-1.0);
        flow_float P_row = -wa*(ga_eff-1.0);
        flow_float z = ua*ua+va*va+wa*wa - P_ro/P_roe;
        flow_float fai = P_ro-ca*ca;

        R[0][0] = 1.0;      R[0][1] = nx;               R[0][2] = ny;              R[0][3] = nz;              R[0][4] = 1.0;
        R[1][0] = Ha-ca*Ua; R[1][1] = z*nx+va*nz-wa*ny; R[1][2] = z*ny+wa*nx-ua*nz;R[1][3] = z*nz+ua*ny-va*nx;R[1][4] = Ha+ca*Ua;
        R[2][0] = ua-ca*nx; R[2][1] = ua*nx;            R[2][2] = ua*ny-nz        ;R[2][3] = ua*nz+ny        ;R[2][4] = ua+ca*nx;
        R[3][0] = va-ca*ny; R[3][1] = va*nx+nz;         R[3][2] = va*ny           ;R[3][3] = va*nz-nx        ;R[3][4] = va+ca*ny;
        R[4][0] = wa-ca*nz; R[4][1] = wa*nx-ny;         R[4][2] = wa*ny+nx        ;R[4][3] = wa*nz;           R[4][4] = wa+ca*nz;

        L[0][0] = 0.5*(P_ro+ca*Ua);            L[0][1] = P_roe*0.5; L[0][2] = (P_rou-ca*nx)*0.5 ; L[0][3] = (P_rov-ca*ny)*0.5 ; L[0][4] =  (P_row-ca*nz)*0.5;
        L[1][0] = -fai*nx-ca*ca*(va*nz-wa*ny); L[1][1] = -P_roe*nx; L[1][2] = -P_rou*nx         ; L[1][3] = -P_rov*nx+ca*ca*nz; L[1][4] = -P_row*nx-ca*ca*ny;
        L[2][0] = -fai*ny-ca*ca*(wa*nx-ua*nz); L[2][1] = -P_roe*ny; L[2][2] = -P_rou*ny-ca*ca*nz; L[2][3] = -P_rov*ny         ; L[2][4] = -P_row*ny+ca*ca*nx;
        L[3][0] = -fai*nz-ca*ca*(ua*ny-va*nx); L[3][1] = -P_roe*nz; L[3][2] = -P_rou*nz+ca*ca*ny; L[3][3] = -P_rov*nz-ca*ca*nx; L[3][4] = -P_row*nz;
        L[4][0] = 0.5*(P_ro-ca*Ua);
        L[4][1] = 0.5*P_roe;
        L[4][2] = 0.5*(P_rou+ca*nx);
        L[4][3] = 0.5*(P_rov+ca*ny);
        L[4][4] = 0.5*(P_row+ca*nz);

        delQ[0] = ro_R - ro_L;
        delQ[1] = (ro_R*Ht_R-P_R) - (ro_L*Ht_L-P_L);
        delQ[2] = ro_R*Ux_R - ro_L*Ux_L;
        delQ[3] = ro_R*Uy_R - ro_L*Uy_L;
        delQ[4] = ro_R*Uz_R - ro_L*Uz_L;

        for (int i=0; i<5; i++) {
            for (int j=0; j<5; j++) {
                L[i][j] = L[i][j]/(ca*ca);
            }
        }

        lam[0] = abs(Ua - ca);
        lam[1] = abs(Ua     );
        lam[2] = abs(Ua     );
        lam[3] = abs(Ua     );
        lam[4] = abs(Ua + ca);

        // van leer's entropy fix
        flow_float eta_vl ;
        if (false) {
            flow_float eta_vl = max(dU/sss-(ca_R-ca_L),0.0);
            if (lam[0] < 2*eta_vl) lam[0] = lam[0]*lam[0]/(4*eta_vl) + eta_vl;

            eta_vl = max(dU/sss,0.0);
            if (lam[1] < 2*eta_vl) {
                lam[1] = lam[1]*lam[1]/(4*eta_vl) + eta_vl;
                lam[2] = lam[1];
                lam[3] = lam[1];
            }
            eta_vl = max(dU/sss+(ca_R-ca_L),0.0);
            if (lam[4] < 2*eta_vl) lam[4] = lam[4]*lam[4]/(4*eta_vl) + eta_vl;
        }
        // harten entropy correction
        if (false) {
            eta_vl = 0.1; // 0.05~0.15
            if (lam[0] <= eta_vl) lam[0] = (lam[0]*lam[0] + eta_vl*eta_vl)/(2*eta_vl);

            if (lam[1] < eta_vl) {
                lam[1] = (lam[1]*lam[1] + eta_vl*eta_vl)/(2*eta_vl);
                lam[2] = lam[1];
                lam[3] = lam[1];
            }
            if (lam[4] <= eta_vl) lam[4] = (lam[4]*lam[4] + eta_vl*eta_vl)/(2*eta_vl);
        }
 

        // harten entropy correction
        if (true) {
            eta_vl = 0.1*(abs(Ua)/ca+1.0); 

            if (lam[0] <= eta_vl) lam[0] = (lam[0]*lam[0] + eta_vl*eta_vl)/(2*eta_vl);

            if (lam[1] < eta_vl) {
                lam[1] = (lam[1]*lam[1] + eta_vl*eta_vl)/(2*eta_vl);
                lam[2] = lam[1];
                lam[3] = lam[1];
            }
            if (lam[4] <= eta_vl) lam[4] = (lam[4]*lam[4] + eta_vl*eta_vl)/(2*eta_vl);
        }
 
        // entropy fix
        //flow_float h_fix = 0.05*(abs(Ua)+ca);

        //flow_float l_temp = lam[0];
        //lam[0] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));

        //l_temp = lam[1];
        //lam[1] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));

        //l_temp = lam[2];
        //lam[2] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));

        //l_temp = lam[3];
        //lam[3] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));

        //l_temp = lam[4];
        //lam[4] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));
        dW[0] = 0.0;
        dW[1] = 0.0;
        dW[2] = 0.0;
        dW[3] = 0.0;
        dW[4] = 0.0;

        for (int i=0; i<5; i++) {
            for (int j=0; j<5; j++) {
                dW[i] += L[i][j]*delQ[j];
            }
        }

        difQ1[0] = 0.0;
        difQ1[1] = 0.0;
        difQ1[2] = 0.0;
        difQ1[3] = 0.0;
        difQ1[4] = 0.0;

        for (int i=0; i<5; i++) {
            difQ1[i] = lam[i]*dW[i];
        }

        difQ2[0] = 0.0;
        difQ2[1] = 0.0;
        difQ2[2] = 0.0;
        difQ2[3] = 0.0;
        difQ2[4] = 0.0;

        for (int i=0; i<5; i++) {
            for (int j=0; j<5; j++) {
                difQ2[i] += R[i][j]*difQ1[j];
            }
        }

        flow_float mdot = 0.5*(ro_L*(Ux_L*nx+Uy_L*ny+Uz_L*nz) + ro_R*(Ux_R*nx+Uy_R*ny+Uz_R*nz))*sss ;
        massflux[ip] = mdot;

        flow_float res_ro_temp   = mdot                                        -0.5*difQ2[0]*sss;
        flow_float res_roe_temp  = mdot*0.5*(Ht_L + Ht_R)                      -0.5*difQ2[1]*sss;
        flow_float res_roUx_temp = mdot*0.5*(Ux_L + Ux_R) +0.5*(P_L + P_R)*sxx -0.5*difQ2[2]*sss;
        flow_float res_roUy_temp = mdot*0.5*(Uy_L + Uy_R) +0.5*(P_L + P_R)*syy -0.5*difQ2[3]*sss;
        flow_float res_roUz_temp = mdot*0.5*(Uz_L + Uz_R) +0.5*(P_L + P_R)*szz -0.5*difQ2[4]*sss;



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

__global__ void convectiveFlux_boundary_d // slau
( 
 flow_float ga,

  // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
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

 // bvar
 flow_float* rob ,
 flow_float* roUxb ,
 flow_float* roUyb ,
 flow_float* roUzb ,
 flow_float* roeb ,
 flow_float* Uxb ,
 flow_float* Uyb ,
 flow_float* Uzb ,
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb ,
 
 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe   
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        flow_float dro;
        flow_float dP;
        flow_float dUx;
        flow_float dUy;
        flow_float dUz;
        flow_float dU;

        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        geom_float nx = sxx/sss;
        geom_float ny = syy/sss;
        geom_float nz = szz/sss;

        flow_float ccx_0 = ccx[ic];
        flow_float ccy_0 = ccy[ic];
        flow_float ccz_0 = ccz[ic];

        flow_float pcx_1 = pcx[ip];
        flow_float pcy_1 = pcy[ip];
        flow_float pcz_1 = pcz[ip];

        flow_float dcc_x = pcx_1 - ccx_0;
        flow_float dcc_y = pcy_1 - ccy_0;
        flow_float dcc_z = pcz_1 - ccz_0;
        flow_float dcc   = sqrt(dcc_x*dcc_x +dcc_y*dcc_y +dcc_z*dcc_z) ;

        flow_float ro_L = ro[ic];
        flow_float Ux_L = Ux[ic];
        flow_float Uy_L = Uy[ic];
        flow_float Uz_L = Uz[ic];
        flow_float P_L  = Ps[ic];
        flow_float Ht_L = Ht[ic];

        flow_float ro_R= rob[ib]; 
        flow_float Ux_R= Uxb[ib];
        flow_float Uy_R= Uyb[ib];
        flow_float Uz_R= Uzb[ib];
        flow_float P_R = Psb[ib];
        //flow_float Ht_R= Htb[ib];
        flow_float Ht_R= (roeb[ib] + P_R)/ro_R;

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

        //flow_float Vn_p = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float Vn_m = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        //flow_float roVn_p = (ro_L*Ux_L*sxx + ro_L*Uy_L*syy + ro_L*Uz_L*szz)/sss;
        //flow_float roVn_m = (ro_R*Ux_R*sxx + ro_R*Uy_R*syy + ro_R*Uz_R*szz)/sss;

        //flow_float VnL = Vn_p;
        //flow_float VnR = Vn_m;

        //flow_float VnL_abs = abs(VnL);
        //flow_float VnR_abs = abs(VnR);

//TODO: change c
        //flow_float c_hat = 0.5*(sonic[ic] + sonic[ig]);

        //flow_float M_p = Vn_p/c_hat;
        //flow_float M_m = Vn_m/c_hat;

        //flow_float ro_del= ro_m - ro_p;
        //flow_float Vn_del= Vn_m - Vn_p;
        //flow_float P_del = P_m  - P_p;

        //flow_float beta_p, beta_m;

        //if (abs(M_p)>=1.0){
        //    beta_p = 0.5*(M_p + abs(M_p))/M_p;
        //} else {
        //    beta_p = 0.25*pow((M_p+1.0),2.0)*(2.0-M_p);
        //}

        //if (abs(M_m)>=1.0){
        //    beta_m = 0.5*(M_m - abs(M_m))/M_m;
        //} else {
        //    beta_m = 0.25*pow((M_m-1.0),2.0)*(2.0+M_m);
        //}

        //flow_float zero = 0.0;
        //flow_float one  = 1.0;
        //flow_float half = 0.5;

        //flow_float g = -max(min(M_p,zero),-one)*min(max(M_m,zero),one);
        //flow_float Vn_hat_abs   = (ro_p*abs(Vn_p) + ro_m*abs(Vn_m))/(ro_p + ro_m);
        //flow_float Vn_hat_p_abs = (1.0-g)*Vn_hat_abs + g*VnL_abs;
        //flow_float Vn_hat_m_abs = (1.0-g)*Vn_hat_abs + g*VnR_abs;

        //flow_float M_hat = min(one, sqrt(half*(u_m*u_m+v_m*v_m+w_m*w_m +u_p*u_p+v_p*v_p+w_p*w_p))/c_hat);
        //flow_float chi = (1.0-M_hat)*(1.0-M_hat);

        flow_float p_tilde = P_R;

        flow_float mdot = sss*(ro_R*Vn_m);

        flow_float res_ro_temp   = mdot;
        flow_float p_tilde_r = p_tilde - d_pRef;   // free-stream 保存: 基準静圧を差し引く (境界面)
        flow_float res_roUx_temp = 0.5*(mdot+abs(mdot))*u_p +0.5*(mdot-abs(mdot))*u_m +p_tilde_r*sxx;
        flow_float res_roUy_temp = 0.5*(mdot+abs(mdot))*v_p +0.5*(mdot-abs(mdot))*v_m +p_tilde_r*syy;
        flow_float res_roUz_temp = 0.5*(mdot+abs(mdot))*w_p +0.5*(mdot-abs(mdot))*w_m +p_tilde_r*szz;
        flow_float res_roe_temp  = 0.5*(mdot+abs(mdot))*h_p +0.5*(mdot-abs(mdot))*h_m ;

        atomicAdd(&res_ro[ic]  , -res_ro_temp);
        atomicAdd(&res_roUx[ic], -res_roUx_temp);
        atomicAdd(&res_roUy[ic], -res_roUy_temp);
        atomicAdd(&res_roUz[ic], -res_roUz_temp);
        atomicAdd(&res_roe[ic] , -res_roe_temp);

    }

    __syncthreads();
}




void convectiveFlux_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    // free-stream 保存: 基準静圧 pRef を device 定数へ転送 (既定 0.0 でビット不変)
    {
        flow_float pRef_h = static_cast<flow_float>(cfg.pRef);
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_pRef, &pRef_h, sizeof(flow_float)));
    }

    // D2a/D4 診断フラグを env から 1 度だけ device へ設定 (既定 off = ビット不変)。
    {
        static bool s_init = false;
        if (!s_init) {
            int c1 = 0, clog = 0; flow_float th = 0.05f, lth = 0.3f, blend = 0.0f;
            if (const char* e = getenv("FORGE_CONTACT_1ST"))       c1   = atoi(e);
            if (const char* e = getenv("FORGE_CONTACT_THRESH"))    th   = (flow_float)atof(e);
            if (const char* e = getenv("FORGE_CONTACT_LOG"))       clog = atoi(e);
            if (const char* e = getenv("FORGE_CONTACT_LOG_THRESH"))lth  = (flow_float)atof(e);
            if (const char* e = getenv("FORGE_CONTACT_BLEND"))     blend= (flow_float)atof(e);
            int fthy = 0;
            if (const char* e = getenv("FORGE_FACE_THERMOY"))      fthy = atoi(e);
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_faceThermoY,      &fthy, sizeof(int)));
            const int sfr = cfg.speciesFaceReconstruction;   // config 由来 (env でない)
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_speciesFaceRecon, &sfr,  sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contact1st,       &c1,   sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactThresh,    &th,   sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactLog,       &clog, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactLogThresh, &lth,  sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactBlend,     &blend,sizeof(flow_float)));
            s_init = true;
        }
    }

    // initialize
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_ro"]  , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roe"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.p_d["massflux"], 0.0, msh.nPlanes*sizeof(flow_float)));

    dim3 dimGrid_normal_halo = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));

    // node-centered 弱形式 (壁のみ): 主対流ループは「内部 + 非壁境界(ゴースト)」を処理し、末尾に並ぶ壁境界 plane
    // (nWallHaloPlanes) を除外する。壁は別途 pressure-only の convectiveFlux_boundary_d (bvar u=0) で扱う
    // (壁ノードは壁上に乗りミラーゴーストの幾何が退化=fx≠0.5 で質量貫通し発散するため)。
    // inlet/outlet/axis は従来どおり主ループ+ゴーストで処理し既存の near-axis corner 修正を保つ。
    // cell モードは全 plane を主ループで処理 (従来どおり)。
    const geom_int convPlaneBound = (cfg.discretization == "node")
                                  ? (msh.nNormal_halo_Planes - msh.nWallHaloPlanes) : msh.nNormal_halo_Planes;

    // -----------------------
    // *** sum over planes ***
    // -----------------------
    // 非平衡凝縮: 二相エネルギー流束補正用の総液相分率 g。off / 未登録なら nullptr (補正なし)。
    // 現状 nCondSpecies=1 (g_0 が総 g)。多成分は総和 g 配列を別途用意する (TODO)。
    flow_float* cond_g = nullptr;
    if (cfg.condensation == 1 && var.nCondSpeciesRegistered >= 1) cond_g = var.c_d["g_0"];

    if (cfg.solver == "SLAU" || cfg.solver == "SLAU2") {
        int slauVariant = (cfg.solver == "SLAU2") ? 2 : 1;
        SLAU_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.convMethod, cfg.limiter, slauVariant,
            cfg.lowMachPrecond, cfg.precondEps,
            cfg.lowMachThornber,

            cfg.gamma,
            cfg.thermalMethod, thermo_species_device_ptr(), cfg.nSpecies,
            species_roY_device_ptr(),
            species_Y_device_ptr(), species_dYdx_device_ptr(), species_dYdy_device_ptr(), species_dYdz_device_ptr(),
            var.c_d["Rmix"],
            cfg.cp, cond_g, var.c_d["T"], cfg.condModel,

            // mesh structure
            msh.nCells,
            msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
            convPlaneBound, msh.normal_halo_planes_d,
            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
            var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],
            var.p_d["massflux"],

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

            var.c_d["limiter_ro"]  ,
            var.c_d["limiter_Ux"]  ,
            var.c_d["limiter_Uy"]  ,
            var.c_d["limiter_Uz"]  ,
            var.c_d["limiter_P"]  ,

            var.c_d["ducros"]  ,

            // gradient
            var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
            var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
            var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
            var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
            var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"]

        ) ;

    } else if (cfg.solver == "HLLE") {
        HLLE_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.convMethod, cfg.limiter,
            cfg.gamma,
            cfg.cp, cond_g, var.c_d["T"], cfg.condModel,

            msh.nCells,
            msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
            convPlaneBound, msh.normal_halo_planes_d,
            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
            var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],
            var.p_d["massflux"],

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

            var.c_d["limiter_ro"]  ,
            var.c_d["limiter_Ux"]  ,
            var.c_d["limiter_Uy"]  ,
            var.c_d["limiter_Uz"]  ,
            var.c_d["limiter_P"]  ,

            var.c_d["ducros"]  ,

            var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
            var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
            var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
            var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
            var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"]

        );

    } else if (cfg.solver == "ROE") {
        //ROE_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        ROE_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.convMethod, cfg.limiter,
            cfg.gamma,
            cfg.thermalMethod, thermo_species_device_ptr(), cfg.nSpecies,
            cfg.cp, cond_g, var.c_d["T"], cfg.condModel,

            // mesh structure
            msh.nCells,
            msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
            convPlaneBound, msh.normal_halo_planes_d,
            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
            var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  
            var.p_d["massflux"],

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

            var.c_d["limiter_ro"]  ,
            var.c_d["limiter_Ux"]  ,
            var.c_d["limiter_Uy"]  ,
            var.c_d["limiter_Uz"]  ,
            var.c_d["limiter_P"]  ,
 
            var.c_d["ducros"]  ,

            // gradient
            var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
            var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
            var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
            var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
            var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"]

        );

    } else {
        std::cerr << "Error: unsupported solver name " << cfg.solver
                  << " (enabled: SLAU, HLLE, ROE)" << std::endl;
        exit(EXIT_FAILURE);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );


    // cell モード: SLAU/ROE/HLLE は主ループが normal_halo_planes_d 経由で全境界 plane をゴースト処理するため、
    // この dedicated 境界カーネルは二重計上になるのでスキップする。
    // node モード (弱形式・壁のみ): 主ループは壁 plane を除外した (convPlaneBound) ので、壁の境界寄与を
    // この convectiveFlux_boundary_d が bvar (u=0 → mdot=0, pressure-only) から担う。非壁 (inlet/outlet/axis)
    // は主ループ+ゴーストのままなので本カーネルからは除外する。
    const bool nodeMode = (cfg.discretization == "node");
    bool skipBoundaryFluxKernel = (!nodeMode)
                               && (cfg.solver == "SLAU"
                                || cfg.solver == "SLAU2"
                                || cfg.solver == "ROE"
                                || cfg.solver == "HLLE");
    auto isWallKind = [](const std::string& k){ return k == "wall" || k == "wall_isothermal"; };

    for (auto& bc : msh.bconds)
    {
        if (bc.bcondKind == "periodic") {
            continue;
        }
        if (skipBoundaryFluxKernel) {
            continue;
        }
        // node モードでは壁 bcond のみ弱形式境界 flux を適用 (非壁は主ループで処理済み)。
        if (nodeMode && !isWallKind(bc.bcondKind)) {
            continue;
        }
        convectiveFlux_boundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
            cfg.gamma,
            // mesh structure
            bc.iPlanes.size(),
            bc.map_bplane_plane_d,  
            bc.map_bplane_cell_d,  
            bc.map_bplane_cell_ghst_d,

            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
            var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

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

            bc.bvar_d["ro"],
            bc.bvar_d["roUx"],
            bc.bvar_d["roUy"],
            bc.bvar_d["roUz"],
            bc.bvar_d["roe"],
            bc.bvar_d["Ux"],
            bc.bvar_d["Uy"],
            bc.bvar_d["Uz"],
            bc.bvar_d["Tt"],
            bc.bvar_d["Pt"],
            bc.bvar_d["Ts"],
            bc.bvar_d["Ps"],
 
            var.c_d["res_ro"] ,
            var.c_d["res_roUx"] ,
            var.c_d["res_roUy"] ,
            var.c_d["res_roUz"] ,
            var.c_d["res_roe"]  
        ) ;
    }


    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );


}