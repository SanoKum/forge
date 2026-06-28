#include "calcGradient_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"


__global__ void WALE_d
( 
 // mesh structure
 geom_int nCells,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* Ht  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* vis_turb , flow_float* wall_dist
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells ) {

        flow_float karman = 0.41;
        //flow_float Cw = 0.5;
        flow_float Cw = 0.325;

        geom_float volume = vol[ic];


        flow_float g11 = dUxdx[ic];
        flow_float g12 = dUxdy[ic];
        flow_float g13 = dUxdz[ic];
        flow_float g21 = dUydx[ic];
        flow_float g22 = dUydy[ic];
        flow_float g23 = dUydz[ic];
        flow_float g31 = dUzdx[ic];
        flow_float g32 = dUzdy[ic];
        flow_float g33 = dUzdz[ic];

        flow_float Sd11 = 0.5*(g11*g11+g11*g11) - 1.0/3.0*g11*g11;
        flow_float Sd12 = 0.5*(g12*g12+g21*g21);
        flow_float Sd13 = 0.5*(g13*g13+g31*g31);
        flow_float Sd21 = 0.5*(g21*g21+g12*g12);
        flow_float Sd22 = 0.5*(g22*g22+g22*g22) - 1.0/3.0*g22*g22;
        flow_float Sd23 = 0.5*(g23*g23+g32*g32);
        flow_float Sd31 = 0.5*(g31*g31+g13*g13);
        flow_float Sd32 = 0.5*(g32*g32+g23*g23);
        flow_float Sd33 = 0.5*(g33*g33+g33*g33) - 1.0/3.0*g33*g33;


        flow_float S11 = 0.5*(g11+g11);
        flow_float S12 = 0.5*(g12+g21);
        flow_float S13 = 0.5*(g13+g31);
        flow_float S21 = 0.5*(g21+g12);
        flow_float S22 = 0.5*(g22+g22);
        flow_float S23 = 0.5*(g23+g32);
        flow_float S31 = 0.5*(g31+g13);
        flow_float S32 = 0.5*(g32+g23);
        flow_float S33 = 0.5*(g33+g33);

        flow_float SdijSdij = Sd11*Sd11 + Sd12*Sd12 + Sd13*Sd13
                            + Sd21*Sd21 + Sd22*Sd22 + Sd23*Sd23
                            + Sd31*Sd31 + Sd32*Sd32 + Sd33*Sd33;

        flow_float SijSij = S11*S11 + S12*S12 + S13*S13
                          + S21*S21 + S22*S22 + S23*S23
                          + S31*S31 + S32*S32 + S33*S33;

        geom_float d = wall_dist[ic];
        flow_float Ls = min(karman*d , Cw*pow(vol[ic],1.0/3.0));

        // WALE は一様流（勾配ゼロ）で分子分母とも 0 になり 0/0=NaN を生む。標準的に分母へ
        // 小さな正則化項を加える（無勾配では vis_turb=0 となり物理的に正しい）。
        const flow_float wale_denom = pow(SijSij, 5.0/2.0) + pow(SdijSdij, 5.0/4.0)
                                    + static_cast<flow_float>(1.0e-30);
        vis_turb[ic] =  ro[ic]*Ls*Ls*pow(SdijSdij, 3.0/2.0)/wale_denom;
    }
}


// SST 渦粘性: mu_t = rho * a1 * k / max(a1 * omega, S * F2)
// 軸対称 (isAxisymmetric=1) では planar 速度勾配に現れないフープひずみ
// S_thetatheta = u_r/r (= axisym_uy_over_r) を S^2 に加える (methods/turbulence/theory.md §7.2)。
__global__ void sst_eddy_viscosity_d(
    geom_int nCells,
    flow_float* ro,
    flow_float* k,
    flow_float* omega,
    flow_float* vis_lam,
    flow_float* wall_dist,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    int isAxisymmetric,
    flow_float* axisym_uy_over_r,
    flow_float* vis_turb)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    constexpr flow_float a1       = static_cast<flow_float>(0.31);
    constexpr flow_float betaStar = static_cast<flow_float>(0.09);
    constexpr flow_float kSmall   = static_cast<flow_float>(1.0e-12);

    const flow_float rho    = max(ro[ic], kSmall);
    const flow_float k_c    = max(k[ic],  static_cast<flow_float>(0.0));
    const flow_float w_c    = max(omega[ic], kSmall);
    const flow_float mu_lam = vis_lam[ic];
    const flow_float y      = max(wall_dist[ic], kSmall);
    const flow_float nu     = mu_lam / rho;

    // ひずみ速度の大きさ S = sqrt(2 Sij Sij)
    const flow_float S11 = dUxdx[ic];
    const flow_float S22 = dUydy[ic];
    const flow_float S33 = dUzdz[ic];
    const flow_float S12 = static_cast<flow_float>(0.5) * (dUxdy[ic] + dUydx[ic]);
    const flow_float S13 = static_cast<flow_float>(0.5) * (dUxdz[ic] + dUzdx[ic]);
    const flow_float S23 = static_cast<flow_float>(0.5) * (dUydz[ic] + dUzdy[ic]);
    flow_float S_sq = static_cast<flow_float>(2.0) *
        (S11*S11 + S22*S22 + S33*S33 +
         static_cast<flow_float>(2.0) * (S12*S12 + S13*S13 + S23*S23));
    // 軸対称: フープひずみ S_thetatheta = u_r/r を加算 (2 * S_tt^2)
    if (isAxisymmetric) {
        const flow_float S_tt = axisym_uy_over_r[ic];
        S_sq += static_cast<flow_float>(2.0) * S_tt * S_tt;
    }
    const flow_float S_mag = sqrt(max(S_sq, static_cast<flow_float>(0.0)));

    // F2 ブレンド関数
    const flow_float arg2_a = static_cast<flow_float>(2.0) * sqrt(k_c) / (betaStar * w_c * y);
    const flow_float arg2_b = static_cast<flow_float>(500.0) * nu / (w_c * y * y);
    const flow_float arg2   = max(arg2_a, arg2_b);
    const flow_float F2     = tanh(arg2 * arg2);

    vis_turb[ic] = rho * a1 * k_c / max(a1 * w_c, S_mag * F2);
}


// ===========================================================================
// SST-DDES length scale (methods/turbulence §8, plan turbulence-iddes-sst.md §4.4)
// ---------------------------------------------------------------------------
// Spalart et al. (2006) シールド関数 f_d で剥離域のみ LES、付着 BL は RANS に保つ。
//   r_d = (ν_t + ν) / (κ² d² √(∂Ui/∂xj·∂Ui/∂xj))      κ=0.41
//   f_d = 1 - tanh((8 r_d)³)
//   l_RANS = √k / (β*^{1/4} ω),  l_LES = C_DES·Δ
//   l_DDES = l_RANS - f_d·max(0, l_RANS - C_DES·Δ)
//
// ★ float32 ガード (plan §4.4): r_d は d² で効くため近壁・よどみ・一様流で飽和しやすい。
//   - 勾配大きさ / 分母に floor を入れ 0 割を防ぐ (一様流 ∇u→0 では r_d→大→f_d→0=RANS が正しい極限)。
//   - r_d を [0, 10] に clamp し f_d が 0/1 へ過度に張り付くのを防ぐ。
//   - r_d / f_d を必ず出力 (NaN が無くても全域張り付きは DES 不全の兆候・l_des だけでは見抜けない)。
//
// 依存順 (plan §4.7・最重要): r_d は ν_t=vis_turb/ρ を読むため、必ず sst_eddy_viscosity_d
// (vis_turb 計算) の後に呼ぶ。本カーネルは CV 抽象 (volume→Δ, wall_dist→d, ∇u, vis_turb) のみに
// 依存し cfg.discretization を参照しないので、node-centered (median-dual) でも無変更で成立する (§5.6)。
__device__ static flow_float compute_rd_ddes(
    flow_float nu_t, flow_float nu, flow_float y, flow_float grad2)
{
    constexpr flow_float kKarman    = static_cast<flow_float>(0.41);
    constexpr flow_float kGradFloor = static_cast<flow_float>(1.0e-30);
    constexpr flow_float kDenomFloor= static_cast<flow_float>(1.0e-30);
    const flow_float gradMag = sqrt(max(grad2, kGradFloor));
    flow_float denom = kKarman * kKarman * y * y * gradMag;
    denom = max(denom, kDenomFloor);
    const flow_float rd = (nu_t + nu) / denom;
    return min(rd, static_cast<flow_float>(10.0));
}

__device__ static flow_float compute_fd_ddes(flow_float rd)
{
    const flow_float arg = static_cast<flow_float>(8.0) * rd;
    return static_cast<flow_float>(1.0) - tanh(arg * arg * arg);
}

__global__ void compute_des_length_d(
    geom_int nCells,
    flow_float* ro,
    flow_float* k,
    flow_float* omega,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float* wall_dist,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* delta_les,
    flow_float C_DES,
    flow_float* l_des,
    flow_float* fd_shield,
    flow_float* rd_des)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    constexpr flow_float kSmall    = static_cast<flow_float>(1.0e-12);
    constexpr flow_float kBetaStar = static_cast<flow_float>(0.09);

    const flow_float rho  = max(ro[ic], kSmall);
    const flow_float k_c  = max(k[ic],  static_cast<flow_float>(0.0));
    const flow_float w_c  = max(omega[ic], kSmall);
    const flow_float nu   = vis_lam[ic] / rho;
    const flow_float nu_t = max(vis_turb[ic], static_cast<flow_float>(0.0)) / rho;
    const flow_float y    = max(wall_dist[ic], kSmall);

    // |∂Ui/∂xj|² = 全 9 成分の二乗和 (= S²+Ω²)。軸対称項は含めない (検証は 3D 直交・plan §4.7)。
    const flow_float g11 = dUxdx[ic], g12 = dUxdy[ic], g13 = dUxdz[ic];
    const flow_float g21 = dUydx[ic], g22 = dUydy[ic], g23 = dUydz[ic];
    const flow_float g31 = dUzdx[ic], g32 = dUzdy[ic], g33 = dUzdz[ic];
    const flow_float grad2 = g11*g11 + g12*g12 + g13*g13
                           + g21*g21 + g22*g22 + g23*g23
                           + g31*g31 + g32*g32 + g33*g33;

    const flow_float rd = compute_rd_ddes(nu_t, nu, y, grad2);
    const flow_float fd = compute_fd_ddes(rd);

    // l_RANS = √k/(β* ω)。これは SST k 消滅項 D_k = β*ρkω = ρk^{3/2}/l_RANS と厳密に整合する
    // 唯一の定義 (Strelets 2001 / Gritskevich 2012)。f_d=0 の完全シールド域で l_DDES=l_RANS となり
    // D_k が標準 SST に一致する (これを満たさないと mode1 はシールド域でも標準 SST に縮退せず MSD を起こす)。
    const flow_float l_rans = sqrt(k_c) / (kBetaStar * w_c);
    const flow_float l_les  = C_DES * delta_les[ic];
    const flow_float l_ddes = l_rans - fd * max(static_cast<flow_float>(0.0), l_rans - l_les);

    l_des[ic]     = max(l_ddes, static_cast<flow_float>(1.0e-20));
    fd_shield[ic] = fd;
    rd_des[ic]    = rd;
}


void turbulent_viscosity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.LESorRANS == 1) { // LES

        if (cfg.LESmodel == 1) {
            // node-centered (median-dual) でも成立: WALE_d は DOF (cell/node) ごとに勾配・体積・wall_dist を
            // 読むだけで cfg.discretization に依存しないため、SST と同じ normalcell グリッドで node も網羅する。
            WALE_d<<<cuda_cfg.dimGrid_normalcell , cuda_cfg.dimBlock>>> (
                // mesh structure
                msh.nCells,
                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

                // basic variables
                var.c_d["ro"] ,
                var.c_d["Ux"] ,
                var.c_d["Uy"] ,
                var.c_d["Uz"] ,
                var.c_d["P"]  , 
                var.c_d["T"]  ,
                var.c_d["Ht"] ,

                // gradient
                var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
                var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
                var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
                var.c_d["vis_turb"] , var.c_d["wall_dist"]
            ) ;
        }

    } else if (cfg.LESorRANS == 2 && cfg.RANSmodel == 1) { // RANS SST
        sst_eddy_viscosity_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells,
            var.c_d["ro"],
            var.c_d["k"],
            var.c_d["omega"],
            var.c_d["vis_lam"],
            var.c_d["wall_dist"],
            var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
            var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
            var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],
            cfg.isAxisymmetric,
            var.c_d["axisym_uy_over_r"],
            var.c_d["vis_turb"]);

        // SST-DES: vis_turb 計算直後に l_DDES とシールド診断を計算する (依存順・plan §4.7)。
        // DESmode==0 (既定) では呼ばない → 既存 SST と完全ビット不変。
        // C_DES は初版で C_DES_kw 固定 (F1 ブレンドは sst_eddy が F1 を持たないため Phase 2・plan §3.2)。
        if (cfg.DESmode > 0) {
            compute_des_length_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
                msh.nCells,
                var.c_d["ro"],
                var.c_d["k"],
                var.c_d["omega"],
                var.c_d["vis_lam"],
                var.c_d["vis_turb"],
                var.c_d["wall_dist"],
                var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
                var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
                var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],
                var.c_d["delta_les"],
                cfg.C_DES_kw,
                var.c_d["l_des"],
                var.c_d["fd_shield"],
                var.c_d["rd_des"]);
        }

    } else {
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["vis_turb"], 0.0, msh.nCells_all*sizeof(flow_float)));
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

}