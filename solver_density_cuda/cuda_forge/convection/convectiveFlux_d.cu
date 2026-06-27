#include <cstdlib>
#include "convectiveFlux_d.cuh"
#include "lowMachPrecond_d.cuh"
#include "speciesTransport_d.cuh"  // species_roY_device_ptr()
#include "condensationProperties_d.cuh"  // n2_latent (二相エネルギー流束の潜熱補正)

// free-stream 保存: 対流流束の圧力項を (p_tilde - d_pRef)*s で組むための基準静圧。
// wrapper で cfg.pRef を cudaMemcpyToSymbol。既定 0.0 で従来挙動 (ビット不変)。
// 非直交メッシュで大きな p*s を float32 加算する際の桁落ち(metric closure 由来の
// 偽運動量源)を抑える。詳細: plans/active/convection-freestream-preserving-flux.md
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
// speciesFaceReconstruction==1: Y_s を ρ/Y 勾配 + min(ψ_ρ,ψ_Y) で face へ再構成し thermo/species 流束で
// 同一 face 組成を使う (proper S2/S3)。wrapper で cfg.speciesFaceReconstruction を設定。
__device__ int g_speciesFaceRecon = 0;
// 再構成 Y が [0,1] を外れた face 数のカウンタ (診断; clamp 前に atomicAdd)。
__device__ unsigned long long g_speciesOvershoot = 0;
// multispeciesRhoYCommonLimiter==1 (opt-in 診断): ρ と全 species に共通リミタ
//   ψ_ρY = min(ψ_ρ, min_s ψ_Y_s) を適用し、ρ_f=ρ(Y_f) の熱力学整合だけを切り分ける。
//   p・速度は従来どおり各自のリミタ。default 0 でビット不変。
__device__ int g_rhoYCommonLim = 0;
// rho-Y 共通リミタ診断カウンタ (すべて g_rhoYCommonLim 時のみ更新)。
__device__ int g_psiRhoY_min_scaled = 2000000;          // min(ψ_ρY)·1e6 (atomicMin)
__device__ unsigned long long g_psiRhoY_lt001 = 0;      // ψ_ρY<0.01 の face-side 数
__device__ unsigned long long g_psiRhoY_lt01  = 0;      // ψ_ρY<0.1  の face-side 数
__device__ unsigned long long g_rhoYMinByRho     = 0;   // min を ρ が決めた face-side 数
__device__ unsigned long long g_rhoYMinBySpecies = 0;   // min を species が決めた face-side 数
__device__ unsigned long long g_rhoYFallback = 0;       // 非実現可能で cell 値へ fallback した face-side 数
__device__ int g_Yface_min_scaled =  2000000;           // min(Y_face)·1e6 (atomicMin)
__device__ int g_Yface_max_scaled = -2000000;           // max(Y_face)·1e6 (atomicMax)

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

#include "convectiveFlux_slau_d.inc.cuh"

inline __device__ flow_float sign_sano(flow_float x)
{
    return x > 0 ? 1 : (x<0 ? -1 : 0);
}

#include "legacy/convectiveFlux_ausm_keep_d.inc.cuh"



#include "convectiveFlux_hlle_d.inc.cuh"



#include "convectiveFlux_roe_d.inc.cuh"

#include "legacy/convectiveFlux_keepslau_d.inc.cuh"

#include "convectiveFlux_boundary_d.inc.cuh"




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
            const int rycl = cfg.multispeciesRhoYCommonLimiter;   // config 由来 (opt-in 診断)
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYCommonLim,    &rycl, sizeof(int)));
            const unsigned long long zero = 0ULL;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_speciesOvershoot, &zero, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contact1st,       &c1,   sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactThresh,    &th,   sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactLog,       &clog, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactLogThresh, &lth,  sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactBlend,     &blend,sizeof(flow_float)));
            s_init = true;
        }
    }

    // rho-Y 共通リミタ診断: 一定間隔で device カウンタを読み出して 1 行印字しリセット (opt-in 時のみ)。
    if (cfg.multispeciesRhoYCommonLimiter == 1) {
        static int s_rycl_call = 0;
        const int interval = 200;
        if ((s_rycl_call % interval) == 0) {
            int psimin=0, ymin=0, ymax=0;
            unsigned long long lt001=0, lt01=0, byrho=0, bysp=0, fb=0, ovs=0;
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&psimin, g_psiRhoY_min_scaled, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&lt001,  g_psiRhoY_lt001, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&lt01,   g_psiRhoY_lt01,  sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&byrho,  g_rhoYMinByRho,  sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&bysp,   g_rhoYMinBySpecies, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&fb,     g_rhoYFallback,  sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&ovs,    g_speciesOvershoot, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&ymin,   g_Yface_min_scaled, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&ymax,   g_Yface_max_scaled, sizeof(int)));
            printf("RHOYLIM call=%d minPsiRhoY=%.4f lt0.01=%llu lt0.1=%llu minBy[rho=%llu sp=%llu] overshoot=%llu fallback=%llu Yface=[%.5f,%.5f]\n",
                   s_rycl_call, (double)psimin*1e-6, lt001, lt01, byrho, bysp, ovs, fb,
                   (double)ymin*1e-6, (double)ymax*1e-6);
            // 次区間用にリセット (min/max は両端へ)。
            const int big=2000000, nbig=-2000000; const unsigned long long z=0ULL;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_psiRhoY_min_scaled, &big, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_Yface_min_scaled,   &big, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_Yface_max_scaled,   &nbig, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_psiRhoY_lt001, &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_psiRhoY_lt01,  &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYMinByRho,  &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYMinBySpecies, &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYFallback,  &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_speciesOvershoot, &z, sizeof(unsigned long long)));
        }
        s_rycl_call++;
    }

    // initialize
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_ro"]  , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roe"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.p_d["massflux"], 0.0, msh.nPlanes*sizeof(flow_float)));

    dim3 dimGrid_normal_halo = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));

    // node-centered 弱形式 (Phase 2): 主対流ループは「内部 + periodic」のみを処理し、末尾に並ぶ全非 periodic 境界 plane
    // (nBoundaryHaloPlanes = wall+inlet+outlet+slip...) を除外する。全境界は別途 convectiveFlux_boundary_d が bvar を
    // R 状態とする弱形式で担う。これは境界ノードが物理境界上に乗る node-centered で ghost を主ループの右状態に食わせる
    // と退化幾何 (d_along_n=0) により出口列/コーナーで近壁 BL が崩壊するため (case/26 で実証)。
    // cell モードは全 plane を主ループでゴースト処理 (従来どおり)。
    const geom_int convPlaneBound = (cfg.discretization == "node")
                                  ? (msh.nNormal_halo_Planes - msh.nBoundaryHaloPlanes) : msh.nNormal_halo_Planes;

    // (検証A) node 弱形式の plane 振り分けを 1 度だけログ: 主ループ面 = 内部(+periodic)、弱形式面 = 全境界半割面。
    if (cfg.discretization == "node") {
        static bool s_planeLogDone = false;
        if (!s_planeLogDone) {
            s_planeLogDone = true;
            const geom_int mainPlanes = convPlaneBound;
            const geom_int periodicInMain = mainPlanes - msh.nNormalPlanes;  // 主ループのうち内部以外 = periodic
            printf("[node weak-form] conv main planes = %ld (internal %ld + periodic %ld), "
                   "boundary weak planes = %ld (wall %ld + non-wall %ld)\n",
                   (long)mainPlanes, (long)msh.nNormalPlanes, (long)periodicInMain,
                   (long)msh.nBoundaryHaloPlanes, (long)msh.nWallHaloPlanes,
                   (long)(msh.nBoundaryHaloPlanes - msh.nWallHaloPlanes));
        }
    }

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
            species_limiterY_device_ptr(),
            (cfg.speciesFaceReconstruction >= 2) ? species_Yface_alloc(msh.nPlanes) : nullptr,
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
    gpuErrchkKernelSync();


    // cell モード: SLAU/ROE/HLLE は主ループが normal_halo_planes_d 経由で全境界 plane をゴースト処理するため、
    // この dedicated 境界カーネルは二重計上になるのでスキップする。
    // node モード (弱形式・Phase 2): 主ループは全非 periodic 境界 plane を除外した (convPlaneBound=内部+periodic) ので、
    // 全境界 (wall/inlet/outlet/slip...) の境界寄与を本 convectiveFlux_boundary_d が bvar を R 状態とする弱形式で担う。
    // (壁: bvar u=0 → mdot=0, pressure-only。出口: bvar Uxb/Psb。slip: Uxb=0+Psb=P[ic]。入口: 固定状態。)
    const bool nodeMode = (cfg.discretization == "node");
    bool skipBoundaryFluxKernel = (!nodeMode)
                               && (cfg.solver == "SLAU"
                                || cfg.solver == "SLAU2"
                                || cfg.solver == "ROE"
                                || cfg.solver == "HLLE");

    for (auto& bc : msh.bconds)
    {
        if (bc.bcondKind == "periodic") {
            continue;
        }
        if (skipBoundaryFluxKernel) {
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
            (nodeMode ? var.p_d["massflux"] : nullptr),   // node 弱形式境界のみ massflux 書き戻し (cell は主ループが書く)

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
    gpuErrchkKernelSync();


}