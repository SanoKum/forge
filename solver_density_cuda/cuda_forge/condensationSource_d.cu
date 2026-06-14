#include "condensationTransport_d.cuh"  // wrapper 宣言 + 必要なクラス
#include "condensationSource_d.cuh"

#include <string>

namespace {

inline bool condensationEnabled(const solverConfig& cfg)
{
    return cfg.condensation == 1 && cfg.nCondSpecies >= 1;
}

// 1 凝縮種 (N2) の相変化ソースを残差へ加える (初期実装: 明示的 + 安定化 clamp/limiter)。
//   res_ro<φ> += S_φ · V       (体積積分ソース、advection 残差へ加算)
// 一温度 T_v=T_d=T。J, r*, dr/dt は現在セル状態から一度だけ評価して freeze する。
// 安定化 (この初期実装の主眼):
//   - J 上限 limiter (overflow 回避)
//   - dr/dt < 0 は 0 に clamp (蒸発は初期実装で扱わない)
//   - r̄ <= r* のとき成長を止める (亜臨界は成長しない)
//   - g,Q0,Q1,Q2 の非負は時間積分の floor で担保 (ソースは加算のみ)
//   - 1 step の Δg・潜熱 ΔT・蒸気枯渇を θ で律速し、全モーメントを同じ θ で縮小 (g=Q3 整合保持)
// src_jac はソース由来を 0 とし、時間項+移流項 (transport_diag) の対角のみで安定化する
//   (ソースの T 依存 dJ/dT, d(drdt)/dT を含む完全線形化は後続)。
__global__ void condensation_source_d(
    geom_int nCells,
    flow_float cp, flow_float gamma, double M,
    double Jmax, double dg_max, double dT_max,
    geom_float* vol, flow_float* dt_local,
    flow_float* T, flow_float* P, flow_float* ro,
    flow_float* rog, flow_float* roQ0, flow_float* roQ1, flow_float* roQ2,
    flow_float* res_rog, flow_float* res_roQ0, flow_float* res_roQ1, flow_float* res_roQ2,
    flow_float* sj_g, flow_float* sj_Q0, flow_float* sj_Q1, flow_float* sj_Q2)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    // ソース由来の src_jac は初期実装では 0 (時間項+移流項の対角のみで安定化)
    sj_Q0[ic] = 0.0; sj_Q1[ic] = 0.0; sj_Q2[ic] = 0.0; sj_g[ic] = 0.0;

    const double rod = (double)ro[ic];
    if (rod <= 1.0e-20) return;

    const double cv   = (double)cp/(double)gamma;
    const double Rgas = ((double)gamma - 1.0)*cv;
    const double Td   = (double)T[ic];
    const double Pd   = (double)P[ic];

    double g = (double)rog[ic]/rod;
    if (g < 0.0) g = 0.0;
    if (g > 0.99) g = 0.99;
    double q0 = (double)roQ0[ic]; if (q0 < 0.0) q0 = 0.0;
    double q1 = (double)roQ1[ic]; if (q1 < 0.0) q1 = 0.0;
    double q2 = (double)roQ2[ic]; if (q2 < 0.0) q2 = 0.0;

    // --- 核生成 J, r* (freeze) ---
    double rho_v = (1.0 - g)*rod; if (rho_v < 0.0) rho_v = 0.0;
    double J, rstar;
    cond_nucleation_N2(Td, Pd, rho_v, Rgas, M, &J, &rstar);
    if (J > Jmax) J = Jmax;          // 上限 limiter
    if (J < 0.0)  J = 0.0;

    // --- 成長 dr/dt (freeze, 亜臨界停止 + 蒸発 clamp) ---
    double r_bar = (q0 > 1.0e-30) ? (q1/q0) : rstar;
    double drdt = 0.0;
    if (q0 > 1.0e-30 && rstar > 0.0 && r_bar > rstar) {
        drdt = cond_growth_N2(Td, Pd, r_bar, rstar, Rgas);
        if (drdt < 0.0) drdt = 0.0;  // 蒸発は扱わない
    }

    // --- ソースベクトル (g=Q3 整合) ---
    const double rho_l = n2_rho_cond(Td);
    double SQ0 = J;
    double SQ1 = J*rstar + q0*drdt;
    double SQ2 = J*rstar*rstar + 2.0*q1*drdt;
    double Sg  = (4.0/3.0)*COND_PI*rho_l*(J*rstar*rstar*rstar + 3.0*q2*drdt);
    if (Sg < 0.0) Sg = 0.0;

    // --- 1 step 律速 θ: Δg, 潜熱 ΔT, 蒸気枯渇 ---
    const double dt = (double)dt_local[ic];
    const double L  = n2_latent(Td);
    double theta = 1.0;
    if (Sg > 0.0 && dt > 0.0) {
        const double dg  = Sg*dt/rod;                       // Δg (= Δ(ρg)/ρ)
        const double dTl = dg*L/(cv + g*Rgas);              // 潜熱による ΔT
        const double avail = (1.0 - g);                     // 利用可能蒸気分率
        if (dg  > dg_max)        theta = fmin(theta, dg_max/dg);
        if (dTl > dT_max)        theta = fmin(theta, dT_max/dTl);
        if (dg  > 0.9*avail)     theta = fmin(theta, 0.9*avail/fmax(dg,1.0e-300));
    }
    SQ0 *= theta; SQ1 *= theta; SQ2 *= theta; Sg *= theta;

    const double v = (double)vol[ic];
    res_roQ0[ic] += (flow_float)(SQ0*v);
    res_roQ1[ic] += (flow_float)(SQ1*v);
    res_roQ2[ic] += (flow_float)(SQ2*v);
    res_rog[ic]  += (flow_float)(Sg *v);
}

}  // namespace

void condensationSource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!condensationEnabled(cfg)) return;

    const double M = condProps_N2().M;  // 現状 N2 のみ
    // 安定化パラメータ (初期実装の既定値。後で config 化可)
    const double Jmax   = 1.0e35;   // 核生成率上限 (overflow 回避)
    const double dg_max = 5.0e-3;   // 1 step あたり最大 Δg
    const double dT_max = 1.0;      // 1 step あたり潜熱による最大 ΔT [K]

    for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
        const std::string i = std::to_string(s);
        condensation_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells,
            cfg.cp, cfg.gamma, M,
            Jmax, dg_max, dT_max,
            var.c_d["volume"], var.c_d["dt_local"],
            var.c_d["T"], var.c_d["P"], var.c_d["ro"],
            var.c_d["rog_"+i], var.c_d["roQ0_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i],
            var.c_d["res_rog_"+i], var.c_d["res_roQ0_"+i], var.c_d["res_roQ1_"+i], var.c_d["res_roQ2_"+i],
            var.c_d["src_jac_g_"+i], var.c_d["src_jac_Q0_"+i], var.c_d["src_jac_Q1_"+i], var.c_d["src_jac_Q2_"+i]);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}
