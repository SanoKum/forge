#include "condensationTransport_d.cuh"  // wrapper 宣言 + 必要なクラス
#include "condensationSource_d.cuh"

#include <string>

namespace {

inline bool condensationEnabled(const solverConfig& cfg)
{
    return cfg.condensation == 1 && cfg.nCondSpecies >= 1;
}

// 凝縮種の蒸気分圧 p_v と蒸気密度 rho_v を返す。
//   pure-condensible (carrier=0, N2 Arthur): 気相=凝縮種。p_v=P(気相圧)、rho_v=(1-g)ρ。
//   carrier+condensible (carrier=1, H2O in N2): p_v=ρ(Y_w-g)R_w T、rho_v=ρ(Y_w-g)。
__host__ __device__ inline void cond_vapor_state(
    int carrier, double rod, double Pd, double Td, double g, double Yw, double Rw,
    double* pv, double* rho_v)
{
    if (carrier) {
        double yv = Yw - g; if (yv < 0.0) yv = 0.0;
        *rho_v = rod*yv;
        *pv    = rod*yv*Rw*Td;
    } else {
        double omg = 1.0 - g; if (omg < 0.0) omg = 0.0;
        *rho_v = rod*omg;
        *pv    = Pd;
    }
}

// 相変化ソース kernel。pure (N2/CPG) と carrier (H2O/TP) の両対応。一温度 T_v=T_d=T。
// J,r*,dr/dt は現在セル状態から freeze。安定化: J 上限、dr/dt<0→0、r̄≤r* 成長停止、g≤g_max、
//   1 step の Δg・潜熱 ΔT・蒸気枯渇を θ で律速 (全モーメント同 θ)。src_jac=潜熱自己抑制 (g)。
__global__ void condensation_source_d(
    geom_int nCells,
    int condModel, int carrier, double Rw, double M,
    int kantrowitz, int growthModel, double gyarC, int twoTemp,
    flow_float cp_cpg, flow_float gamma_cpg,
    double Jmax, double dg_max, double dT_max,
    geom_float* vol, flow_float* dt_local,
    flow_float* T, flow_float* P, flow_float* ro, flow_float* cp_cell, flow_float* Rmix_cell,
    flow_float* roY_w,   // carrier: 凝縮気相種の保存量 ρY_w (pure では nullptr)
    flow_float* rog, flow_float* roQ0, flow_float* roQ1, flow_float* roQ2,
    flow_float* res_rog, flow_float* res_roQ0, flow_float* res_roQ1, flow_float* res_roQ2,
    flow_float* sj_g, flow_float* sj_Q0, flow_float* sj_Q1, flow_float* sj_Q2)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    sj_Q0[ic] = 0.0; sj_Q1[ic] = 0.0; sj_Q2[ic] = 0.0; sj_g[ic] = 0.0;

    const double rod = (double)ro[ic];
    if (rod <= 1.0e-20) return;

    const CondSpeciesProps cprops = (condModel == 1) ? condProps_H2O() : condProps_N2();
    const double Td = (double)T[ic];
    const double Pd = (double)P[ic];

    // 気相熱容量 (ΔT/src_jac 用)。carrier(TP)=cp_cell、pure(CPG)=cp_cpg。
    const double cpg = (carrier && cp_cell) ? (double)cp_cell[ic] : (double)cp_cpg;
    const double Rg  = (carrier && Rmix_cell) ? (double)Rmix_cell[ic] : ((double)gamma_cpg-1.0)*(double)cp_cpg/(double)gamma_cpg;
    const double cvg = (cpg - Rg > 1.0e-3) ? (cpg - Rg) : 1.0e-3;

    double g = (double)rog[ic]/rod; if (g < 0.0) g = 0.0;
    const double Yw = (carrier && roY_w) ? (double)roY_w[ic]/rod : 1.0;
    const double gmax = carrier ? Yw : 0.99;
    if (g > gmax) g = (gmax > 0.0 ? gmax : 0.0);
    double q0 = (double)roQ0[ic]; if (q0 < 0.0) q0 = 0.0;
    double q1 = (double)roQ1[ic]; if (q1 < 0.0) q1 = 0.0;
    double q2 = (double)roQ2[ic]; if (q2 < 0.0) q2 = 0.0;

    // 蒸気状態
    double pv, rho_v;
    cond_vapor_state(carrier, rod, Pd, Td, g, Yw, Rw, &pv, &rho_v);

    // Kantrowitz 用の気相比熱比、Gyarmathy Kn 用の全圧 (carrier=Pd, pure=pv)
    const double gamma_gas = cpg/cvg;
    const double p_gas = carrier ? Pd : pv;

    // 核生成・成長 (freeze, 亜臨界停止・蒸発 clamp)
    double J, rstar;
    cond_nucleation(cprops, Td, pv, rho_v, &J, &rstar, kantrowitz, gamma_gas);
    if (J > Jmax) J = Jmax;
    if (J < 0.0)  J = 0.0;
    double r_bar = (q0 > 1.0e-30) ? (q1/q0) : rstar;
    double drdt = 0.0;
    if (q0 > 1.0e-30 && rstar > 0.0 && r_bar > rstar) {
        drdt = cond_growth(cprops, Td, pv, r_bar, rstar, growthModel, p_gas, gyarC, twoTemp);
        if (drdt < 0.0) drdt = 0.0;
    }
    const double rho_l = cond_rho_cond(cprops, Td);
    const double r_nuc = COND_RNUC_FAC*rstar;   // わずかに超臨界で生み成長を起動 (r*ちょうどだと(1-r*/r)=0で停止)
    double SQ0 = J;
    double SQ1 = J*r_nuc + q0*drdt;
    double SQ2 = J*r_nuc*r_nuc + 2.0*q1*drdt;
    double Sg  = (4.0/3.0)*COND_PI*rho_l*(J*r_nuc*r_nuc*r_nuc + 3.0*q2*drdt);
    if (Sg < 0.0) Sg = 0.0;

    // θ 律速: Δg, 潜熱 ΔT, 蒸気枯渇 (carrier は利用可能蒸気=Yw-g)
    const double dt = (double)dt_local[ic];
    const double L  = cond_latent(cprops, Td);
    double theta = 1.0;
    if (Sg > 0.0 && dt > 0.0) {
        const double dg  = Sg*dt/rod;
        const double dTl = dg*L/cvg;                       // 潜熱による ΔT (気相 cv で割る)
        const double avail = carrier ? (Yw - g) : (1.0 - g);
        if (dg  > dg_max)    theta = fmin(theta, dg_max/dg);
        if (dTl > dT_max)    theta = fmin(theta, dT_max/dTl);
        if (avail > 0.0 && dg > 0.9*avail) theta = fmin(theta, 0.9*avail/fmax(dg,1.0e-300));
        else if (avail <= 0.0) theta = 0.0;
    }

    // src_jac: 潜熱自己抑制 (g↑→T↑→psat↑→S↓→S_g↓)。∂T/∂(ρg)=(L-R_w T)/(ρ c_vg)。∂S_g/∂T 数値。
    {
        const double dTp = 0.1;
        double pvp, rvp;
        cond_vapor_state(carrier, rod, Pd, Td+dTp, g, Yw, Rw, &pvp, &rvp);
        double a0,a1,a2,ag;
        cond_source_vector(cprops, Td+dTp, pvp, rvp, q0, q1, q2, &a0,&a1,&a2,&ag,
                           kantrowitz, growthModel, gamma_gas, p_gas, gyarC, twoTemp);
        if (ag < 0.0) ag = 0.0;
        const double dSgdT  = (ag - Sg)/dTp;
        const double dTdrog = (L - (carrier?Rw:Rg)*Td)/(rod*cvg);
        double sjg = -theta*dSgdT*dTdrog;
        if (sjg < 0.0) sjg = 0.0;
        sj_g[ic] = (flow_float)sjg;
        if (q0 > 1.0e-30) {
            const double dq1 = (q1 > 0.0 ? 0.01*q1 : 1.0e-3);
            double b0,b1,b2,bg;
            cond_source_vector(cprops, Td, pv, rho_v, q0, q1+dq1, q2, &b0,&b1,&b2,&bg,
                               kantrowitz, growthModel, gamma_gas, p_gas, gyarC, twoTemp);
            double sjq1 = -theta*(b1 - SQ1)/dq1;
            if (sjq1 < 0.0) sjq1 = 0.0;
            sj_Q1[ic] = (flow_float)sjq1;
        }
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

    const int carrier = (cfg.condGasSpecies >= 0) ? 1 : 0;
    const CondSpeciesProps cprops = (cfg.condModel == 1) ? condProps_H2O() : condProps_N2();
    const double M  = cprops.M;
    const double Rw = cprops.R;
    const double Jmax   = 1.0e35;
    const double dg_max = 5.0e-3;
    const double dT_max = 1.0;

    for (int s = 0; s < var.nCondSpeciesRegistered; ++s) {
        const std::string i = std::to_string(s);
        flow_float* roY_w = nullptr;
        if (carrier) roY_w = var.c_d["roY" + std::to_string(cfg.condGasSpecies)];
        // CPG (thermalMethod!=2) では per-cell cp/Rmix 配列は未充填。nullptr で cfg.cp/γ フォールバック
        // (= carrier N2 の cp/R)。TP のみ per-cell 配列を渡す。
        flow_float* cp_cell   = (cfg.thermalMethod == 2) ? var.c_d["cp"]   : nullptr;
        flow_float* Rmix_cell = (cfg.thermalMethod == 2) ? var.c_d["Rmix"] : nullptr;
        condensation_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells,
            cfg.condModel, carrier, Rw, M,
            cfg.condKantrowitz, cfg.condGrowthModel, cfg.condGyarmathyC, cfg.condTwoTemp,
            cfg.cp, cfg.gamma,
            Jmax, dg_max, dT_max,
            var.c_d["volume"], var.c_d["dt_local"],
            var.c_d["T"], var.c_d["P"], var.c_d["ro"], cp_cell, Rmix_cell,
            roY_w,
            var.c_d["rog_"+i], var.c_d["roQ0_"+i], var.c_d["roQ1_"+i], var.c_d["roQ2_"+i],
            var.c_d["res_rog_"+i], var.c_d["res_roQ0_"+i], var.c_d["res_roQ1_"+i], var.c_d["res_roQ2_"+i],
            var.c_d["src_jac_g_"+i], var.c_d["src_jac_Q0_"+i], var.c_d["src_jac_Q1_"+i], var.c_d["src_jac_Q2_"+i]);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
