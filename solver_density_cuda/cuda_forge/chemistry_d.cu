// =============================================================================
// chemistry_d.cu — 有限速度化学ソース項 (Phase 1: 陽ソース + 対角 point-implicit Jacobian)
//   methods/chemistry.md 実装 §2。反応表は chemistry_mech_io.hpp で読み device へ 1 度アップロード。
// =============================================================================
#include <iostream>
#include <vector>
#include <cstring>
#include "chemistry_d.cuh"
#include "chemistry_mech_io.hpp"
#include "chemistrySource_d.cuh"
#include "speciesTransport_d.cuh"

static ReactionTable* g_rt_dev = nullptr;
static ReactionTable  g_rt_host;
static bool g_chem_ready = false;
// jacobianMode==2 用: 種ブロック Jacobian 残差行列 R [nCells × ns × ns] と反応熱のエネルギー対角 [nCells]
static flow_float* g_jac_dev = nullptr;      // R_sk = J_total_sk + d_s δ_sk (d_s = max(0,−J_ss) は src_jac_Y に入れる)
static flow_float* g_jacroe_dev = nullptr;   // max(0, −∂Q̇/∂(ρe)) [1/s]
static flow_float* g_cq_dev = nullptr;       // g_k = ∂Q̇/∂(ρY_k) = −Σ_s c_s J_total_sk [nCells × ns] (反応熱の陰的注入用)
static geom_int g_jac_cap = 0;
static int g_jac_ns = 0;
static bool g_block_active = false;

flow_float* chemistry_jac_device_ptr()    { return g_block_active ? g_jac_dev : nullptr; }
flow_float* chemistry_jacroe_device_ptr() { return g_block_active ? g_jacroe_dev : nullptr; }
flow_float* chemistry_cq_device_ptr()     { return g_block_active ? g_cq_dev : nullptr; }
bool chemistry_block_active()             { return g_block_active; }

bool chemistryEnabled(const solverConfig& cfg)
{
    return cfg.chemEnabled != 0 && cfg.thermalMethod == 2 && cfg.nSpecies >= 2;
}

bool chemistry_strang_active(const solverConfig& cfg)
{
    return chemistryEnabled(cfg) && cfg.chemStrang == 1 && cfg.unsteady == 1 && cfg.dualTime == 0
        && (cfg.timeIntegration == 1 || cfg.timeIntegration == 3);
}

// -----------------------------------------------------------------------------
// Strang 分離: セル内 ODE  d(ρY_s)/dt = ω_s,  d(ρe)/dt = Q̇  を dt_half だけ進める。
//   backward Euler を sub-cycle (h ≤ 0.5 τ_c, τ_c=1/max|J_ss|; 最大 CHEM_STRANG_MAXSUB 回) し、各 sub-step は
//   線形陰的 Euler (I − h J_tot) δ = h ω を 2 回 Newton 反復。エネルギーは離散恒等式 Δ(ρe) = −Σ_s c_s Δ(ρY_s) で更新
//   (絶対 datum では 0)。ρ・運動量は不変。T は e から Newton 反転 (thermo_T_from_e)。
// -----------------------------------------------------------------------------
#define CHEM_STRANG_MAXSUB 64
__global__ void chemistry_strang_d(
    geom_int nCells, int nSpecies, double dt_half, double Tmax, double Tfreeze,
    const SpeciesThermo* sp, const ReactionTable* rt,
    const flow_float* ro, const flow_float* roUx, const flow_float* roUy, const flow_float* roUz,
    flow_float* const* roY_dev, flow_float* roe, flow_float* chemQdot, flow_float* chemTau)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const int n = nSpecies;
    const double rho = (double)ro[ic];
    if (!(rho > 0.0)) return;
    const double ux = (double)roUx[ic]/rho, uy = (double)roUy[ic]/rho, uz = (double)roUz[ic]/rho;
    const double ek = 0.5*(ux*ux + uy*uy + uz*uz);
    double e = (double)roe[ic]/rho - ek;
    double U[THERMO_MAX_SPECIES], Y[THERMO_MAX_SPECIES], om[THERMO_MAX_SPECIES], dOdT[THERMO_MAX_SPECIES];
    double J[THERMO_MAX_SPECIES*THERMO_MAX_SPECIES], M[THERMO_MAX_SPECIES*THERMO_MAX_SPECIES], rhs[THERMO_MAX_SPECIES], dTdU[THERMO_MAX_SPECIES];
    double usum = 0.0;
    for (int s = 0; s < n; ++s) { double u = (double)roY_dev[s][ic]; if (u < 0.0) u = 0.0; U[s] = u; usum += u; }
    if (usum <= 0.0) return;
    for (int s = 0; s < n; ++s) { U[s] *= rho/usum; Y[s] = U[s]/rho; }
    double T = thermo_T_from_e(sp, n, Y, e, 1000.0, 200.0, 6000.0);
    if (T < Tfreeze) { chemQdot[ic] = 0.0; return; }
    const double U0e = e;
    double t = 0.0; int nsub = 0;
    while (t < dt_half && nsub < CHEM_STRANG_MAXSUB) {
        // Jacobian at current state
        double Tc = (T > Tmax) ? Tmax : ((T < 200.0) ? 200.0 : T);
        double Qd;
        chem_source(sp, rt, rho, Tc, Y, om, &Qd, 2, J, dOdT);
        const double cv = thermo_cp_mix(sp, n, Y, Tc) - thermo_R_mix(sp, n, Y);
        const double inv_rhocv = 1.0/(rho*((cv > 1.0e-3) ? cv : 1.0e-3));
        double jmax = 0.0;
        for (int k = 0; k < n; ++k) { const double e_k = thermo_h_mass(sp[k], Tc) - thermo_R_species(sp[k])*Tc; dTdU[k] = -e_k*inv_rhocv; }
        for (int s = 0; s < n; ++s) { const double jss = J[s*n+s] + dOdT[s]*dTdU[s]; if (fabs(jss) > jmax) jmax = fabs(jss); }
        double h = (jmax > 0.0) ? 0.5/jmax : dt_half;
        if (nsub == CHEM_STRANG_MAXSUB - 1 || h > dt_half - t) h = dt_half - t;
        if (h < 1.0e-30) break;
        // 線形陰的 Euler + 1 回の Newton 補正: (I − h J_tot) δ = h ω(U^k) − (U^k − U^n)
        double Un[THERMO_MAX_SPECIES]; for (int s = 0; s < n; ++s) Un[s] = U[s];
        for (int it = 0; it < 2; ++it) {
            if (it > 0) {
                double Tk = thermo_T_from_e(sp, n, Y, e, T, 200.0, 6000.0); Tc = (Tk > Tmax) ? Tmax : ((Tk < 200.0) ? 200.0 : Tk); T = Tk;
                chem_source(sp, rt, rho, Tc, Y, om, &Qd, 2, J, dOdT);
                const double cv2 = thermo_cp_mix(sp, n, Y, Tc) - thermo_R_mix(sp, n, Y);
                const double inv2 = 1.0/(rho*((cv2 > 1.0e-3) ? cv2 : 1.0e-3));
                for (int k = 0; k < n; ++k) { const double e_k = thermo_h_mass(sp[k], Tc) - thermo_R_species(sp[k])*Tc; dTdU[k] = -e_k*inv2; }
            }
            for (int s = 0; s < n; ++s) {
                for (int k = 0; k < n; ++k) M[s*n+k] = -h*(J[s*n+k] + dOdT[s]*dTdU[k]);
                M[s*n+s] += 1.0;
                rhs[s] = h*om[s] - (U[s] - Un[s]);
            }
            // LU (部分ピボット)
            for (int k = 0; k < n; ++k) {
                int p = k; double amax = fabs(M[k*n+k]);
                for (int i = k+1; i < n; ++i) if (fabs(M[i*n+k]) > amax) { amax = fabs(M[i*n+k]); p = i; }
                if (p != k) { for (int j = 0; j < n; ++j) { const double tt = M[k*n+j]; M[k*n+j] = M[p*n+j]; M[p*n+j] = tt; } const double tt = rhs[k]; rhs[k] = rhs[p]; rhs[p] = tt; }
                const double d = (fabs(M[k*n+k]) > 1.0e-300) ? M[k*n+k] : 1.0e-300;
                for (int i = k+1; i < n; ++i) { const double f = M[i*n+k]/d; if (f == 0.0) continue; for (int j = k; j < n; ++j) M[i*n+j] -= f*M[k*n+j]; rhs[i] -= f*rhs[k]; }
            }
            for (int i = n-1; i >= 0; --i) { double tt = rhs[i]; for (int j = i+1; j < n; ++j) tt -= M[i*n+j]*rhs[j]; rhs[i] = tt/M[i*n+i]; }
            // 更新 + エネルギー (Δe = −Σ c_s ΔU_s/ρ) + 実現可能性
            double de = 0.0, us = 0.0;
            for (int s = 0; s < n; ++s) { double Us = U[s] + rhs[s]; if (Us < 0.0) Us = 0.0; de -= (sp[s].h_datum/sp[s].MW)*(Us - U[s])/rho; U[s] = Us; us += Us; }
            e += de;
            for (int s = 0; s < n; ++s) { U[s] *= rho/us; Y[s] = U[s]/rho; }
        }
        T = thermo_T_from_e(sp, n, Y, e, T, 200.0, 6000.0);
        t += h; ++nsub;
    }
    for (int s = 0; s < n; ++s) roY_dev[s][ic] = (flow_float)U[s];
    roe[ic] = (flow_float)(rho*(e + ek));
    chemQdot[ic] = (flow_float)(rho*(e - U0e)/dt_half);   // 平均反応熱 [W/m3]
    chemTau[ic]  = (flow_float)((nsub > 0) ? dt_half/nsub : 1.0e30);   // 平均 sub-step (診断)
}

void chemistryStrangHalfStep_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, double dt_half)
{
    if (!g_chem_ready || !chemistry_strang_active(cfg)) return;
    flow_float** roY = species_roY_device_ptr();
    if (roY == nullptr) return;
    chemistry_strang_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells, var.nSpeciesRegistered, dt_half, cfg.chemTmaxReaction, cfg.chemFreezeBelowT,
        thermo_species_device_ptr(), g_rt_dev,
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"],
        roY, var.c_d["roe"], var.c_d["chemQdot"], var.c_d["chemTau"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void chemistry_init(solverConfig& cfg)
{
    g_chem_ready = false;
    if (cfg.chemEnabled == 0) return;
    if (cfg.thermalMethod != 2 || cfg.nSpecies < 2) {
        std::cerr << "[chemistry] chemistry.enabled=1 requires thermalMethod=2 and >=2 species" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::vector<std::string> eqs;
    try {
        chem_io::loadMechanism(cfg.chemMechanismFile, cfg.speciesNames, g_rt_host, &eqs);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl; std::exit(EXIT_FAILURE);
    }
    if (g_rt_dev) { cudaFree(g_rt_dev); g_rt_dev = nullptr; }
    gpuErrchk( cudaMalloc((void**)&g_rt_dev, sizeof(ReactionTable)) );
    gpuErrchk( cudaMemcpy(g_rt_dev, &g_rt_host, sizeof(ReactionTable), cudaMemcpyHostToDevice) );
    g_chem_ready = true;

    std::cout << "[chemistry] mechanism '" << cfg.chemMechanismFile << "': " << g_rt_host.nReac
              << " reactions, " << g_rt_host.nSpecies << " species (flow order)."
              << " Tmax=" << cfg.chemTmaxReaction << "K freezeBelowT=" << cfg.chemFreezeBelowT
              << "K jacobianMode=" << cfg.chemJacobianMode << std::endl;
    for (size_t r = 0; r < eqs.size(); ++r) {
        std::cout << "  (" << r+1 << ") " << eqs[r] << "  A=" << g_rt_host.A[r] << " b=" << g_rt_host.b[r]
                  << " Ea=" << g_rt_host.Ea[r] << " J/mol" << (g_rt_host.thirdBody[r] ? " [M]" : "")
                  << (g_rt_host.reversible[r] ? "" : " (irrev)") << "\n";
    }
    // datum の確認 (反応熱は h_datum から出る)
    const SpeciesThermo* sp = thermo_species_host();
    std::cout << "[chemistry] species datum c_s=h_abs(Tref) [J/kg]:";
    for (int s = 0; s < g_rt_host.nSpecies; ++s) std::cout << " " << cfg.speciesNames[s] << "=" << sp[s].h_datum/sp[s].MW;
    std::cout << std::endl;
    if (cfg.thermoHrefTemp <= 0.0) {
        std::cout << "[chemistry] WARNING: thermoHrefTemp=0 (absolute datum). Reaction heat enters through absolute enthalpies;"
                     " block-DPLUR is known to be unstable with absolute datum (methods/thermophysics.md)." << std::endl;
    }
}

const ReactionTable* chemistry_table_host() { return g_chem_ready ? &g_rt_host : nullptr; }

// -----------------------------------------------------------------------------
// セルごとに ω_s, Q̇, 対角 Jacobian を評価し残差へ加算する。
//   res_roY[s] += V ω_s,  src_jac_Y[s] += max(0, −∂ω_s/∂(ρY_s)),  res_roe += V Q̇
//   chemQdot = Q̇ [W/m3], chemTau = 1/max_s|∂ω_s/∂(ρY_s)| [s] (診断)
// -----------------------------------------------------------------------------
__global__ void chemistry_source_d(
    geom_int nCells, int nSpecies, int jacMode,
    double Tmax, double Tfreeze,
    const SpeciesThermo* sp, const ReactionTable* rt,
    const flow_float* vol,      // 周期 node では部分体積 (seam 二重計上防止, bodyForce_d と同じ)
    const flow_float* ro, const flow_float* T,
    flow_float* const* roY_dev,
    flow_float* const* res_roY_dev,
    flow_float* const* src_jac_dev,
    flow_float* res_roe,
    flow_float* chemQdot, flow_float* chemTau,
    flow_float* chem_jac,      // jacMode==2: [nCells*ns*ns] (nullptr 可)
    flow_float* chem_jacroe,   // jacMode==2: [nCells] (nullptr 可)
    flow_float* chem_cq)       // jacMode==2: [nCells*ns] ∂Q̇/∂(ρY_k) (nullptr 可)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    const double rho = (double)ro[ic];
    double Tc = (double)T[ic];
    chemQdot[ic] = 0.0; chemTau[ic] = 1.0e30;
    if (chem_jacroe) chem_jacroe[ic] = 0.0;
    if (chem_jac) for (int i = 0; i < nSpecies*nSpecies; ++i) chem_jac[(size_t)ic*nSpecies*nSpecies + i] = 0.0;
    if (chem_cq)  for (int k = 0; k < nSpecies; ++k) chem_cq[(size_t)ic*nSpecies + k] = 0.0;
    if (!(rho > 0.0) || !(Tc > 0.0) || Tc < Tfreeze) return;
    if (Tc > Tmax) Tc = Tmax;
    if (Tc < 200.0) Tc = 200.0;

    double Y[THERMO_MAX_SPECIES], omega[THERMO_MAX_SPECIES], dOdT[THERMO_MAX_SPECIES];
    double J[THERMO_MAX_SPECIES*THERMO_MAX_SPECIES];
    double ysum = 0.0;
    for (int s = 0; s < nSpecies; ++s) { double y = (double)roY_dev[s][ic] / rho; if (y < 0.0) y = 0.0; Y[s] = y; ysum += y; }
    if (ysum > 0.0) for (int s = 0; s < nSpecies; ++s) Y[s] /= ysum;

    double Qdot = 0.0;
    const int mode = (jacMode >= 2 && chem_jac != nullptr) ? 2 : (jacMode >= 1 ? 1 : 0);
    chem_source(sp, rt, rho, Tc, Y, omega, &Qdot, mode, J, (mode > 0) ? dOdT : nullptr);

    // 温度結合: ∂T/∂(ρY_k) = −e_k/(ρ c_v),  ∂T/∂(ρe) = 1/(ρ c_v)  (methods/chemistry.md §3)
    double dTdU[THERMO_MAX_SPECIES]; double inv_rhocv = 0.0;
    if (mode > 0) {
        const double cv = thermo_cp_mix(sp, nSpecies, Y, Tc) - thermo_R_mix(sp, nSpecies, Y);
        const double cvf = (cv > 1.0e-3) ? cv : 1.0e-3;
        inv_rhocv = 1.0 / (rho * cvf);
        for (int k = 0; k < nSpecies; ++k) {
            const double e_k = thermo_h_mass(sp[k], Tc) - thermo_R_species(sp[k]) * Tc;
            dTdU[k] = -e_k * inv_rhocv;
        }
    }

    const double V = (double)vol[ic];
    double jmax = 0.0;
    for (int s = 0; s < nSpecies; ++s) {
        res_roY_dev[s][ic] += (flow_float)(V * omega[s]);
        if (mode > 0) {
            // 対角 (温度結合込み): J_ss + ∂ω_s/∂T·∂T/∂(ρY_s)
            const double jss = J[s*nSpecies+s] + dOdT[s] * dTdU[s];
            const double d_s = (jss < 0.0) ? -jss : 0.0;
            src_jac_dev[s][ic] += (flow_float)d_s;
            if (fabs(jss) > jmax) jmax = fabs(jss);
            if (mode == 2) {
                // R_sk = J_total_sk + d_s δ_sk  (対角の消費部分は src_jac_Y に、残りをブロックへ)
                const double cs = sp[s].h_datum / sp[s].MW;
                for (int k = 0; k < nSpecies; ++k) {
                    const double jsk = J[s*nSpecies+k] + dOdT[s] * dTdU[k];
                    chem_jac[((size_t)ic*nSpecies + s)*nSpecies + k] = (flow_float)(jsk + ((s == k) ? d_s : 0.0));
                    if (chem_cq) chem_cq[(size_t)ic*nSpecies + k] += (flow_float)(-cs * jsk);   // ∂Q̇/∂(ρY_k)
                }
            }
        }
    }
    res_roe[ic] += (flow_float)(V * Qdot);
    chemQdot[ic] = (flow_float)Qdot;
    chemTau[ic]  = (flow_float)((jmax > 0.0) ? 1.0/jmax : 1.0e30);
    if (mode == 2 && chem_jacroe) {
        // ∂Q̇/∂(ρe) = −Σ_s c_s ∂ω_s/∂T /(ρ c_v)。安定化側 (負) のみ対角へ (Patankar)。
        double dQdT = 0.0;
        for (int s = 0; s < nSpecies; ++s) dQdT -= (sp[s].h_datum / sp[s].MW) * dOdT[s];
        const double dQde = dQdT * inv_rhocv;
        chem_jacroe[ic] = (flow_float)((dQde < 0.0) ? -dQde : 0.0);
    }
}

void chemistrySource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_chem_ready || !chemistryEnabled(cfg)) return;
    if (chemistry_strang_active(cfg)) return;   // Strang 分離: ソース項は RK に入れない (advanceExplicitRK の半ステップで進める)
    // jacobianMode==2: 種ブロック Jacobian バッファを遅延確保 (nCells × ns × ns)
    g_block_active = false;
    if (cfg.chemJacobianMode >= 2) {
        const int ns = var.nSpeciesRegistered;
        if (g_jac_cap < msh.nCells || g_jac_ns != ns) {
            if (g_jac_dev) cudaFree(g_jac_dev);
            if (g_jacroe_dev) cudaFree(g_jacroe_dev);
            gpuErrchk( cudaMalloc((void**)&g_jac_dev, (size_t)msh.nCells*ns*ns*sizeof(flow_float)) );
            gpuErrchk( cudaMalloc((void**)&g_jacroe_dev, (size_t)msh.nCells_all*sizeof(flow_float)) );
            if (g_cq_dev) cudaFree(g_cq_dev);
            gpuErrchk( cudaMalloc((void**)&g_cq_dev, (size_t)msh.nCells*ns*sizeof(flow_float)) );
            gpuErrchk( cudaMemset(g_jacroe_dev, 0, (size_t)msh.nCells_all*sizeof(flow_float)) );
            g_jac_cap = msh.nCells; g_jac_ns = ns;
        }
        g_block_active = true;
    }
    flow_float** roY = species_roY_device_ptr();
    flow_float** res = species_resroY_device_ptr();
    flow_float** sj  = species_srcjac_device_ptr();
    if (roY == nullptr || res == nullptr || sj == nullptr) return;

    chemistry_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells, var.nSpeciesRegistered, cfg.chemJacobianMode,
        cfg.chemTmaxReaction, cfg.chemFreezeBelowT,
        thermo_species_device_ptr(), g_rt_dev,
        // 周期 node の volume は合併体積なので体積ソースには部分体積を使う (periodicNodeGather で seam 2 倍/コーナー 4 倍になる)。
        (msh.volumePartial_d != nullptr) ? msh.volumePartial_d : var.c_d["volume"],
        var.c_d["ro"], var.c_d["T"],
        roY, res, sj,
        var.c_d["res_roe"], var.c_d["chemQdot"], var.c_d["chemTau"],
        g_block_active ? g_jac_dev : nullptr, g_block_active ? g_jacroe_dev : nullptr, g_block_active ? g_cq_dev : nullptr);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
