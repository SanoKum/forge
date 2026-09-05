// =============================================================================
// 混合分率インフラ (Phase A of methods/chemistry_cmc.md): Bilger ξ̃ 診断、分散 ξ''² 輸送、χ̃ 診断。
// =============================================================================
#include "cmc_d.cuh"
#include "scalarTransport_d.cuh"
#include "chemistrySource_d.cuh"
#include "chemistry_d.cuh"
#include <vector>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <stdexcept>

namespace {

// Bilger 係数: β = Σ_s c_s Y_s、ξ = (β - β_O)/(β_F - β_O)
struct BilgerCoef {
    int nSpecies;
    double c[THERMO_MAX_SPECIES];
    double betaF, betaO;
};
__constant__ BilgerCoef g_bilger;
static BilgerCoef       g_bilger_host;
static flow_float**     g_Y_dev = nullptr;   // 原始質量分率 Y{s} の device ポインタ配列
static bool             g_ready = false;

__global__ void cmc_mixfrac_d(geom_int nCells, flow_float** Y, flow_float* xi)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    double beta = 0.0;
    for (int s = 0; s < g_bilger.nSpecies; ++s) beta += g_bilger.c[s] * (double)Y[s][ic];
    double z = (beta - g_bilger.betaO) / (g_bilger.betaF - g_bilger.betaO);
    z = fmin(fmax(z, 0.0), 1.0);
    xi[ic] = (flow_float)z;
}

// 1 スカラの Green-Gauss 勾配 (ransTransport の k/ω 版と同じ ghostless 規約)
__global__ void cmc_gradient_face_d(geom_int nPlanes, geom_int nCells, int isNode, geom_int* plane_cells,
                                    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz,
                                    flow_float* phi, flow_float* dpdx, flow_float* dpdy, flow_float* dpdz)
{
    const geom_int ip = blockDim.x * blockIdx.x + threadIdx.x;
    if (ip >= nPlanes) return;
    const geom_int ic0 = plane_cells[2 * ip + 0];
    const geom_int ic1 = plane_cells[2 * ip + 1];
    const geom_float f = fx[ip];
    flow_float pf;
    if (isNode != 0 && ic1 >= nCells) pf = phi[ic0];
    else                              pf = f * phi[ic0] + (1.0 - f) * phi[ic1];
    const geom_float sxx = sx[ip], syy = sy[ip], szz = sz[ip];
    atomicAdd(&dpdx[ic0],  sxx * pf); atomicAdd(&dpdy[ic0],  syy * pf); atomicAdd(&dpdz[ic0],  szz * pf);
    atomicAdd(&dpdx[ic1], -sxx * pf); atomicAdd(&dpdy[ic1], -syy * pf); atomicAdd(&dpdz[ic1], -szz * pf);
}
__global__ void cmc_gradient_div_vol_d(geom_int nCells, geom_float* vol, flow_float* dpdx, flow_float* dpdy, flow_float* dpdz)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const geom_float v = vol[ic];
    dpdx[ic] /= v; dpdy[ic] /= v; dpdz[ic] /= v;
}

// 分散のソース: 生成 2 μ_t/Sc_t |∇ξ̃|²、散逸 ρχ̃ = cChi β* ω (ρξ''²) (src_jac に対角)。χ̃ を診断出力。
__global__ void cmc_variance_source_d(geom_int nCells, geom_float* vol, flow_float* ro, flow_float* vis_turb, flow_float* omega,
                                      flow_float* dXidx, flow_float* dXidy, flow_float* dXidz,
                                      flow_float* roXiVar, flow_float* res_roXiVar, flow_float* src_jac,
                                      flow_float* chi, double Sc_t, double cChi, double betaStar)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const double g2   = (double)dXidx[ic]*dXidx[ic] + (double)dXidy[ic]*dXidy[ic] + (double)dXidz[ic]*dXidz[ic];
    const double prod = 2.0 * (double)vis_turb[ic] / Sc_t * g2;                 // [kg/m3/s]
    const double rate = cChi * betaStar * fmax((double)omega[ic], 0.0);          // [1/s]
    const double rxv  = fmax((double)roXiVar[ic], 0.0);
    const double diss = rate * rxv;                                              // [kg/m3/s]
    res_roXiVar[ic] += (flow_float)((prod - diss) * (double)vol[ic]);
    src_jac[ic]     += (flow_float)rate;
    chi[ic] = (flow_float)(rate * rxv / fmax((double)ro[ic], 1.0e-30));         // χ̃ [1/s]
}

__global__ void cmc_primitive_d(geom_int nCells_all, flow_float* ro, flow_float* xi, flow_float* roXiVar, flow_float* xiVar)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells_all) return;
    const double r  = fmax((double)ro[ic], 1.0e-30);
    const double z  = (double)xi[ic];
    const double vmax = r * fmax(z * (1.0 - z), 0.0);                            // 実現可能性: ξ''² ≤ ξ̃(1-ξ̃)
    double rv = fmin(fmax((double)roXiVar[ic], 0.0), vmax);
    roXiVar[ic] = (flow_float)rv;
    xiVar[ic]   = (flow_float)(rv / r);
}

__global__ void cmc_dirichlet_zero_d(geom_int nb, int isNode, geom_int* bplane_cell, geom_int* bplane_cell_ghst,
                                     flow_float* rophi, flow_float* phi)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib >= nb) return;
    const geom_int ig = bplane_cell_ghst[ib];
    rophi[ig] = 0.0; phi[ig] = 0.0;
    if (isNode != 0) { const geom_int ic = bplane_cell[ib]; rophi[ic] = 0.0; phi[ic] = 0.0; }   // 境界ノードへピン
}
__global__ void cmc_neumann_d(geom_int nb, geom_int* bplane_cell, geom_int* bplane_cell_ghst, flow_float* rophi, flow_float* phi)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib >= nb) return;
    const geom_int ic = bplane_cell[ib], ig = bplane_cell_ghst[ib];
    rophi[ig] = rophi[ic]; phi[ig] = phi[ic];
}

ScalarTransportDesc buildVarianceDesc(variables& var, const solverConfig& cfg)
{
    ScalarTransportDesc d{};
    d.phi = var.c_d.at("xiVar");
    d.dphidx = var.c_d.at("dXiVardx"); d.dphidy = var.c_d.at("dXiVardy"); d.dphidz = var.c_d.at("dXiVardz");
    d.rho_phi = var.c_d.at("roXiVar"); d.rho_phi_N = var.c_d.at("roXiVarN"); d.rho_phi_M = var.c_d.at("roXiVarM");
    d.res_rho_phi = var.c_d.at("res_roXiVar"); d.res_rho_phi_m = var.c_d.at("res_roXiVar_m");
    d.src_jac = var.c_d.at("src_jac_xiVar"); d.transport_diag = var.c_d.at("transport_diag_xiVar");
    d.floor = static_cast<flow_float>(0.0);
    d.sigma = static_cast<flow_float>(1.0 / cfg.Sc_t);   // 有効拡散 = μ + μ_t/Sc_t
    d.diffusion = 1;
    return d;
}

void gradientOf(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                const char* phi, const char* gx, const char* gy, const char* gz)
{
    flow_float* grad_vol = (cfg.isAxisymmetric == 1) ? var.c_d.at("A_planar") : var.c_d.at("volume");
    flow_float* grad_sx  = (cfg.isAxisymmetric == 1) ? var.p_d.at("sx_planar") : var.p_d.at("sx");
    flow_float* grad_sy  = (cfg.isAxisymmetric == 1) ? var.p_d.at("sy_planar") : var.p_d.at("sy");
    flow_float* grad_sz  = (cfg.isAxisymmetric == 1) ? var.p_d.at("sz_planar") : var.p_d.at("sz");
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at(gx), 0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at(gy), 0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at(gz), 0, msh.nCells_all * sizeof(flow_float)));
    cmc_gradient_face_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
        msh.nPlanes, msh.nCells, (cfg.discretization == "node") ? 1 : 0, msh.map_plane_cells_d,
        var.p_d.at("fx"), grad_sx, grad_sy, grad_sz, var.c_d.at(phi), var.c_d.at(gx), var.c_d.at(gy), var.c_d.at(gz));
    cmc_gradient_div_vol_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells, grad_vol, var.c_d.at(gx), var.c_d.at(gy), var.c_d.at(gz));
}

} // namespace

void cmcInit_d(solverConfig& cfg, variables& var)
{
    g_ready = false;
    if (!cmcMixfracEnabled(cfg, var)) return;
    const ReactionTable* rt = chemistry_table_host();
    if (rt == nullptr || rt->nElem == 0)
        throw std::runtime_error("[cmc] mixfrac needs element composition: mechanismFile must have species[].composition");
    if (rt->elemH < 0 && rt->elemC < 0)
        throw std::runtime_error("[cmc] Bilger mixture fraction needs H or C in the mechanism elements");

    // 種の分子量 (元素組成から) と Bilger 係数 c_s = Σ_e coef_e · atoms[s][e] W_e / W_s
    const int ns = rt->nSpecies;
    std::vector<double> W(ns, 0.0);
    for (int s = 0; s < ns; ++s) for (int e = 0; e < rt->nElem; ++e) W[s] += rt->atoms[s][e] * rt->elemW[e];
    auto coefOf = [&](int e) -> double {
        if (e == rt->elemC) return 2.0 / rt->elemW[e];
        if (e == rt->elemH) return 0.5 / rt->elemW[e];
        if (e == rt->elemO) return -1.0 / rt->elemW[e];
        return 0.0;
    };
    BilgerCoef& b = g_bilger_host;
    b.nSpecies = ns;
    for (int s = 0; s < ns; ++s) {
        double c = 0.0;
        if (W[s] > 0.0) for (int e = 0; e < rt->nElem; ++e) c += coefOf(e) * rt->atoms[s][e] * rt->elemW[e] / W[s];
        b.c[s] = c;
        if (W[s] <= 0.0) std::printf("[cmc] WARNING: species %s has no composition in the mechanism (treated as inert for xi)\n", cfg.speciesNames[s].c_str());
    }
    // 燃料/酸化剤流の β (モル分率 → 質量分率)
    auto betaOfStream = [&](const std::map<std::string, double>& X, const char* label) -> double {
        double mmix = 0.0;
        for (const auto& kv : X) {
            int s = -1; for (int i = 0; i < ns; ++i) if (cfg.speciesNames[i] == kv.first) s = i;
            if (s < 0) throw std::runtime_error(std::string("[cmc] mixfrac.") + label + ": species '" + kv.first + "' not in physProp.species");
            mmix += kv.second * W[s];
        }
        double beta = 0.0;
        for (const auto& kv : X) {
            int s = 0; for (int i = 0; i < ns; ++i) if (cfg.speciesNames[i] == kv.first) s = i;
            beta += b.c[s] * (kv.second * W[s] / mmix);
        }
        return beta;
    };
    b.betaF = betaOfStream(cfg.chemMixfracFuelX, "fuelX");
    b.betaO = betaOfStream(cfg.chemMixfracOxidX, "oxidizerX");
    if (std::fabs(b.betaF - b.betaO) < 1.0e-12) throw std::runtime_error("[cmc] fuelX and oxidizerX give the same Bilger beta");
    gpuErrchk( cudaMemcpyToSymbol(g_bilger, &b, sizeof(BilgerCoef)) );

    std::vector<flow_float*> hY(ns);
    for (int s = 0; s < ns; ++s) hY[s] = var.c_d.at("Y" + std::to_string(s));
    if (g_Y_dev) { cudaFree(g_Y_dev); g_Y_dev = nullptr; }
    gpuErrchk( cudaMalloc((void**)&g_Y_dev, ns * sizeof(flow_float*)) );
    gpuErrchk( cudaMemcpy(g_Y_dev, hY.data(), ns * sizeof(flow_float*), cudaMemcpyHostToDevice) );
    g_ready = true;
    std::printf("cmcInit_d: Bilger beta_F=%.6g beta_O=%.6g (nElem=%d, cChi=%.3g)\n", b.betaF, b.betaO, rt->nElem, cfg.chemCchi);
}

void cmcMixfrac_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    cmc_mixfrac_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, g_Y_dev, var.c_d.at("xi"));
    gradientOf(cfg, cuda_cfg, msh, var, "xi", "dXidx", "dXidy", "dXidz");
    gradientOf(cfg, cuda_cfg, msh, var, "xiVar", "dXiVardx", "dXiVardy", "dXiVardz");
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void cmcVarianceTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at("res_roXiVar"), 0, msh.nCells * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at("transport_diag_xiVar"), 0, msh.nCells * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d.at("src_jac_xiVar"), 0, msh.nCells * sizeof(flow_float)));
    const ScalarTransportDesc desc = buildVarianceDesc(var, cfg);
    scalarTransportResidual_d(cfg, cuda_cfg, msh, var, desc);
    const double betaStar = 0.09;
    cmc_variance_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        (msh.volumePartial_d != nullptr) ? msh.volumePartial_d : var.c_d.at("volume"),
        var.c_d.at("ro"), var.c_d.at("vis_turb"), var.c_d.at("omega"),
        var.c_d.at("dXidx"), var.c_d.at("dXidy"), var.c_d.at("dXidz"),
        var.c_d.at("roXiVar"), var.c_d.at("res_roXiVar"), var.c_d.at("src_jac_xiVar"), var.c_d.at("chi"),
        (double)cfg.Sc_t, cfg.chemCchi, betaStar);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void cmcVarianceTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    const ScalarTransportDesc desc = buildVarianceDesc(var, cfg);
    scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void cmcPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    // xi は最新の Y から (species primitive の後に呼ぶこと)
    cmc_mixfrac_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, g_Y_dev, var.c_d.at("xi"));
    cmc_primitive_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d.at("ro"), var.c_d.at("xi"), var.c_d.at("roXiVar"), var.c_d.at("xiVar"));
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void applyCmcBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_ready) return;
    const int isNode = (cfg.discretization == "node") ? 1 : 0;
    for (auto& bc : msh.bconds) {
        if (bc.iPlanes.empty()) continue;
        const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
        const bool isInlet = bc.bcondKind.rfind("inlet_", 0) == 0;
        if (isInlet) {
            cmc_dirichlet_zero_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(nb, isNode, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                                                                                 var.c_d.at("roXiVar"), var.c_d.at("xiVar"));
        } else {
            cmc_neumann_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(nb, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                                                                          var.c_d.at("roXiVar"), var.c_d.at("xiVar"));
        }
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void cmcUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!g_ready) return;
    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    gpuErrchk( cudaMemcpy(var.c_d.at("roXiVarN"), var.c_d.at("roXiVar"), bytes, cudaMemcpyDeviceToDevice) );
    gpuErrchk( cudaMemcpy(var.c_d.at("roXiVarM"), var.c_d.at("roXiVar"), bytes, cudaMemcpyDeviceToDevice) );
}

void cmcUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!g_ready) return;
    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    gpuErrchk( cudaMemcpy(var.c_d.at("roXiVarM"), var.c_d.at("roXiVar"), bytes, cudaMemcpyDeviceToDevice) );
}

// =============================================================================
// Phase B/C: 条件付きスカラー Q(η) — 混合線初期化・スライス輸送・η 拡散 (AMC) + 化学 (点陰解)・β-PDF 平均ソース
// =============================================================================
#include "thermo_d.cuh"
#include "speciesTransport_d.cuh"

namespace {

#define CMC_MAX_ETA 129

struct CmcParams {
    int nEta, nSpecies, chemOn, couple, kSt;
    double pdfFloor, Tmax, Tfreeze, dtScale, relax;
    double eta[CMC_MAX_ETA];       // η 格子
    double w[CMC_MAX_ETA];         // 台形重み (Σ=1)
    double G[CMC_MAX_ETA];         // AMC 形状 exp(-2 erfinv(2η-1)^2) (端は 0)
    double YF[THERMO_MAX_SPECIES], YO[THERMO_MAX_SPECIES];   // 燃料/酸化剤流の質量分率
    double hF, hO;                 // 燃料/酸化剤流の顕エンタルピー [J/kg]
};
__constant__ CmcParams g_cmc;
static CmcParams g_cmc_host;
static bool g_q_ready = false;
static int  g_nSlice = 0;
static geom_int g_nCellsAll = 0, g_nCells = 0;
static flow_float *g_Q = nullptr, *g_roQ = nullptr, *g_resQ = nullptr, *g_diagQ = nullptr, *g_sjQ = nullptr, *g_resmQ = nullptr;
static flow_float *g_gdum = nullptr;                       // 勾配ダミー (拡散カーネルは 1 次で勾配を使わない)
static flow_float *g_omega = nullptr, *g_qdot = nullptr, *g_jac = nullptr;   // PDF 平均ソース
static flow_float *g_ypdf = nullptr, *g_hpdf = nullptr, *g_tau = nullptr;   // couple 2: PDF 積分 Ỹ_s, h̃, 緩和時間 τ
static int g_cmcStep = 0;

__device__ __forceinline__ size_t sl(int var, int k) { return (size_t)(var * g_cmc.nEta + k); }

// 混合線: Y(η) = η Y_F + (1-η) Y_O, h(η) = η h_F + (1-η) h_O
__global__ void cmc_q_init_d(geom_int nCells_all, const flow_float* ro, flow_float* Q, flow_float* roQ)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells_all) return;
    const int ns = g_cmc.nSpecies, ne = g_cmc.nEta;
    const double r = (double)ro[ic];
    for (int k = 0; k < ne; ++k) {
        const double e = g_cmc.eta[k];
        for (int s = 0; s < ns; ++s) {
            const double y = e * g_cmc.YF[s] + (1.0 - e) * g_cmc.YO[s];
            Q[sl(s, k) * nCells_all + ic] = (flow_float)y; roQ[sl(s, k) * nCells_all + ic] = (flow_float)(r * y);
        }
        const double h = e * g_cmc.hF + (1.0 - e) * g_cmc.hO;
        Q[sl(ns, k) * nCells_all + ic] = (flow_float)h; roQ[sl(ns, k) * nCells_all + ic] = (flow_float)(r * h);
    }
}

// 燃料/酸化剤流の顕エンタルピーを device thermo で評価 (1 スレッド)
__global__ void cmc_stream_enthalpy_d(const SpeciesThermo* sp, int ns, double TF, double TO, double* out)
{
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    out[0] = thermo_h_mix(sp, ns, g_cmc.YF, TF);
    out[1] = thermo_h_mix(sp, ns, g_cmc.YO, TO);
}

// 境界: 入口は混合線 (ghost + node ピン)、他は ghost = owner。grid: (bplane × slice)
__global__ void cmc_q_boundary_d(geom_int nb, int nSlice, geom_int nCells_all, int isNode, int isInlet,
                                 const geom_int* bplane_cell, const geom_int* bplane_cell_ghst,
                                 const flow_float* ro, flow_float* Q, flow_float* roQ)
{
    const size_t t = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    if (t >= (size_t)nb * nSlice) return;
    const int s = (int)(t / nb); const geom_int ib = (geom_int)(t % nb);
    const geom_int ic = bplane_cell[ib], ig = bplane_cell_ghst[ib];
    const size_t off = (size_t)s * nCells_all;
    if (isInlet) {
        const int var = s / g_cmc.nEta, k = s % g_cmc.nEta; const double e = g_cmc.eta[k];
        const double q = (var < g_cmc.nSpecies) ? (e * g_cmc.YF[var] + (1.0 - e) * g_cmc.YO[var]) : (e * g_cmc.hF + (1.0 - e) * g_cmc.hO);
        Q[off + ig] = (flow_float)q; roQ[off + ig] = (flow_float)((double)ro[ig] * q);
        if (isNode) { Q[off + ic] = (flow_float)q; roQ[off + ic] = (flow_float)((double)ro[ic] * q); }
    } else {
        Q[off + ig] = Q[off + ic]; roQ[off + ig] = roQ[off + ic];
    }
}

// 保存量 roQ = ρ̄ Q を現在の ρ̄ で作り直す (輸送 point-implicit の直前)。前ステップの ρ̄ で作った roQ を更新後の ρ̄ で割ると
// 平均密度の変化分だけ Q が偽に動く (凍結検証で条件付き T が 1045 K を超えた指紋) ので、更新の前後で同じ ρ̄ を使う。
__global__ void cmc_q_resync_d(geom_int nCells_all, int nSlice, const flow_float* ro, const flow_float* Q, flow_float* roQ)
{
    const size_t t = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    if (t >= (size_t)nCells_all * nSlice) return;
    const int s = (int)(t / nCells_all); const geom_int ic = (geom_int)(t % nCells_all);
    const size_t i = (size_t)s * nCells_all + ic;
    roQ[i] = ro[ic] * Q[i];
}

// 非保存補正: 収束していない平均流では Σ_f ṁ_f ≠ 0 なので、保存形の移流残差は一様な Q でも −Q Σṁ_f = Q·res_ro が残る。
// 流れ側の密度更新 (block-DPLUR) と点陰解のスカラー更新は同じ Δρ を与えないため、これが毎ステップ Q を (ρ+Δρ_s)/(ρ+Δρ_f) 倍
// 偽に動かす (凍結検証で条件付き T が ±40 K ずれた指紋)。res_Q -= Q·res_ro で「Q の輸送方程式」(非保存形) に揃える。
__global__ void cmc_q_nonconservative_fix_d(geom_int nCells, int nSlice, geom_int nCells_all, const flow_float* res_ro, const flow_float* Q, flow_float* resQ)
{
    const size_t t = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    if (t >= (size_t)nCells * nSlice) return;
    const int s = (int)(t / nCells); const geom_int ic = (geom_int)(t % nCells);
    const size_t i = (size_t)s * nCells_all + ic;
    resQ[i] -= Q[i] * res_ro[ic];
}

// 原始量 Q = roQ/ρ (輸送後)
__global__ void cmc_q_primitive_d(geom_int nCells, int nSlice, geom_int nCells_all, const flow_float* ro, const flow_float* roQ, flow_float* Q)
{
    const size_t t = (size_t)blockDim.x * blockIdx.x + threadIdx.x;
    if (t >= (size_t)nCells * nSlice) return;
    const int s = (int)(t / nCells); const geom_int ic = (geom_int)(t % nCells);
    const size_t i = (size_t)s * nCells_all + ic;
    Q[i] = roQ[i] / fmaxf(ro[ic], 1.0e-30f);
}

// β-PDF の離散重み Ω_k (Σ=1)。分散が微小ならデルタ (隣接 2 点へ線形配分)。
__device__ void cmc_pdf_weights(double xi, double var, double* Om)
{
    const int ne = g_cmc.nEta;
    xi = fmin(fmax(xi, 0.0), 1.0);
    const double vmax = xi * (1.0 - xi);
    for (int k = 0; k < ne; ++k) Om[k] = 0.0;
    if (var < 1.0e-4 * vmax + 1.0e-10 || vmax < 1.0e-12) {
        int k = 0; while (k < ne - 2 && g_cmc.eta[k + 1] <= xi) ++k;
        const double d = g_cmc.eta[k + 1] - g_cmc.eta[k];
        const double f = fmin(fmax((xi - g_cmc.eta[k]) / d, 0.0), 1.0);
        Om[k] = 1.0 - f; Om[k + 1] = f;
        return;
    }
    double gam = vmax / fmin(var, 0.999 * vmax) - 1.0;
    gam = fmax(gam, 0.05);
    const double a = xi * gam, b = (1.0 - xi) * gam;
    double lnOm[CMC_MAX_ETA]; double lmax = -1.0e300;
    for (int k = 1; k < ne - 1; ++k) {
        const double e = g_cmc.eta[k];
        lnOm[k] = log(g_cmc.w[k]) + (a - 1.0) * log(e) + (b - 1.0) * log(1.0 - e);
        if (lnOm[k] > lmax) lmax = lnOm[k];
    }
    // 端ビン: ∫_0^{η1/2} η^{a-1} dη = (η1/2)^a / a (特異性の解析処理), 同様に η=1 側
    lnOm[0]      = a * log(0.5 * g_cmc.eta[1]) - log(a);
    lnOm[ne - 1] = b * log(0.5 * (1.0 - g_cmc.eta[ne - 2])) - log(b);
    if (lnOm[0] > lmax) lmax = lnOm[0];
    if (lnOm[ne - 1] > lmax) lmax = lnOm[ne - 1];
    double sum = 0.0;
    for (int k = 0; k < ne; ++k) { Om[k] = exp(lnOm[k] - lmax); sum += Om[k]; }
    for (int k = 0; k < ne; ++k) Om[k] /= sum;
}

// 小行列の Gauss 消去 (部分ピボット): A x = b を解く (A は n×n 行優先, 破壊)
__device__ bool cmc_solve_dense(int n, double* A, double* b)
{
    for (int c = 0; c < n; ++c) {
        int p = c; double pm = fabs(A[c*n + c]);
        for (int r = c + 1; r < n; ++r) if (fabs(A[r*n + c]) > pm) { pm = fabs(A[r*n + c]); p = r; }
        if (pm < 1.0e-300) return false;
        if (p != c) { for (int j = 0; j < n; ++j) { double t = A[c*n+j]; A[c*n+j] = A[p*n+j]; A[p*n+j] = t; } double t = b[c]; b[c] = b[p]; b[p] = t; }
        for (int r = c + 1; r < n; ++r) {
            const double f = A[r*n + c] / A[c*n + c];
            if (f == 0.0) continue;
            for (int j = c; j < n; ++j) A[r*n + j] -= f * A[c*n + j];
            b[r] -= f * b[c];
        }
    }
    for (int c = n - 1; c >= 0; --c) {
        double s = b[c];
        for (int j = c + 1; j < n; ++j) s -= A[c*n + j] * b[j];
        b[c] = s / A[c*n + c];
    }
    return true;
}

// node 毎: η 拡散 (陰的, AMC) → 各 η 点で化学 (点陰解) → PDF 平均 (ω̄, Q̇̄, J̄, Ỹ_pdf)
__global__ void cmc_step_d(geom_int nCells, geom_int nCells_all, const SpeciesThermo* sp, const ReactionTable* rt,
                           const flow_float* ro, const flow_float* P, const flow_float* xi, const flow_float* xiVar, const flow_float* chi,
                           const flow_float* dt_local, flow_float* const* roY_dev,
                           flow_float* Q, flow_float* roQ,
                           flow_float* omegaBar, flow_float* qdotBar, flow_float* jacBar, flow_float* cmc_dY, flow_float* cmc_TQmax,
                           flow_float* ypdfOut, flow_float* hpdfOut, flow_float* tauOut, flow_float* cmc_TQst)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const int ns = g_cmc.nSpecies, ne = g_cmc.nEta;
    const double rho = (double)ro[ic], p = (double)P[ic];
    const double dtau = (double)dt_local[ic] * g_cmc.dtScale;

    double Om[CMC_MAX_ETA], chiEta[CMC_MAX_ETA];
    cmc_pdf_weights((double)xi[ic], (double)xiVar[ic], Om);
    // AMC: <χ|η> = χ̃ G(η) / Σ Ω_k G_k
    double gnorm = 0.0;
    for (int k = 0; k < ne; ++k) gnorm += Om[k] * g_cmc.G[k];
    const double chim = fmax((double)chi[ic], 0.0);
    for (int k = 0; k < ne; ++k) chiEta[k] = (gnorm > 1.0e-30) ? chim * g_cmc.G[k] / gnorm : chim * g_cmc.G[k];

    // ---- η 拡散 (陰的三重対角, 端 Dirichlet) を各変数に ----
    double a[CMC_MAX_ETA], bb[CMC_MAX_ETA], cc[CMC_MAX_ETA], d[CMC_MAX_ETA];
    for (int var = 0; var <= ns; ++var) {
        for (int k = 1; k < ne - 1; ++k) {
            const double hm = g_cmc.eta[k] - g_cmc.eta[k-1], hp = g_cmc.eta[k+1] - g_cmc.eta[k], hc = 0.5 * (hm + hp);
            const double D = 0.5 * chiEta[k] * dtau;
            a[k] = -D / (hm * hc); cc[k] = -D / (hp * hc); bb[k] = 1.0 - a[k] - cc[k];
            d[k] = (double)Q[sl(var, k) * nCells_all + ic];
        }
        const double q0 = (double)Q[sl(var, 0) * nCells_all + ic], qN = (double)Q[sl(var, ne-1) * nCells_all + ic];
        d[1] -= a[1] * q0; d[ne-2] -= cc[ne-2] * qN;
        // Thomas
        for (int k = 2; k < ne - 1; ++k) { const double m = a[k] / bb[k-1]; bb[k] -= m * cc[k-1]; d[k] -= m * d[k-1]; }
        d[ne-2] /= bb[ne-2];
        for (int k = ne - 3; k >= 1; --k) d[k] = (d[k] - cc[k] * d[k+1]) / bb[k];
        for (int k = 1; k < ne - 1; ++k) Q[sl(var, k) * nCells_all + ic] = (flow_float)d[k];
    }

    // ---- 化学 (各 η 点, 点陰解) + PDF 平均 ----
    double obar[THERMO_MAX_SPECIES], ypdf[THERMO_MAX_SPECIES];
    double jbar[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES];
    for (int s = 0; s < ns; ++s) { obar[s] = 0.0; ypdf[s] = 0.0; }
    for (int i = 0; i < ns * ns; ++i) jbar[i] = 0.0;
    double qbar = 0.0, TQmax = 0.0, hpdf = 0.0, TQst = 0.0;
    double Y[THERMO_MAX_SPECIES], omega[THERMO_MAX_SPECIES], dOdT[THERMO_MAX_SPECIES], J[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES];
    double A[THERMO_MAX_SPECIES * THERMO_MAX_SPECIES], rhs[THERMO_MAX_SPECIES];
    for (int k = 0; k < ne; ++k) {
        double ysum = 0.0;
        for (int s = 0; s < ns; ++s) { double y = (double)Q[sl(s, k) * nCells_all + ic]; if (y < 0.0) y = 0.0; Y[s] = y; ysum += y; }
        if (ysum > 0.0) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
        double h = (double)Q[sl(ns, k) * nCells_all + ic];
        double T = thermo_T_from_h(sp, ns, Y, h, fmin(fmax(300.0 + h / 1200.0, 200.0), 4000.0), 200.0, g_cmc.Tmax);
        if (T > TQmax) TQmax = T;
        if (k == g_cmc.kSt) TQst = T;
        for (int s = 0; s < ns; ++s) ypdf[s] += Om[k] * Y[s];
        hpdf += Om[k] * h;
        const bool interior = (k > 0 && k < ne - 1);
        if (!g_cmc.chemOn || !interior || Om[k] < g_cmc.pdfFloor || T < g_cmc.Tfreeze) continue;
        const double rhoEta = p / (thermo_R_mix(sp, ns, Y) * T);
        double Qd = 0.0;
        chem_source(sp, rt, rhoEta, fmin(T, g_cmc.Tmax), Y, omega, &Qd, 2, J, dOdT);
        // PDF 平均 (更新前の状態で評価; 1 step の遅れは擬似時間収束で消える)
        for (int s = 0; s < ns; ++s) obar[s] += Om[k] * omega[s];
        for (int i = 0; i < ns * ns; ++i) jbar[i] += Om[k] * J[i];
        qbar += Om[k] * Qd;
        // 点陰解: (I - Δτ J) δY = Δτ ω/ρ_η
        for (int s = 0; s < ns; ++s) { for (int q = 0; q < ns; ++q) A[s*ns + q] = -dtau * J[s*ns + q]; A[s*ns + s] += 1.0; rhs[s] = dtau * omega[s] / rhoEta; }
        if (cmc_solve_dense(ns, A, rhs)) {
            ysum = 0.0;
            for (int s = 0; s < ns; ++s) { Y[s] = fmax(Y[s] + rhs[s], 0.0); ysum += Y[s]; }
            if (ysum > 0.0) for (int s = 0; s < ns; ++s) Y[s] /= ysum;
            h += dtau * Qd / rhoEta;
            for (int s = 0; s < ns; ++s) Q[sl(s, k) * nCells_all + ic] = (flow_float)Y[s];
            Q[sl(ns, k) * nCells_all + ic] = (flow_float)h;
        }
    }
    // 保存量 roQ = ρ̄ Q (輸送の重み), 診断, PDF 平均ソース
    for (int s = 0; s <= ns; ++s) for (int k = 0; k < ne; ++k) { const size_t i = sl(s, k) * nCells_all + ic; roQ[i] = (flow_float)(rho * (double)Q[i]); }
    double dy = 0.0;
    for (int s = 0; s < ns; ++s) { const double ym = (double)roY_dev[s][ic] / fmax(rho, 1.0e-30); dy = fmax(dy, fabs(ypdf[s] - ym)); }
    cmc_dY[ic] = (flow_float)dy; cmc_TQmax[ic] = (flow_float)TQmax; cmc_TQst[ic] = (flow_float)TQst;
    for (int s = 0; s < ns; ++s) ypdfOut[(size_t)s * nCells + ic] = (flow_float)ypdf[s];
    hpdfOut[ic] = (flow_float)hpdf;
    tauOut[ic]  = (flow_float)fmax((double)dt_local[ic] * g_cmc.relax, 1.0e-12);
    for (int s = 0; s < ns; ++s) omegaBar[(size_t)s * nCells + ic] = (flow_float)obar[s];
    qdotBar[ic] = (flow_float)qbar;
    for (int i = 0; i < ns * ns; ++i) jacBar[(size_t)ic * ns * ns + i] = (flow_float)jbar[i];
}

ScalarTransportDesc sliceDesc(int s, const solverConfig& cfg)
{
    const size_t off = (size_t)s * g_nCellsAll;
    ScalarTransportDesc d{};
    d.phi = g_Q + off; d.dphidx = g_gdum; d.dphidy = g_gdum; d.dphidz = g_gdum;
    d.rho_phi = g_roQ + off; d.rho_phi_N = g_roQ + off; d.rho_phi_M = g_roQ + off;
    d.res_rho_phi = g_resQ + off; d.res_rho_phi_m = g_resmQ + off;
    d.src_jac = g_sjQ + off; d.transport_diag = g_diagQ + off;
    d.floor = static_cast<flow_float>(-1.0e30);   // Q は負のエンタルピーもあり得る (フロアなし)
    d.sigma = static_cast<flow_float>(1.0 / cfg.Sc_t);
    d.diffusion = 1;
    return d;
}

} // namespace

bool cmc_coupling_active() { return g_q_ready && g_cmc_host.couple != 0; }
int  cmc_coupling_mode()   { return g_q_ready ? g_cmc_host.couple : 0; }
const flow_float* cmc_ypdf_device_ptr() { return g_ypdf; }
const flow_float* cmc_hpdf_device_ptr() { return g_hpdf; }
const flow_float* cmc_tau_device_ptr()  { return g_tau; }
const flow_float* cmc_omega_device_ptr() { return g_omega; }
const flow_float* cmc_qdot_device_ptr()  { return g_qdot; }
const flow_float* cmc_jac_device_ptr()   { return g_jac; }

void cmcQInit_d(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    g_q_ready = false;
    if (cfg.chemCmc == 0 || !g_ready) return;
    const int ns = var.nSpeciesRegistered, ne = cfg.chemCmcNEta;
    if (ne > CMC_MAX_ETA) throw std::runtime_error("[cmc] nEta > CMC_MAX_ETA");
    CmcParams& c = g_cmc_host;
    std::memset(&c, 0, sizeof(CmcParams));
    c.nEta = ne; c.nSpecies = ns; c.chemOn = cfg.chemCmcChem; c.couple = cfg.chemCmcCouple;
    c.pdfFloor = cfg.chemCmcPdfFloor; c.Tmax = cfg.chemTmaxReaction; c.Tfreeze = cfg.chemFreezeBelowT; c.dtScale = cfg.chemCmcDtScale;
    c.relax = cfg.chemCmcRelax;
    for (int k = 0; k < ne; ++k) { const double s = (double)k / (double)(ne - 1); c.eta[k] = std::pow(s, cfg.chemCmcEtaPow); }
    c.eta[0] = 0.0; c.eta[ne - 1] = 1.0;
    { int kb = 0; double db = 1.0; for (int k = 0; k < ne; ++k) { const double dd = std::fabs(c.eta[k] - cfg.chemCmcXiSt); if (dd < db) { db = dd; kb = k; } } c.kSt = kb; }
    for (int k = 0; k < ne; ++k) {
        const double lo = (k == 0) ? 0.0 : 0.5 * (c.eta[k] + c.eta[k-1]);
        const double hi = (k == ne - 1) ? 1.0 : 0.5 * (c.eta[k] + c.eta[k+1]);
        c.w[k] = hi - lo;
        c.G[k] = 0.0;   // 端は 0、内点は下で AMC 形状
    }
    // AMC 形状 G(η) = exp(-2 [erfinv(2η-1)]^2): erfinv は host に無いので Newton で求める
    for (int k = 1; k < ne - 1; ++k) {
        const double y = 2.0 * c.eta[k] - 1.0; double x = 0.0;
        for (int it = 0; it < 60; ++it) { const double f = std::erf(x) - y; const double df = 2.0 / std::sqrt(M_PI) * std::exp(-x * x); x -= f / df; }
        c.G[k] = std::exp(-2.0 * x * x);
    }
    // 燃料/酸化剤流の質量分率 (fuelX/oxidizerX → Y)
    const ReactionTable* rt = chemistry_table_host();
    std::vector<double> W(ns, 0.0);
    for (int s = 0; s < ns; ++s) for (int e = 0; e < rt->nElem; ++e) W[s] += rt->atoms[s][e] * rt->elemW[e];
    auto toY = [&](const std::map<std::string, double>& X, double* Yout) {
        double mm = 0.0;
        for (const auto& kv : X) { int s = -1; for (int i = 0; i < ns; ++i) if (cfg.speciesNames[i] == kv.first) s = i; mm += kv.second * W[s]; }
        for (int s = 0; s < ns; ++s) Yout[s] = 0.0;
        for (const auto& kv : X) { int s = 0; for (int i = 0; i < ns; ++i) if (cfg.speciesNames[i] == kv.first) s = i; Yout[s] = kv.second * W[s] / mm; }
    };
    toY(cfg.chemMixfracFuelX, c.YF); toY(cfg.chemMixfracOxidX, c.YO);
    gpuErrchk( cudaMemcpyToSymbol(g_cmc, &c, sizeof(CmcParams)) );
    // 流れのエンタルピー (device thermo)
    double* hdev = nullptr; double hh[2];
    gpuErrchk( cudaMalloc((void**)&hdev, 2 * sizeof(double)) );
    cmc_stream_enthalpy_d<<<1, 32>>>(thermo_species_device_ptr(), ns, cfg.chemCmcTfuel, cfg.chemCmcTox, hdev);
    gpuErrchk( cudaMemcpy(hh, hdev, 2 * sizeof(double), cudaMemcpyDeviceToHost) ); cudaFree(hdev);
    c.hF = hh[0]; c.hO = hh[1];
    gpuErrchk( cudaMemcpyToSymbol(g_cmc, &c, sizeof(CmcParams)) );

    g_nSlice = (ns + 1) * ne; g_nCellsAll = msh.nCells_all; g_nCells = msh.nCells;
    const size_t bytes = (size_t)g_nSlice * g_nCellsAll * sizeof(flow_float);
    for (flow_float** pp : {&g_Q, &g_roQ, &g_resQ, &g_diagQ, &g_sjQ, &g_resmQ}) { if (*pp) cudaFree(*pp); gpuErrchk( cudaMalloc((void**)pp, bytes) ); gpuErrchk( cudaMemset(*pp, 0, bytes) ); }
    if (g_gdum) cudaFree(g_gdum); gpuErrchk( cudaMalloc((void**)&g_gdum, (size_t)g_nCellsAll * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_gdum, 0, (size_t)g_nCellsAll * sizeof(flow_float)) );
    if (g_omega) cudaFree(g_omega); gpuErrchk( cudaMalloc((void**)&g_omega, (size_t)ns * msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_omega, 0, (size_t)ns * msh.nCells * sizeof(flow_float)) );
    if (g_qdot) cudaFree(g_qdot); gpuErrchk( cudaMalloc((void**)&g_qdot, (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_qdot, 0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_ypdf) cudaFree(g_ypdf); gpuErrchk( cudaMalloc((void**)&g_ypdf, (size_t)ns * msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_ypdf, 0, (size_t)ns * msh.nCells * sizeof(flow_float)) );
    if (g_hpdf) cudaFree(g_hpdf); gpuErrchk( cudaMalloc((void**)&g_hpdf, (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_hpdf, 0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_tau)  cudaFree(g_tau);  gpuErrchk( cudaMalloc((void**)&g_tau,  (size_t)msh.nCells * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_tau,  0, (size_t)msh.nCells * sizeof(flow_float)) );
    if (g_jac) cudaFree(g_jac); gpuErrchk( cudaMalloc((void**)&g_jac, (size_t)msh.nCells * ns * ns * sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_jac, 0, (size_t)msh.nCells * ns * ns * sizeof(flow_float)) );

    cmc_q_init_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d.at("ro"), g_Q, g_roQ);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    g_q_ready = true; g_cmcStep = 0;
    std::printf("cmcQInit_d: nEta=%d nSlice=%d (%.1f MB x6), hF=%.4g hO=%.4g J/kg, couple=%d chem=%d\n",
                ne, g_nSlice, bytes / 1.0e6, c.hF, c.hO, c.couple, c.chemOn);
}

void cmcQTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_q_ready) return;
    const size_t bytes = (size_t)g_nSlice * g_nCellsAll * sizeof(flow_float);
    CHECK_CUDA_ERROR(cudaMemset(g_resQ, 0, bytes)); CHECK_CUDA_ERROR(cudaMemset(g_diagQ, 0, bytes)); CHECK_CUDA_ERROR(cudaMemset(g_sjQ, 0, bytes));
    for (int s = 0; s < g_nSlice; ++s) scalarTransportResidual_d(cfg, cuda_cfg, msh, var, sliceDesc(s, cfg));
    { const size_t n = (size_t)g_nCells * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
      cmc_q_nonconservative_fix_d<<<grid, blk>>>(g_nCells, g_nSlice, g_nCellsAll, var.c_d.at("res_ro"), g_Q, g_resQ); }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}

void applyCmcQBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_q_ready) return;
    const int isNode = (cfg.discretization == "node") ? 1 : 0;
    for (auto& bc : msh.bconds) {
        if (bc.iPlanes.empty()) continue;
        const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
        const bool isInlet = bc.bcondKind.rfind("inlet_", 0) == 0;
        const size_t n = (size_t)nb * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
        cmc_q_boundary_d<<<grid, blk>>>(nb, g_nSlice, g_nCellsAll, isNode, isInlet ? 1 : 0, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                                         var.c_d.at("ro"), g_Q, g_roQ);
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}

void cmcQUpdate_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_q_ready) return;
    // (0) roQ を現在の ρ̄ で作り直す (ρ̄ は直前の流れ更新で変わっている)
    { const size_t n = (size_t)g_nCellsAll * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
      cmc_q_resync_d<<<grid, blk>>>(g_nCellsAll, g_nSlice, var.c_d.at("ro"), g_Q, g_roQ); }
    // (1) 物理空間の輸送 point-implicit (各スライス; N=M=現在値のエイリアス)
    for (int s = 0; s < g_nSlice; ++s) scalarTimeIntegration_d(0, cfg, cuda_cfg, msh, var, sliceDesc(s, cfg));
    { const size_t n = (size_t)g_nCells * g_nSlice; const int blk = 256; const unsigned grid = (unsigned)((n + blk - 1) / blk);
      cmc_q_primitive_d<<<grid, blk>>>(g_nCells, g_nSlice, g_nCellsAll, var.c_d.at("ro"), g_roQ, g_Q); }
    // (2) η 拡散 + 化学 + PDF 平均 (interval 毎)
    ++g_cmcStep;
    if (cfg.chemCmcInterval <= 1 || (g_cmcStep % cfg.chemCmcInterval) == 0) {
        cmc_step_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            g_nCells, g_nCellsAll, thermo_species_device_ptr(), chemistry_table_device(),
            var.c_d.at("ro"), var.c_d.at("P"), var.c_d.at("xi"), var.c_d.at("xiVar"), var.c_d.at("chi"),
            var.c_d.at("dt_local"), species_roY_device_ptr(),
            g_Q, g_roQ, g_omega, g_qdot, g_jac, var.c_d.at("cmc_dY"), var.c_d.at("cmc_TQmax"),
            g_ypdf, g_hpdf, g_tau, var.c_d.at("cmc_TQst"));
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}
