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
