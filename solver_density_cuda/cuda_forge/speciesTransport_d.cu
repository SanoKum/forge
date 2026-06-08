#include "speciesTransport_d.cuh"

#include "scalarTransport_d.cuh"

#include <string>
#include <vector>

namespace {

// device の roY ポインタ配列 (flow_float*[nSpecies]) と化学種数。speciesInit_d で 1 度だけ構築。
flow_float** g_roY_dev = nullptr;
int          g_nSpecies = 0;

constexpr flow_float kSmall = static_cast<flow_float>(1.0e-30);

inline bool speciesEnabled(const variables& var)
{
    return var.nSpeciesRegistered >= 2;
}

// 化学種 s 用のスカラ輸送記述子を構築する (RANS buildScalarDescs と同形)。
// floor=0 (Y_s>=0)、sigma=0 (M2 は移流のみ。拡散は M4)。
ScalarTransportDesc buildSpeciesDesc(variables& var, int s)
{
    const std::string i = std::to_string(s);
    return ScalarTransportDesc{
        var.c_d["Y"+i], var.c_d["roY"+i], var.c_d["roY"+i+"N"], var.c_d["roY"+i+"M"],
        var.c_d["res_roY"+i], var.c_d["res_roY"+i+"_m"],
        var.c_d["src_jac_Y"+i], var.c_d["transport_diag_Y"+i],
        static_cast<flow_float>(0.0), static_cast<flow_float>(0.0)
    };
}

__global__ void species_primitive_d(
    geom_int nCells_all,
    flow_float* ro,
    flow_float* roY,
    flow_float* Y)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells_all) {
        Y[ic] = roY[ic] / max(ro[ic], kSmall);
    }
}

// Neumann (zero-gradient) ghost 充填: roY[ig]=roY[ic], Y[ig]=Y[ic]。
__global__ void species_neumann_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* roY,
    flow_float* Y)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];
        roY[ig] = roY[ic];
        Y[ig]   = Y[ic];
    }
}

// 実現可能性 + 再正規化: 各 ρY_s>=0 にクランプ後、Σ_s ρY_s = ρ となるよう再スケール (ΣY_s=1)。
__global__ void species_renormalize_d(
    geom_int nCells,
    int nSpecies,
    flow_float** roY,
    flow_float* ro)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        double sum = 0.0;
        for (int s = 0; s < nSpecies; s++) {
            flow_float v = roY[s][ic];
            if (v < 0.0) v = 0.0;
            roY[s][ic] = v;
            sum += (double)v;
        }
        const double factor = (double)ro[ic] / (sum > (double)kSmall ? sum : (double)kSmall);
        for (int s = 0; s < nSpecies; s++) {
            roY[s][ic] = (flow_float)((double)roY[s][ic] * factor);
        }
    }
}

}  // namespace

void speciesInit_d(solverConfig& cfg, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) {
        g_roY_dev = nullptr;
        g_nSpecies = 0;
        return;
    }
    g_nSpecies = var.nSpeciesRegistered;

    std::vector<flow_float*> h(g_nSpecies);
    for (int s = 0; s < g_nSpecies; s++) {
        h[s] = var.c_d["roY"+std::to_string(s)];
    }
    gpuErrchk( cudaMalloc((void**)&g_roY_dev, g_nSpecies*sizeof(flow_float*)) );
    gpuErrchk( cudaMemcpy(g_roY_dev, h.data(), g_nSpecies*sizeof(flow_float*), cudaMemcpyHostToDevice) );

    std::cout << "speciesInit_d: built device roY[] for nSpecies=" << g_nSpecies << "\n";
}

flow_float** species_roY_device_ptr()
{
    return g_roY_dev;
}

void speciesPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        species_primitive_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            msh.nCells_all,
            var.c_d["ro"],
            var.c_d["roY"+i],
            var.c_d["Y"+i]);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var)
{
    (void)msh;
    if (!speciesEnabled(var)) return;
    if (bc.iPlanes.empty()) return;

    const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        // M2: 全境界種別を Neumann で扱う (slip 閉領域・内部 contact の検証に十分)。
        species_neumann_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
            nb,
            bc.map_bplane_cell_d,
            bc.map_bplane_cell_ghst_d,
            var.c_d["roY"+i],
            var.c_d["Y"+i]);
    }
}

void applySpeciesBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (auto& bc : msh.bconds) {
        speciesBoundary_d_wrapper(cfg, cuda_cfg, bc, msh, var);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roY"+i], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["transport_diag_Y"+i], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["src_jac_Y"+i], 0, msh.nCells * sizeof(flow_float)));
    }

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const ScalarTransportDesc desc = buildSpeciesDesc(var, s);
        scalarTransportResidual_d(cfg, cuda_cfg, msh, var, desc);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const ScalarTransportDesc desc = buildSpeciesDesc(var, s);
        scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesRenormalize_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var) || g_roY_dev == nullptr) return;

    species_renormalize_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        g_nSpecies,
        g_roY_dev,
        var.c_d["ro"]);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!speciesEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"N"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"M"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
    }
}

void speciesUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!speciesEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"M"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
    }
}
