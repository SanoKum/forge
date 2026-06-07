#include "ransScalarBoundary_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"

namespace {

constexpr flow_float kSstBeta1 = static_cast<flow_float>(0.075);
constexpr flow_float kOmegaMin = static_cast<flow_float>(1.0e-10);
constexpr flow_float kSmall = static_cast<flow_float>(1.0e-12);

__global__ void rans_wall_scalar_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* ro,
    flow_float* vis_lam,
    flow_float* wall_dist,
    flow_float* k,
    flow_float* omega,
    flow_float* kb,
    flow_float* omegab)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;

    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];

        const flow_float rho_c = max(ro[ic], kSmall);
        const flow_float nu_c = vis_lam[ic] / rho_c;
        const flow_float y_w = max(wall_dist[ic], kSmall);

        kb[ib] = static_cast<flow_float>(0.0);
        omegab[ib] = max(static_cast<flow_float>(60.0) * nu_c / (kSstBeta1 * y_w * y_w), kOmegaMin);

        k[ig] = static_cast<flow_float>(2.0) * kb[ib] - k[ic];
        omega[ig] = static_cast<flow_float>(2.0) * omegab[ib] - omega[ic];
    }
}

__global__ void rans_dirichlet_scalar_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* k,
    flow_float* omega,
    flow_float* kb,
    flow_float* omegab)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;

    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];

        k[ig] = static_cast<flow_float>(2.0) * kb[ib] - k[ic];
        omega[ig] = static_cast<flow_float>(2.0) * omegab[ib] - omega[ic];
    }
}

__global__ void rans_neumann_scalar_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* k,
    flow_float* omega,
    flow_float* kb,
    flow_float* omegab)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;

    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];

        k[ig] = k[ic];
        omega[ig] = omega[ic];
        kb[ib] = k[ic];
        omegab[ib] = omega[ic];
    }
}

}

void ransScalarBoundary_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var)
{
    (void)msh;

    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1)) {
        return;
    }

    if (bc.iPlanes.empty()) {
        return;
    }

    if (bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal") {
        rans_wall_scalar_boundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
            static_cast<geom_int>(bc.iPlanes.size()),
            bc.map_bplane_cell_d,
            bc.map_bplane_cell_ghst_d,
            var.c_d["ro"],
            var.c_d["vis_lam"],
            var.c_d["wall_dist"],
            var.c_d["k"],
            var.c_d["omega"],
            bc.bvar_d["kb"],
            bc.bvar_d["omegab"]);
        return;
    }

    if (bc.bcondKind == "inlet_uniformVelocity" ||
        bc.bcondKind == "inlet_fluctVelocity" ||
        bc.bcondKind == "inlet_Pressure" ||
        bc.bcondKind == "inlet_Pressure_dir") {
        rans_dirichlet_scalar_boundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
            static_cast<geom_int>(bc.iPlanes.size()),
            bc.map_bplane_cell_d,
            bc.map_bplane_cell_ghst_d,
            var.c_d["k"],
            var.c_d["omega"],
            bc.bvar_d["k"],
            bc.bvar_d["omega"]);
        return;
    }

    if (bc.bcondKind == "slip" ||
        bc.bcondKind == "axis" ||
        bc.bcondKind == "outlet_statPress" ||
        bc.bcondKind == "outflow" ||
        bc.bcondKind == "periodic") {
        rans_neumann_scalar_boundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
            static_cast<geom_int>(bc.iPlanes.size()),
            bc.map_bplane_cell_d,
            bc.map_bplane_cell_ghst_d,
            var.c_d["k"],
            var.c_d["omega"],
            bc.bvar_d["kb"],
            bc.bvar_d["omegab"]);
        return;
    }
}