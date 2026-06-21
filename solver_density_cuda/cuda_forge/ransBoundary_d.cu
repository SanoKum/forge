#include "ransBoundary_d.cuh"
#include "ransWallFunction_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"

namespace {

constexpr flow_float kSstBeta1  = static_cast<flow_float>(0.075);
constexpr flow_float kSstBetaSt = static_cast<flow_float>(0.09);   // β*
constexpr flow_float kKappa     = static_cast<flow_float>(0.41);   // von Karman 定数
constexpr flow_float kOmegaMin  = static_cast<flow_float>(1.0e-10);
constexpr flow_float kSmall = static_cast<flow_float>(1.0e-12);

// SST 壁面 ω 境界。
// wallTreatment==0: low-Re 壁解像型 ω_w = 60ν/(β₁ y²) (docs/turbulence §6.1)。
// wallTreatment==1: Menter automatic ブレンド ω_w = √(ω_vis² + ω_log²) (§6.5 (b))。
//   ω_vis = 6ν/(β₁ y²), ω_log = u_τ/(√β* κ y)。u_τ は事前に utau に格納済み。
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
    flow_float* omegab,
    int wallTreatment,
    flow_float* utau)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;

    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];

        const flow_float rho_c = max(ro[ic], kSmall);
        const flow_float nu_c = vis_lam[ic] / rho_c;
        const flow_float y_w = max(wall_dist[ic], kSmall);

        // k は mode に関わらず k_w = 0 (Dirichlet)。教科書的には wall-function では k を
        // zero-gradient にするが、それは近壁 P_k も log則 τ_w から計算し第一セル k を平衡値
        // u_τ²/√β* に固定する整合とセットで成立する。forge は P_k が解像勾配のままなので
        // (plan §10 将来課題)、k を zero-gradient にすると誤った解像生産で近壁 k が暴走し
        // μ_t 過大 → u_τ 過大 → Cf 過大になる (検証: y⁺30 で 0.89→1.80, y⁺10 で 1.10→1.68)。
        // P_k の wall-function 化までは k=0 Dirichlet の方が良い (docs/turbulence §6.5 (b'))。
        kb[ib] = static_cast<flow_float>(0.0);

        if (wallTreatment == 1) {
            const flow_float omega_vis = static_cast<flow_float>(6.0) * nu_c / (kSstBeta1 * y_w * y_w);
            const flow_float omega_log = utau[ib] / (sqrt(kSstBetaSt) * kKappa * y_w);
            omegab[ib] = max(sqrt(omega_vis * omega_vis + omega_log * omega_log), kOmegaMin);
        } else {
            omegab[ib] = max(static_cast<flow_float>(60.0) * nu_c / (kSstBeta1 * y_w * y_w), kOmegaMin);
        }

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

void ransBoundary_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var)
{
    (void)msh;

    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1)) {
        return;
    }

    if (bc.iPlanes.empty()) {
        return;
    }

    if (bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal") {
        // automatic wall treatment (wallTreatmentSST==1) では先に摩擦速度 u_τ を Reichardt 逆解きで
        // 求め bc.bvar_d["utau"] に格納する。ω 壁 BC と viscousFlux 壁せん断が同じ u_τ を共有する。
        computeWallFrictionSST_d_wrapper(cfg , cuda_cfg , bc , msh , var);

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
            bc.bvar_d["omegab"],
            cfg.wallTreatmentSST,
            bc.bvar_d["utau"]);
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