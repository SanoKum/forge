#include "axisymmetricSource_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"


// 軸対称 (B 流儀) — 半径方向運動量の圧力・粘性ソース項を residual に加算する。
//
// 円柱座標で書いた半径運動量保存則を r 倍してから (r,z) 平面で積分すると、
// 圧力勾配項のうち r 倍では落ちないテルムとして
//
//     R_{rho u_r} = ∫∫ (p / r) · r dr dz = p_cell · A_planar
//
// が残る。粘性は $-\tau_{\theta\theta}$ を同じ面積で積分した source として加える
// (詳しくは docs/axisymmetric/theory.md)。
//
// res_roUy には -∑(F·S) が累積されており、時間進行は
//   roUy^{n+1} = roUy^n + coef_DT * res_roUy * dt / V
// で書かれる (timeIntegration_d.cu)。よってソース項は res_roUy に
// "+= (p - tau_theta_theta) * A_planar" として加算する。
__global__ void axisymmetricSource_d
(
    geom_int nCells,
    flow_float* P,
    flow_float* A_planar,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float* axisym_divU,
    flow_float* axisym_uy_over_r,
    flow_float* res_roUy
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float mu_total = vis_lam[ic] + vis_turb[ic];
        const flow_float tau_theta_theta =
            (flow_float)2.0 * mu_total * axisym_uy_over_r[ic]
            - (flow_float)(2.0 / 3.0) * mu_total * axisym_divU[ic];
        res_roUy[ic] += (P[ic] - tau_theta_theta) * A_planar[ic];
    }
}

// 軸対称ソース項の前計算。半径方向ひずみ u_r/r と発散 divU = ∂xUx + ∂yUy + u_r/r を
// セル毎に求める。これらは turbulent_viscosity (S_θθ) / ransSource / axisymmetricSource が
// 消費するため、それらより前に実行する必要がある。
// 有効半径 r_eff は revolved 体積 (per-radian, V = r·A_planar) を planar 面積で割って得る。
__global__ void axisymmetricGeomTerms_d
(
    geom_int nCells,
    flow_float* Uy,
    flow_float* volume,
    flow_float* A_planar,
    flow_float* dUxdx,
    flow_float* dUydy,
    flow_float* axisym_uy_over_r,
    flow_float* axisym_divU
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float area = A_planar[ic];
        const flow_float r_eff = (area > (flow_float)0.0) ? volume[ic] / area : (flow_float)0.0;

        if (r_eff > (flow_float)0.0) {
            axisym_uy_over_r[ic] = Uy[ic] / r_eff;
        } else {
            axisym_uy_over_r[ic] = (flow_float)0.0;
        }

        axisym_divU[ic] = dUxdx[ic] + dUydy[ic] + axisym_uy_over_r[ic];
    }
}


void axisymmetricSource_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isAxisymmetric != 1) return;

    axisymmetricSource_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["P"],
        var.c_d["A_planar"],
        var.c_d["vis_lam"],
        var.c_d["vis_turb"],
        var.c_d["axisym_divU"],
        var.c_d["axisym_uy_over_r"],
        var.c_d["res_roUy"]
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void axisymmetricGeomTerms_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isAxisymmetric != 1) return;

    axisymmetricGeomTerms_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["Uy"],
        var.c_d["volume"],
        var.c_d["A_planar"],
        var.c_d["dUxdx"],
        var.c_d["dUydy"],
        var.c_d["axisym_uy_over_r"],
        var.c_d["axisym_divU"]
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}
