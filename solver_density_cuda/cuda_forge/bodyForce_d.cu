#include "bodyForce_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"

// 一様体積力ソース。設計は bodyForce_d.cuh を参照。

namespace {

__global__ void bodyForce_d
(
    geom_int nCells,
    flow_float fx, flow_float fy, flow_float fz,
    flow_float* volume,
    flow_float* Ux, flow_float* Uy, flow_float* Uz,
    flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe
)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float V = volume[ic];
    res_roUx[ic] += fx*V;
    res_roUy[ic] += fy*V;
    res_roUz[ic] += fz*V;
    // 体積力の仕事 (エネルギー整合)。これを落とすと駆動運動量だけ増えてエネルギー収支が破れる。
    res_roe[ic]  += (fx*Ux[ic] + fy*Uy[ic] + fz*Uz[ic])*V;
}

} // namespace

void bodyForce_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.bodyForceX == 0.0 && cfg.bodyForceY == 0.0 && cfg.bodyForceZ == 0.0) return;
    bodyForce_d<<<cuda_cfg.dimGrid_normalcell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        cfg.bodyForceX, cfg.bodyForceY, cfg.bodyForceZ,
        var.c_d["volume"],
        var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"],
        var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
