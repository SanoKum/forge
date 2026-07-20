#include "bodyForce_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

// 一様体積力ソース。設計は bodyForce_d.cuh を参照。

namespace {

// (roUx, volume) → double 積。float32 ビルドでも積算は double で行う (1M CV の桁落ち回避)。
struct WeightedMomentum {
    __host__ __device__ double operator()(const thrust::tuple<flow_float, flow_float>& t) const {
        return static_cast<double>(thrust::get<0>(t)) * static_cast<double>(thrust::get<1>(t));
    }
};

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

void bodyForceCtrlUpdate(solverConfig& cfg , mesh& msh , variables& var , int iStep)
{
    if (cfg.bodyForceCtrl != 1) return;

    // file-static 状態 (単一ソルバインスタンス前提)。V は周期 seam の合併体積二重計上を含むが、
    // M も同一レンジで測るため目標 M_t = target·V で二重計上は自己整合的に相殺する。
    static bool   initialized = false;
    static double Vtot = 0.0;
    static double Mprev = 0.0;
    static std::ofstream hist;

    const geom_int n = msh.nCells;
    thrust::device_ptr<flow_float> roUx(var.c_d["roUx"]);
    thrust::device_ptr<flow_float> vol(var.c_d["volume"]);

    const double M = thrust::transform_reduce(
        thrust::make_zip_iterator(thrust::make_tuple(roUx, vol)),
        thrust::make_zip_iterator(thrust::make_tuple(roUx + n, vol + n)),
        WeightedMomentum(), 0.0, thrust::plus<double>());

    if (!initialized) {
        Vtot = thrust::reduce(vol, vol + n, 0.0, thrust::plus<double>());
        Mprev = M;                    // 初回は P 項のみで立ち上げる
        hist.open("bodyforce_history.csv");
        hist << "step,time,fx,M,M_target\n";
        initialized = true;
    }

    // 場が NaN 化していたら制御を止めて現状維持 (NaN を体積力経由で全域へ撒かない)。
    if (!std::isfinite(M)) {
        if (iStep % 100 == 0) {
            std::cerr << "[bodyForceCtrl] WARNING: non-finite momentum integral, skipping update\n";
        }
        return;
    }

    const double Mt = static_cast<double>(cfg.bodyForceCtrlTarget) * Vtot;
    const double dt = static_cast<double>(cfg.dt);
    const double gamma = static_cast<double>(cfg.bodyForceCtrlRelax);
    const double dfx = gamma * (Mt - 2.0 * M + Mprev) / (Vtot * std::max(dt, 1.0e-30));
    cfg.bodyForceX = static_cast<flow_float>(static_cast<double>(cfg.bodyForceX) + dfx);
    Mprev = M;

    if (iStep < 1000 || iStep % 10 == 0) {
        hist << iStep << "," << cfg.totalTime << "," << cfg.bodyForceX << ","
             << M << "," << Mt << "\n";
        if (iStep % 100 == 0) hist.flush();
    }
    if (iStep % 100 == 0) {
        std::cout << "[bodyForceCtrl] fx=" << cfg.bodyForceX
                  << "  M/Mt=" << M / Mt << "\n";
    }
}
