#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// RANS (SST k-ω) のスカラ輸送オーケストレーション。
// 汎用コア scalarTransport_d.* の輸送・時間積分を k/ω に適用し、
// k/ω 勾配 (源項で使用) を計算する。
void ransTransport_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
void ransTimeIntegration_d_wrapper(int loop , solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
void ransGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
