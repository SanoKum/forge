#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// 一様体積力ソース (config top-level `bodyForce: [fx, fy, fz]`, 既定 [0,0,0]=off)。
// 周期チャネル駆動 (wmles plan §5-7) 等に使う:
//   res_roU_i += f_i · V,   res_roe += (f · u) · V
// 平衡では f_x·δ = τ_w が成り立ち、u_τ = √(f_x δ/ρ) を直接指定できる。
// 軸対称は未対応 (config 読込で拒否)。陰解法では明示ソース (Jacobian 無し)。
void bodyForce_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
