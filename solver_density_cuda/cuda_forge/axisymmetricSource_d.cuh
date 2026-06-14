#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

void axisymmetricSource_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
void axisymmetricGeomTerms_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// node-centered 軸対称: 軸上ノード (R=0) で半径方向運動量 roUy=0 (対称条件) を課す。
// 軸上 CV が特異点になり半径方向圧力ソースで偽の Uy が駆動されるのを防ぐ。cell モードや非軸対称では no-op。
void enforceAxisSymmetry_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// SU2 流 (MARKER_SYM): node-centered 軸対称で軸上 CV の半径方向運動量「残差」を 0 にし、
// roUy=0 (対称条件) を保つ。assembleResidual の最後 (全 flux+source 積算後) に呼ぶ。
void zeroAxisRadialResidual_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// node-centered 壁 Dirichlet: 壁ノード (wall_flag) の速度を厳密に 0 に固定する (壁ゴースト撤廃の代替)。
// enforceWallNoSlip: 一度きり (IC 確定後) に state を u=0 へ初期化 (KE を roe から除去)。
// zeroWallMomentumResidual: 毎反復 assembleResidual 末尾で壁ノードの運動量残差を 0 に射影 (陰解法整合)。
// cell モードや wall_flag 未構築では no-op。
void enforceWallNoSlip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
void zeroWallMomentumResidual_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
