#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

void axisymmetricSource_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
void axisymmetricGeomTerms_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// axisymMethod==1 (SU2 流 planar+ソース): 非粘性+粘性の 1/y 軸対称ソースを res_ro..res_roe に加算。
// AuxVar (μv/y 系) の GG 勾配も内部で計算。method!=1 では no-op。SST の 1/y 項は ransSource 側。
void axisymmetricSourceSU2_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// node-centered 軸対称: 軸上ノード (R=0) で半径方向運動量 roUy=0 (対称条件) を課す。
// 軸上 CV が特異点になり半径方向圧力ソースで偽の Uy が駆動されるのを防ぐ。cell モードや非軸対称では no-op。
void enforceAxisSymmetry_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// SU2 流 (MARKER_SYM): node-centered 軸対称で軸上 CV の半径方向運動量「残差」を 0 にし、
// roUy=0 (対称条件) を保つ。assembleResidual の最後 (全 flux+source 積算後) に呼ぶ。
void zeroAxisRadialResidual_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// nodeAxisDirichlet=1: 軸上ノード状態を radial 代表点 (axis_rep) からの対称 Dirichlet で置換する
// (∂q/∂r=0 の 1 次、u_r=0)。assembleResidual 冒頭 (dependentVariables の前) に呼ぶ。
// 軸半 CV の真空過膨張 (plan boundary-node-nozzle-wall-outlet-stability §2.6) の根治。
void enforceAxisMirror_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// nodeAxisDirichlet=1: 軸上ノードの全保存量残差 (+RANS 時 roK/roOmega) を 0 化する。
// assembleResidual 末尾で zeroAxisRadialResidual の直後に呼ぶ (rms 汚染防止 + explicit 状態固定)。
void zeroAxisAllResiduals_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// node-centered 壁 Dirichlet: 壁ノード (wall_flag) の速度を厳密に 0 に固定する (壁ゴースト撤廃の代替)。
