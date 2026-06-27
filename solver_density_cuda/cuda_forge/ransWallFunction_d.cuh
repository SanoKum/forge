#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// SST automatic (enhanced / y⁺ 非依存) wall treatment 用の摩擦速度算出 (methods/turbulence §6.5)。
// 壁面ごとに Reichardt 普遍速度則を逆解きして u_τ を求め、bc.bvar_d["utau"]/["ypls"] を埋める。
// cfg.wallTreatmentSST==1 かつ wall / wall_isothermal の場合のみ起動する。
// ここで算出した u_τ を ω 壁面 BC (ransBoundary) と運動量壁せん断 (viscousFlux) が共有する。
void computeWallFrictionSST_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var);

// wall-function 生産 wf_pk を全セル -1 (inactive) に初期化する (bc ループ前に 1 回呼ぶ)。
void initWallFunctionPk_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
