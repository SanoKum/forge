#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// WMLES 代数壁応力モデル (Reichardt + Kader)。理論: methods/turbulence/theory.md §10、
// 計画: plans/active/turbulence-wmles-wall-stress.md。
// 壁境界面ごとにマッチング点 (壁から 1 つ内側の解点) の瞬時 LES 状態から τ_w / q_w を算出し、
//   cell: viscousFlux_wall_d の壁面粘性流束を置換 (wallTreatment==2, twall_*/qwall bvar 経由)
//   node: Tau_Wall / Qw_Wall (per-CV) に格納し viscousFlux_d の W↔I AddTauWall / AddQWall で適用
// する。SST automatic wall treatment (wallTreatmentSST) とはゲートを分離 (LESorRANS!=2 のみ有効)。

// この bcond で WMLES 壁モデルが有効か (kind が wall 系 かつ ints: wallModelLES==1 かつ 非 RANS)
bool wmlesActiveForBcond(const solverConfig& cfg, const bcond& bc);

// いずれかの壁 bcond で WMLES が有効か (viscousFlux の Tau_Wall/Qw_Wall 受け渡しゲート用)
bool wmlesAnyActive(const solverConfig& cfg, const mesh& msh);

// node モードで WMLES が有効か
bool wmlesNodeActive(const solverConfig& cfg, const mesh& msh);

// node モードで等温 WMLES 壁があるか (壁ノード T ピン + res_roe ゼロ化ゲート用)
bool wmlesNodeIsothermalActive(const solverConfig& cfg, const mesh& msh);

// 全 wall bcond の WMLES 壁モデルを実行する (applyBconds の後・粘性流束の前に毎 step 呼ぶ)。
// node では Tau_Wall/Qw_Wall の -1 初期化と等温壁ノードの T ピンも行う。
void applyWmlesWallModel(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 素の node 等温壁 (非 WMLES) の壁ノード温度ピンが有効か (node && nodeWallDirichlet && wall_isothermal あり)
bool nodeIsothermalPinActive(const solverConfig& cfg, const mesh& msh);

// 素の node 等温壁の壁ノード温度状態ピン (nodeWallDirichlet の熱版。applyWmlesWallModel と同位相で呼ぶ)
void applyNodeIsothermalWallPin(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 素の node 等温壁の壁ノード res_roe ゼロ化 (zeroWallMomentumResidual と同位相で呼ぶ)
void zeroNodeIsothermalEnergyResidual(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
