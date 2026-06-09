#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// 多成分化学種輸送 (M2)。汎用スカラ輸送コア scalarTransport_d (ScalarTransportDesc) を
// 化学種ごとに再利用し、保存質量分率 ρY_s を流れと共に移流する (M2 は移流のみ。拡散は M4)。
// 単成分 (var.nSpeciesRegistered < 2) のときは全て no-op で、M1 と同一経路を保つ。

// device の roY ポインタ配列 (flow_float*[nSpecies]) を 1 度だけ構築する。
// allocVariables 後 (c_d["roY{s}"] が確定後) に呼ぶこと。
void speciesInit_d(solverConfig& cfg, variables& var);

// dependentVariables_d_wrapper へ渡す device roY 配列ポインタ。
// 単成分時は nullptr (混合則 thermo は Y={1} に縮退)。
flow_float** species_roY_device_ptr();

// 原始質量分率 Y_s = ρY_s/ρ を全セル (ghost 含む) について更新する。
void speciesPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 化学種 ghost を Neumann (zero-gradient) で埋める。bcond ごとに呼ぶ applySpeciesBoundaries から使用。
void speciesBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var);

// 全 bcond をループして化学種 ghost を埋める (applyRansScalarBoundaries と同形)。
void applySpeciesBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 化学種移流残差を組み立てる (res_roY{s} / transport_diag をゼロ初期化してから集計)。
// 粘性 (viscMethod!=0) かつ nSpecies>=2 のとき M4 の Fick 拡散 + ΣJ=0 補正 + エンタルピー拡散も加える。
void speciesTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 化学種の時間積分 (scalarTimeIntegration_d を化学種ごとに呼ぶ)。
void speciesTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 化学種の実現可能性・再正規化: ρY_s>=0 にクランプし Σ_s ρY_s = ρ となるよう再スケール (ΣY_s=1)。
void speciesRenormalize_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// RK ステップ/ステージ始点の保存 (roY{s}N / roY{s}M)。NS の updateVariablesOuter/Inner に対応。
void speciesUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void speciesUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
