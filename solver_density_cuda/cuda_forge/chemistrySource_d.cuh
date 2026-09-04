#pragma once
// 有限速度化学ソース項の host 側インタフェース (chemistry_d.cu)。理論・設計: methods/chemistry.md。
#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"
#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"
#include "cuda_forge/chemistry_d.cuh"

// chemistry.enabled==1 かつ TP 多成分のとき true。
bool chemistryEnabled(const solverConfig& cfg);

// 反応機構ファイルを読み ReactionTable を device へアップロードする。thermo_init_db の後に 1 度呼ぶ。
// chemistry.enabled==0 では何もしない。
void chemistry_init(solverConfig& cfg);

// jacobianMode==2 (種ブロック point-implicit) 用アクセサ。ブロックが無効 (mode<2 / chemistry off) なら nullptr / false。
//   chemistry_jac_device_ptr(): R [nCells*ns*ns] (R_sk = J_total_sk + d_s δ_sk, d_s は src_jac_Y に入れた対角消費部分)
//   chemistry_jacroe_device_ptr(): max(0,−∂Q̇/∂(ρe)) [nCells] → block-DPLUR の (5,5) 対角へ
flow_float* chemistry_jac_device_ptr();
flow_float* chemistry_jacroe_device_ptr();
//   chemistry_cq_device_ptr(): ∂Q̇/∂(ρY_k) [nCells*ns]。案C 予測子 δ(ρY)* から反応熱を線形化して陰的に注入する用
//   (res_roe += V Σ_k ∂Q̇/∂(ρY_k) δ(ρY_k)*)。speciesEOSCrossPredictInject が使う。
flow_float* chemistry_cq_device_ptr();
bool chemistry_block_active();

// host 側の反応表 (未初期化なら nullptr)。
const ReactionTable* chemistry_table_host();

// assembleResidual 内 (speciesTransport の直後) で呼ぶ: res_roY{s} += Vω_s, res_roe += VQ̇,
// src_jac_Y{s} += max(0,−∂ω_s/∂(ρY_s))。chemistry off / 単成分では no-op。
void chemistrySource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
