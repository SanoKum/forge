#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5)。
// setPeriodicPartner + buildPeriodicNodeGroups で構築した周期ノード group (periodicRoot) に対し、
// 保存量残差 res_* を group 全員で足し合わせ、全員へ同じ和を書き戻す (gather + broadcast)。
// 合併体積 (buildPeriodicNodeGroups で vol を group 合算) と組み合わせると、両側部分 CV が
// 同 res・同 vol で更新され「1 つの CV」として bit 一致同期する。assembleResidual の末尾
// (全 flux/source 積算 + 壁/軸射影の後) に毎反復呼ぶ。cell モード / 非周期では no-op。
void periodicNodeGather_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
