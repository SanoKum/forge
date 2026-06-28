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

// block-DPLUR の各 sweep 後に呼ぶ (§4.5.7)。周期 group の補正 dq を root(master) から member(slave) へ
// ミラー (slave=master) し、同一視ノードが多値化 (drift) して発散するのを防ぐ。master/slave は別 plane を
// 持つため別 dq が出るが、両者は同じ DOF なので master 解を共有させる。blockDPLUR=1 は dq_block_old_*、
// =0 は dq_*_old を対象にする。cell/非周期では no-op。
void periodicMirrorDq_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 勾配の periodic gather (§4.5 拡張)。calcGradient 後に呼び、合併体積を利用して boundary periodic node の
// Green-Gauss 勾配を「和→broadcast」で厳密合併に直す (片側勾配 → 両側)。2次再構成・粘性の精度に効く。
// 前提: calcGradient_b_d で periodic 半割面を除外しておく。非軸対称限定。cell/非周期では no-op。
void periodicGradientGather_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
