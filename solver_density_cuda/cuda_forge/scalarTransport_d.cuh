#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// 汎用スカラ輸送コア。特定の物理 (RANS / species など) には依存せず、
// 1 変数ぶんの保存量 ρφ について移流・拡散の残差組み立てと時間積分を提供する。
// RANS k/ω への適用は ransTransport_d.* が ScalarTransportDesc を構築して呼び出す。
struct ScalarTransportDesc {
    flow_float* phi;
    flow_float* rho_phi;
    flow_float* rho_phi_N;
    flow_float* rho_phi_M;
    flow_float* res_rho_phi;
    flow_float* res_rho_phi_m;
    flow_float* src_jac;        // 源項消散ヤコビアン対角 (∂D/∂(ρφ)≥0)。陽解法 RK の point-implicit 減衰に使用
    flow_float* transport_diag; // 輸送(移流+拡散)ヤコビアン対角 [m³/s]。advection/diffusion kernel で集計
    flow_float floor;           // realizability 下限。陰解法 point-implicit と整合
    flow_float sigma;           // 拡散係数スケール (有効粘性 = vis_lam + sigma·vis_turb)
};

// 1 変数ぶんの移流 + (任意) 拡散残差を組み立てる。
// res_rho_phi / transport_diag は呼び出し側でゼロ初期化済みであること。同期は行わない。
void scalarTransportResidual_d(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                               const ScalarTransportDesc& desc);

// 1 変数ぶんの時間積分 (RK4 / point-implicit RK)。同期は行わない。
void scalarTimeIntegration_d(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                             const ScalarTransportDesc& desc);
