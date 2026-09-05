#pragma once
// =============================================================================
// 混合分率インフラ / CMC (Conditional Moment Closure) — methods/chemistry_cmc.md
//   Phase A (本ファイル): Bilger 混合分率 ξ̃ の診断、分散 ξ''² の輸送 (生成 2μ_t/Sc_t|∇ξ̃|², 散逸 ρχ̃)、
//   χ̃ = cChi β* ω ξ''² の診断。physProp.chemistry.mixfrac.enabled==1 のときだけ動く (既定 off で全経路ビット不変)。
//   Phase B/C (条件付きスカラー Q(η)・β-PDF ソース結合) は同じモジュールに追加する。
// =============================================================================
#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"
#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// Bilger 係数を機構の元素組成と config の燃料/酸化剤モル分率から作る。chemistry_init の後、allocVariables の後に 1 度。
void cmcInit_d(solverConfig& cfg, variables& var);

// ξ̃ = Bilger(Ỹ) と ∇ξ̃、∇ξ''² を更新する (assembleResidual 冒頭、species primitive の後)。
void cmcMixfrac_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 分散 roXiVar の移流+拡散残差 (scalarTransport) と生成/散逸ソース (散逸は src_jac で点陰解)、χ̃ 診断。
void cmcVarianceTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 時間積分 (陽解法 RK / 陰解法 point-implicit) — scalarTimeIntegration_d に委譲。
void cmcVarianceTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 原始量 xiVar = roXiVar/ρ と実現可能性クランプ (0 ≤ ξ''² ≤ ξ̃(1-ξ̃))。
void cmcPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 境界: 入口は分散 0 (Dirichlet, cell は ghost へ)、他は zero-gradient。
void applyCmcBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// N/M 始点の更新 (凝縮モーメントと同じ規約)。
void cmcUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void cmcUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

inline bool cmcMixfracEnabled(const solverConfig& cfg, const variables& var)
{
    return cfg.chemMixfrac != 0 && var.mixfracRegistered != 0;
}
