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

// ---------------------------------------------------------------------------
// Phase B/C: 条件付きスカラー Q(η) (physProp.chemistry.cmc)。methods/chemistry_cmc.md 実装 §1。
//   Q は [slice][node] の packed 配列 (slice = var*nEta + k, var = 0..n_s-1 が Y_α, n_s が h)。
//   物理空間の移流+拡散は各スライスを scalarTransport_d に渡して流用、η 拡散 (AMC) と化学 (各 η 点の点陰解) は node 毎のカーネル。
//   結合 (couple==1): PDF 平均の ω̄, Q̇̄, J̄ を chemistry_source_d に渡す (tci 2 相当)。
// ---------------------------------------------------------------------------
void cmcQInit_d(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);   // 確保 + 混合線初期化 (restart 読込後)
void cmcQTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);   // assembleResidual
void cmcQUpdate_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);      // 陰解法 inner
void applyCmcQBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void cmcQWrite(const std::string& path);   // Q(η) を生バイナリで保存 (run dir の cmc_Q.bin; restart は cmc.restartQ)
bool cmc_coupling_active();
const flow_float* cmc_omega_device_ptr();   // [ns][nCells] PDF 平均 ω̄_s [kg/m3/s]
const flow_float* cmc_qdot_device_ptr();    // [nCells] Q̇̄ [W/m3]
const flow_float* cmc_jac_device_ptr();     // [nCells*ns*ns] PDF 平均 ∂ω_s/∂(ρY_k) [1/s]
int  cmc_coupling_mode();                    // 0: なし, 1: PDF 平均ソース (ω̄), 2: PDF 積分状態への緩和 ρ̄(Ỹ_pdf−Ỹ)/τ
const flow_float* cmc_ypdf_device_ptr();    // [ns][nCells] Ỹ_pdf = Σ Ω_k Q_s(η_k)
const flow_float* cmc_hpdf_device_ptr();    // [nCells] h̃_pdf [J/kg]
const flow_float* cmc_tau_device_ptr();     // [nCells] 緩和時間 τ = relax·Δτ_local [s]
const flow_float* cmc_qcap_device_ptr();    // [nCells] couple 5: この step の発熱 (リミッタ後) [J/m3]
const flow_float* cmc_gate_device_ptr();    // [nCells] couple 6/7: 燃焼領域ゲート g∈[0,1]
