#pragma once

// Weiss-Smith 型低マッハ前処理の共有 device ヘルパ。
//
// 低マッハ域でフラックス散逸を音速 c でスケールすると圧力-速度カップリングが O(M) に縮退し、
// 低マッハの自励振動 (チェッカーボード) を招く。これを是正するため、基準速度 Ur と
// 前処理音速 c' を提供する。
//
// 使用箇所: convectiveFlux_d.cu の SLAU_d のみ — 質量流束の圧力散逸項 chi/c_hat·ΔP の c_hat → c'。
//   (setDT のスペクトル半径前処理・block DPLUR の固有値前処理も試したが、いずれも block DPLUR の
//    対角優位性を崩して有害と判明し不採用。詳細は計画
//    .github/plans/time_integration-lowmach-preconditioning.md §9、docs/{convection,time_integration}/theory.md)
//
// 設計・理論・検証は docs/convection/theory.md「低マッハ前処理」節を参照。

#include "flowFormat.hpp"

// 基準速度 Ur = min(c, max(|u|, eps·c))。
// 停留点フロア eps·c で Ur→0 を防ぎ、M>=1 では |u|>=c により Ur=c に飽和する。
__device__ __forceinline__
flow_float lowMachUr(flow_float c, flow_float velMag, flow_float eps)
{
    return min(c, max(velMag, eps * c));
}

// 前処理音速 c' = 0.5*sqrt((1-beta)^2 Un^2 + 4 Ur^2),  beta = (Ur/c)^2。
//   beta=1 (M>=1) で c'=c に厳密復帰、低マッハで c'→Ur=O(|u|)。
// Un は界面法線速度 (符号は二乗で消えるため任意)、velMag は速度の大きさ |u|。
__device__ __forceinline__
flow_float lowMachCprime(flow_float c, flow_float velMag, flow_float Un, flow_float eps)
{
    const flow_float Ur   = lowMachUr(c, velMag, eps);
    const flow_float beta = (Ur / c) * (Ur / c);
    const flow_float omb  = static_cast<flow_float>(1.0) - beta;
    return static_cast<flow_float>(0.5)
         * sqrt(omb * omb * Un * Un + static_cast<flow_float>(4.0) * Ur * Ur);
}

// =============================================================================
// Phase 4 (完全 Γ⁻¹A 前処理, lowMachPrecond=2) 用の土台ヘルパ。
// 計画 .github/plans/time_integration-lowmach-preconditioning.md §5 Phase 4、
// 導出は docs/time_integration/theory.md「低マッハ前処理固有系 (Weiss-Smith)」節。
// （Phase 1 の lowMachPrecond=1 経路は上の lowMachCprime のみを使い、以下には依存しない。）
// =============================================================================

// 前処理係数 β = (Ur/c)^2 ∈ (0,1]。β=1 (M>=1) で無前処理。
__device__ __forceinline__
flow_float lowMachBeta(flow_float c, flow_float velMag, flow_float eps)
{
    const flow_float Ur = lowMachUr(c, velMag, eps);
    return (Ur / c) * (Ur / c);
}

// Weiss-Smith 保存変数前処理のランク 1 共通形 I + coeff·g·rᵀ を double[5][5] に組む。
//   g = (1, u, v, w, H)ᵀ,  rᵀ = (γ-1)(ek, -u, -v, -w, 1) = ∂p/∂Q_c,  H = c²/(γ-1)+ek,  ek = ½|u|²。
//   rᵀg = c² (検証済)。Γ_c は coeff=(1-β)/(βc²)、Γ_c⁻¹ は coeff=-(1-β)/c²。
//   いずれも β=1 で coeff=0 → 単位行列に厳密復帰。条件数 ~1/β を扱うため double で組む。
__device__ __forceinline__
void lowMachGammaRank1(flow_float gamma, flow_float c,
                       flow_float u, flow_float v, flow_float w,
                       double coeff, double mat[5][5])
{
    const double gm1 = static_cast<double>(gamma) - 1.0;
    const double ek  = 0.5 * (static_cast<double>(u) * u
                            + static_cast<double>(v) * v
                            + static_cast<double>(w) * w);
    const double H   = static_cast<double>(c) * static_cast<double>(c) / gm1 + ek;
    const double g[5] = {1.0, static_cast<double>(u), static_cast<double>(v),
                         static_cast<double>(w), H};
    const double r[5] = {gm1 * ek, -gm1 * static_cast<double>(u), -gm1 * static_cast<double>(v),
                         -gm1 * static_cast<double>(w), gm1};
    #pragma unroll
    for (int i = 0; i < 5; ++i) {
        #pragma unroll
        for (int j = 0; j < 5; ++j) {
            mat[i][j] = (i == j ? 1.0 : 0.0) + coeff * g[i] * r[j];
        }
    }
}

// 前処理行列 Γ_c = I + ((1-β)/(βc²)) g rᵀ。擬似時間項 Γ_c·V/Δτ' の構築に使う。
__device__ __forceinline__
void lowMachGammaC(flow_float gamma, flow_float c,
                   flow_float u, flow_float v, flow_float w,
                   flow_float beta, double Gc[5][5])
{
    const double b = static_cast<double>(beta);
    const double coeff = (1.0 - b) / (b * static_cast<double>(c) * static_cast<double>(c));
    lowMachGammaRank1(gamma, c, u, v, w, coeff, Gc);
}

// 逆行列 Γ_c⁻¹ = I - ((1-β)/c²) g rᵀ (Sherman-Morrison)。前処理ヤコビアン P=Γ_c⁻¹A の構築に使う。
__device__ __forceinline__
void lowMachGammaCinv(flow_float gamma, flow_float c,
                      flow_float u, flow_float v, flow_float w,
                      flow_float beta, double Gcinv[5][5])
{
    const double coeff = -(1.0 - static_cast<double>(beta))
                       / (static_cast<double>(c) * static_cast<double>(c));
    lowMachGammaRank1(gamma, c, u, v, w, coeff, Gcinv);
}
