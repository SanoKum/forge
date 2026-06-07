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
