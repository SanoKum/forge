#pragma once

// =============================================================================
// species_eos_coupling_d.cuh
//   多成分 TP 陰解法の案C (block-triangular roe↔roY coupling) の中核となる
//   EOS クロス応答の解析評価。
//
//   組成変化を flow 密度固定の接空間成分 z_s = ρ δY_s (Σ_s z_s = 0) で表すとき、
//   保存内部エネルギー e 一定条件下での温度・圧力の一次応答は (数値微分なし):
//
//     δT_Y = -(1/(ρ c_{v,mix})) Σ_s e_s(T) z_s
//     δp_Y =   T Σ_s R_s z_s  -  (R_mix/c_{v,mix}) Σ_s e_s(T) z_s
//
//   ここで R_mix=Σ Y_s R_s, c_{v,mix}=c_{p,mix}(T)-R_mix, e_s(T)=h_s(T)-R_s T。
//   e_s(T) は保存エネルギー (roe) と同一の NASA-9 datum (thermo_h_mass) を使う。
//
//   これは ∂p/∂(ρY_s) の全列を陽に行列保存せず、species 仮更新方向 z に対する
//   Jacobian-vector product δp_Y = Σ_s (∂p/∂(ρY_s)) z_s を直接評価したもの。
//   本番計算はこの解析式のみを使い、有限差分は検証ハーネス専用とする。
//
//   同一コードを __host__ __device__ で提供し、device カーネルと host 単体検証で
//   ビット同一の評価を保証する (thermo_d.cuh と同じ THERMO_HD 規約)。
// =============================================================================

#include "thermo_d.cuh"

// 組成接空間方向 z_s = ρ δY_s (Σ z_s = 0 を仮定) に対する EOS クロス応答。
//   sp  : 化学種 thermo 配列, n: 化学種数, Y: 質量分率 (Σ Y=1), T: 温度 [K], ro: 密度 [kg/m^3]
//   z   : ρ δY_s [kg/m^3] (length n)
//   dpY : δp_Y [Pa] (nullptr 可), dTY: δT_Y [K] (nullptr 可)
THERMO_HD void species_eos_cross_response(
    const SpeciesThermo* sp, int n, const double* Y, double T, double ro,
    const double* z, double* dpY, double* dTY)
{
    const double R_mix  = thermo_R_mix(sp, n, Y);
    const double cp_mix = thermo_cp_mix(sp, n, Y, T);
    const double cv_mix = cp_mix - R_mix;
    // thermo_T_from_e と同じ 0 割フロア (cv が小さい領域での暴走防止)。
    const double cvf = (cv_mix > 1.0e-2 * R_mix ? cv_mix : 1.0e-2 * R_mix);

    double sumRz = 0.0;   // Σ_s R_s z_s
    double sumEz = 0.0;   // Σ_s e_s(T) z_s
    for (int s = 0; s < n; ++s) {
        const double R_s = thermo_R_species(sp[s]);
        const double e_s = thermo_h_mass(sp[s], T) - R_s * T;   // e_s(T) = h_s(T) - R_s T
        sumRz += R_s * z[s];
        sumEz += e_s * z[s];
    }
    if (dTY) *dTY = -sumEz / (ro * cvf);
    if (dpY) *dpY =  T * sumRz - (R_mix / cvf) * sumEz;
}
