#pragma once

#include "condensationProperties_d.cuh"

// 非平衡凝縮: 核生成 (CNT+Iland) と成長 (Goodheart) のソース項 device 関数 (Phase 2)。
// docs/condensation/theory.md 2-3 節。内部 double。現状 N2 のみ (n2_* 物性)。
//
// モーメント保存ベクトル後半の相変化ソース (論文 Eq.1):
//   d(ρQ0)/dt = J
//   d(ρQ1)/dt = J r*   + ρQ0 dr/dt
//   d(ρQ2)/dt = J r*^2 + 2 ρQ1 dr/dt
//   d(ρg)/dt  = (4/3)π ρ_l (J r*^3 + 3 ρQ2 dr/dt)
// J: 核生成率 [1/(m^3 s)], r*: 臨界半径 [m], dr/dt: 平均成長率 [m/s] (平均半径 r̄=Q1/Q0 で評価)。
// 一温度: T_v=T_d=T (核生成・成長とも気相 T)。

#define COND_KB 1.380649e-23      // Boltzmann [J/K]
#define COND_NA 6.02214076e23     // Avogadro [1/mol]
#define COND_PI 3.141592653589793

// 核生成 (CNT × Iland 経験補正)。J [1/(m^3 s)], r* [m] を返す。過飽和でなければ 0。
//   J_CNT = √(2σ/(π m^3))·(ρ_v²/ρ_l)·exp(-ΔG*/(k_B T)),  ΔG*=(4/3)π r*²σ,
//   r* = 2σ/(ρ_l R T ln S),  S=p_v/p_sat(T),  J = J_CNT·exp(A+B/T) [A=-55,B=4270].
__host__ __device__ inline void cond_nucleation_N2(
    double T, double p_v, double rho_v, double R, double M,
    double* J_out, double* rstar_out)
{
    const double psat = n2_psat(T);
    const double S = p_v / (psat > 1.0e-300 ? psat : 1.0e-300);
    if (S <= 1.0 || rho_v <= 0.0) { *J_out = 0.0; *rstar_out = 0.0; return; }

    const double lnS   = log(S);
    const double sigma = n2_sigma(T);
    const double rho_l = n2_rho_cond(T);
    const double rstar = 2.0*sigma / (rho_l*R*T*lnS);
    const double dG    = (4.0/3.0)*COND_PI*rstar*rstar*sigma;
    const double m     = M / COND_NA;                       // 分子質量 [kg]
    const double Kcnt  = sqrt(2.0*sigma/(COND_PI*m*m*m)) * (rho_v*rho_v/rho_l);
    const double expo  = -dG/(COND_KB*T) + (-55.0 + 4270.0/T); // CNT 障壁 + Iland 補正
    // exp 引数の下限ガード (アンダーフロー時 J=0)
    const double J = (expo > -700.0) ? Kcnt*exp(expo) : 0.0;
    *J_out = J;
    *rstar_out = rstar;
}

// 成長率 dr/dt [m/s] (Goodheart, 平均半径 r_bar で評価)。一温度 T_d=T。
//   dr/dt = (1/ρ_l) f_FS(Kn) (k R T²/L²) ln(p/p_sat(T)),
//   f_FS = (1+2Kn)/(r(1+3.42Kn+5.32Kn²))·(1-r*/r),  Kn=λ/2r。
__host__ __device__ inline double cond_growth_N2(
    double T, double p_v, double r_bar, double rstar, double R)
{
    if (r_bar <= 0.0) return 0.0;
    const double psat    = n2_psat(T);
    const double driving = log(p_v / (psat > 1.0e-300 ? psat : 1.0e-300)); // ln(p/psat) (>0 で成長)
    const double L   = n2_latent(T);
    const double k   = n2_kgas(T);
    const double lam = cond_mean_free_path(T, p_v, R);
    const double Kn  = lam/(2.0*r_bar);
    const double fFS = (1.0+2.0*Kn)/(r_bar*(1.0+3.42*Kn+5.32*Kn*Kn)) * (1.0 - rstar/r_bar);
    const double rho_l = n2_rho_cond(T);
    return (1.0/rho_l)*fFS*(k*R*T*T/(L*L))*driving;
}

// 相変化ソースベクトル S=(S_Q0,S_Q1,S_Q2,S_g)。モーメントは保存量 (ρQn)。N2。
__host__ __device__ inline void cond_source_vector_N2(
    double T, double p_v, double rho, double g, double R, double M,
    double roQ0, double roQ1, double roQ2,
    double* SQ0, double* SQ1, double* SQ2, double* Sg)
{
    double rho_v = (1.0 - g)*rho;
    if (rho_v < 0.0) rho_v = 0.0;
    double J, rstar;
    cond_nucleation_N2(T, p_v, rho_v, R, M, &J, &rstar);
    const double r_bar = (roQ0 > 1.0e-30) ? (roQ1/roQ0) : rstar;
    const double drdt  = (roQ0 > 1.0e-30) ? cond_growth_N2(T, p_v, r_bar, rstar, R) : 0.0;
    const double rho_l = n2_rho_cond(T);
    *SQ0 = J;
    *SQ1 = J*rstar + roQ0*drdt;
    *SQ2 = J*rstar*rstar + 2.0*roQ1*drdt;
    *Sg  = (4.0/3.0)*COND_PI*rho_l*(J*rstar*rstar*rstar + 3.0*roQ2*drdt);
}
