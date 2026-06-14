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

// 核生成 (CNT × 種ごと補正)。J [1/(m^3 s)], r* [m]。過飽和でなければ 0。cprops で N2/H2O を切替。
//   J_CNT = √(2σ/(π m^3))·(ρ_v²/ρ_l)·exp(-ΔG*/(k_B T)),  ΔG*=(4/3)π r*²σ,  r*=2σ/(ρ_l R T ln S)。
//   N2: ×exp(A+B/T) [Iland]。H2O: 等温 CNT (キャリア N2 がクラスタを熱平衡化し Kantrowitz 非等温抑制は
//        ~1 に量子化されるため、希薄水-N2 では等温近似。過抑制を避ける; 必要なら carrier-Kantrowitz を後段)。
// 核生成 (CNT × 種ごと補正 × 任意 Kantrowitz 非等温補正)。
//   kantrowitz!=0 のとき J を 1/(1+θ) 倍する (Feder/Kantrowitz 非等温補正):
//     θ = 2(γ-1)/(γ+1) · b(b-1/2),  b = L/(R_v T)  (γ=熱を運ぶ気相の比熱比 gamma_gas)。
//   純蒸気では θ が大きく J を桁で抑える。キャリア気体 (N2) 中ではキャリアが潜熱を奪い θ→小 となるが、
//   ここでは感度評価用に純蒸気形 (gamma_gas=carrier γ) をそのまま掛ける on/off スイッチとして実装。
__host__ __device__ inline void cond_nucleation(
    const CondSpeciesProps& cp, double T, double p_v, double rho_v,
    double* J_out, double* rstar_out,
    int kantrowitz = 0, double gamma_gas = 1.4)
{
    const double psat = cond_psat(cp, T);
    const double S = p_v / (psat > 1.0e-300 ? psat : 1.0e-300);
    if (S <= 1.0 || rho_v <= 0.0) { *J_out = 0.0; *rstar_out = 0.0; return; }
    const double R     = cp.R;
    const double lnS   = log(S);
    const double sigma = cond_sigma(cp, T);
    const double rho_l = cond_rho_cond(cp, T);
    const double rstar = 2.0*sigma / (rho_l*R*T*lnS);
    const double dG    = (4.0/3.0)*COND_PI*rstar*rstar*sigma;
    const double m     = cp.M / COND_NA;
    const double Kcnt  = sqrt(2.0*sigma/(COND_PI*m*m*m)) * (rho_v*rho_v/rho_l);
    double corr = 1.0;
    if (cp.model == COND_MODEL_N2) corr = exp(-55.0 + 4270.0/T); // Iland 経験補正
    // H2O: corr=1 (carrier-thermalized 等温 CNT)
    if (kantrowitz) {
        const double b = cond_latent(cp, T) / (R*T);
        const double theta = (2.0*(gamma_gas-1.0)/(gamma_gas+1.0)) * b * (b - 0.5);
        corr /= (1.0 + (theta > 0.0 ? theta : 0.0));   // 非等温抑制 (θ<0 は掛けない)
    }
    const double expo = -dG/(COND_KB*T);
    const double J = (expo > -700.0) ? Kcnt*exp(expo)*corr : 0.0;
    *J_out = J;
    *rstar_out = rstar;
}

// 成長率 dr/dt [m/s]。一温度 T_d=T。
//   既定 (growthModel=0): N2=Goodheart, H2O=Hertz-Knudsen(自由分子, 質量輸送律速)。
//   growthModel=1 (Gyarmathy): 熱伝導律速の Gyarmathy(1982) 式 (極超音速ノズル凝縮で標準)。
//     dr/dt = (k R T²/L²) ln(S) /ρ_l · (1-r*/r)/(r(1+3.18 Kn)),  Kn=l/(2r) (キャリア気体の l)。
//     N2 Goodheart と前因子 (k R T² ln S /ρ_l L²) を共有し、Knudsen 内挿のみ 1/(1+3.18Kn) に差し替えた形。
//   p_gas: Kn 用の気相全圧 (<0 で p_v にフォールバック=純蒸気)。k_gas は n2_kgas(T) (キャリア N2 熱伝導率)。
__host__ __device__ inline double cond_growth(
    const CondSpeciesProps& cp, double T, double p_v, double r_bar, double rstar,
    int growthModel = 0, double p_gas = -1.0, double gyarC = 3.18)
{
    if (r_bar <= 0.0) return 0.0;
    const double R    = cp.R;
    const double psat = cond_psat(cp, T);
    const double rho_l= cond_rho_cond(cp, T);
    const double pK   = (p_gas > 0.0) ? p_gas : p_v;   // Kn 用は全圧 (carrier)、純蒸気では p_v
    if (growthModel == 1) {
        // Gyarmathy 熱伝導律速 (種共通)。過飽和 → 過冷却 ΔT_sub=(R T²/L)ln S を駆動力に。
        // Knudsen 補正 1/(1+gyarC·Kn): gyarC は Gyarmathy 標準 3.18 (config condGyarmathyC で可変)。
        if (p_v <= psat) return 0.0;
        const double driving = log(p_v / (psat > 1.0e-300 ? psat : 1.0e-300));
        const double L   = cond_latent(cp, T);
        const double k   = n2_kgas(T);                 // キャリア (N2) 熱伝導率
        const double lam = cond_mean_free_path(T, pK, R);
        const double Kn  = lam/(2.0*r_bar);
        const double fac = (1.0 - rstar/r_bar) / (r_bar*(1.0 + gyarC*Kn));
        return (1.0/rho_l)*fac*(k*R*T*T/(L*L))*driving;
    }
    if (cp.model == COND_MODEL_H2O) {
        const double sigma = cond_sigma(cp, T);
        const double pd = psat*exp(2.0*sigma/(rho_l*R*T*r_bar)); // Kelvin 効果
        const double alpha = 1.0; // 質量適応係数
        return (alpha/rho_l)*(p_v - pd)/sqrt(2.0*COND_PI*R*T);
    } else {
        const double driving = log(p_v / (psat > 1.0e-300 ? psat : 1.0e-300));
        const double L   = cond_latent(cp, T);
        const double k   = n2_kgas(T);
        const double lam = cond_mean_free_path(T, pK, R);
        const double Kn  = lam/(2.0*r_bar);
        const double fFS = (1.0+2.0*Kn)/(r_bar*(1.0+3.42*Kn+5.32*Kn*Kn)) * (1.0 - rstar/r_bar);
        return (1.0/rho_l)*fFS*(k*R*T*T/(L*L))*driving;
    }
}

// 相変化ソースベクトル S=(S_Q0,S_Q1,S_Q2,S_g)。モーメントは保存量 (ρQn)。
//   p_v: 凝縮種の蒸気分圧 (pure=P, carrier=ρ(Y-g)R_w T)。rho_v: 蒸気密度。cprops で種を切替。
__host__ __device__ inline void cond_source_vector(
    const CondSpeciesProps& cp, double T, double p_v, double rho_v,
    double roQ0, double roQ1, double roQ2,
    double* SQ0, double* SQ1, double* SQ2, double* Sg,
    int kantrowitz = 0, int growthModel = 0, double gamma_gas = 1.4, double p_gas = -1.0,
    double gyarC = 3.18)
{
    if (rho_v < 0.0) rho_v = 0.0;
    double J, rstar;
    cond_nucleation(cp, T, p_v, rho_v, &J, &rstar, kantrowitz, gamma_gas);
    const double r_bar = (roQ0 > 1.0e-30) ? (roQ1/roQ0) : rstar;
    double drdt = (roQ0 > 1.0e-30) ? cond_growth(cp, T, p_v, r_bar, rstar, growthModel, p_gas, gyarC) : 0.0;
    const double rho_l = cond_rho_cond(cp, T);
    *SQ0 = J;
    *SQ1 = J*rstar + roQ0*drdt;
    *SQ2 = J*rstar*rstar + 2.0*roQ1*drdt;
    *Sg  = (4.0/3.0)*COND_PI*rho_l*(J*rstar*rstar*rstar + 3.0*roQ2*drdt);
}
