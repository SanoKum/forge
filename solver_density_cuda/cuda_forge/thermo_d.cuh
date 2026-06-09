#pragma once

// =============================================================================
// thermo_d.cuh
//   多成分 thermally-perfect gas の熱力学コア。
//   - 比熱 cp(T) は NASA-9 係数多項式 (CEA 準拠, McBride-Gordon 2002, 2 温度域)。
//   - 混合則は ideal-gas mixing (質量重み cp/h、モル重み M)。
//   - エネルギーは NASA 絶対エンタルピー基準 (生成エンタルピー込み)。
//
//   設計メモ:
//   - flow_float は float (FP32)。NASA-9 多項式和・Newton 反転・混合和は
//     FP32 では破綻するため、ここでの内部計算は全て double で行う。
//   - 同じ関数を __host__ __device__ で提供し、初期条件 (host) と
//     カーネル (device) で同一の熱力学を保証する。
//   - 化学種データ (SpeciesThermo 配列) は device global memory に置き、
//     ポインタをカーネル引数で渡す (translation unit をまたぐ __constant__ は
//     -rdc を要するため避ける)。アクセスは全スレッド同係数読みで L2/読取専用
//     キャッシュに乗る。
// =============================================================================

#include "flowFormat.hpp"

#ifndef __CUDACC__
  #include <cmath>
  #define THERMO_HD inline
#else
  #define THERMO_HD __host__ __device__ inline
#endif

// 普遍気体定数 [J/(mol·K)]
#define THERMO_RU 8.314462618
// コンパイル時の最大化学種数 (DB 容量の上限)
#define THERMO_MAX_SPECIES 16

// 1 化学種の NASA-9 係数 (2 温度域) と Lennard-Jones パラメータ。
// 係数配列 a[0..8] の規約:
//   cp/R = a0/T^2 + a1/T + a2 + a3 T + a4 T^2 + a5 T^3 + a6 T^4
//   h/RT = -a0/T^2 + a1 ln(T)/T + a2 + a3 T/2 + a4 T^2/3 + a5 T^3/4 + a6 T^4/5 + a7/T
//   s/R  = -a0/(2T^2) - a1/T + a2 ln(T) + a3 T + a4 T^2/2 + a5 T^3/3 + a6 T^4/4 + a8
struct SpeciesThermo {
    double MW;        // 分子量 [kg/mol]
    double sigma_LJ;  // Lennard-Jones 衝突直径 [Angstrom] (kinetic theory 輸送用)
    double eps_kB;    // Lennard-Jones ポテンシャル深さ ε/kB [K]
    double Tlo;       // 低温端 [K]
    double Tmid;      // 温度域境界 [K]
    double Thi;       // 高温端 [K]
    double low[9];    // Tlo  <= T <  Tmid の係数
    double high[9];   // Tmid <= T <= Thi  の係数
};

// -----------------------------------------------------------------------------
// 単一化学種の NASA-9 評価 (全て double)。範囲外は端でクランプし、
// エンタルピーは線形外挿 (h(T) ~ h(Tc) + cp(Tc)(T-Tc)) して衝撃波での暴走を防ぐ。
// -----------------------------------------------------------------------------

// 係数選択 (クランプ済み温度 Tc に対応する係数配列を返す)
THERMO_HD const double* thermo_pick_coeffs(const SpeciesThermo& sp, double Tc)
{
    return (Tc < sp.Tmid) ? sp.low : sp.high;
}

// cp_molar [J/(mol·K)] (クランプ温度で評価)
THERMO_HD double thermo_cp_molar_clamped(const SpeciesThermo& sp, double Tc)
{
    const double* a = thermo_pick_coeffs(sp, Tc);
    const double Ti  = 1.0/Tc;
    const double Ti2 = Ti*Ti;
    return THERMO_RU * ( a[0]*Ti2 + a[1]*Ti + a[2]
                       + a[3]*Tc + a[4]*Tc*Tc + a[5]*Tc*Tc*Tc + a[6]*Tc*Tc*Tc*Tc );
}

// h_molar [J/mol] (クランプ温度で評価, 絶対エンタルピー)
THERMO_HD double thermo_h_molar_clamped(const SpeciesThermo& sp, double Tc)
{
    const double* a = thermo_pick_coeffs(sp, Tc);
    const double Ti  = 1.0/Tc;
    const double Ti2 = Ti*Ti;
    const double lnT = log(Tc);
    // h/(R T) = -a0/T^2 + a1 ln(T)/T + a2 + a3 T/2 + a4 T^2/3 + a5 T^3/4 + a6 T^4/5 + a7/T
    const double hRT = -a[0]*Ti2 + a[1]*lnT*Ti + a[2]
                     + a[3]*Tc/2.0 + a[4]*Tc*Tc/3.0 + a[5]*Tc*Tc*Tc/4.0
                     + a[6]*Tc*Tc*Tc*Tc/5.0 + a[7]*Ti;
    return THERMO_RU * Tc * hRT;
}

// 範囲クランプ + 外挿込みの cp_molar [J/(mol·K)]
THERMO_HD double thermo_cp_molar(const SpeciesThermo& sp, double T)
{
    double Tc = T;
    if (Tc < sp.Tlo) Tc = sp.Tlo;
    if (Tc > sp.Thi) Tc = sp.Thi;
    return thermo_cp_molar_clamped(sp, Tc);
}

// 範囲クランプ + 線形外挿込みの h_molar [J/mol]
THERMO_HD double thermo_h_molar(const SpeciesThermo& sp, double T)
{
    if (T < sp.Tlo) {
        const double h0  = thermo_h_molar_clamped(sp, sp.Tlo);
        const double cp0 = thermo_cp_molar_clamped(sp, sp.Tlo);
        return h0 + cp0*(T - sp.Tlo);
    }
    if (T > sp.Thi) {
        const double h1  = thermo_h_molar_clamped(sp, sp.Thi);
        const double cp1 = thermo_cp_molar_clamped(sp, sp.Thi);
        return h1 + cp1*(T - sp.Thi);
    }
    return thermo_h_molar_clamped(sp, T);
}

// 単一化学種の比気体定数 R_s = Ru/MW [J/(kg·K)]
THERMO_HD double thermo_R_species(const SpeciesThermo& sp)
{
    return THERMO_RU / sp.MW;
}

// 単一化学種 質量基準: cp_s [J/(kg·K)], h_s [J/kg]
THERMO_HD double thermo_cp_mass(const SpeciesThermo& sp, double T)
{ return thermo_cp_molar(sp, T) / sp.MW; }

THERMO_HD double thermo_h_mass(const SpeciesThermo& sp, double T)
{ return thermo_h_molar(sp, T) / sp.MW; }

// -----------------------------------------------------------------------------
// 混合則 (質量分率 Y[0..n-1] を与える, ideal-gas mixing)
// -----------------------------------------------------------------------------

// 混合分子量の逆数 1/M_mix = Σ Y_s/MW_s [mol/kg]
THERMO_HD double thermo_inv_Mmix(const SpeciesThermo* sp, int n, const double* Y)
{
    double s = 0.0;
    for (int i=0;i<n;i++) s += Y[i]/sp[i].MW;
    return s;
}

// 混合比気体定数 R_mix = Ru·Σ Y_s/MW_s [J/(kg·K)]
THERMO_HD double thermo_R_mix(const SpeciesThermo* sp, int n, const double* Y)
{
    return THERMO_RU * thermo_inv_Mmix(sp, n, Y);
}

// 混合 cp [J/(kg·K)] (質量重み)
THERMO_HD double thermo_cp_mix(const SpeciesThermo* sp, int n, const double* Y, double T)
{
    double cp = 0.0;
    for (int i=0;i<n;i++) cp += Y[i]*thermo_cp_mass(sp[i], T);
    return cp;
}

// 混合絶対エンタルピー [J/kg] (質量重み)
THERMO_HD double thermo_h_mix(const SpeciesThermo* sp, int n, const double* Y, double T)
{
    double h = 0.0;
    for (int i=0;i<n;i++) h += Y[i]*thermo_h_mass(sp[i], T);
    return h;
}

// 混合内部エネルギー e_mix(T) = h_mix(T) - R_mix·T [J/kg]
THERMO_HD double thermo_e_mix(const SpeciesThermo* sp, int n, const double* Y, double T)
{
    return thermo_h_mix(sp, n, Y, T) - thermo_R_mix(sp, n, Y)*T;
}

// 混合比熱比 gamma_mix = cp/(cp-R)
THERMO_HD double thermo_gamma_mix(const SpeciesThermo* sp, int n, const double* Y, double T)
{
    const double cp = thermo_cp_mix(sp, n, Y, T);
    const double R  = thermo_R_mix(sp, n, Y);
    const double cv = cp - R;
    return cp / (cv > 1.0e-6 ? cv : 1.0e-6);
}

// -----------------------------------------------------------------------------
// 輸送係数 (kinetic theory, M3)。Chapman-Enskog 単成分 + Wilke/Mason-Saxena 混合。
//   - 衝突積分 Ω(2,2)*, Ω(1,1)* は Neufeld et al. (1972) の閉形式近似 (0.3<=T*<=100)。
//   - 単位は CGS 換算定数で Pa·s / W·m⁻¹·K⁻¹ / m²·s⁻¹ に整える。内部 double。
// -----------------------------------------------------------------------------

// 換算衝突積分 Ω*(2,2) (粘性・熱伝導)。T* = T/(ε/kB)。
THERMO_HD double thermo_omega22(double Tstar)
{
    if (Tstar < 0.3)   Tstar = 0.3;
    if (Tstar > 100.0) Tstar = 100.0;
    return 1.16145*pow(Tstar, -0.14874)
         + 0.52487*exp(-0.77320*Tstar)
         + 2.16178*exp(-2.43787*Tstar);
}

// 換算衝突積分 Ω*(1,1) (拡散)。M4 の二元拡散係数で使用。
THERMO_HD double thermo_omega11(double Tstar)
{
    if (Tstar < 0.3)   Tstar = 0.3;
    if (Tstar > 100.0) Tstar = 100.0;
    return 1.06036*pow(Tstar, -0.15610)
         + 0.19300*exp(-0.47635*Tstar)
         + 1.03587*exp(-1.52996*Tstar)
         + 1.76474*exp(-3.89411*Tstar);
}

// 単成分 Chapman-Enskog 粘性 [Pa·s]。
//   μ = 2.6693e-6 √(M[g/mol]·T) / (σ[Å]²·Ω(2,2)*)
THERMO_HD double thermo_mu_species(const SpeciesThermo& sp, double T)
{
    const double Tstar = T / (sp.eps_kB > 1.0e-30 ? sp.eps_kB : 1.0e-30);
    const double Mg     = sp.MW * 1000.0;                 // kg/mol -> g/mol
    const double sig2   = sp.sigma_LJ * sp.sigma_LJ;      // Å²
    const double om     = thermo_omega22(Tstar);
    return 2.6693e-6 * sqrt(Mg * T) / (sig2 * om);
}

// 単成分 熱伝導 [W/(m·K)] (modified Eucken: λ = μ(cp_mass + 1.25 R_s))。
THERMO_HD double thermo_lambda_species(const SpeciesThermo& sp, double T)
{
    const double mu = thermo_mu_species(sp, T);
    const double cp = thermo_cp_mass(sp, T);
    const double R  = thermo_R_species(sp);
    return mu * (cp + 1.25*R);
}

// 質量分率 Y -> モル分率 X (X[i]=(Y_i/MW_i)/Σ_j Y_j/MW_j)。
THERMO_HD void thermo_X_from_Y(const SpeciesThermo* sp, int n, const double* Y, double* X)
{
    double s = 0.0;
    for (int i=0;i<n;i++) { X[i] = Y[i]/sp[i].MW; s += X[i]; }
    const double inv = 1.0/(s > 1.0e-300 ? s : 1.0e-300);
    for (int i=0;i<n;i++) X[i] *= inv;
}

// Wilke の φ_ij = [1+√(μ_i/μ_j)(M_j/M_i)^¼]² / √(8(1+M_i/M_j))。
THERMO_HD double thermo_wilke_phi(double mu_i, double mu_j, double M_i, double M_j)
{
    const double r  = sqrt(mu_i/mu_j) * pow(M_j/M_i, 0.25);
    const double num = (1.0 + r)*(1.0 + r);
    const double den = sqrt(8.0*(1.0 + M_i/M_j));
    return num/den;
}

// Wilke 混合粘性 [Pa·s]。X はモル分率。
THERMO_HD double thermo_mu_mix(const SpeciesThermo* sp, int n, const double* X, double T)
{
    if (n == 1) return thermo_mu_species(sp[0], T);
    double mu_s[THERMO_MAX_SPECIES];
    for (int i=0;i<n;i++) mu_s[i] = thermo_mu_species(sp[i], T);
    double mu = 0.0;
    for (int i=0;i<n;i++) {
        if (X[i] <= 0.0) continue;
        double denom = 0.0;
        for (int j=0;j<n;j++)
            denom += X[j]*thermo_wilke_phi(mu_s[i], mu_s[j], sp[i].MW, sp[j].MW);
        if (denom > 1.0e-300) mu += X[i]*mu_s[i]/denom;
    }
    return mu;
}

// Mason-Saxena (Wassiljewa, Wilke の φ 流用) 混合熱伝導 [W/(m·K)]。
THERMO_HD double thermo_lambda_mix(const SpeciesThermo* sp, int n, const double* X, double T)
{
    if (n == 1) return thermo_lambda_species(sp[0], T);
    double mu_s[THERMO_MAX_SPECIES];
    double la_s[THERMO_MAX_SPECIES];
    for (int i=0;i<n;i++) { mu_s[i] = thermo_mu_species(sp[i], T); la_s[i] = thermo_lambda_species(sp[i], T); }
    double la = 0.0;
    for (int i=0;i<n;i++) {
        if (X[i] <= 0.0) continue;
        double denom = 0.0;
        for (int j=0;j<n;j++)
            denom += X[j]*thermo_wilke_phi(mu_s[i], mu_s[j], sp[i].MW, sp[j].MW);
        if (denom > 1.0e-300) la += X[i]*la_s[i]/denom;
    }
    return la;
}

// 二元拡散係数 D_ij [m²/s] (Chapman-Enskog)。P [Pa]。M4 の混合平均拡散で使用。
//   D_ij = 1.8583e-7 √(T³(1/M_i+1/M_j)[1/(g/mol)]) / (P[atm]·σ_ij²[Å²]·Ω(1,1)*)  [cm²/s] -> m²/s
THERMO_HD double thermo_Dbinary(const SpeciesThermo& a, const SpeciesThermo& b, double T, double P)
{
    const double Mi = a.MW*1000.0, Mj = b.MW*1000.0;       // g/mol
    const double sig = 0.5*(a.sigma_LJ + b.sigma_LJ);       // Å
    const double eps = sqrt(a.eps_kB * b.eps_kB);           // K
    const double Tstar = T / (eps > 1.0e-30 ? eps : 1.0e-30);
    const double om = thermo_omega11(Tstar);
    const double Patm = P / 101325.0;
    const double Dcm2 = 1.8583e-3 * sqrt(T*T*T*(1.0/Mi + 1.0/Mj))
                        / (Patm * sig*sig * om);             // cm²/s
    return Dcm2 * 1.0e-4;                                    // -> m²/s
}

// 混合平均拡散係数 D_i [m²/s] (M4): D_i = (1-X_i)/Σ_{j≠i} X_j/D_ij。
THERMO_HD double thermo_Dmix_species(const SpeciesThermo* sp, int n, const double* X,
                                     int i, double T, double P)
{
    if (n == 1) return 0.0;
    double denom = 0.0;
    for (int j=0;j<n;j++) {
        if (j==i) continue;
        const double Dij = thermo_Dbinary(sp[i], sp[j], T, P);
        denom += X[j]/(Dij > 1.0e-30 ? Dij : 1.0e-30);
    }
    if (denom < 1.0e-300) return thermo_Dbinary(sp[i], sp[i], T, P);
    return (1.0 - X[i])/denom;
}

// -----------------------------------------------------------------------------
// 温度反転: 比内部エネルギー e [J/kg] と組成 Y から T を Newton 反転で求める。
//   f(T) = e_mix(T) - e = 0, f'(T) = cv_mix(T) = cp_mix(T) - R_mix (厳密微分)
//   初期値 T_guess (前ステップ T) でウォームスタート, 毎反復クランプ。
// -----------------------------------------------------------------------------
THERMO_HD double thermo_T_from_e(const SpeciesThermo* sp, int n, const double* Y,
                                 double e, double T_guess,
                                 double T_min, double T_max)
{
    const double R = thermo_R_mix(sp, n, Y);
    double T = T_guess;
    if (!(T > T_min)) T = T_min;      // NaN/未初期化ガード
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it=0; it<20; ++it) {
        const double e_T = thermo_h_mix(sp, n, Y, T) - R*T;   // e_mix(T)
        const double cv  = thermo_cp_mix(sp, n, Y, T) - R;    // de/dT
        const double cvf = (cv > 1.0e-2*R ? cv : 1.0e-2*R);   // 0 割回避フロア
        double dT = (e_T - e)/cvf;
        // 1 反復あたりのステップを制限 (発散防止)
        if (dT >  0.5*T) dT =  0.5*T;
        if (dT < -0.5*T) dT = -0.5*T;
        T -= dT;
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*T) break;
    }
    return T;
}

// 温度反転 (エンタルピー版): 比エンタルピー h [J/kg] と組成 Y から T を求める。
//   f(T) = h_mix(T) - h = 0, f'(T) = cp_mix(T) (厳密微分)。Roe 平均状態の有効 γ 用。
THERMO_HD double thermo_T_from_h(const SpeciesThermo* sp, int n, const double* Y,
                                 double h, double T_guess,
                                 double T_min, double T_max)
{
    double T = T_guess;
    if (!(T > T_min)) T = T_min;
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it=0; it<20; ++it) {
        const double h_T = thermo_h_mix(sp, n, Y, T);
        const double cp  = thermo_cp_mix(sp, n, Y, T);
        const double cpf = (cp > 1.0e-3 ? cp : 1.0e-3);
        double dT = (h_T - h)/cpf;
        if (dT >  0.5*T) dT =  0.5*T;
        if (dT < -0.5*T) dT = -0.5*T;
        T -= dT;
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*T) break;
    }
    return T;
}

// 標準状態エントロピー S^o/R (圧力項を含まない). Roe/等エントロピー BC 用。
THERMO_HD double thermo_s_molar_clamped(const SpeciesThermo& sp, double Tc)
{
    const double* a = thermo_pick_coeffs(sp, Tc);
    const double Ti = 1.0/Tc; const double Ti2 = Ti*Ti; const double lnT = log(Tc);
    const double sR = -a[0]*Ti2/2.0 - a[1]*Ti + a[2]*lnT + a[3]*Tc
                    + a[4]*Tc*Tc/2.0 + a[5]*Tc*Tc*Tc/3.0 + a[6]*Tc*Tc*Tc*Tc/4.0 + a[8];
    return THERMO_RU * sR;
}
THERMO_HD double thermo_s0_mass(const SpeciesThermo& sp, double T)
{
    double Tc = T;
    if (Tc < sp.Tlo) Tc = sp.Tlo;
    if (Tc > sp.Thi) Tc = sp.Thi;
    return thermo_s_molar_clamped(sp, Tc) / sp.MW;   // [J/(kg·K)]
}

// 単成分 (M1) の等エントロピー total→static 反転。
//   全温 Tt・全圧 Pt と局所マッハ M から静的 (Ts, Ps, ρ, |u|) を求める。
//     h(Tt) = h(Ts) + ½u²,  u = M·a(Ts),  a=√(γ(Ts)R Ts)
//     Ps = Pt·exp(-(s0(Tt)-s0(Ts))/R)
THERMO_HD void thermo_isentropic_from_total_single(
    const SpeciesThermo* sp, double Pt, double Tt, double mach,
    double* Ts_out, double* Ps_out, double* ro_out, double* umag_out)
{
    const double R  = thermo_R_species(sp[0]);
    const double h0 = thermo_h_mass(sp[0], Tt);
    double Ts = Tt/(1.0 + 0.2*mach*mach);   // CPG 近似の初期値
    if (Ts < 10.0) Ts = 10.0;
    #pragma unroll 1
    for (int it=0; it<25; ++it) {
        const double cp = thermo_cp_mass(sp[0], Ts);
        const double g  = cp/((cp-R) > 1.0e-6 ? (cp-R) : 1.0e-6);
        const double a2 = g*R*Ts;
        const double f  = thermo_h_mass(sp[0], Ts) + 0.5*mach*mach*a2 - h0;
        const double df = cp + 0.5*mach*mach*g*R;   // 近似微分 (dγ/dT 無視)
        double dT = f/df;
        if (dT >  0.4*Ts) dT =  0.4*Ts;
        if (dT < -0.4*Ts) dT = -0.4*Ts;
        Ts -= dT;
        if (Ts < 10.0) Ts = 10.0;
        if (Ts > Tt)   Ts = Tt;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*Ts) break;
    }
    const double cp = thermo_cp_mass(sp[0], Ts);
    const double g  = cp/((cp-R) > 1.0e-6 ? (cp-R) : 1.0e-6);
    const double a  = sqrt(g*R*Ts);
    *umag_out = mach*a;
    *Ps_out   = Pt*exp(-(thermo_s0_mass(sp[0],Tt) - thermo_s0_mass(sp[0],Ts))/R);
    *Ts_out   = Ts;
    *ro_out   = (*Ps_out)/(R*Ts);
}

// 単成分の等エントロピー total→static 反転 (静圧 Ps 指定版)。
//   全温 Tt・全圧 Pt と静圧 Ps から静温 Ts・密度 ρ・速度 |u| を求める。
//     s0(Ts) = s0(Tt) + R·ln(Ps/Pt)   (等エントロピー: s0(T)-R ln P 一定)
//     |u| = √(2 max(h(Tt)-h(Ts), 0)),  ρ = Ps/(R Ts)
THERMO_HD void thermo_isentropic_from_total_Ps_single(
    const SpeciesThermo* sp, double Pt, double Tt, double Ps,
    double* Ts_out, double* ro_out, double* umag_out)
{
    const double R  = thermo_R_species(sp[0]);
    const double s_target = thermo_s0_mass(sp[0], Tt) + R*log((Ps > 1.0e-30 ? Ps : 1.0e-30)
                                                              /(Pt > 1.0e-30 ? Pt : 1.0e-30));
    double Ts = Tt;   // Ps<=Pt なら Ts<=Tt。Tt から下降側へ Newton。
    if (Ts < 10.0) Ts = 10.0;
    #pragma unroll 1
    for (int it=0; it<30; ++it) {
        const double f  = thermo_s0_mass(sp[0], Ts) - s_target;   // s0(Ts)-target
        const double df = thermo_cp_mass(sp[0], Ts)/Ts;            // ds0/dT = cp/T
        double dT = f/(df > 1.0e-12 ? df : 1.0e-12);
        if (dT >  0.4*Ts) dT =  0.4*Ts;
        if (dT < -0.4*Ts) dT = -0.4*Ts;
        Ts -= dT;
        if (Ts < 10.0) Ts = 10.0;
        if (Ts > Tt)   Ts = Tt;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*Ts) break;
    }
    const double dh = thermo_h_mass(sp[0], Tt) - thermo_h_mass(sp[0], Ts);
    *umag_out = sqrt(dh > 0.0 ? 2.0*dh : 0.0);
    *Ts_out   = Ts;
    *ro_out   = Ps/(R*Ts);
}

// =============================================================================
// device 側化学種データへのアクセス (thermo_d.cu が所有)
//   - thermo_init_db: solverConfig から DB を構築し device へアップロード
//   - thermo_species_device_ptr / thermo_num_species: カーネル引数用に取得
//   - thermo_species_host: host 計算 (初期条件) 用の CPU 側配列
// =============================================================================
class solverConfig;

void                 thermo_init_db(solverConfig& cfg);
const SpeciesThermo* thermo_species_device_ptr();   // device global memory
int                  thermo_num_species();
const SpeciesThermo* thermo_species_host();          // host array (length = thermo_num_species)
