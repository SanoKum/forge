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
