#pragma once

#include "thermo_d.cuh"
#include "condensationProperties_d.cuh"

// 一温度 二相 温度反転 (carrier + condensible, thermally perfect; Phase 3 H2O)。
// 気相 = 全化学種 (蒸気とみなした NASA-9 混合)。凝縮種 (H2O) の液相分率 g。
//   混合内部エネルギー e_mix(T) = e_gas^全蒸気(Y,T) + g(R_w T - L_w(T))
//     (液 e_l = e_v + R_w T - L_w、凝縮ぶん g の補正)。
// e_mix(T)=e_in を Newton で解く。Y は総質量分率 (N2 + 総水)。g≤Y_H2O は呼び出し側で保証。
//   Rw = R_u/M_H2O (凝縮種の比気体定数)、cprops = 凝縮種物性 (model=H2O)。
__host__ __device__ inline double cond_T_from_e_carrier(
    const SpeciesThermo* sp, int nSp, const double* Y,
    double e_in, double g, double Rw, const CondSpeciesProps& cprops,
    double T_guess, double T_min, double T_max)
{
    const double Rmix = thermo_R_mix(sp, nSp, Y);
    double T = T_guess;
    if (!(T > T_min)) T = T_min;
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it = 0; it < 30; ++it) {
        double cp_T, h_T;
        thermo_cph_mix(sp, nSp, Y, T, &cp_T, &h_T);          // 全蒸気混合 cp,h
        const double e_allvap = h_T - Rmix*T;                 // 全蒸気の内部エネルギー
        const double cv = cp_T - Rmix;
        const double L  = cond_latent(cprops, T);
        const double dL = (cond_latent(cprops, T+0.1) - cond_latent(cprops, T-0.1)) / 0.2;
        const double e_mix = e_allvap + g*(Rw*T - L);
        const double demix = cv + g*(Rw - dL);
        const double demf  = (demix > 1.0e-2*cv) ? demix : 1.0e-2*cv;
        double dT = (e_mix - e_in)/demf;
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

// 非平衡凝縮: 一温度 二相 EOS の温度反転 (Phase 2)。methods/condensation.md 5 節。
//
// 【一温度近似】気相温度 T_v と液滴温度 T_d を分けず T_v=T_d=T とする (初期実装)。
// 将来は e=(1-g)e_v(T_v)+g e_l(T_d) と Hill T_d 式の 2 変数局所 Newton へ拡張する設計
// (本関数は T 1 変数 Newton; 拡張時は (T_v,T_d) の 2x2 Newton に置換)。
//
// 【混合内部エネルギー】論文形 e = (1-g)e_v(T) + g e_l(T)、液相 e_l(T)=e_v(T)-L(T) より
//   e = e_v(T) - g L(T)
// となる。気相 e_v(T) は forge の thermally-perfect thermo (NASA-9, thermo_cph_mix) で評価し
// **定比熱近似は使わない**。g は液相質量分率 (総和)。L は凝縮種の潜熱 (condensationProperties)。
//
// 【温度反転】保存量から e_in = roe/ρ - ½|u|² を作り、
//   G(T) = e_v(T) - g L(T) - e_in = 0,   G'(T) = c_v(T) - g L'(T)
// を Newton で解く。g=0 のとき G=e_v(T)-e_in で thermo_T_from_e と一致 → 単相 TP に厳密縮約
// (呼び出し側は g≈0 で従来 thermo_T_from_e を使い厳密同一を保証する)。
//
// 【圧力】実在気体補正なし。液滴は圧力を持たないとして p = (1-g)ρ R_v T。
//
// 現状 N2 のみ (潜熱 n2_latent)。多成分凝縮では Σ_sp g_sp L_sp(T) に一般化する (将来)。

// 一温度 二相 温度反転 (calorically perfect gas; 初期実装)。論文 Eq.16-18 (g0=0):
//   液相 e_l = e_v + R_v T - L(T) = c_p T - L  (∵ h_v-h_l=L, h_v=e_v+R_vT, h_l≈e_l)
//   e = (1-g)e_v + g e_l = (c_v + g R_v) T - g L(T)   (e_v=c_v T)
//   ⇒ T = (e_in + g L(T))/(c_v + g R_v)  (= Eq.18, Cv0=Cvv=c_v)
// **L の温度依存は入れる** (n2_latent(T))。g=0 で T=e_in/cv (従来 CPG と一致)。cv=cp/γ, R=(γ-1)cv。
__host__ __device__ inline double cond_T_from_e_cpg(
    double e_in, double g_tot, double cv, double R, double T_guess,
    const CondSpeciesProps& cprops)
{
    const double a = cv + g_tot*R;   // 実効熱容量 (T に対し一定)
    double T = T_guess;
    if (!(T > 1.0)) T = 1.0;
    #pragma unroll 1
    for (int it = 0; it < 30; ++it) {
        const double L  = cond_latent(cprops, T);
        const double dL = (cond_latent(cprops, T + 0.1) - cond_latent(cprops, T - 0.1)) / 0.2;
        const double G  = a*T - g_tot*L - e_in;
        const double Gp = a - g_tot*dL;
        const double Gpf = (Gp > 1.0e-2*a) ? Gp : 1.0e-2*a;
        double dT = G / Gpf;
        if (dT >  0.5*T) dT =  0.5*T;
        if (dT < -0.5*T) dT = -0.5*T;
        T -= dT;
        if (T < 1.0)    T = 1.0;
        if (T > 6000.0) T = 6000.0;
        if (dT < 0.0) dT = -dT;
        if (dT < 1.0e-3 + 1.0e-6*T) break;
    }
    return T;
}

// 一温度 二相 温度反転 (thermally perfect; 後続拡張)。e_in [J/kg], g_tot=総液相質量分率, T_guess。
__host__ __device__ inline double cond_T_from_e_onetemp(
    const SpeciesThermo* sp, int nSp, const double* Y,
    double e_in, double g_tot,
    double T_guess, double T_min, double T_max)
{
    const double R = thermo_R_mix(sp, nSp, Y);
    double T = T_guess;
    if (!(T > T_min)) T = T_min;
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it = 0; it < 30; ++it) {
        double cp_T, h_T;
        thermo_cph_mix(sp, nSp, Y, T, &cp_T, &h_T);   // 気相 cp,h (TP)
        const double e_v = h_T - R*T;                  // 気相内部エネルギー
        const double cv  = cp_T - R;                   // de_v/dT
        const double L   = n2_latent(T);               // 潜熱 (N2)
        const double dL  = (n2_latent(T + 0.1) - n2_latent(T - 0.1)) / 0.2; // dL/dT (数値)
        // e = (1-g)e_v + g e_l, e_l = e_v + R T - L ⇒ e = e_v + g R T - g L
        const double G   = e_v + g_tot*R*T - g_tot*L - e_in;
        const double Gp  = cv + g_tot*R - g_tot*dL;
        const double Gpf = (Gp > 1.0e-2*R) ? Gp : 1.0e-2*R; // 0 割回避フロア
        double dT = G / Gpf;
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
