#pragma once

#include "thermo_d.cuh"
#include "condensationProperties_d.cuh"

// 非平衡凝縮: 一温度 二相 EOS の温度反転 (Phase 2)。docs/condensation/theory.md 5 節。
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

// 一温度 二相 温度反転 (calorically perfect gas; 初期実装)。
// 気相 e_v(T)=c_v T (cv 一定)、液相 e_l=e_v-L(T) より e = c_v T - g L(T) = e_in を Newton で解く。
// **L の温度依存は入れる** (n2_latent(T))。論文は完全理想気体 (cp 一定) なのでこの形が基本。
// g=0 で T=e_in/cv (従来 CPG と一致)。cv = cp/γ。
__host__ __device__ inline double cond_T_from_e_cpg(
    double e_in, double g_tot, double cv, double T_guess)
{
    double T = T_guess;
    if (!(T > 1.0)) T = 1.0;
    #pragma unroll 1
    for (int it = 0; it < 30; ++it) {
        const double L  = n2_latent(T);
        const double dL = (n2_latent(T + 0.1) - n2_latent(T - 0.1)) / 0.2;
        const double G  = cv*T - g_tot*L - e_in;
        const double Gp = cv - g_tot*dL;
        const double Gpf = (Gp > 1.0e-2*cv) ? Gp : 1.0e-2*cv;
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
        const double G   = e_v - g_tot*L - e_in;
        const double Gp  = cv - g_tot*dL;
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
