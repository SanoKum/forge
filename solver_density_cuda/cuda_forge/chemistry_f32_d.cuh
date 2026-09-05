#pragma once
// =============================================================================
// fp32 版の反応速度・Jacobian 評価 (CMC の条件付き空間 Q(η) 専用)。methods/chemistry_cmc.md 実装 §3。
//   RTX 3060 級 (fp64 = fp32 の 1/64) では chem_source (double) が CMC 化学カーネルの FP64 パイプを 87% 占め律速だった
//   (ncu 2026-09-05)。式は chemistry_d.cuh の chem_source と同じで、精度だけ float (速度定数は対数空間で評価するので
//   範囲は足りる; ω の平衡近傍の桁落ちは 1e-7·max(q_f,q_b) で点陰解が減衰)。点陰解の 9×9 Gauss は double のまま
//   (剛性で条件数が大きい; コストは無視できる)。平均場の化学 (chemistry_d.cu) は double のまま変えない。
//   係数表は float 構造体を __constant__ に置く (全 thread が同じ反応 r を読むので broadcast が効き、fp64→fp32 変換も不要)。
// =============================================================================
#include "cuda_forge/thermo_d.cuh"
#include "cuda_forge/chemistry_d.cuh"

struct SpeciesThermoF {
    float MW, Tlo, Tmid, Thi, low[9], high[9], h_datum;
};
struct ReactionTableF {
    int nReac, nSpecies;
    float lnA[CHEM_MAX_REACTIONS], lnA0[CHEM_MAX_REACTIONS], b[CHEM_MAX_REACTIONS], Ea[CHEM_MAX_REACTIONS];
    int nR[CHEM_MAX_REACTIONS], nP[CHEM_MAX_REACTIONS];
    int rIdx[CHEM_MAX_REACTIONS][CHEM_MAX_SIDE], pIdx[CHEM_MAX_REACTIONS][CHEM_MAX_SIDE];
    int rNu[CHEM_MAX_REACTIONS][CHEM_MAX_SIDE], pNu[CHEM_MAX_REACTIONS][CHEM_MAX_SIDE];
    int thirdBody[CHEM_MAX_REACTIONS], reversible[CHEM_MAX_REACTIONS], falloff[CHEM_MAX_REACTIONS];
    float dnu[CHEM_MAX_REACTIONS];
    float eff[CHEM_MAX_REACTIONS][THERMO_MAX_SPECIES];
    float b0[CHEM_MAX_REACTIONS], Ea0[CHEM_MAX_REACTIONS];
    float troeA[CHEM_MAX_REACTIONS], troeT3[CHEM_MAX_REACTIONS], troeT1[CHEM_MAX_REACTIONS], troeT2[CHEM_MAX_REACTIONS];
};

inline void chemf_build_tables(const SpeciesThermo* sp, int n, const ReactionTable* rt, SpeciesThermoF* spf, ReactionTableF* rtf)
{
    for (int s = 0; s < n; ++s) {
        spf[s].MW = (float)sp[s].MW; spf[s].Tlo = (float)sp[s].Tlo; spf[s].Tmid = (float)sp[s].Tmid; spf[s].Thi = (float)sp[s].Thi;
        for (int i = 0; i < 9; ++i) { spf[s].low[i] = (float)sp[s].low[i]; spf[s].high[i] = (float)sp[s].high[i]; }
        spf[s].h_datum = (float)sp[s].h_datum;
    }
    std::memset(rtf, 0, sizeof(ReactionTableF));
    rtf->nReac = rt->nReac; rtf->nSpecies = rt->nSpecies;
    for (int r = 0; r < rt->nReac; ++r) {
        rtf->lnA[r] = (float)rt->lnA[r]; rtf->lnA0[r] = (float)rt->lnA0[r]; rtf->b[r] = (float)rt->b[r]; rtf->Ea[r] = (float)rt->Ea[r];
        rtf->nR[r] = rt->nR[r]; rtf->nP[r] = rt->nP[r];
        for (int j = 0; j < CHEM_MAX_SIDE; ++j) { rtf->rIdx[r][j] = rt->rIdx[r][j]; rtf->pIdx[r][j] = rt->pIdx[r][j]; rtf->rNu[r][j] = rt->rNu[r][j]; rtf->pNu[r][j] = rt->pNu[r][j]; }
        rtf->thirdBody[r] = rt->thirdBody[r]; rtf->reversible[r] = rt->reversible[r]; rtf->falloff[r] = rt->falloff[r];
        rtf->dnu[r] = (float)rt->dnu[r];
        for (int s = 0; s < THERMO_MAX_SPECIES; ++s) rtf->eff[r][s] = (float)rt->eff[r][s];
        rtf->b0[r] = (float)rt->b0[r]; rtf->Ea0[r] = (float)rt->Ea0[r];
        rtf->troeA[r] = (float)rt->troeA[r]; rtf->troeT3[r] = (float)rt->troeT3[r]; rtf->troeT1[r] = (float)rt->troeT1[r]; rtf->troeT2[r] = (float)rt->troeT2[r];
    }
}

#define CHEMF_RU 8.314462618f

__device__ __forceinline__ const float* chemf_pick(const SpeciesThermoF& sp, float Tc) { return (Tc < sp.Tmid) ? sp.low : sp.high; }

// cp, h (絶対基準 = 係数どおり), s°/R を 1 スイープ。範囲外はクランプ + 線形外挿 (thermo_d.cuh と同じ規約)。
__device__ __forceinline__ void chemf_cph_molar(const SpeciesThermoF& sp, float T, float Ti, float lnT, float* cp, float* h)
{
    float Tc = T, Tci = Ti, lnTc = lnT; bool ext = false;
    if (T < sp.Tlo) { Tc = sp.Tlo; Tci = 1.0f / Tc; lnTc = logf(Tc); ext = true; }
    else if (T > sp.Thi) { Tc = sp.Thi; Tci = 1.0f / Tc; lnTc = logf(Tc); ext = true; }
    const float* a = chemf_pick(sp, Tc);
    const float Ti2 = Tci * Tci, T2 = Tc * Tc, T3 = T2 * Tc, T4 = T3 * Tc;
    const float cpc = CHEMF_RU * (a[0]*Ti2 + a[1]*Tci + a[2] + a[3]*Tc + a[4]*T2 + a[5]*T3 + a[6]*T4);
    const float hc  = CHEMF_RU * Tc * (-a[0]*Ti2 + a[1]*lnTc*Tci + a[2] + a[3]*Tc*0.5f + a[4]*T2*(1.0f/3.0f) + a[5]*T3*0.25f + a[6]*T4*0.2f + a[7]*Tci);
    *cp = cpc; *h = ext ? hc + cpc * (T - Tc) : hc;
}
__device__ __forceinline__ float chemf_sR(const SpeciesThermoF& sp, float T, float Ti, float lnT)
{
    float Tc = T, Tci = Ti, lnTc = lnT; bool ext = false;
    if (T < sp.Tlo) { Tc = sp.Tlo; Tci = 1.0f / Tc; lnTc = logf(Tc); ext = true; }
    else if (T > sp.Thi) { Tc = sp.Thi; Tci = 1.0f / Tc; lnTc = logf(Tc); ext = true; }
    const float* a = chemf_pick(sp, Tc);
    const float Ti2 = Tci * Tci, T2 = Tc * Tc, T3 = T2 * Tc, T4 = T3 * Tc;
    float sR = -a[0]*Ti2*0.5f - a[1]*Tci + a[2]*lnTc + a[3]*Tc + a[4]*T2*0.5f + a[5]*T3*(1.0f/3.0f) + a[6]*T4*0.25f + a[8];
    if (ext) { const float cpR = a[0]*Ti2 + a[1]*Tci + a[2] + a[3]*Tc + a[4]*T2 + a[5]*T3 + a[6]*T4; sR += cpR * (lnT - lnTc); }   // 定 cp 外挿
    return sR;
}

__device__ __forceinline__ float chemf_R_mix(const SpeciesThermoF* sp, int n, const float* Y)
{ float s = 0.0f; for (int i = 0; i < n; ++i) s += Y[i] / sp[i].MW; return CHEMF_RU * s; }

// T(h, Y) Newton 反転 (thermo_T_from_h と同じ反復・収束判定)
__device__ __forceinline__ float chemf_T_from_h(const SpeciesThermoF* sp, int n, const float* Y, float h, float T_guess, float T_min, float T_max)
{
    float T = T_guess;
    if (!(T > T_min)) T = T_min;
    if (T > T_max) T = T_max;
    #pragma unroll 1
    for (int it = 0; it < 20; ++it) {
        const float Ti = 1.0f / T, lnT = logf(T);
        float cp = 0.0f, h_T = 0.0f;
        for (int i = 0; i < n; ++i) { float cpi, hi; chemf_cph_molar(sp[i], T, Ti, lnT, &cpi, &hi); cp += Y[i] * (cpi / sp[i].MW); h_T += Y[i] * (hi / sp[i].MW); }
        const float cpf = (cp > 1.0e-3f ? cp : 1.0e-3f);
        float dT = (h_T - h) / cpf;
        if (dT >  0.5f*T) dT =  0.5f*T;
        if (dT < -0.5f*T) dT = -0.5f*T;
        T -= dT;
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        if (dT < 0.0f) dT = -dT;
        if (dT < 1.0e-3f + 1.0e-6f*T) break;
    }
    return T;
}

__device__ __forceinline__ float chemf_ln_kf_falloff(const ReactionTableF* rt, int r, float T, float lnT, float RuT, float M)
{
    const float lnkinf = rt->lnA[r]  + rt->b[r]*lnT  - rt->Ea[r]/RuT;
    const float lnk0   = rt->lnA0[r] + rt->b0[r]*lnT - rt->Ea0[r]/RuT;
    const float Mp = (M > 1.0e-30f) ? M : 1.0e-30f;
    const float lnPr = lnk0 + logf(Mp) - lnkinf;
    const float Pr = expf(lnPr);
    float lnF = 0.0f;
    if (rt->falloff[r] == 2) {
        float Fc = (1.0f - rt->troeA[r]) * expf(-T / rt->troeT3[r]) + rt->troeA[r] * expf(-T / rt->troeT1[r]);
        if (rt->troeT2[r] > 0.0f) Fc += expf(-rt->troeT2[r] / T);
        if (Fc < 1.0e-30f) Fc = 1.0e-30f;
        const float lFc = log10f(Fc);
        const float c = -0.4f - 0.67f*lFc, nn = 0.75f - 1.27f*lFc, d = 0.14f;
        const float lPr = lnPr * 0.4342944819f;
        const float x = (lPr + c) / (nn - d*(lPr + c));
        lnF = 2.302585093f * lFc / (1.0f + x*x);
    }
    return lnkinf + logf(Pr / (1.0f + Pr)) + lnF;
}

// 生成率 ω_s [kg/m3/s], 反応熱 Qdot [W/m3], 全 n×n Jacobian J[s*n+k] = ∂ω_s/∂(ρY_k) at fixed T [1/s] (chem_source jacMode 2 と同じ式)
__device__ void chemf_source(const SpeciesThermoF* sp, const ReactionTableF* rt, float ro, float T, const float* Y, float* omega, float* Qdot, float* J)
{
    const int n = rt->nSpecies;
    float X[THERMO_MAX_SPECIES], g[THERMO_MAX_SPECIES];
    const float RuT = CHEMF_RU * T, invT = 1.0f / T, lnT = logf(T);
    for (int s = 0; s < n; ++s) { const float y = (Y[s] > 0.0f) ? Y[s] : 0.0f; X[s] = ro * y / sp[s].MW; omega[s] = 0.0f; }
    for (int i = 0; i < n*n; ++i) J[i] = 0.0f;
    for (int s = 0; s < n; ++s) { float cp, h; chemf_cph_molar(sp[s], T, invT, lnT, &cp, &h); g[s] = chemf_sR(sp[s], T, invT, lnT) - (h + sp[s].h_datum) / RuT; }
    const float lncstd = logf((float)CHEM_P_STD / RuT);
    for (int r = 0; r < rt->nReac; ++r) {
        float M = 1.0f;
        if (rt->thirdBody[r] || rt->falloff[r]) { M = 0.0f; for (int s = 0; s < n; ++s) M += rt->eff[r][s] * X[s]; }
        float lnkf, dkf_dM_over_kf = 0.0f;
        if (rt->falloff[r]) {
            lnkf = chemf_ln_kf_falloff(rt, r, T, lnT, RuT, M);
            const float dM = 1.0e-3f * ((M > 1.0e-30f) ? M : 1.0e-30f);   // float なので相対刻みは 1e-3
            const float mp = chemf_ln_kf_falloff(rt, r, T, lnT, RuT, M + dM);
            const float mm = chemf_ln_kf_falloff(rt, r, T, lnT, RuT, (M - dM > 0.0f) ? M - dM : 0.0f);
            dkf_dM_over_kf = (mp - mm) / (((M - dM > 0.0f) ? 2.0f*dM : M + dM));
        } else {
            lnkf = rt->lnA[r] + rt->b[r]*lnT - rt->Ea[r]/RuT;
        }
        const float kf = expf(lnkf);
        float kb = 0.0f;
        if (rt->reversible[r]) {
            float sumg = 0.0f;
            for (int j = 0; j < rt->nP[r]; ++j) sumg += rt->pNu[r][j]*g[rt->pIdx[r][j]];
            for (int j = 0; j < rt->nR[r]; ++j) sumg -= rt->rNu[r][j]*g[rt->rIdx[r][j]];
            kb = expf(lnkf - (sumg + rt->dnu[r]*lncstd));
        }
        float Pf = 1.0f, Pb = 1.0f, dPf[CHEM_MAX_SIDE], dPb[CHEM_MAX_SIDE];
        for (int j = 0; j < rt->nR[r]; ++j) { const float x = X[rt->rIdx[r][j]]; const int nu = rt->rNu[r][j]; float xp = 1.0f; for (int k = 0; k < nu; ++k) xp *= x; Pf *= xp; }
        for (int j = 0; j < rt->nP[r]; ++j) { const float x = X[rt->pIdx[r][j]]; const int nu = rt->pNu[r][j]; float xp = 1.0f; for (int k = 0; k < nu; ++k) xp *= x; Pb *= xp; }
        for (int j = 0; j < rt->nR[r]; ++j) {
            float prod = 1.0f;
            for (int i = 0; i < rt->nR[r]; ++i) { const float x = X[rt->rIdx[r][i]]; int nu = rt->rNu[r][i]; if (i == j) { nu -= 1; prod *= (float)rt->rNu[r][i]; } for (int k = 0; k < nu; ++k) prod *= x; }
            dPf[j] = prod;
        }
        for (int j = 0; j < rt->nP[r]; ++j) {
            float prod = 1.0f;
            for (int i = 0; i < rt->nP[r]; ++i) { const float x = X[rt->pIdx[r][i]]; int nu = rt->pNu[r][i]; if (i == j) { nu -= 1; prod *= (float)rt->pNu[r][i]; } for (int k = 0; k < nu; ++k) prod *= x; }
            dPb[j] = prod;
        }
        const float Mq = rt->thirdBody[r] ? M : 1.0f;
        const float qf = kf*Pf, qb = kb*Pb, q = Mq * (qf - qb);
        for (int j = 0; j < rt->nR[r]; ++j) { const int s = rt->rIdx[r][j]; omega[s] -= sp[s].MW * rt->rNu[r][j] * q; }
        for (int j = 0; j < rt->nP[r]; ++j) { const int s = rt->pIdx[r][j]; omega[s] += sp[s].MW * rt->pNu[r][j] * q; }
        auto accumulate = [&](int k, float dq_dXk) {
            const float invWk = 1.0f / sp[k].MW;
            for (int j = 0; j < rt->nR[r]; ++j) { const int s = rt->rIdx[r][j]; J[s*n+k] -= sp[s].MW * rt->rNu[r][j] * dq_dXk * invWk; }
            for (int j = 0; j < rt->nP[r]; ++j) { const int s = rt->pIdx[r][j]; J[s*n+k] += sp[s].MW * rt->pNu[r][j] * dq_dXk * invWk; }
        };
        for (int j = 0; j < rt->nR[r]; ++j) accumulate(rt->rIdx[r][j],  Mq * kf * dPf[j]);
        for (int j = 0; j < rt->nP[r]; ++j) accumulate(rt->pIdx[r][j], -Mq * kb * dPb[j]);
        if (rt->thirdBody[r]) { for (int k = 0; k < n; ++k) if (rt->eff[r][k] != 0.0f) accumulate(k, rt->eff[r][k] * (qf - qb)); }
        else if (rt->falloff[r]) { for (int k = 0; k < n; ++k) if (rt->eff[r][k] != 0.0f) accumulate(k, rt->eff[r][k] * dkf_dM_over_kf * (qf - qb)); }
    }
    float Q = 0.0f;
    for (int s = 0; s < n; ++s) Q -= (sp[s].h_datum / sp[s].MW) * omega[s];
    *Qdot = Q;
}
