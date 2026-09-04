#pragma once
// =============================================================================
// chemistry_d.cuh — 有限速度化学 (finite-rate chemistry) の反応表と生成率評価
//
//   理論: methods/chemistry.md (Arrhenius・三体・平衡定数 K_c・反応熱・解析 Jacobian)。
//   thermo_d.cuh と同じ THERMO_HD 規約で host/device 共通 (ホスト単体テスト tools/test_chemistry.cpp
//   がビット同一の評価を検証する)。内部計算は全て double (Phase 1。FP32 化は後続)。
//
//   反応表 ReactionTable は固定長配列 (CHEM_MAX_REACTIONS × CHEM_MAX_SIDE) の POD で、
//   host で組んで device global memory へ 1 度コピーし、カーネルはポインタ経由で読む。
// =============================================================================

#include "thermo_d.cuh"

#define CHEM_MAX_REACTIONS 96
#define CHEM_MAX_SIDE 4          // 片側の異なる化学種数の上限 (H2 系は 3 以下)
#define CHEM_P_STD 1.0e5         // NASA-9 標準状態エントロピーの基準圧 [Pa] (1 bar)

struct ReactionTable {
    int nReac;
    int nSpecies;                          // 流れ側 (cfg.speciesNames) の化学種数
    double A [CHEM_MAX_REACTIONS];         // SI: (m^3/mol)^(order-1)/s
    double lnA[CHEM_MAX_REACTIONS];        // log(A) (読込時に計算。カーネルの log を省く)
    double lnA0[CHEM_MAX_REACTIONS];       // log(A0) (falloff 低圧極限)
    double b [CHEM_MAX_REACTIONS];         // 温度指数
    double Ea[CHEM_MAX_REACTIONS];         // 活性化エネルギー [J/mol]
    int    nR[CHEM_MAX_REACTIONS];         // 反応物側の化学種数
    int    nP[CHEM_MAX_REACTIONS];         // 生成物側の化学種数
    int    rIdx[CHEM_MAX_REACTIONS][CHEM_MAX_SIDE];   // 反応物 化学種 index (流れ側)
    int    pIdx[CHEM_MAX_REACTIONS][CHEM_MAX_SIDE];
    int    rNu [CHEM_MAX_REACTIONS][CHEM_MAX_SIDE];   // 化学量論係数 (整数)
    int    pNu [CHEM_MAX_REACTIONS][CHEM_MAX_SIDE];
    int    thirdBody [CHEM_MAX_REACTIONS];            // 1: 三体反応 (+M)
    int    reversible[CHEM_MAX_REACTIONS];            // 1: 逆反応を K_c から
    double dnu[CHEM_MAX_REACTIONS];                   // Σ(ν''-ν') (三体 M を除く)
    double eff[CHEM_MAX_REACTIONS][THERMO_MAX_SPECIES];   // 三体効率 (既定 1)
    // falloff (+M): k = k∞ Pr/(1+Pr) F, Pr = k0 [M]/k∞。A/b/Ea が高圧極限 k∞、A0/b0/Ea0 が低圧極限 k0。
    int    falloff[CHEM_MAX_REACTIONS];               // 0: なし, 1: Lindemann (F=1), 2: Troe
    double A0[CHEM_MAX_REACTIONS], b0[CHEM_MAX_REACTIONS], Ea0[CHEM_MAX_REACTIONS];
    double troeA[CHEM_MAX_REACTIONS], troeT3[CHEM_MAX_REACTIONS], troeT1[CHEM_MAX_REACTIONS], troeT2[CHEM_MAX_REACTIONS]; // T2<=0 で項なし
};

// falloff の ln k_f(T, [M]) (Lindemann / Troe)。微分は呼び出し側で有限差分 (滑らかなスカラ関数、相対刻み 1e-4)。
THERMO_HD double chem_ln_kf_falloff(const ReactionTable* rt, int r, double T, double lnT, double RuT, double M)
{
    const double lnkinf = rt->lnA[r]  + rt->b[r]*lnT  - rt->Ea[r]/RuT;
    const double lnk0   = rt->lnA0[r] + rt->b0[r]*lnT - rt->Ea0[r]/RuT;
    const double Mp = (M > 1.0e-300) ? M : 1.0e-300;
    const double lnPr = lnk0 + log(Mp) - lnkinf;
    const double Pr = exp(lnPr);
    double lnF = 0.0;
    if (rt->falloff[r] == 2) {
        double Fc = (1.0 - rt->troeA[r]) * exp(-T / rt->troeT3[r]) + rt->troeA[r] * exp(-T / rt->troeT1[r]);
        if (rt->troeT2[r] > 0.0) Fc += exp(-rt->troeT2[r] / T);
        if (Fc < 1.0e-300) Fc = 1.0e-300;
        const double lFc = log10(Fc);
        const double c = -0.4 - 0.67*lFc, n = 0.75 - 1.27*lFc, d = 0.14;
        const double lPr = lnPr / 2.302585092994046;
        const double x = (lPr + c) / (n - d*(lPr + c));
        lnF = 2.302585092994046 * lFc / (1.0 + x*x);
    }
    return lnkinf + log(Pr / (1.0 + Pr)) + lnF;
}

// -----------------------------------------------------------------------------
// 化学種ごとの無次元ギブス関数 g_s = S°_s/Ru − H_s/(Ru T) (絶対エンタルピー基準)。
//   sensible datum (thermoHrefTemp) を焼き込んだ係数では H が h_datum だけずれているので、
//   H_abs = thermo_h_molar + h_datum で戻す。S° は datum 無関係。
// -----------------------------------------------------------------------------
THERMO_HD void chem_species_gibbs(const SpeciesThermo* sp, int n, double T,
                                  double* g, double* H_abs)
{
    const double RuT = THERMO_RU * T;
    for (int s = 0; s < n; ++s) {
        const double H = thermo_h_molar(sp[s], T) + sp[s].h_datum;   // [J/mol] 絶対基準
        const double S = (T < sp[s].Tlo || T > sp[s].Thi)
                       ? thermo_s0_mass(sp[s], T) * sp[s].MW          // 域外は cp 一定外挿と整合
                       : thermo_s_molar_clamped(sp[s], T);           // [J/(mol K)]
        g[s] = S / THERMO_RU - H / RuT;
        if (H_abs) H_abs[s] = H;
    }
}

// -----------------------------------------------------------------------------
// 生成率・反応熱・Jacobian を 1 セル分評価する。
//   ro [kg/m3], T [K], Y[n] 質量分率 (Σ=1, 各 ≥0)
//   omega[n]  : ω_s [kg/m3/s]
//   Qdot      : 反応熱 −Σ_s c_s ω_s [W/m3] (c_s = h_datum_s/W_s。絶対 datum なら 0)
//   jacMode   : 0 = Jacobian なし, 1 = 対角のみ J[s*n+s], 2 = 全 n×n (J[s*n+k] = ∂ω_s/∂(ρY_k) at fixed T [1/s])
//   dOmegadT  : ∂ω_s/∂T at fixed [X] [kg/m3/s/K] (nullptr 可)
//   T_eval は [Tlo_eval, Tmax] にクランプして速度式を評価する (低温の exp(−Ea/RT) アンダーフロー・
//   高温の暴走防止)。反応の凍結 (T < Tfreeze) は呼び出し側で ω=0 とする。
// -----------------------------------------------------------------------------
THERMO_HD void chem_source(const SpeciesThermo* sp, const ReactionTable* rt,
                           double ro, double T, const double* Y,
                           double* omega, double* Qdot,
                           int jacMode, double* J, double* dOmegadT)
{
    const int n = rt->nSpecies;
    double X[THERMO_MAX_SPECIES];          // 濃度 [mol/m3]
    double g[THERMO_MAX_SPECIES];          // 無次元ギブス
    double H[THERMO_MAX_SPECIES];          // 絶対 H [J/mol]
    for (int s = 0; s < n; ++s) {
        const double y = (Y[s] > 0.0) ? Y[s] : 0.0;
        X[s] = ro * y / sp[s].MW;
        omega[s] = 0.0;
        if (dOmegadT) dOmegadT[s] = 0.0;
    }
    if (jacMode == 1) { for (int s = 0; s < n; ++s) J[s*n+s] = 0.0; }
    else if (jacMode == 2) { for (int i = 0; i < n*n; ++i) J[i] = 0.0; }

    chem_species_gibbs(sp, n, T, g, H);
    const double RuT   = THERMO_RU * T;
    const double invT  = 1.0 / T;
    const double lnT   = log(T);
    const double cstd  = CHEM_P_STD / RuT;     // (p°/RuT) [mol/m3]
    const double lncstd = log(cstd);

    for (int r = 0; r < rt->nReac; ++r) {
        // ---- 三体濃度 [M] (三体反応・falloff で使用) ----
        double M = 1.0;
        if (rt->thirdBody[r] || rt->falloff[r]) {
            M = 0.0;
            for (int s = 0; s < n; ++s) M += rt->eff[r][s] * X[s];
        }
        // ---- 正反応速度定数 kf = A T^b exp(-Ea/RuT) (falloff は k∞ Pr/(1+Pr) F) ----
        double lnkf, dlnkf_dT, dkf_dM_over_kf = 0.0;
        if (rt->falloff[r]) {
            lnkf = chem_ln_kf_falloff(rt, r, T, lnT, RuT, M);
            const double dT = 1.0e-4 * T, dM = 1.0e-4 * ((M > 1.0e-30) ? M : 1.0e-30);
            const double lp = chem_ln_kf_falloff(rt, r, T + dT, log(T + dT), THERMO_RU*(T + dT), M);
            const double lm = chem_ln_kf_falloff(rt, r, T - dT, log(T - dT), THERMO_RU*(T - dT), M);
            dlnkf_dT = (lp - lm) / (2.0*dT);
            const double mp = chem_ln_kf_falloff(rt, r, T, lnT, RuT, M + dM);
            const double mm = chem_ln_kf_falloff(rt, r, T, lnT, RuT, (M - dM > 0.0) ? M - dM : 0.0);
            dkf_dM_over_kf = (mp - mm) / (((M - dM > 0.0) ? 2.0*dM : M + dM));   // ∂ln kf/∂[M]
        } else {
            lnkf = rt->lnA[r] + rt->b[r]*lnT - rt->Ea[r]/RuT;
            dlnkf_dT = rt->b[r]*invT + rt->Ea[r]/(RuT*T);
        }
        const double kf = exp(lnkf);

        // ---- 逆反応速度定数 kb = kf / Kc ----
        double kb = 0.0, dlnkb_dT = 0.0;
        if (rt->reversible[r]) {
            double sumg = 0.0, dH = 0.0;
            for (int j = 0; j < rt->nP[r]; ++j) { sumg += rt->pNu[r][j]*g[rt->pIdx[r][j]]; dH += rt->pNu[r][j]*H[rt->pIdx[r][j]]; }
            for (int j = 0; j < rt->nR[r]; ++j) { sumg -= rt->rNu[r][j]*g[rt->rIdx[r][j]]; dH -= rt->rNu[r][j]*H[rt->rIdx[r][j]]; }
            // ln Kc = ln Kp + dnu ln(p°/RuT), ln Kp = Σν g
            const double lnKc = sumg + rt->dnu[r]*lncstd;
            kb = exp(lnkf - lnKc);
            // d ln Kc/dT = ΔH/(Ru T²) − dnu/T  (van 't Hoff)
            const double dlnKc_dT = dH/(RuT*T) - rt->dnu[r]*invT;
            dlnkb_dT = dlnkf_dT - dlnKc_dT;
        }

        // ---- 質量作用則 Pf = Π X^ν', Pb = Π X^ν'' と各反応物/生成物に対する微分 ----
        double Pf = 1.0, Pb = 1.0;
        double dPf[CHEM_MAX_SIDE], dPb[CHEM_MAX_SIDE];
        for (int j = 0; j < rt->nR[r]; ++j) {
            const double x = X[rt->rIdx[r][j]]; const int nu = rt->rNu[r][j];
            double xp = 1.0; for (int k = 0; k < nu; ++k) xp *= x;
            Pf *= xp;
        }
        for (int j = 0; j < rt->nP[r]; ++j) {
            const double x = X[rt->pIdx[r][j]]; const int nu = rt->pNu[r][j];
            double xp = 1.0; for (int k = 0; k < nu; ++k) xp *= x;
            Pb *= xp;
        }
        if (jacMode > 0) {
            for (int j = 0; j < rt->nR[r]; ++j) {
                // ∂Pf/∂X_j = ν_j X_j^(ν_j−1) Π_{i≠j} X_i^ν_i
                double prod = 1.0;
                for (int i = 0; i < rt->nR[r]; ++i) {
                    const double x = X[rt->rIdx[r][i]]; int nu = rt->rNu[r][i];
                    if (i == j) { nu -= 1; prod *= rt->rNu[r][i]; }
                    for (int k = 0; k < nu; ++k) prod *= x;
                }
                dPf[j] = prod;
            }
            for (int j = 0; j < rt->nP[r]; ++j) {
                double prod = 1.0;
                for (int i = 0; i < rt->nP[r]; ++i) {
                    const double x = X[rt->pIdx[r][i]]; int nu = rt->pNu[r][i];
                    if (i == j) { nu -= 1; prod *= rt->pNu[r][i]; }
                    for (int k = 0; k < nu; ++k) prod *= x;
                }
                dPb[j] = prod;
            }
        }

        // falloff では [M] は Pr 経由で kf に入るので進行率の乗数にしない
        const double Mq = rt->thirdBody[r] ? M : 1.0;
        const double qf = kf*Pf, qb = kb*Pb;
        const double q  = Mq * (qf - qb);                // [mol/m3/s]

        // ---- ω_s 加算 ----
        for (int j = 0; j < rt->nR[r]; ++j) { const int s = rt->rIdx[r][j]; omega[s] -= sp[s].MW * rt->rNu[r][j] * q; }
        for (int j = 0; j < rt->nP[r]; ++j) { const int s = rt->pIdx[r][j]; omega[s] += sp[s].MW * rt->pNu[r][j] * q; }

        // ---- ∂ω_s/∂T (固定 [X]) ----
        if (dOmegadT) {
            const double dq_dT = Mq * (qf*dlnkf_dT - qb*dlnkb_dT);
            for (int j = 0; j < rt->nR[r]; ++j) { const int s = rt->rIdx[r][j]; dOmegadT[s] -= sp[s].MW * rt->rNu[r][j] * dq_dT; }
            for (int j = 0; j < rt->nP[r]; ++j) { const int s = rt->pIdx[r][j]; dOmegadT[s] += sp[s].MW * rt->pNu[r][j] * dq_dT; }
        }

        // ---- Jacobian ∂ω_s/∂(ρY_k) = W_s (ν''_s−ν'_s) ∂q/∂X_k / W_k ----
        if (jacMode > 0) {
            // ∂q/∂X_k を、この反応に関与する k (反応物・生成物・三体効率 ≠0 の全種) について評価
            // 三体項: ∂q/∂X_k += eff_k (qf − qb)
            // 濃度項: ∂q/∂X_k += M (kf ∂Pf/∂X_k − kb ∂Pb/∂X_k)
            // 寄与先 s は反応物・生成物のみ。
            auto accumulate = [&](int k, double dq_dXk) {
                const double invWk = 1.0 / sp[k].MW;
                for (int j = 0; j < rt->nR[r]; ++j) {
                    const int s = rt->rIdx[r][j];
                    if (jacMode == 1 && s != k) continue;
                    J[s*n+k] -= sp[s].MW * rt->rNu[r][j] * dq_dXk * invWk;
                }
                for (int j = 0; j < rt->nP[r]; ++j) {
                    const int s = rt->pIdx[r][j];
                    if (jacMode == 1 && s != k) continue;
                    J[s*n+k] += sp[s].MW * rt->pNu[r][j] * dq_dXk * invWk;
                }
            };
            for (int j = 0; j < rt->nR[r]; ++j) accumulate(rt->rIdx[r][j],  Mq * kf * dPf[j]);
            for (int j = 0; j < rt->nP[r]; ++j) accumulate(rt->pIdx[r][j], -Mq * kb * dPb[j]);
            if (rt->thirdBody[r]) {
                for (int k = 0; k < n; ++k) {
                    if (rt->eff[r][k] != 0.0) accumulate(k, rt->eff[r][k] * (qf - qb));
                }
            } else if (rt->falloff[r]) {
                // ∂q/∂X_k += (∂kf/∂[M]) eff_k (Pf − Pb/Kc) = eff_k (∂ln kf/∂[M]) (qf − qb)  (kb = kf/Kc)
                for (int k = 0; k < n; ++k) {
                    if (rt->eff[r][k] != 0.0) accumulate(k, rt->eff[r][k] * dkf_dM_over_kf * (qf - qb));
                }
            }
        }
    }

    // ---- 反応熱 (sensible datum の残差項) ----
    double Q = 0.0;
    for (int s = 0; s < n; ++s) Q -= (sp[s].h_datum / sp[s].MW) * omega[s];
    *Qdot = Q;
}
