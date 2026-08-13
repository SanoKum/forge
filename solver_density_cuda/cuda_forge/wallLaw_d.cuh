#pragma once

#include "flowFormat.hpp"

// 壁法則デバイス関数群 (Reichardt / Kader / u_τ Newton)。
// SST automatic wall treatment (ransWallFunction_d.cu) と WMLES 代数壁応力モデル
// (wmlesWallModel_d.cu) が共有する。理論: methods/turbulence/theory.md §6.5(a), §10。
// Reichardt 2 関数は ransWallFunction_d.cu のファイルローカル実装をそのまま昇格したもので、
// 演算列は不変 (SST 経路のビット不変を維持)。

constexpr flow_float kWallLawKappa = static_cast<flow_float>(0.41);   // von Karman 定数
constexpr flow_float kWallLawSmall = static_cast<flow_float>(1.0e-12);

// Reichardt 則 u⁺(y⁺)
__device__ __forceinline__ flow_float wallLaw_reichardt_uplus(flow_float yp)
{
    const flow_float lg = log(static_cast<flow_float>(1.0) + kWallLawKappa * yp) / kWallLawKappa;
    const flow_float e1 = exp(-yp / static_cast<flow_float>(11.0));
    const flow_float e2 = exp(-yp / static_cast<flow_float>(3.0));
    const flow_float tr = static_cast<flow_float>(7.8) *
        (static_cast<flow_float>(1.0) - e1 - (yp / static_cast<flow_float>(11.0)) * e2);
    return lg + tr;
}

// du⁺/dy⁺
__device__ __forceinline__ flow_float wallLaw_reichardt_duplus_dyp(flow_float yp)
{
    const flow_float dlg = static_cast<flow_float>(1.0) / (static_cast<flow_float>(1.0) + kWallLawKappa * yp);
    const flow_float e1  = exp(-yp / static_cast<flow_float>(11.0));
    const flow_float e2  = exp(-yp / static_cast<flow_float>(3.0));
    // d/dy⁺ [1 - e1 - (y⁺/11) e2] = e1/11 - e2/11 + (y⁺/33) e2
    const flow_float dtr = static_cast<flow_float>(7.8) *
        (e1 / static_cast<flow_float>(11.0) - e2 / static_cast<flow_float>(11.0)
         + (yp / static_cast<flow_float>(33.0)) * e2);
    return dlg + dtr;
}

// ---------------------------------------------------------------------------
// Spalding 則 (診断実験専用) — plans/active/turbulence-reichardt-gap-residual.md §5.2 ②
//
// Spalding は **陰形式** y⁺ = H(u⁺) で与えられる:
//   H(u⁺) = u⁺ + e^{-κB}[ e^{κu⁺} − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6 ]
// したがって u⁺(y⁺) を内側 Newton で解いてから外側で u_τ を解く**二重 Newton は不要**。
//   F(u_τ) = H(U_t/u_τ) − u_τ y/ν = 0
// を直接 Newton で解けばよい。勾配も解析的に得られる:
//   du⁺/dy⁺ = 1/H'(u⁺)
//
// 定数は case/26 の後処理ツール (tools/wall_law_inverse.py) と同一: κ=0.41, B=5.0。
// **この経路は診断実験専用**であり、既定 (selector=0) では一切呼ばれない。
constexpr flow_float kWallLawSpaldingB = static_cast<flow_float>(5.0);

// y⁺ = H(u⁺)
__device__ __forceinline__ flow_float wallLaw_spalding_yplus(flow_float up)
{
    const flow_float e = kWallLawKappa * up;
    const flow_float ex = exp(e);
    return up + exp(-kWallLawKappa * kWallLawSpaldingB) *
        (ex - static_cast<flow_float>(1.0) - e
         - e * e / static_cast<flow_float>(2.0)
         - e * e * e / static_cast<flow_float>(6.0));
}

// H'(u⁺) = dy⁺/du⁺
__device__ __forceinline__ flow_float wallLaw_spalding_dyplus_dup(flow_float up)
{
    const flow_float e = kWallLawKappa * up;
    const flow_float ex = exp(e);
    return static_cast<flow_float>(1.0) + exp(-kWallLawKappa * kWallLawSpaldingB) * kWallLawKappa *
        (ex - static_cast<flow_float>(1.0) - e - e * e / static_cast<flow_float>(2.0));
}

// du⁺/dy⁺ = 1/H'(u⁺)。**整合実験 (§5.2 ③) 用**。診断実験 (②) では使わない
// (g は Reichardt 固定が仕様)。
__device__ __forceinline__ flow_float wallLaw_spalding_duplus_dyp_from_up(flow_float up)
{
    return static_cast<flow_float>(1.0) / max(wallLaw_spalding_dyplus_dup(up), kWallLawSmall);
}

// F(u_τ) = H(U_t/u_τ) − u_τ y/ν = 0 の直接 Newton。
//   dF/du_τ = H'(u⁺)·(−U_t/u_τ²) − y/ν  < 0 (単調)
__device__ __forceinline__ flow_float wallLaw_spalding_solve_utau(
    flow_float Ut, flow_float y, flow_float nu, int maxIt)
{
    flow_float utau = sqrt(max(nu * Ut / max(y, kWallLawSmall), kWallLawSmall));  // 層流推定
    for (int it = 0; it < maxIt; ++it) {
        utau = max(utau, kWallLawSmall);
        const flow_float up = Ut / utau;
        const flow_float F  = wallLaw_spalding_yplus(up) - utau * y / nu;
        const flow_float dF = wallLaw_spalding_dyplus_dup(up) * (-Ut / (utau * utau)) - y / nu;
        if (fabs(dF) < kWallLawSmall) break;
        const flow_float du = -F / dF;
        utau += du;
        if (utau <= static_cast<flow_float>(0.0))
            utau = sqrt(max(nu * Ut / max(y, kWallLawSmall), kWallLawSmall));
        if (fabs(du) < static_cast<flow_float>(1.0e-10) * max(utau, kWallLawSmall)) break;
    }
    return max(utau, static_cast<flow_float>(0.0));
}

// WMLES 用 u_τ ソルバ (theory §10.2): F(u_τ) = u_τ·f(y⁺(u_τ)) − u_∥ = 0 の Newton。
// SST 壁関数の残差形 U_t/u_τ − u⁺ と根は同一だが、この形は u_τ→0 で正則で warm start 向き。
// F' = f + y⁺ f' > 0 (f, f' > 0) なので単調・ゼロ割なし。
//   utau: in/out。入力は前 step の warm start 値 (≤0 / NaN なら層流推定に置き直す)。
//   戻り値: 収束フラグ (false = maxIt で未収束。呼び出し側で層流フォールバック + カウンタ)。
__device__ __forceinline__ bool wallLaw_solve_utau(
    flow_float upar, flow_float d, flow_float rho_w, flow_float mu_w,
    flow_float tol, int maxIt, flow_float& utau)
{
    const flow_float nu = mu_w / rho_w;
    const flow_float utauLam = sqrt(mu_w * upar / (rho_w * d));   // 層流推定 (初期値・復帰値)
    if (!(utau > static_cast<flow_float>(0.0)) || isnan(utau)) utau = utauLam;
    for (int it = 0; it < maxIt; ++it) {
        const flow_float yp = utau * d / nu;
        const flow_float f  = wallLaw_reichardt_uplus(yp);
        const flow_float fp = wallLaw_reichardt_duplus_dyp(yp);
        const flow_float F  = utau * f - upar;
        const flow_float dF = f + yp * fp;
        const flow_float du = -F / max(dF, kWallLawSmall);
        utau += du;
        if (utau <= static_cast<flow_float>(0.0)) utau = utauLam;   // 下限クリップ: 層流推定に置き直し
        if (fabs(du) < tol * max(utau, kWallLawSmall)) return true;
    }
    return false;
}

// Kader 温度壁法則 T⁺(y⁺; Pr) (Kader 1981 の原式, theory §10.3):
//   T⁺ = Pr·y⁺·e^{−Γ} + [2.12·ln(1+y⁺) + β(Pr)]·e^{−1/Γ}
//   β(Pr) = (3.85·Pr^{1/3} − 1.3)² + 2.12·ln(Pr),  Γ = 0.01(Pr·y⁺)⁴/(1+5Pr³y⁺)
//   Γ→0 (y⁺→0) では exp(-1/Γ)→0 で T⁺→Pr·y⁺ (伝導低層)、大 y⁺ では 2.12·ln y⁺+β に漸近。
//   【修正 2026-08-11】旧実装は対数部を Pr_t·(u⁺+β) としていた — Jayatilleke 型の運動量
//   アナロジー形 (Pr_t(u⁺+P_jat)) と Kader の β (2.12·ln(1+y⁺) と対で較正) の混用で、
//   log 層 T⁺ を +30% 級過大評価し等温壁 q_w が系統的に過小になっていた
//   (case/26 平板 y+30: Colburn/low-Re 基準比 −13%)。原式に戻す。u⁺/Pr_t は不要になった。
__device__ __forceinline__ flow_float wallLaw_kader_tplus(
    flow_float Pr, flow_float yp)
{
    const flow_float pry   = Pr * yp;
    const flow_float gam   = static_cast<flow_float>(0.01) * pry * pry * pry * pry
                           / (static_cast<flow_float>(1.0)
                              + static_cast<flow_float>(5.0) * Pr * Pr * Pr * yp);
    const flow_float pfn0  = static_cast<flow_float>(3.85) * cbrt(Pr) - static_cast<flow_float>(1.3);
    const flow_float beta  = pfn0 * pfn0 + static_cast<flow_float>(2.12) * log(Pr);
    // Γ=0 のとき -1/Γ = -inf → exp = 0 (IEEE)。NaN 防止に下限だけガード。
    const flow_float einv  = exp(-static_cast<flow_float>(1.0) / max(gam, static_cast<flow_float>(1.0e-30)));
    return pry * exp(-gam)
         + (static_cast<flow_float>(2.12) * log(static_cast<flow_float>(1.0) + yp) + beta) * einv;
}
