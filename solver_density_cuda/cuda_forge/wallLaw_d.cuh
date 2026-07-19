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

// Kader 温度壁法則 T⁺(y⁺; Pr, Pr_t) (theory §10.3)。u⁺ は収束済み u_τ の Reichardt 値を渡す。
//   Γ→0 (y⁺→0) では exp(-1/Γ)→0 で T⁺→Pr·y⁺ (伝導低層)、大 y⁺ では対数層漸近。
__device__ __forceinline__ flow_float wallLaw_kader_tplus(
    flow_float Pr, flow_float Prt, flow_float yp, flow_float uplus)
{
    const flow_float pry   = Pr * yp;
    const flow_float gam   = static_cast<flow_float>(0.01) * pry * pry * pry * pry
                           / (static_cast<flow_float>(1.0)
                              + static_cast<flow_float>(5.0) * Pr * Pr * Pr * yp);
    const flow_float pfn0  = static_cast<flow_float>(3.85) * cbrt(Pr) - static_cast<flow_float>(1.3);
    const flow_float pfn   = pfn0 * pfn0 + static_cast<flow_float>(2.12) * log(Pr);
    // Γ=0 のとき -1/Γ = -inf → exp = 0 (IEEE)。NaN 防止に下限だけガード。
    const flow_float einv  = exp(-static_cast<flow_float>(1.0) / max(gam, static_cast<flow_float>(1.0e-30)));
    return pry * exp(-gam) + Prt * (uplus + pfn) * einv;
}
