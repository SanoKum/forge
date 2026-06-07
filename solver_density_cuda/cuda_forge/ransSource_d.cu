#include "ransSource_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"

namespace {

constexpr flow_float kBetaStar  = static_cast<flow_float>(0.09);
constexpr flow_float kAlpha1    = static_cast<flow_float>(5.0 / 9.0);
constexpr flow_float kAlpha2    = static_cast<flow_float>(0.44);
constexpr flow_float kBeta1     = static_cast<flow_float>(0.075);
constexpr flow_float kBeta2     = static_cast<flow_float>(0.0828);
constexpr flow_float kSigmaW2   = static_cast<flow_float>(0.856);
constexpr flow_float kA1        = static_cast<flow_float>(0.31);
constexpr flow_float kSmall     = static_cast<flow_float>(1.0e-12);

// SST source terms (production, destruction, cross-diffusion) for k and omega.
// vis_turb にはこの段階での渦粘性が入る（Stage 4 以前はゼロまたは簡易値）。
// 生産項の mu_t が 0 の場合は内部で rho*a1*k/omega で代替計算する。
__global__ void rans_sst_source_d(
    geom_int nCells,
    geom_float* vol,
    flow_float* ro,
    flow_float* k,
    flow_float* omega,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float* wall_dist,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* dKdx,     flow_float* dKdy,     flow_float* dKdz,
    flow_float* dOmegadx, flow_float* dOmegady, flow_float* dOmegadz,
    int isAxisymmetric,
    flow_float* axisym_uy_over_r,
    int dilatationCorrection,
    flow_float* res_roK,
    flow_float* res_roOmega,
    // SST 陰解法 (point-implicit) 用: 消散項ヤコビアン対角（β* ω, 2 β ω）
    flow_float* src_jac_k,
    flow_float* src_jac_omega)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    const flow_float rho    = max(ro[ic], kSmall);
    const flow_float k_c    = max(k[ic],  static_cast<flow_float>(0.0));
    const flow_float w_c    = max(omega[ic], kSmall);
    const flow_float mu_lam = vis_lam[ic];
    const flow_float mu_t   = max(vis_turb[ic], static_cast<flow_float>(0.0));
    const flow_float y      = max(wall_dist[ic], kSmall);
    const geom_float v      = vol[ic];

    // 渦粘性がゼロ（SST eddy viscosity 未実装段階）では k/omega から代替
    const flow_float mu_t_eff = (mu_t > kSmall) ? mu_t : rho * kA1 * k_c / w_c;

    // ひずみ速度の大きさ squared: S^2 = 2 Sij Sij
    const flow_float S11 = dUxdx[ic];
    const flow_float S22 = dUydy[ic];
    const flow_float S33 = dUzdz[ic];
    const flow_float S12 = static_cast<flow_float>(0.5) * (dUxdy[ic] + dUydx[ic]);
    const flow_float S13 = static_cast<flow_float>(0.5) * (dUxdz[ic] + dUzdx[ic]);
    const flow_float S23 = static_cast<flow_float>(0.5) * (dUydz[ic] + dUzdy[ic]);
    flow_float S_sq = static_cast<flow_float>(2.0) *
        (S11*S11 + S22*S22 + S33*S33 +
         static_cast<flow_float>(2.0) * (S12*S12 + S13*S13 + S23*S23));
    // 軸対称: フープひずみ S_thetatheta = u_r/r を加算 (docs/turbulence/theory.md §7.2)
    if (isAxisymmetric) {
        const flow_float S_tt = axisym_uy_over_r[ic];
        S_sq += static_cast<flow_float>(2.0) * S_tt * S_tt;
    }

    // 速度の発散 divU = ∂xUx + ∂yUy + ∂zUz (+ 軸対称 u_r/r) = axisym_divU と一致
    flow_float divU = S11 + S22 + S33;
    if (isAxisymmetric) {
        divU += axisym_uy_over_r[ic];
    }

    // (A) deviatoric トレース除去 (dilatationCorrection >= 1, docs/turbulence/theory.md §7.3)
    //   S^2 = 2 Sij Sij から 2/3 (divU)^2 を引く (= 2 * 1/3 (divU)^2)
    if (dilatationCorrection >= 1) {
        S_sq -= static_cast<flow_float>(2.0 / 3.0) * divU * divU;
        S_sq = max(S_sq, static_cast<flow_float>(0.0));
    }

    // 交差拡散項（F1 ブレンドに使用）
    const flow_float grad_k_dot_w =
        dKdx[ic] * dOmegadx[ic] + dKdy[ic] * dOmegady[ic] + dKdz[ic] * dOmegadz[ic];
    const flow_float CD_kw = max(
        static_cast<flow_float>(2.0) * rho * kSigmaW2 / w_c * grad_k_dot_w,
        static_cast<flow_float>(1.0e-10));

    // F1 ブレンド関数
    const flow_float nu     = mu_lam / rho;
    const flow_float arg1_a = sqrt(k_c) / (kBetaStar * w_c * y);
    const flow_float arg1_b = static_cast<flow_float>(500.0) * nu / (w_c * y * y);
    const flow_float arg1_c = static_cast<flow_float>(4.0) * rho * kSigmaW2 * k_c / (CD_kw * y * y);
    const flow_float arg1   = min(max(arg1_a, arg1_b), arg1_c);
    const flow_float F1     = tanh(arg1 * arg1 * arg1 * arg1);

    // ブレンドされた係数
    const flow_float alpha = F1 * kAlpha1 + (static_cast<flow_float>(1.0) - F1) * kAlpha2;
    const flow_float beta  = F1 * kBeta1  + (static_cast<flow_float>(1.0) - F1) * kBeta2;

    // k 生産項（10 beta* rho k omega でリミット）
    flow_float Pk = min(mu_t_eff * S_sq,
        static_cast<flow_float>(10.0) * kBetaStar * rho * k_c * w_c);

    // (B) 等方項 -2/3 rho k divU を k 生産に加算 (dilatationCorrection >= 2, theory.md §7.3)
    //   膨張(divU>0)でシンク, 圧縮(divU<0)でソース。負生産を防ぐため Pk>=0 でクリップ。
    if (dilatationCorrection >= 2) {
        Pk -= static_cast<flow_float>(2.0 / 3.0) * rho * k_c * divU;
        Pk = max(Pk, static_cast<flow_float>(0.0));
    }

    // k 消滅項
    const flow_float Dk = kBetaStar * rho * k_c * w_c;

    // omega 生産項
    const flow_float Pw = alpha * rho * S_sq;

    // omega 消滅項
    const flow_float Dw = beta * rho * w_c * w_c;

    // omega 交差拡散ソース（k-epsilon 側でのみ有効: 1 - F1）
    const flow_float CDw = (static_cast<flow_float>(1.0) - F1) *
        static_cast<flow_float>(2.0) * rho * kSigmaW2 / w_c * grad_k_dot_w;

    // 残差に加算（符号: ソース項は保存変数を増加させる正方向）
    atomicAdd(&res_roK[ic],     (Pk - Dk) * v);
    atomicAdd(&res_roOmega[ic], (Pw - Dw + CDw) * v);

    // SST 陰解法用の消散ヤコビアン対角（point-implicit で D に加える正の量）。
    //   Dk = β* ρ k ω = β* (ρk) ω  →  ∂Dk/∂(ρk) = β* ω
    //   Dω = β ρ ω²  = β (ρω)²/ρ   →  ∂Dω/∂(ρω) = 2 β ω
    // 生産項は limiter 付きで bounded のため lagged（陰化しない）。
    src_jac_k[ic]     = kBetaStar * w_c;
    src_jac_omega[ic] = static_cast<flow_float>(2.0) * beta * w_c;
}

}

void ransSource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1)) {
        return;
    }

    rans_sst_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["volume"],
        var.c_d["ro"],
        var.c_d["k"],
        var.c_d["omega"],
        var.c_d["vis_lam"],
        var.c_d["vis_turb"],
        var.c_d["wall_dist"],
        var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
        var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
        var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],
        var.c_d["dKdx"],     var.c_d["dKdy"],     var.c_d["dKdz"],
        var.c_d["dOmegadx"], var.c_d["dOmegady"], var.c_d["dOmegadz"],
        cfg.isAxisymmetric,
        var.c_d["axisym_uy_over_r"],
        cfg.dilatationCorrection,
        var.c_d["res_roK"],
        var.c_d["res_roOmega"],
        var.c_d["src_jac_k"],
        var.c_d["src_jac_omega"]);

    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
}
