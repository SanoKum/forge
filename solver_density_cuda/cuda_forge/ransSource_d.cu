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
    int katoLaunder,
    int wallTreatment,
    flow_float* wf_pk,
    flow_float* res_roK,
    flow_float* res_roOmega,
    // SST 陰解法 (point-implicit) 用: 消散項ヤコビアン対角（β* ω, 2 β ω）
    flow_float* src_jac_k,
    flow_float* src_jac_omega,
    // SST-DES (docs/turbulence §8): DESmode>0 で k 消滅項を ρk^{3/2}/l_des に切替える。
    int DESmode,
    flow_float* l_des)
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

    // Kato-Launder 補正 (docs/turbulence/theory.md §7.5): 生産の片方の S を渦度 Omega に
    // 置換し、非回転の強加速 (よどみ点・ノズル喉中心線) での偽生産を抑える。
    //   S_prod = S^2 (標準, katoLaunder==0)  /  S*Omega (katoLaunder==1)
    // せん断層 S~Omega では従来同等、非回転 Omega->0 では生産->0。Omega は無旋回でも
    // フープ補正不要 (反対称部、ひずみ S^2 の +2 S_tt^2 と非対称)。
    flow_float S_prod = S_sq;
    if (katoLaunder == 1) {
        const flow_float wx = dUzdy[ic] - dUydz[ic];   // omega = rot u
        const flow_float wy = dUxdz[ic] - dUzdx[ic];
        const flow_float wz = dUydx[ic] - dUxdy[ic];
        const flow_float Om_sq = wx*wx + wy*wy + wz*wz; // = 2 Omega_ij Omega_ij = |omega|^2
        S_prod = sqrt(S_sq * Om_sq);                    // = sqrt(S_sq)*sqrt(Om_sq) = S*Omega
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

    // k 生産項（10 beta* rho k omega でリミット）。S_prod は katoLaunder で S^2 / S*Omega。
    flow_float Pk = min(mu_t_eff * S_prod,
        static_cast<flow_float>(10.0) * kBetaStar * rho * k_c * w_c);

    // (B) 等方項 -2/3 rho k divU を k 生産に加算 (dilatationCorrection >= 2, theory.md §7.3)
    //   膨張(divU>0)でシンク, 圧縮(divU<0)でソース。負生産を防ぐため Pk>=0 でクリップ。
    if (dilatationCorrection >= 2) {
        Pk -= static_cast<flow_float>(2.0 / 3.0) * rho * k_c * divU;
        Pk = max(Pk, static_cast<flow_float>(0.0));
    }

    // automatic wall treatment (docs/turbulence §6.5(d)): wall-adjacent セル (wf_pk>=0) では
    // 解像勾配ベースの P_k を wall-function 生産 P_k=ρu_τ⁴/ν·g(1-g) に置換する。ω はピン留め済み
    // (ransBoundary) なので k は平衡値 u_τ²/√β* に自律収束する。非 wall セルは wf_pk<0 で不変。
    if (wallTreatment == 1 && wf_pk[ic] >= static_cast<flow_float>(0.0)) {
        Pk = wf_pk[ic];
    }

    // k 消滅項。標準 SST は D_k = β* ρ k ω = ρ k^{3/2}/l_RANS。
    // SST-DES (DESmode>0) では l_RANS を l_DDES に置換: D_k = ρ k^{3/2}/l_des (docs/turbulence §8)。
    // 剥離域 (l_des < l_RANS) で消滅が強まり渦粘性が下がる。ω 方程式・渦粘性式は不変。
    // DESmode==0 では従来式そのまま (係数・順序・キャスト一致) = ビット不変。
    flow_float Dk;
    flow_float jac_k;  // 陰解法ヤコビアン対角 ∂Dk/∂(ρk)
    if (DESmode > 0) {
        const flow_float l_hyb = max(l_des[ic], static_cast<flow_float>(1.0e-20));
        Dk    = rho * pow(k_c, static_cast<flow_float>(1.5)) / l_hyb;
        jac_k = static_cast<flow_float>(1.5) * sqrt(k_c) / l_hyb;
    } else {
        Dk    = kBetaStar * rho * k_c * w_c;
        jac_k = kBetaStar * w_c;
    }

    // omega 生産項
    const flow_float Pw = alpha * rho * S_prod;

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
    // 生産項は limiter 付きで bounded のため lagged（陰化しない）。生産振幅を正の対角に足す
    // under-relaxation も試したが、安定 cfl_pseudo を一切上げず（律速は k/ω 輸送項の陽的扱い）
    // 不採用。輸送項の point-implicit 化は scalarTransport_d.cu（transport_diag）を参照。
    src_jac_k[ic]     = jac_k;  // DESmode==0: β* ω (従来式と一致), DESmode>0: 1.5 √k / l_des
    src_jac_omega[ic] = static_cast<flow_float>(2.0) * beta * w_c;

    // automatic wall treatment: wall-adjacent セルは ω をピン留め (Dirichlet, ransBoundary) する。
    // Dirichlet セルの残差は 0 が正しい。res_roOmega を 0 にし (transport+source 込みを上書き)、
    // 無駄な update と rms_roOmega の汚染 (ピン留めの巨大残差) を防ぐ。k は解くので残差は残す。
    // 注: 非軸対称ケース前提 (axisymmetricSource が後段で res_roOmega を足す場合は別途零化が要る)。
    if (wallTreatment == 1 && wf_pk[ic] >= static_cast<flow_float>(0.0)) {
        res_roOmega[ic]   = static_cast<flow_float>(0.0);
        src_jac_omega[ic] = static_cast<flow_float>(0.0);
    }
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
        cfg.katoLaunder,
        cfg.wallTreatmentSST,
        var.c_d["wf_pk"],
        var.c_d["res_roK"],
        var.c_d["res_roOmega"],
        var.c_d["src_jac_k"],
        var.c_d["src_jac_omega"],
        cfg.DESmode,
        var.c_d["l_des"]);

    gpuErrchk(cudaPeekAtLastError());
    gpuErrchkKernelSync();
}
