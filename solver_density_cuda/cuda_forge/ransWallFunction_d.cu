#include "ransWallFunction_d.cuh"

#include <cstdio>
#include <cstdlib>

#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/wallLaw_d.cuh"

// SST automatic wall treatment の摩擦速度 u_τ を Reichardt 普遍速度則の逆解きで求める。
// 理論は methods/turbulence/theory.md §6.5 (a)。粘性低層・バッファ・対数を 1 式で繋ぐため
// 第一セルの y⁺ 位置に依存せず妥当な u_τ を返す。
// Reichardt u⁺(y⁺)/du⁺dy⁺ は wallLaw_d.cuh へ昇格 (WMLES 壁モデルと共有、演算列は不変)。
namespace {

constexpr flow_float kKappa    = kWallLawKappa;                  // von Karman 定数
constexpr flow_float kSmall    = kWallLawSmall;
constexpr int        kNewtonIt = 5;                              // Newton 反復回数
constexpr flow_float kSstBeta1  = static_cast<flow_float>(0.075); // SST β₁ (ω_vis 用)
constexpr flow_float kSstBetaSt = static_cast<flow_float>(0.09);  // SST β* (ω_log・k_wf 用)
constexpr flow_float kOmegaMin  = static_cast<flow_float>(1.0e-10);

// wf_pk・Tau_Wall・roK_wf を全セル -1 (inactive) に初期化する。
// Qw_Wall も -1 化する (sstEnergyWallFunction の AddQWall マーカ。§6.5(g)。WMLES と SST は
// 共存しない (wmlesActiveForBcond は RANS で false) ため二重初期化にならない)。
// Taw_Wall_Flux も -1 化する (§6.5(f) 断熱壁の流束モデル置換マーカ)。
__global__ void init_wf_pk_d(geom_int nCells, flow_float* wf_pk, flow_float* Tau_Wall, flow_float* roK_wf,
                             flow_float* roOmega_wf, flow_float* Taw_diag, flow_float* Qw_Wall,
                             flow_float* Taw_Wall_Flux, flow_float* Taw_Rep_Id)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        wf_pk[ic] = static_cast<flow_float>(-1.0);
        Tau_Wall[ic] = static_cast<flow_float>(-1.0);
        roK_wf[ic] = static_cast<flow_float>(-1.0);
        roOmega_wf[ic] = static_cast<flow_float>(-1.0);
        Taw_diag[ic] = static_cast<flow_float>(-1.0);
        Qw_Wall[ic] = static_cast<flow_float>(-1.0);
        Taw_Wall_Flux[ic] = static_cast<flow_float>(-1.0);
        Taw_Rep_Id[ic] = static_cast<flow_float>(-1.0);
    }
}

__global__ void compute_wall_friction_sst_d(
    geom_int nb,
    geom_int* bplane_plane,
    geom_int* bplane_cell,
    geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    flow_float* ro,
    flow_float* vis_lam,
    flow_float* wall_dist,
    flow_float* Ux, flow_float* Uy, flow_float* Uz,
    flow_float* utau_b,
    flow_float* ypls_b,
    flow_float* wf_pk,
    // node-centered (isNode=1): bplane_cell は壁ノード W (u=0,Dirichlet)。代表内部点を選ぶための
    // CV 座標と CSR 隣接。SU2 Normal_Neighbor 流。cell モード (isNode=0) では未使用 (nullptr 可)。
    int isNode,
    geom_int nNormalPlanes,
    geom_int* cell_planes_index, geom_int* cell_planes, geom_int* plane_cells,
    geom_int* wall_flag,
    geom_float* ccx, geom_float* ccy, geom_float* ccz,
    // node: 壁ノードに τ_w=ρu_τ² を格納 (viscousFlux_d の AddTauWall 再スケール用)。cell では書かない (-1 維持)。
    flow_float* Tau_Wall,
    // node: 第一内層ノードの k Dirichlet 値 (保存量 ρ·k_wf, SU2 SetTurbVars_WF 流) を格納。cell では書かない (-1 維持)。
    int enableKwf, flow_float* roK_wf, int enableOmgWf, flow_float* roOmega_wf,
    // 診断: Crocco 型 断熱回復温度 T_aw = T_rep + r·U_t²/(2cp) (r=Pr_lam^{1/3}, SU2 壁関数の熱的閉包と同式)。
    // 出力のみで壁状態には未適用。ic (cell=第一セル / node=壁ノード) に格納。
    flow_float* T_arr, flow_float* cp_arr, flow_float recovery, flow_float* Taw_diag,
    // §6.5(f) 流束置換の対象エッジ制限用: 壁ノード ic の代表点 (Normal_Neighbor) の CV id を
    // flow_float に exact 格納 (nCells < 2^24 で厳密, 既存 Tau_Wall 等と同じ流儀。非対象/未解決は
    // -1 のまま)。viscousFlux_d は Taw_Wall_Flux[W] による端点置換を、この代表点との内部辺
    // **のみ**に限定する (全内部辺に適用すると、T_aw が物理的意味を持たない接線方向隣接ノードとの
    // 温度差まで流束化し発散した実測に基づく制限, 2026-08-11)。nullptr 可。
    flow_float* Taw_Rep_Id,
    // エネルギー壁関数 (§6.5(g), sstEnergyWallFunction==1 × wall_isothermal のみ energyWf=1):
    //   Kader T⁺ で q_w = ρ cp u_τ (T_w − T_aw,rep)/T⁺ を計算し bvar qwall と (node) Qw_Wall へ書く。
    //   Tw_b は等温壁の指定壁温 (bvar Ts, 固定入力)。thermCond_arr は Pr 評価用。
    //   energyWf=0 では qwall_b/Qw_Wall に一切触れない (ビット不変)。
    int energyWf, flow_float* Tw_b, flow_float* thermCond_arr, flow_float Prt,
    flow_float* qwall_b, flow_float* Qw_Wall)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib >= nb) return;

    const geom_int ip = bplane_plane[ib];
    const geom_int ic = bplane_cell[ib];

    // 壁単位法線 (外向き)
    const flow_float sss = max(ss[ip], kSmall);
    const flow_float nx = sx[ip] / sss;
    const flow_float ny = sy[ip] / sss;
    const flow_float nz = sz[ip] / sss;

    // 代表点 irep と壁法線距離 y。
    //   cell: irep=ic (壁隣接内部セル), y=wall_dist[ic]  … 従来挙動 (ビット不変)。
    //   node: ic=壁ノード W は u=0/wall_dist≈0 で退化するため、W の入射内部双対面から壁内向き法線 -n
    //         との cos 最大の内部ノード I を SU2 Normal_Neighbor として選び irep=I, y=(x_I-x_W)·(-n)。
    geom_int irep = ic;
    flow_float y  = max(wall_dist[ic], kSmall);
    if (isNode != 0) {
        flow_float best_cos = static_cast<flow_float>(-2.0);
        geom_int   bestI = -1;
        flow_float bestDn = kSmall;
        const flow_float xw = ccx[ic], yw = ccy[ic], zw = ccz[ic];
        for (geom_int j = cell_planes_index[ic]; j < cell_planes_index[ic + 1]; ++j) {
            const geom_int ipn = cell_planes[j];
            if (ipn >= nNormalPlanes) continue;                 // 内部双対面のみ
            const geom_int a = plane_cells[2 * ipn + 0];
            const geom_int b = plane_cells[2 * ipn + 1];
            const geom_int cand = (a == ic) ? b : a;
            if (wall_flag[cand] != 0) continue;                 // 内部ノードのみ
            const flow_float dx = ccx[cand] - xw;
            const flow_float dy = ccy[cand] - yw;
            const flow_float dz = ccz[cand] - zw;
            const flow_float dist = sqrt(dx * dx + dy * dy + dz * dz);
            if (dist <= kSmall) continue;
            const flow_float dn_in = -(dx * nx + dy * ny + dz * nz);  // 壁内向き距離 (>0 が内部側)
            if (dn_in <= static_cast<flow_float>(0.0)) continue;
            const flow_float cosv = dn_in / dist;               // 内向き法線との cos
            if (cosv > best_cos) { best_cos = cosv; bestI = cand; bestDn = dn_in; }
        }
        if (bestI < 0) {                                        // 代表点なし: 退避
            utau_b[ib] = static_cast<flow_float>(0.0);
            ypls_b[ib] = static_cast<flow_float>(0.0);
            wf_pk[ic]  = static_cast<flow_float>(0.0);
            if (energyWf != 0) { qwall_b[ib] = static_cast<flow_float>(0.0); Qw_Wall[ic] = static_cast<flow_float>(0.0); }
            return;
        }
        irep = bestI;
        y    = max(bestDn, kSmall);
    }

    const flow_float rho = max(ro[irep], kSmall);
    const flow_float nu  = vis_lam[irep] / rho;

    // 壁接線速度の大きさ U_t = |U_c - (U_c·n)n| (代表点 irep の速度)
    const flow_float uc = Ux[irep];
    const flow_float vc = Uy[irep];
    const flow_float wc = Uz[irep];
    const flow_float un = uc * nx + vc * ny + wc * nz;
    const flow_float utx = uc - un * nx;
    const flow_float uty = vc - un * ny;
    const flow_float utz = wc - un * nz;
    const flow_float Ut  = sqrt(utx * utx + uty * uty + utz * utz);

    if (Ut <= kSmall) {
        // 淀み域: せん断なし → u_τ=0 (ω 壁 BC は ω_vis に縮退), 生産なし
        utau_b[ib] = static_cast<flow_float>(0.0);
        ypls_b[ib] = static_cast<flow_float>(0.0);
        wf_pk[ic]  = static_cast<flow_float>(0.0);
        Taw_diag[ic] = T_arr[irep];   // 淀み: 運動加熱なし
        if (isNode != 0 && Taw_Rep_Id != nullptr) Taw_Rep_Id[ic] = static_cast<flow_float>(irep);
        // エネルギー壁関数: 淀みは純伝導へ退避 (§6.5(g)。T⁺ 評価が u_τ=0 で退化するため)。
        if (energyWf != 0) {
            const flow_float qw = thermCond_arr[irep] * (Tw_b[ib] - T_arr[irep]) / y;
            qwall_b[ib] = qw;
            if (isNode != 0) Qw_Wall[ic] = qw;
        }
        // node: 第一内層ノード irep も wall-function 生産 (=0) で覆う。壁ノードだけだと第一内層に標準
        // 解像生産 P_k=μ_t S² が残り淀みコーナーで k 暴走するため (methods/turbulence §6.5(e))。
        if (isNode != 0) {
            wf_pk[irep] = static_cast<flow_float>(0.0); Tau_Wall[ic] = static_cast<flow_float>(0.0);
            // 淀み域 u_τ=0 → k_wf=0。node k Dirichlet (SU2 SetTurbVars_WF 流, wallTreatment==1 のみ)。
            if (enableKwf == 1) {
                roK_wf[irep] = static_cast<flow_float>(0.0);
                // 淀み: ω は粘性則値 (u_τ=0 → ω_log=0)。y は代表内点距離。
                if (enableOmgWf == 1) {
                    const flow_float ov = static_cast<flow_float>(6.0) * nu / (kSstBeta1 * y * y);
                    roOmega_wf[irep] = rho * max(ov, kOmegaMin);
                }
            }
        }
        return;
    }

    // Newton 法で f(u_τ) = U_t/u_τ - u⁺(y⁺(u_τ)) = 0 を解く。初期値は粘性則。
    flow_float utau = sqrt(nu * Ut / y);
    for (int it = 0; it < kNewtonIt; ++it) {
        utau = max(utau, kSmall);
        const flow_float yp = utau * y / nu;
        const flow_float f  = Ut / utau - wallLaw_reichardt_uplus(yp);
        const flow_float df = -Ut / (utau * utau) - wallLaw_reichardt_duplus_dyp(yp) * (y / nu);
        if (fabs(df) < kSmall) break;
        utau -= f / df;
    }
    utau = max(utau, static_cast<flow_float>(0.0));

    // wall-function 生産 P_k = (τ_w - τ_visc)·∂U/∂y = ρu_τ⁴/ν · g(1-g), g = du⁺/dy⁺(y⁺₁)
    // (methods/turbulence §6.5(d))。粘性低層 g→1 で P_k→0 (壁解像極限を保つ)、対数層 g→1/(κy⁺) で
    // P_k→ρu_τ³/(κy)。これと ω ピン留めで k が平衡値 u_τ²/√β* に収束し runaway を断つ。
    const flow_float yp1 = utau * y / nu;
    const flow_float g   = wallLaw_reichardt_duplus_dyp(yp1);
    const flow_float pk_wf = max(rho * utau * utau * utau * utau / nu * g * (static_cast<flow_float>(1.0) - g),
                    static_cast<flow_float>(0.0));
    wf_pk[ic] = pk_wf;
    // node: 第一内層ノード irep も wall-function 生産で覆う (壁ノードのみだと第一内層に標準解像生産が
    // 残り k 暴走。methods/turbulence §6.5(e))。ω ピン/decouple は壁ノードのまま (ここでは触らない)。
    if (isNode != 0) {
        wf_pk[irep] = pk_wf;
        // node k Dirichlet (SU2 SetTurbVars_WF 流, enableKwf==1 のみ): 第一内層ノードの k を
        // k_wf = ω_w·μ_t,wall/ρ = ω_w·ν·(1/g - 1) に固定し、近壁 k 蓄積 (再付着 μ_t ピーク) を断つ。
        // ω_w は Menter ブレンド (壁ノードのピンと同式)、μ_t,wall は壁モデル渦粘性 (Reichardt g 由来)。
        if (enableKwf == 1) {
            const flow_float omega_vis = static_cast<flow_float>(6.0) * nu / (kSstBeta1 * y * y);
            const flow_float omega_log = utau / (sqrt(kSstBetaSt) * kKappa * y);
            const flow_float omega_w   = max(sqrt(omega_vis * omega_vis + omega_log * omega_log), kOmegaMin);
            const flow_float mut_wall   = nu * max(static_cast<flow_float>(1.0) / max(g, kSmall) - static_cast<flow_float>(1.0), static_cast<flow_float>(0.0));
            const flow_float k_wf       = omega_w * mut_wall;  // = ω_w·μ_t,wall/ρ (mut_wall は運動学渦粘性 ν_t,wall)
            roK_wf[irep] = rho * max(k_wf, static_cast<flow_float>(0.0));
            // ω も第一内層ノードへピン (SU2 SetTurbVars_WF は第一点の k と ω の両方を設定する)。
            // 壁ノードのみのピンでは第一内点の ω がせん断駆動 P_ω=γρS² で暴騰し νt=ρk/ω が
            // 立たない (case/38 チャネルで実測: ω→3e5, νt→2.5μ)。ν_t(irep)=k_wf/ω_w=ν_t,wall と
            // なり mixing-length 級の近壁渦粘性が構造的に保証される。
            if (enableOmgWf == 1) roOmega_wf[irep] = rho * omega_w;
        }
    }

    utau_b[ib] = utau;
    ypls_b[ib] = utau * y / nu;

    // 断熱回復温度モデル (診断): T_aw = T_rep + r·U_t²/(2cp)。SU2 壁関数はこれを壁温として
    // 更新するが forge は未適用 (壁温 −230K の主因, notes/sessions/wall-temperature-three-way-analysis.md)。
    Taw_diag[ic] = T_arr[irep] + recovery * Ut * Ut / (static_cast<flow_float>(2.0) * max(cp_arr[irep], kSmall));
    if (isNode != 0 && Taw_Rep_Id != nullptr) Taw_Rep_Id[ic] = static_cast<flow_float>(irep);

    // node: 壁ノードに τ_w=ρu_τ² を格納 → viscousFlux_d が W-I エッジの接線応力を τ_w に再スケール
    // (SU2 AddTauWall)。代表点 ρ を使う (u_τ と整合)。cell では Tau_Wall を書かず -1 維持 (再スケール無効)。
    if (isNode != 0) Tau_Wall[ic] = rho * utau * utau;

    // エネルギー壁関数 (§6.5(g)): q_w = ρ cp u_τ (T_w − T_aw,rep)/T⁺(y⁺; Pr, Pr_t)。
    //   u_τ/y⁺/u⁺ は上の Reichardt Newton の収束値を再利用 (運動量経路はビット不変)。
    //   駆動温度差は回復温度 T_aw,rep (Taw_diag と同式)。q_w>0 = 壁→流体。
    //   低 y⁺ は純伝導式へ退避 (Kader の層流極限 T⁺→Pr·y⁺ と連続。WMLES §10.3 と同ガード)。
    if (energyWf != 0) {
        const flow_float cp_rep  = max(cp_arr[irep], kSmall);
        const flow_float lam_rep = max(thermCond_arr[irep], kSmall);
        const flow_float Pr      = vis_lam[irep] * cp_rep / lam_rep;
        const flow_float Taw     = T_arr[irep] + recovery * Ut * Ut / (static_cast<flow_float>(2.0) * cp_rep);
        const flow_float Tw      = Tw_b[ib];
        const flow_float ypq     = utau * y / nu;
        flow_float qw;
        if (ypq < static_cast<flow_float>(1.0e-1)) {
            qw = lam_rep * (Tw - Taw) / y;
        } else {
            const flow_float Tp = wallLaw_kader_tplus(Pr, ypq);
            qw = rho * cp_rep * utau * (Tw - Taw) / max(Tp, kSmall);
        }
        qwall_b[ib] = qw;
        if (isNode != 0) Qw_Wall[ic] = qw;   // AddQWall (W-I 熱流束置換, viscousFlux_d)
    }
}

}

void initWallFunctionPk_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1 && cfg.wallTreatmentSST == 1)) {
        return;
    }
    init_wf_pk_d<<<cuda_cfg.dimGrid_normalcell , cuda_cfg.dimBlock>>>(msh.nCells, var.c_d["wf_pk"], var.c_d["Tau_Wall"], var.c_d["roK_wf"], var.c_d["roOmega_wf"], var.c_d["Taw_diag"], var.c_d["Qw_Wall"], var.c_d["Taw_Wall_Flux"], var.c_d["Taw_Rep_Id"]);
}

namespace {
// node k Dirichlet の状態ピン: roK_wf>=0 の第一内層ノードで roK (と primitive k) を固定する。
// 従来この状態固定は point-implicit 更新カーネル (update_d.cu) にしかなく、explicit RK3 では
// res_roK ゼロ化だけが効いて k が初期値のまま凍結するバグがあった (2026-07-20, case/38 IDDES
// チャネルで発見: k=IC のまま νt が立たず近壁混合欠損)。applyBconds 位相 (両積分経路共通) で
// 状態を直接ピンして解消する。implicit の update 時ピンは残す (同値・冪等)。
__global__ void apply_roKwf_pin_d(geom_int nCells, flow_float* roK_wf, flow_float* roOmega_wf,
                                  flow_float* ro, flow_float* roK, flow_float* k,
                                  flow_float* roOmega, flow_float* omega)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float roi = max(ro[ic], static_cast<flow_float>(1.0e-30));
    const flow_float rw = roK_wf[ic];
    if (rw >= static_cast<flow_float>(0.0)) {
        roK[ic] = rw;
        k[ic]   = rw / roi;
    }
    const flow_float rww = roOmega_wf[ic];
    if (rww >= static_cast<flow_float>(0.0)) {
        roOmega[ic] = rww;
        omega[ic]   = rww / roi;
    }
}

// ============================================================================
// 【実験専用パッチ・恒久実装禁止】νt 床の判別実験
// (notes/sessions/channel-wmles-nutfloor-session-prompt.md, case/38 run_0016)。
// 解像できない帯 d < C·h_max (env FORGE_NUTFLOOR_DMAX [m] で指定, 未設定=完全不活性)
// の未ピン内部ノードにも壁法則整合の k/ω を書き、既存 wf ピン経路 (状態ピン +
// Bradshaw 迂回 νt=ρk/ω) に乗せる。k = ω_w·κ·d·u_τ により νt = ρκd·u_τ
// (mixing-length) を厳密に満たす。u_τ は env FORGE_NUTFLOOR_UTAU (チャネルは
// 体積力から既知 3.85)。第一内層 (wf 由来ピン) は上書きしない。
__global__ void apply_nutfloor_ext_d(geom_int nCells, geom_int* wall_flag,
                                     flow_float* wall_dist, flow_float* ro, flow_float* vis_lam,
                                     flow_float dmax, flow_float utau,
                                     flow_float* roK_wf, flow_float* roOmega_wf)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    if (wall_flag[ic] != 0) return;                                  // 壁ノード対象外
    if (roK_wf[ic] >= static_cast<flow_float>(0.0)) return;          // wf ピン済みは保持
    const flow_float d = wall_dist[ic];
    if (d <= static_cast<flow_float>(0.0) || d >= dmax) return;
    const flow_float rho = max(ro[ic], kSmall);
    const flow_float nu  = vis_lam[ic] / rho;
    const flow_float omega_vis = static_cast<flow_float>(6.0) * nu / (kSstBeta1 * d * d);
    const flow_float omega_log = utau / (sqrt(kSstBetaSt) * kKappa * d);
    const flow_float omega_w   = max(sqrt(omega_vis * omega_vis + omega_log * omega_log), kOmegaMin);
    const flow_float k_fl      = omega_w * kKappa * d * utau;        // νt=k/ω=κ·d·u_τ
    roK_wf[ic]     = rho * k_fl;
    roOmega_wf[ic] = rho * omega_w;
}
} // namespace

void applyNodeKwfStatePin_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1 && cfg.wallTreatmentSST == 1)) return;
    if (!(cfg.discretization == "node" && cfg.nodeKwfDirichlet == 1)) return;
    // 【実験専用・恒久実装禁止】νt 床 (env ゲート)。未設定なら旧経路と完全同一。
    static const flow_float nutfloorDmax = [] {
        const char* s = std::getenv("FORGE_NUTFLOOR_DMAX");
        return s ? static_cast<flow_float>(std::atof(s)) : static_cast<flow_float>(0.0);
    }();
    static const flow_float nutfloorUtau = [] {
        const char* s = std::getenv("FORGE_NUTFLOOR_UTAU");
        return s ? static_cast<flow_float>(std::atof(s)) : static_cast<flow_float>(3.85);
    }();
    static bool nutfloorLogged = false;
    if (!nutfloorLogged) {
        nutfloorLogged = true;
        if (nutfloorDmax > static_cast<flow_float>(0.0)) {
            printf("[EXPERIMENT] nut-floor active: dmax=%g m, utau=%g m/s (FORGE_NUTFLOOR_*)\n",
                   static_cast<double>(nutfloorDmax), static_cast<double>(nutfloorUtau));
        }
    }
    if (nutfloorDmax > static_cast<flow_float>(0.0) &&
        cfg.nodeOmegaWfDirichlet == 1 && msh.wall_flag_d != nullptr) {
        apply_nutfloor_ext_d<<<cuda_cfg.dimGrid_normalcell , cuda_cfg.dimBlock>>>(
            msh.nCells, msh.wall_flag_d, var.c_d["wall_dist"], var.c_d["ro"], var.c_d["vis_lam"],
            nutfloorDmax, nutfloorUtau, var.c_d["roK_wf"], var.c_d["roOmega_wf"]);
    }
    apply_roKwf_pin_d<<<cuda_cfg.dimGrid_normalcell , cuda_cfg.dimBlock>>>(
        msh.nCells, var.c_d["roK_wf"], var.c_d["roOmega_wf"], var.c_d["ro"],
        var.c_d["roK"], var.c_d["k"], var.c_d["roOmega"], var.c_d["omega"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void computeWallFrictionSST_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var)
{
    (void)msh;

    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1 && cfg.wallTreatmentSST == 1)) {
        return;
    }
    if (!(bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal")) {
        return;
    }
    if (bc.iPlanes.empty()) {
        return;
    }

    const int isNode = (cfg.discretization == "node" && msh.wall_flag_d != nullptr) ? 1 : 0;
    compute_wall_friction_sst_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
        static_cast<geom_int>(bc.iPlanes.size()),
        bc.map_bplane_plane_d,
        bc.map_bplane_cell_d,
        var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        var.c_d["ro"],
        var.c_d["vis_lam"],
        var.c_d["wall_dist"],
        var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"],
        bc.bvar_d["utau"],
        bc.bvar_d["ypls"],
        var.c_d["wf_pk"],
        // node: SU2 Normal_Neighbor 代表内部点選択用 (cell では未使用)
        isNode,
        msh.nNormalPlanes,
        msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.map_plane_cells_d,
        msh.wall_flag_d,
        var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.c_d["Tau_Wall"],
        // node k Dirichlet ゲート: wallTreatmentSST==1 かつ nodeKwfDirichlet==1 のときだけ roK_wf を書く (既定 OFF)。
        (cfg.wallTreatmentSST == 1 && cfg.nodeKwfDirichlet == 1) ? 1 : 0, var.c_d["roK_wf"],
        (cfg.wallTreatmentSST == 1 && cfg.nodeOmegaWfDirichlet == 1) ? 1 : 0, var.c_d["roOmega_wf"],
        var.c_d["T"], var.c_d["cp"],
        (flow_float)pow((double)cfg.prandtlLam, 1.0/3.0), var.c_d["Taw_diag"],
        var.c_d["Taw_Rep_Id"],
        // エネルギー壁関数 (§6.5(g)): 等温壁のみ。Tw は bvar Ts (固定入力)。
        (cfg.sstEnergyWallFunction == 1 && bc.bcondKind == "wall_isothermal") ? 1 : 0,
        bc.bvar_d["Ts"], var.c_d["thermCond"], cfg.turbulentPrandtl,
        bc.bvar_d["qwall"], var.c_d["Qw_Wall"]);
}
