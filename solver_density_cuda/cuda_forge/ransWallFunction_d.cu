#include "ransWallFunction_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"

// SST automatic wall treatment の摩擦速度 u_τ を Reichardt 普遍速度則の逆解きで求める。
// 理論は docs/turbulence/theory.md §6.5 (a)。粘性低層・バッファ・対数を 1 式で繋ぐため
// 第一セルの y⁺ 位置に依存せず妥当な u_τ を返す。
namespace {

constexpr flow_float kKappa    = static_cast<flow_float>(0.41);  // von Karman 定数
constexpr flow_float kSmall    = static_cast<flow_float>(1.0e-12);
constexpr int        kNewtonIt = 5;                              // Newton 反復回数

// Reichardt 則 u⁺(y⁺)
__device__ inline flow_float reichardt_uplus(flow_float yp)
{
    const flow_float lg = log(static_cast<flow_float>(1.0) + kKappa * yp) / kKappa;
    const flow_float e1 = exp(-yp / static_cast<flow_float>(11.0));
    const flow_float e2 = exp(-yp / static_cast<flow_float>(3.0));
    const flow_float tr = static_cast<flow_float>(7.8) *
        (static_cast<flow_float>(1.0) - e1 - (yp / static_cast<flow_float>(11.0)) * e2);
    return lg + tr;
}

// du⁺/dy⁺
__device__ inline flow_float reichardt_duplus_dyp(flow_float yp)
{
    const flow_float dlg = static_cast<flow_float>(1.0) / (static_cast<flow_float>(1.0) + kKappa * yp);
    const flow_float e1  = exp(-yp / static_cast<flow_float>(11.0));
    const flow_float e2  = exp(-yp / static_cast<flow_float>(3.0));
    // d/dy⁺ [1 - e1 - (y⁺/11) e2] = e1/11 - e2/11 + (y⁺/33) e2
    const flow_float dtr = static_cast<flow_float>(7.8) *
        (e1 / static_cast<flow_float>(11.0) - e2 / static_cast<flow_float>(11.0)
         + (yp / static_cast<flow_float>(33.0)) * e2);
    return dlg + dtr;
}

// wf_pk を全セル -1 (inactive) に初期化する。
__global__ void init_wf_pk_d(geom_int nCells, flow_float* wf_pk)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) wf_pk[ic] = static_cast<flow_float>(-1.0);
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
    geom_float* ccx, geom_float* ccy, geom_float* ccz)
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
        return;
    }

    // Newton 法で f(u_τ) = U_t/u_τ - u⁺(y⁺(u_τ)) = 0 を解く。初期値は粘性則。
    flow_float utau = sqrt(nu * Ut / y);
    for (int it = 0; it < kNewtonIt; ++it) {
        utau = max(utau, kSmall);
        const flow_float yp = utau * y / nu;
        const flow_float f  = Ut / utau - reichardt_uplus(yp);
        const flow_float df = -Ut / (utau * utau) - reichardt_duplus_dyp(yp) * (y / nu);
        if (fabs(df) < kSmall) break;
        utau -= f / df;
    }
    utau = max(utau, static_cast<flow_float>(0.0));

    // wall-function 生産 P_k = (τ_w - τ_visc)·∂U/∂y = ρu_τ⁴/ν · g(1-g), g = du⁺/dy⁺(y⁺₁)
    // (docs/turbulence §6.5(d))。粘性低層 g→1 で P_k→0 (壁解像極限を保つ)、対数層 g→1/(κy⁺) で
    // P_k→ρu_τ³/(κy)。これと ω ピン留めで k が平衡値 u_τ²/√β* に収束し runaway を断つ。
    const flow_float yp1 = utau * y / nu;
    const flow_float g   = reichardt_duplus_dyp(yp1);
    wf_pk[ic] = max(rho * utau * utau * utau * utau / nu * g * (static_cast<flow_float>(1.0) - g),
                    static_cast<flow_float>(0.0));

    utau_b[ib] = utau;
    ypls_b[ib] = utau * y / nu;
}

}

void initWallFunctionPk_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1 && cfg.wallTreatmentSST == 1)) {
        return;
    }
    init_wf_pk_d<<<cuda_cfg.dimGrid_normalcell , cuda_cfg.dimBlock>>>(msh.nCells, var.c_d["wf_pk"]);
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
        var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"]);
}
