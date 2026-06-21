#include "timeIntegration_d.cuh"
#include "lowMachPrecond_d.cuh"   // Phase 4: β (lowMachBeta)・c' (lowMachCprime) device ヘルパ
#include "cuda_forge/eos_jacobian_d.cuh"  // 一般EOS固有系 (eos_split_jacobian_general_closed)。precond 経路で使用
#include "cuda_forge/block_dplur_jacobian_d.cuh"  // block_dplur::accumulate_split_jacobian_cf (共有ヘッダ)
#include <cstdlib>

// 診断 (env FORGE_AXIS_DIAG_ALPHA, 既定 0=不変): 近軸 (r→0) の半径方向音響モード安定化。
// roUy 対角に α·A_planar·c を加える。revolved 軸面積 (r_f·S→0) が落とす半径音響スペクトル半径を
// planar 面積 A_pl で補うもの。near-axis radial-momentum 不安定 (case/28 TP) の最小修正診断。
__device__ float g_axisDiagAlpha = 0.0f;

namespace block_dplur {

template<typename T>
__device__ __forceinline__ void zero5(T* vec)
{
    #pragma unroll
    for (int i = 0; i < 5; ++i) {
        vec[i] = 0.0;
    }
}

template<typename T>
__device__ __forceinline__ void add_identity_scaled(T mat[5][5], T scale)
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        mat[row][row] += scale;
    }
}

// accumulate_split_jacobian_cf は cuda_forge/block_dplur_jacobian_d.cuh へ移設 (host/device 共有・Level3
// 単体テスト tools/test_eos_jacobian.cpp から呼ぶため)。block_dplur:: で参照する。

__device__ __forceinline__ void multiply_add_5x5_vec(
    const flow_float mat[5][5],
    const flow_float vec[5],
    flow_float out[5]
)
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        flow_float sum = 0.0;
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            sum += mat[row][col] * vec[col];
        }
        out[row] += sum;
    }
}

__device__ __forceinline__ void copy5(const flow_float* src, flow_float* dst)
{
    #pragma unroll
    for (int i = 0; i < 5; ++i) {
        dst[i] = src[i];
    }
}

template<typename T>
__device__ __forceinline__ void zero5x5(T mat[5][5])
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            mat[row][col] = 0.0;
        }
    }
}

__device__ __forceinline__ void add_scaled_5x5(flow_float dst[5][5], const flow_float src[5][5], flow_float scale)
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            dst[row][col] += scale * src[row][col];
        }
    }
}

__device__ __forceinline__ void build_abs_jacobian(
    flow_float gamma,
    flow_float nx,
    flow_float ny,
    flow_float nz,
    flow_float u,
    flow_float v,
    flow_float w,
    flow_float H,
    flow_float c,
    flow_float abs_jac[5][5]
)
{
    zero5x5(abs_jac);

    const flow_float sonic = max(c, static_cast<flow_float>(1.0e-8));
    const flow_float ek = 0.5 * (u * u + v * v + w * w);
    const flow_float U = u * nx + v * ny + w * nz;
    const flow_float chi = (gamma - 1.0) / sonic;
    const flow_float inv_sqrt2 = static_cast<flow_float>(0.7071067811865475244);

    flow_float lambda[5] = {
        fabs(U + sonic),
        fabs(U),
        fabs(U),
        fabs(U),
        fabs(U - sonic)
    };

    flow_float R[5][5];
    flow_float L[5][5];

    R[0][0] = inv_sqrt2 / sonic;             R[0][1] = ny / sonic;               R[0][2] = nz / sonic;               R[0][3] = nx / sonic;               R[0][4] = inv_sqrt2 / sonic;
    R[1][0] = (u / sonic + nx) * inv_sqrt2; R[1][1] = u * ny / sonic + nz;      R[1][2] = u * nz / sonic - ny;      R[1][3] = u * nx / sonic;           R[1][4] = (u / sonic - nx) * inv_sqrt2;
    R[2][0] = (v / sonic + ny) * inv_sqrt2; R[2][1] = v * ny / sonic;           R[2][2] = v * nz / sonic + nx;      R[2][3] = v * nx / sonic - nz;      R[2][4] = (v / sonic - ny) * inv_sqrt2;
    R[3][0] = (w / sonic + nz) * inv_sqrt2; R[3][1] = w * ny / sonic - nx;      R[3][2] = w * nz / sonic;           R[3][3] = w * nx / sonic + ny;      R[3][4] = (w / sonic - nz) * inv_sqrt2;
    R[4][0] = (ek / sonic + 1.0 / chi + U) * inv_sqrt2;
    R[4][1] = ek * ny / sonic + nz * u - nx * w;
    R[4][2] = ek * nz / sonic + nx * v - ny * u;
    R[4][3] = ek * nx / sonic + ny * w - nz * v;
    R[4][4] = (ek / sonic + 1.0 / chi - U) * inv_sqrt2;

    L[0][0] = (chi * ek - U) * inv_sqrt2;   L[0][1] = (-chi * u + nx) * inv_sqrt2; L[0][2] = (-chi * v + ny) * inv_sqrt2; L[0][3] = (-chi * w + nz) * inv_sqrt2; L[0][4] = chi * inv_sqrt2;
    L[1][0] = ny * (-chi * ek + sonic) - nz * u + nx * w;
    L[1][1] = ny * chi * u + nz;            L[1][2] = ny * chi * v;              L[1][3] = ny * chi * w - nx;         L[1][4] = -ny * chi;
    L[2][0] = nz * (-chi * ek + sonic) - nx * v + ny * u;
    L[2][1] = nz * chi * u - ny;            L[2][2] = nz * chi * v + nx;         L[2][3] = nz * chi * w;              L[2][4] = -nz * chi;
    L[3][0] = nx * (-chi * ek + sonic) - ny * w + nz * v;
    L[3][1] = nx * chi * u;                 L[3][2] = nx * chi * v - nz;         L[3][3] = nx * chi * w + ny;         L[3][4] = -nx * chi;
    L[4][0] = (chi * ek + U) * inv_sqrt2;   L[4][1] = (-chi * u - nx) * inv_sqrt2; L[4][2] = (-chi * v - ny) * inv_sqrt2; L[4][3] = (-chi * w - nz) * inv_sqrt2; L[4][4] = chi * inv_sqrt2;

    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            flow_float sum = 0.0;
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                sum += R[row][k] * lambda[k] * L[k][col];
            }
            abs_jac[row][col] = sum;
        }
    }
}

// LU-SGS の通量分割 Jacobian を同時構築する。
//   a_plus = A^+ = R Λ^+ L,  k_off = -A^- = R(-Λ^-)L = ½(|A|-A)。
// 検証済みの一般EOS閉形式 eos_split_jacobian_general_closed (eos_jacobian_d.cuh, Level1〜3 検証済) を流用し
// CPG/TP を統一 (旧来の CPG 専用べた打ち R/L=RL≠I の近似を廃止)。H は実全エンタルピー Ht[ic]、
// κ=γ−1、χ_eos=c²−κh。CPG では χ_eos≈0 で従来の CPG 固有系に簡約 (収束先は不変・厳密 Jacobian で僅かに向上)。
// double で組み立て float へ格納 (precond カーネルは元々 double 相当)。標準経路 accumulate_split_jacobian_cf は
// 別実装 (CPG ビット不変) のまま。
__device__ __forceinline__ void build_jacobian_split(
    flow_float gamma,
    flow_float nx, flow_float ny, flow_float nz,
    flow_float u, flow_float v, flow_float w,
    flow_float H, flow_float c,
    flow_float a_plus[5][5], flow_float k_off[5][5]
)
{
    const double sonic = (double)max(c, static_cast<flow_float>(1.0e-8));
    const double ek    = 0.5*((double)u*u + (double)v*v + (double)w*w);
    const double kappa = (double)gamma - 1.0;
    const double chi   = sonic*sonic - kappa*((double)H - ek);   // χ_eos = c²−κh (CPG で ≈0)
    double Ap[5][5], Am[5][5];
    eos_split_jacobian_general_closed((double)u,(double)v,(double)w, (double)nx,(double)ny,(double)nz,
                                      sonic, (double)H, kappa, chi, Ap, Am);   // Ap=A⁺, Am=A⁻
    #pragma unroll
    for (int i=0;i<5;++i)
        #pragma unroll
        for (int j=0;j<5;++j){ a_plus[i][j]=(flow_float)Ap[i][j]; k_off[i][j]=(flow_float)(-Am[i][j]); }
}

template<typename T>
__device__ __forceinline__ bool solve_5x5(T mat[5][5], T rhs[5], T sol[5])
{
    #pragma unroll
    for (int col = 0; col < 5; ++col) {
        int pivot = col;
        T pivot_abs = fabs(mat[col][col]);
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const T candidate = fabs(mat[row][col]);
            if (candidate > pivot_abs) {
                pivot = row;
                pivot_abs = candidate;
            }
        }

        if (pivot_abs < static_cast<T>(1.0e-20)) {
            zero5(sol);
            return false;
        }

        if (pivot != col) {
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                const T tmp = mat[col][k];
                mat[col][k] = mat[pivot][k];
                mat[pivot][k] = tmp;
            }
            const T rhs_tmp = rhs[col];
            rhs[col] = rhs[pivot];
            rhs[pivot] = rhs_tmp;
        }

        const T inv_pivot = static_cast<T>(1.0) / mat[col][col];
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const T factor = mat[row][col] * inv_pivot;
            mat[row][col] = 0.0;
            #pragma unroll
            for (int k = col + 1; k < 5; ++k) {
                mat[row][k] -= factor * mat[col][k];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    for (int row = 4; row >= 0; --row) {
        T sum = rhs[row];
        #pragma unroll
        for (int col = row + 1; col < 5; ++col) {
            sum -= mat[row][col] * sol[col];
        }
        sol[row] = sum / mat[row][row];
    }

    return true;
}

// 5×5 を 2 つの RHS について同時に解く (部分ピボット Gauss 消去を 1 回・float)。
// Phase 4 (lowMachPrecond=2) の Sherman-Morrison 解法で D0⁻¹b と D0⁻¹g を同じ分解で得るのに使う。
// D0 は物理ブロックで良条件 (既存 0/1 カーネルが float で解いているのと同形) ゆえ float で十分。
__device__ __forceinline__ bool solve_5x5_2rhs(flow_float mat[5][5],
                                               flow_float r1[5], flow_float r2[5],
                                               flow_float s1[5], flow_float s2[5])
{
    #pragma unroll
    for (int col = 0; col < 5; ++col) {
        int pivot = col;
        flow_float pivot_abs = fabs(mat[col][col]);
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const flow_float candidate = fabs(mat[row][col]);
            if (candidate > pivot_abs) { pivot = row; pivot_abs = candidate; }
        }

        if (pivot_abs < static_cast<flow_float>(1.0e-20)) {
            zero5(s1); zero5(s2);
            return false;
        }

        if (pivot != col) {
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                const flow_float tmp = mat[col][k]; mat[col][k] = mat[pivot][k]; mat[pivot][k] = tmp;
            }
            flow_float t = r1[col]; r1[col] = r1[pivot]; r1[pivot] = t;
            t = r2[col]; r2[col] = r2[pivot]; r2[pivot] = t;
        }

        const flow_float inv_pivot = static_cast<flow_float>(1.0) / mat[col][col];
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const flow_float factor = mat[row][col] * inv_pivot;
            mat[row][col] = 0.0;
            #pragma unroll
            for (int k = col + 1; k < 5; ++k) mat[row][k] -= factor * mat[col][k];
            r1[row] -= factor * r1[col];
            r2[row] -= factor * r2[col];
        }
    }

    for (int row = 4; row >= 0; --row) {
        flow_float sum1 = r1[row], sum2 = r2[row];
        #pragma unroll
        for (int col = row + 1; col < 5; ++col) { sum1 -= mat[row][col] * s1[col]; sum2 -= mat[row][col] * s2[col]; }
        const flow_float inv = static_cast<flow_float>(1.0) / mat[row][row];
        s1[row] = sum1 * inv;
        s2[row] = sum2 * inv;
    }

    return true;
}

__device__ __forceinline__ void load_block_vec(
    geom_int ic,
    flow_float* v0,
    flow_float* v1,
    flow_float* v2,
    flow_float* v3,
    flow_float* v4,
    flow_float out[5]
)
{
    out[0] = v0[ic];
    out[1] = v1[ic];
    out[2] = v2[ic];
    out[3] = v3[ic];
    out[4] = v4[ic];
}

__device__ __forceinline__ void store_block_vec(
    geom_int ic,
    const flow_float in[5],
    flow_float* v0,
    flow_float* v1,
    flow_float* v2,
    flow_float* v3,
    flow_float* v4
)
{
    v0[ic] = in[0];
    v1[ic] = in[1];
    v2[ic] = in[2];
    v3[ic] = in[3];
    v4[ic] = in[4];
}

}

__global__ void runge_kutta_exp_4th_d
// see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
// N: previous outer step , M: previous inner loop
( 
 int loop, 
 flow_float coef_DT,
 flow_float coef_Res,

 flow_float dt ,
 flow_float* dt_local ,

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol ,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe,

 flow_float* res_ro_m,
 flow_float* res_roUx_m,
 flow_float* res_roUy_m,
 flow_float* res_roUz_m,
 flow_float* res_roe_m

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    geom_float v = vol[ic];

    if (ic < nCells) {
        flow_float dt_l = dt_local[ic];

        if (loop == 0) {
            res_ro_m[ic]   = 0.0;
            res_roUx_m[ic] = 0.0;
            res_roUy_m[ic] = 0.0;
            res_roUz_m[ic] = 0.0;
            res_roe_m[ic]  = 0.0;
        }
        // N: previous outer step , M: previous inner loop
        res_ro_m[ic]   += coef_Res*res_ro[ic]*dt_l/v;
        res_roUx_m[ic] += coef_Res*res_roUx[ic]*dt_l/v;
        res_roUy_m[ic] += coef_Res*res_roUy[ic]*dt_l/v;
        res_roUz_m[ic] += coef_Res*res_roUz[ic]*dt_l/v;
        res_roe_m[ic]  += coef_Res*res_roe[ic]*dt_l/v;

        if (loop < 3) {
            ro[ic]   = roN[ic]   +coef_DT*res_ro[ic]*dt_l/v;
            roUx[ic] = roUxN[ic] +coef_DT*res_roUx[ic]*dt_l/v;
            roUy[ic] = roUyN[ic] +coef_DT*res_roUy[ic]*dt_l/v;
            roUz[ic] = roUzN[ic] +coef_DT*res_roUz[ic]*dt_l/v;
            roe[ic]  = roeN[ic]  +coef_DT*res_roe[ic]*dt_l/v;
        } else {
            res_ro[ic]   = res_ro_m[ic]*v/dt_l ;
            res_roUx[ic] = res_roUx_m[ic]*v/dt_l ;
            res_roUy[ic] = res_roUy_m[ic]*v/dt_l ;
            res_roUz[ic] = res_roUz_m[ic]*v/dt_l ;
            res_roe[ic]  = res_roe_m[ic]*v/dt_l ;

            ro[ic]   = roN[ic]  + res_ro_m[ic] ;
            roUx[ic] = roUxN[ic]+ res_roUx_m[ic] ;
            roUy[ic] = roUyN[ic]+ res_roUy_m[ic] ;
            roUz[ic] = roUzN[ic]+ res_roUz_m[ic] ;
            roe[ic]  = roeN[ic] + res_roe_m[ic] ;
        }
    }
}

__global__ void runge_kutta_exp_d
// see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
// N: previous outer step , M: previous inner loop
( 
 int loop,
 flow_float coef_N,
 flow_float coef_M,
 flow_float coef_Res,

 flow_float dt ,
 flow_float* dt_local ,

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol ,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        geom_float v = vol[ic];
        flow_float dt_l = dt_local[ic];
        // N: previous outer step , M: previous inner loop
        ro[ic]   = coef_N*roN[ic]   + coef_M*roM[ic]   + coef_Res*res_ro[ic]*dt_l/v;
        roUx[ic] = coef_N*roUxN[ic] + coef_M*roUxM[ic] + coef_Res*res_roUx[ic]*dt_l/v;
        roUy[ic] = coef_N*roUyN[ic] + coef_M*roUyM[ic] + coef_Res*res_roUy[ic]*dt_l/v;
        roUz[ic] = coef_N*roUzN[ic] + coef_M*roUzM[ic] + coef_Res*res_roUz[ic]*dt_l/v;
        roe[ic]  = coef_N*roeN[ic]  + coef_M*roeM[ic]  + coef_Res*res_roe[ic]*dt_l/v;
    }
}
__global__ void implicit_defect_correction_d
(
 int loop,
 flow_float dt,
 flow_float* dt_local,
 flow_float implicit_relax,
 flow_float* gamma_arr,   // per-cell γ (TP: γ_mix(T), CPG: cfg.gamma)。軸対称ソース Jacobian 用

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol,
 geom_int* plane_cells,
 geom_int* cell_planes_index,
 geom_int* cell_planes,
 geom_float* ccx,
 geom_float* ccy,
 geom_float* ccz,
 geom_float* sx,
 geom_float* sy,
 geom_float* sz,
 geom_float* ss,

 // variables
 flow_float* ro,
 flow_float* roUx,
 flow_float* roUy,
 flow_float* roUz,
 flow_float* roe,

 flow_float* roN,
 flow_float* roUxN,
 flow_float* roUyN,
 flow_float* roUzN,
 flow_float* roeN,

 flow_float laminar_visc,
 flow_float* vis_turb,
 flow_float* sonic,
 flow_float* Ux,
 flow_float* Uy,
 flow_float* Uz,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe,

 flow_float* corr_ro_old,
 flow_float* corr_roUx_old,
 flow_float* corr_roUy_old,
 flow_float* corr_roUz_old,
 flow_float* corr_roe_old,

 flow_float* corr_ro_new,
 flow_float* corr_roUx_new,
 flow_float* corr_roUy_new,
 flow_float* corr_roUz_new,
 flow_float* corr_roe_new,

 // 軸対称ソース項 Jacobian 用 (CPG/TP 共通)
 int isAxisymmetric,
 flow_float* A_planar,
 // dual-time 物理時間項の対角係数 a/Δt（定常は 0）
 flow_float unsteady_diag
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        geom_float v = vol[ic];
        flow_float dt_l = dt_local[ic];
        const flow_float density = max(ro[ic], static_cast<flow_float>(1.0e-30));
        const flow_float velocity_x = Ux[ic];
        const flow_float velocity_y = Uy[ic];
        const flow_float velocity_z = Uz[ic];
        const flow_float local_sonic = max(sonic[ic], static_cast<flow_float>(0.0));
        const flow_float nu_eff = (laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0))) / density;

        if (loop == 0) {
            corr_ro_old[ic] = 0.0;
            corr_roUx_old[ic] = 0.0;
            corr_roUy_old[ic] = 0.0;
            corr_roUz_old[ic] = 0.0;
            corr_roe_old[ic] = 0.0;
        }

        flow_float diag_face_sum = 0.0;
        flow_float neighbor_ro = 0.0;
        flow_float neighbor_roUx = 0.0;
        flow_float neighbor_roUy = 0.0;
        flow_float neighbor_roUz = 0.0;
        flow_float neighbor_roe = 0.0;
        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end = cell_planes_index[ic + 1];
        for (geom_int plane_offset = plane_begin; plane_offset < plane_end; ++plane_offset) {
            const geom_int ip = cell_planes[plane_offset];
            const flow_float face_area = max(ss[ip], static_cast<flow_float>(1.0e-30));
            const flow_float advective_radius = fabs(
                velocity_x * sx[ip] + velocity_y * sy[ip] + velocity_z * sz[ip]
            ) / face_area + local_sonic;

            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;
            const flow_float dcc_x = ccx[other_ic] - ccx[ic];
            const flow_float dcc_y = ccy[other_ic] - ccy[ic];
            const flow_float dcc_z = ccz[other_ic] - ccz[ic];
            const flow_float dcc = max(
                sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float dcc_dot_s = max(
                fabs(dcc_x * sx[ip] + dcc_y * sy[ip] + dcc_z * sz[ip]),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float delta = max(dcc * face_area * face_area / dcc_dot_s, static_cast<flow_float>(1.0e-30));
            // 粘性対角は residual の粘性流束 Jacobian (2ν·ss²/dcc_dot_s = 2ν·delta/dcc) と整合させる。
            // 旧 face_area·(2ν/delta) は面積 delta を長さ²扱いし ≈2ν に潰れ、軸対称近軸で r 重みを失い、
            // ゼロ面積(対称)面にもスプリアス項を載せていた (residual は ip<nNormalPlanes で軸面を除外)。
            const flow_float viscous_diag = static_cast<flow_float>(2.0) * nu_eff * delta / dcc;
            const flow_float face_coeff = face_area * advective_radius + viscous_diag;
            const flow_float offdiag_coeff = static_cast<flow_float>(0.5) * face_coeff;
            diag_face_sum += face_coeff;

            if (other_ic < nCells) {
                neighbor_ro += offdiag_coeff * corr_ro_old[other_ic];
                neighbor_roUx += offdiag_coeff * corr_roUx_old[other_ic];
                neighbor_roUy += offdiag_coeff * corr_roUy_old[other_ic];
                neighbor_roUz += offdiag_coeff * corr_roUz_old[other_ic];
                neighbor_roe += offdiag_coeff * corr_roe_old[other_ic];
            }
        }

        const flow_float diag = max(
            static_cast<flow_float>(v / max(dt_l, static_cast<flow_float>(1.0e-30)))
                + static_cast<flow_float>(v) * unsteady_diag + diag_face_sum,
            static_cast<flow_float>(1.0e-30)
        );
        const flow_float inv_diag = 1.0 / diag;

        // 軸対称フープ源 (res_roUy += (P-τθθ)·A_planar, axisymmetricSource_d.cu) の Jacobian 対角成分を
        // roUy 方程式の対角に陰化する。これが無いと近軸の剛フープ源 (∝1/r) が陽 (lagged) 扱いになり、
        // block-DPLUR が安定な CFL でも scalar が発散する (block は diag_block[2][2] で陰化済。
        // 切り分けは case/29 README / plan time_integration-scalar-dplur-axisym-source.md)。
        // block 版 diag_block[2][2] と同形: A_pl·((γ-1)u_y + 2μ/(ρ r_eff))。
        // γ は per-cell gamma_arr[ic] (TP=γ_mix(T) / CPG=cfg.gamma) を使い thermally perfect でも整合。
        // scalar 対角の正値性 (対角優位) を保つため非負側のみ加える (defect-correction の不動点は不変)。
        flow_float diag_roUy = diag;
        if (isAxisymmetric == 1) {
            const flow_float A_pl = max(A_planar[ic], static_cast<flow_float>(1.0e-30));
            const flow_float r_eff = max(static_cast<flow_float>(v) / A_pl, static_cast<flow_float>(1.0e-30));
            const flow_float g1 = gamma_arr[ic] - static_cast<flow_float>(1.0);
            const flow_float mu_total = laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0));
            const flow_float hoop = static_cast<flow_float>(2.0) * mu_total / (density * r_eff);
            const flow_float src_diag = A_pl * (g1 * velocity_y + hoop);
            diag_roUy = diag + max(src_diag, static_cast<flow_float>(0.0));
        }
        const flow_float inv_diag_roUy = 1.0 / diag_roUy;

        const flow_float jacobi_ro = (res_ro[ic] + neighbor_ro) * inv_diag;
        const flow_float jacobi_roUx = (res_roUx[ic] + neighbor_roUx) * inv_diag;
        const flow_float jacobi_roUy = (res_roUy[ic] + neighbor_roUy) * inv_diag_roUy;
        const flow_float jacobi_roUz = (res_roUz[ic] + neighbor_roUz) * inv_diag;
        const flow_float jacobi_roe = (res_roe[ic] + neighbor_roe) * inv_diag;

        corr_ro_new[ic] = implicit_relax * jacobi_ro;
        corr_roUx_new[ic] = implicit_relax * jacobi_roUx;
        corr_roUy_new[ic] = implicit_relax * jacobi_roUy;
        corr_roUz_new[ic] = implicit_relax * jacobi_roUz;
        corr_roe_new[ic] = implicit_relax * jacobi_roe;
    }
}

// 5×5 行列を複数保持しレジスタ消費が大きいため、block 上限を超えないよう __launch_bounds__ で
// 1 block あたりスレッド数を 128 に制限する（起動時の "too many resources" を回避）。
#define BLOCK_DPLUR_THREADS 128
// block-DPLUR の閉形式 FVS 版。線形 solve の内部精度を ST (float 既定 / double で軸対称近軸を根治) で
// テンプレート化。残差/状態 (flow_float=float) を ST へキャストして取り込み、R/L を作らず閉形式で
// diag/nbr を畳み、ST で in-place 5×5 solve、補正を float dq_new へ書戻す (混合精度 iterative refinement)。
// 詳細: .github/plans/precision-mixed-axisym.md。
template<typename ST>
__global__ void __launch_bounds__(BLOCK_DPLUR_THREADS) implicit_defect_correction_block_d
(
 int loop,
 flow_float dt,
 flow_float* dt_local,
 flow_float implicit_relax,
 flow_float* gamma_arr,   // per-cell γ (TP: γ_mix(T), CPG: cfg.gamma)。frozen-coefficient Jacobian 用
 int thermallyPerfect,    // 1: TP 固有系 (実 Ht・χ_eos=c²−κh, κ=γ−1), 0: CPG 閉形式 (従来・ビット不変)

 geom_int nCells_all , geom_int nCells,
 geom_float* vol,
 geom_int* plane_cells,
 geom_int* cell_planes_index,
 geom_int* cell_planes,
 geom_float* ccx,
 geom_float* ccy,
 geom_float* ccz,
 geom_float* sx,
 geom_float* sy,
 geom_float* sz,
 geom_float* ss,

 flow_float* ro,
 flow_float* roUx,
 flow_float* roUy,
 flow_float* roUz,
 flow_float* roe,

 flow_float laminar_visc,
 flow_float* vis_turb,
 flow_float* sonic,
 flow_float* Ux,
 flow_float* Uy,
 flow_float* Uz,
 flow_float* Ht,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe,

 flow_float* dq_old_0,
 flow_float* dq_old_1,
 flow_float* dq_old_2,
 flow_float* dq_old_3,
 flow_float* dq_old_4,

 flow_float* dq_new_0,
 flow_float* dq_new_1,
 flow_float* dq_new_2,
 flow_float* dq_new_3,
 flow_float* dq_new_4,

 flow_float* rhs_0,
 flow_float* rhs_1,
 flow_float* rhs_2,
 flow_float* rhs_3,
 flow_float* rhs_4,

 flow_float* diag_00, flow_float* diag_01, flow_float* diag_02, flow_float* diag_03, flow_float* diag_04,
 flow_float* diag_10, flow_float* diag_11, flow_float* diag_12, flow_float* diag_13, flow_float* diag_14,
 flow_float* diag_20, flow_float* diag_21, flow_float* diag_22, flow_float* diag_23, flow_float* diag_24,
 flow_float* diag_30, flow_float* diag_31, flow_float* diag_32, flow_float* diag_33, flow_float* diag_34,
 flow_float* diag_40, flow_float* diag_41, flow_float* diag_42, flow_float* diag_43, flow_float* diag_44,

 // 軸対称ソースヤコビアン用（isAxisymmetric==1 のときのみ使用）
 int isAxisymmetric,
 flow_float* A_planar,

 // dual-time 物理時間項の対角係数 a/Δt（定常は 0）
 flow_float unsteady_diag,

 // node-centered 軸対称: 軸上 CV で半径方向運動量 (roUy, index2) 行を decouple する (nullptr 可)。
 // SU2 流の対称面を Jacobian 内で課す = solve の外で状態を手術せず一貫して dq_roUy=0 を得る。
 geom_int* axis_flag,

 // node-centered 壁 no-slip: 壁ノードで運動量3行 (index1=roUx,2=roUy,3=roUz) を decouple する (nullptr 可)。
 // SU2 `DeleteValsRowi` 相当。残差射影だけでは block-DPLUR が壁運動量を連成したまま dq≠0 を返し速度 drift
 // するのを防ぐ。連続(行0)・エネルギー(行4)は保持。docs/discretization/implementation.md §7.2.1。
 geom_int* wall_flag
)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        // 状態/残差 (flow_float=float) を solve 精度 ST へキャストして取り込む (混合精度)。
        const ST gamma = static_cast<ST>(gamma_arr[ic]);   // 局所 γ
        const ST v = static_cast<ST>(vol[ic]);
        const ST dt_l = static_cast<ST>(dt_local[ic]);
        const ST density = max(static_cast<ST>(ro[ic]), static_cast<ST>(1.0e-30));
        const ST velocity_x = static_cast<ST>(Ux[ic]);
        const ST velocity_y = static_cast<ST>(Uy[ic]);
        const ST velocity_z = static_cast<ST>(Uz[ic]);
        const ST local_sonic = max(static_cast<ST>(sonic[ic]), static_cast<ST>(1.0e-8));
        const ST local_Ht = static_cast<ST>(Ht[ic]);   // 一般EOS固有系のエネルギー成分 (TP)
        const ST nu_eff = (static_cast<ST>(laminar_visc) + max(static_cast<ST>(vis_turb[ic]), static_cast<ST>(0.0))) / density;

        if (loop == 0) {
            dq_old_0[ic] = 0.0; dq_old_1[ic] = 0.0; dq_old_2[ic] = 0.0; dq_old_3[ic] = 0.0; dq_old_4[ic] = 0.0;
        }

        ST diag_block[5][5];
        block_dplur::zero5x5(diag_block);
        block_dplur::add_identity_scaled(diag_block, static_cast<ST>(v / max(dt_l, static_cast<ST>(1.0e-30))));
        // dual-time: 物理時間項 a·V/Δt を対角へ（定常は unsteady_diag==0）。
        block_dplur::add_identity_scaled(diag_block, v * static_cast<ST>(unsteady_diag));

        ST rhs[5] = {
            static_cast<ST>(res_ro[ic]),
            static_cast<ST>(res_roUx[ic]),
            static_cast<ST>(res_roUy[ic]),
            static_cast<ST>(res_roUz[ic]),
            static_cast<ST>(res_roe[ic])
        };
        ST neighbor_accum[5];
        block_dplur::zero5(neighbor_accum);

        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end = cell_planes_index[ic + 1];
        for (geom_int plane_offset = plane_begin; plane_offset < plane_end; ++plane_offset) {
            const geom_int ip = cell_planes[plane_offset];
            const ST face_area = max(static_cast<ST>(ss[ip]), static_cast<ST>(1.0e-30));
            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;

            // 格納法線 (sx,sy,sz) は ic0→ic1。ic が neighbor 側 (ic1) のとき符号反転。
            const ST nsign = (ic0 == ic) ? static_cast<ST>(1.0) : static_cast<ST>(-1.0);
            const ST nx = nsign * static_cast<ST>(sx[ip]) / face_area;
            const ST ny = nsign * static_cast<ST>(sy[ip]) / face_area;
            const ST nz = nsign * static_cast<ST>(sz[ip]) / face_area;

            // 閉形式 FVS (R/L を作らず rank-2 外積で A⁺S を diag・k_off·sdq を nbr へ)。
            const bool has_nbr = (other_ic < nCells);
            ST sdq[5] = {static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0)};
            if (has_nbr) {
                sdq[0] = face_area * static_cast<ST>(dq_old_0[other_ic]);
                sdq[1] = face_area * static_cast<ST>(dq_old_1[other_ic]);
                sdq[2] = face_area * static_cast<ST>(dq_old_2[other_ic]);
                sdq[3] = face_area * static_cast<ST>(dq_old_3[other_ic]);
                sdq[4] = face_area * static_cast<ST>(dq_old_4[other_ic]);
            }
            block_dplur::accumulate_split_jacobian_cf<ST>(
                gamma, nx, ny, nz, velocity_x, velocity_y, velocity_z,
                local_sonic, local_Ht, thermallyPerfect != 0,
                face_area, has_nbr, sdq, diag_block, neighbor_accum
            );

            const ST dcc_x = static_cast<ST>(ccx[other_ic]) - static_cast<ST>(ccx[ic]);
            const ST dcc_y = static_cast<ST>(ccy[other_ic]) - static_cast<ST>(ccy[ic]);
            const ST dcc_z = static_cast<ST>(ccz[other_ic]) - static_cast<ST>(ccz[ic]);
            const ST dcc = max(sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z), static_cast<ST>(1.0e-30));
            const ST dcc_dot_s = max(
                fabs(dcc_x * static_cast<ST>(sx[ip]) + dcc_y * static_cast<ST>(sy[ip]) + dcc_z * static_cast<ST>(sz[ip])),
                static_cast<ST>(1.0e-30)
            );
            const ST delta = max(dcc * face_area * face_area / dcc_dot_s, static_cast<ST>(1.0e-30));
            // 粘性対角は residual の粘性流束 Jacobian (2ν·ss²/dcc_dot_s = 2ν·delta/dcc) と整合させる
            // (旧 face_area·(2ν/delta) は ≈2ν に潰れ近軸で r 重み喪失・ゼロ面積面にスプリアス。詳細は site1 コメント)。
            const ST viscous_diag = static_cast<ST>(2.0) * nu_eff * delta / dcc;
            block_dplur::add_identity_scaled(diag_block, viscous_diag);
        }

        #pragma unroll
        for (int i = 0; i < 5; ++i) {
            rhs[i] += neighbor_accum[i];
        }

        // 軸対称ソース項のヤコビアンを対角ブロックに加える（roUy 行 = index 2）。詳細は実装ドキュメント参照。
        if (isAxisymmetric == 1) {
            const ST A_pl = static_cast<ST>(A_planar[ic]);
            const ST r_eff = max(v / max(A_pl, static_cast<ST>(1.0e-30)), static_cast<ST>(1.0e-30));
            const ST g1 = gamma - static_cast<ST>(1.0);
            const ST q2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
            const ST mu_total = static_cast<ST>(laminar_visc) + max(static_cast<ST>(vis_turb[ic]), static_cast<ST>(0.0));
            const ST hoop = static_cast<ST>(2.0) * mu_total / (density * r_eff);
            diag_block[2][0] += -A_pl * (static_cast<ST>(0.5)*g1*q2 + hoop * velocity_y);
            diag_block[2][1] += A_pl * (g1 * velocity_x);
            diag_block[2][2] += A_pl * (g1 * velocity_y + hoop);
            diag_block[2][3] += A_pl * (g1 * velocity_z);
            diag_block[2][4] += -A_pl * g1;
            // 診断: 近軸半径音響スペクトル半径 α·A_pl·c を roUy 対角に補う (FORGE_AXIS_DIAG_ALPHA>0 のみ)。
            diag_block[2][2] += static_cast<ST>(g_axisDiagAlpha) * A_pl * local_sonic;
        }

        // SU2 流の軸対称対称面 (MARKER_SYM) を Jacobian 内で課す: 軸上 CV で roUy 行 (index 2) を単位行に
        // 置換し rhs[2]=0 とする → solve が一貫して dq_roUy=0 を返す (他方程式は dq_roUy=0 固定で解かれる)。
        // solve の外で状態を手術する方式は block-DPLUR と非整合で Mach~1000 に発散したが、本方式は整合的。
        if (axis_flag != nullptr && axis_flag[ic] == 1) {
            for (int jj = 0; jj < 5; ++jj) diag_block[2][jj] = static_cast<ST>(0.0);
            diag_block[2][2] = static_cast<ST>(1.0);
            rhs[2] = static_cast<ST>(0.0);
        }

        // SU2 `DeleteValsRowi` 相当の壁 no-slip Dirichlet: 壁ノードで運動量3行 (index 1,2,3) を単位行に
        // 置換し rhs=0 → solve が一貫して dq_roUx=dq_roUy=dq_roUz=0 を返す。連続(0)・エネルギー(4)行は
        // 保持され ρ,ρe は保存式で発展、圧力は EOS が復元 (CPG/TP 共通)。残差射影だけでは block-DPLUR が
        // 壁運動量を連成し dq≠0 を返して壁速度が drift する問題を Jacobian 整合で根治する。
        if (wall_flag != nullptr && wall_flag[ic] == 1) {
            for (int row = 1; row <= 3; ++row) {
                for (int jj = 0; jj < 5; ++jj) diag_block[row][jj] = static_cast<ST>(0.0);
                diag_block[row][row] = static_cast<ST>(1.0);
                rhs[row] = static_cast<ST>(0.0);
            }
        }

        // diag_block を破壊して in-place で解く (solve_mat コピー排除)。
        ST correction[5] = {static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0)};
        const bool ok = block_dplur::solve_5x5(diag_block, rhs, correction);
        if (!ok) {
            block_dplur::zero5(correction);
        }

        const ST relax = static_cast<ST>(implicit_relax);
        #pragma unroll
        for (int i = 0; i < 5; ++i) {
            correction[i] *= relax;
        }

        // 古典 DPLUR: float dq_new へ書戻し。Q への commit は applyBlockImplicitCorrection。
        dq_new_0[ic] = static_cast<flow_float>(correction[0]);
        dq_new_1[ic] = static_cast<flow_float>(correction[1]);
        dq_new_2[ic] = static_cast<flow_float>(correction[2]);
        dq_new_3[ic] = static_cast<flow_float>(correction[3]);
        dq_new_4[ic] = static_cast<flow_float>(correction[4]);
        // rhs_** は読者が無いが従来通り残置（診断用）。
        rhs_0[ic] = static_cast<flow_float>(rhs[0]);
        rhs_1[ic] = static_cast<flow_float>(rhs[1]);
        rhs_2[ic] = static_cast<flow_float>(rhs[2]);
        rhs_3[ic] = static_cast<flow_float>(rhs[3]);
        rhs_4[ic] = static_cast<flow_float>(rhs[4]);
    }
}

// =============================================================================
// Phase 4 (a): 完全 Γ⁻¹A 低マッハ前処理の block DPLUR (lowMachPrecond=2)。
// 既存 implicit_defect_correction_block_d (lowMachPrecond 0/1) とは別カーネルにし、
// 0/1 経路のレジスタ・ビットを一切変えない。
// 保存形は (Γ_c·V/Δτ' + A_c)ΔQ = -R で、**前処理は擬似時間項 Γ_c のみ**。フラックス A_c は
// 既存と同じ物理厳密 FVS (a_plus=A_c⁺・k_off=-A_c⁻) をそのまま使う (非前処理が正しい)。収束は
// (Δτ'/V)Γ_c⁻¹A_c の固有値 λ' で一様に前処理され、スカラー Δτ'=cell/ρ' で効く。
//   - 擬似時間項: Γ_c·V/Δτ' (Δτ' は setDT 側で前処理スペクトル半径から拡大した dt_local)。
//   - フラックス: 物理 a_plus を対角・k_off を近傍 (既存 block と同一)。
//   - 物理 BDF 項 a·V/Δt·I は非前処理 (dual-time 所有)。
// **Sherman-Morrison 解法**: Γ_c=I+α g rᵀ がランク1なので D=D0+γ g rᵀ (D0=物理ブロック・良条件)。
//   D0 を float で 2 RHS 同時 (solve_5x5_2rhs) に解き、悪条件 ~1/β は分母スカラーのみ double に隔離。
//   FP64 を回避して 0/1 カーネルに近い速度。β=1 (超音速) で Γ_c=I・Δτ'=Δτ・フラックス同一ゆえ現行と解一致。
// 理論: docs/time_integration/theory.md「低マッハ前処理固有系」、計画 §5 Phase 4。
__global__ void __launch_bounds__(BLOCK_DPLUR_THREADS) implicit_defect_correction_block_precond_d
(
 int loop,
 flow_float dt,
 flow_float* dt_local,
 flow_float implicit_relax,
 flow_float* gamma_arr,   // per-cell γ (TP: γ_mix(T), CPG: cfg.gamma)
 flow_float precondEps,

 geom_int nCells_all , geom_int nCells,
 geom_float* vol,
 geom_int* plane_cells,
 geom_int* cell_planes_index,
 geom_int* cell_planes,
 geom_float* ccx, geom_float* ccy, geom_float* ccz,
 geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,

 flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz, flow_float* roe,

 flow_float laminar_visc,
 flow_float* vis_turb,
 flow_float* sonic,
 flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* Ht,

 flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe,

 flow_float* dq_old_0, flow_float* dq_old_1, flow_float* dq_old_2, flow_float* dq_old_3, flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,

 int isAxisymmetric,
 flow_float* A_planar,
 flow_float unsteady_diag
)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        const flow_float gamma = gamma_arr[ic];   // 局所 γ (TP の γ_mix。CPG は cfg.gamma で不変)
        const geom_float v = vol[ic];
        const flow_float dt_l = dt_local[ic];
        const flow_float density = max(ro[ic], static_cast<flow_float>(1.0e-30));
        const flow_float vx = Ux[ic];
        const flow_float vy = Uy[ic];
        const flow_float vz = Uz[ic];
        const flow_float local_sonic = max(sonic[ic], static_cast<flow_float>(1.0e-8));
        const flow_float local_enthalpy = max(Ht[ic], static_cast<flow_float>(1.0e-8));
        const flow_float nu_eff = (laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0))) / density;
        const flow_float velMag = sqrt(vx*vx + vy*vy + vz*vz);
        const flow_float beta = lowMachBeta(local_sonic, velMag, precondEps);

        if (loop == 0) {
            dq_old_0[ic] = 0.0; dq_old_1[ic] = 0.0; dq_old_2[ic] = 0.0; dq_old_3[ic] = 0.0; dq_old_4[ic] = 0.0;
        }

        // 対角ブロックは Γ_c=I+α g rᵀ がランク1ゆえ D = D0 + γ g rᵀ と書ける:
        //   D0 = V/Δτ'·I + a·V/Δt·I + Σ A_c⁺ S + 粘性 + 軸対称  (物理ブロック・良条件・float 可)
        //   γ g rᵀ = (V/Δτ')·α g rᵀ                             (Γ_c 前処理寄与、悪条件 ~1/β の源)
        // Sherman-Morrison: x = y - [γ(rᵀy)/(1+γ(rᵀz))] z,  y=D0⁻¹b, z=D0⁻¹g。
        // D0 を float で 2 RHS 同時に解き、悪条件は分母スカラー(double)に隔離 → FP64 を回避 (RTX 等で高速)。
        flow_float D0[5][5];
        block_dplur::zero5x5(D0);
        const flow_float v_over_dtau = static_cast<flow_float>(v / max(dt_l, static_cast<flow_float>(1.0e-30)));
        block_dplur::add_identity_scaled(D0, v_over_dtau);
        block_dplur::add_identity_scaled(D0, static_cast<flow_float>(v) * unsteady_diag);  // dual-time BDF (非前処理)

        flow_float b[5] = { res_ro[ic], res_roUx[ic], res_roUy[ic], res_roUz[ic], res_roe[ic] };
        flow_float nbr[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end = cell_planes_index[ic + 1];
        for (geom_int plane_offset = plane_begin; plane_offset < plane_end; ++plane_offset) {
            const geom_int ip = cell_planes[plane_offset];
            const flow_float face_area = max(ss[ip], static_cast<flow_float>(1.0e-30));
            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;
            const flow_float nsign = (ic0 == ic) ? static_cast<flow_float>(1.0) : static_cast<flow_float>(-1.0);
            const flow_float nx = nsign * sx[ip] / face_area;
            const flow_float ny = nsign * sy[ip] / face_area;
            const flow_float nz = nsign * sz[ip] / face_area;

            // フラックスは物理の厳密 FVS。前処理は時間項のみ (保存形 (Γ_c V/Δτ'+A_c)ΔQ=-R)。
            flow_float a_plus[5][5];
            flow_float k_off[5][5];
            block_dplur::build_jacobian_split(gamma, nx, ny, nz, vx, vy, vz,
                                              local_enthalpy, local_sonic, a_plus, k_off);
            block_dplur::add_scaled_5x5(D0, a_plus, face_area);

            const flow_float dcc_x = ccx[other_ic] - ccx[ic];
            const flow_float dcc_y = ccy[other_ic] - ccy[ic];
            const flow_float dcc_z = ccz[other_ic] - ccz[ic];
            const flow_float dcc = max(sqrt(dcc_x*dcc_x + dcc_y*dcc_y + dcc_z*dcc_z), static_cast<flow_float>(1.0e-30));
            const flow_float dcc_dot_s = max(fabs(dcc_x*sx[ip] + dcc_y*sy[ip] + dcc_z*sz[ip]), static_cast<flow_float>(1.0e-30));
            const flow_float delta = max(dcc * face_area * face_area / dcc_dot_s, static_cast<flow_float>(1.0e-30));
            // 粘性対角は residual の粘性流束 Jacobian (2ν·ss²/dcc_dot_s = 2ν·delta/dcc) と整合させる
            // (旧 face_area·(2ν/delta) は ≈2ν に潰れ近軸で r 重み喪失・ゼロ面積面にスプリアス。詳細は site1 コメント)。
            const flow_float viscous_diag = static_cast<flow_float>(2.0) * nu_eff * delta / dcc;
            block_dplur::add_identity_scaled(D0, viscous_diag);

            // 近傍 += k_off S ΔQ_nbr (= -A_c⁻ S ΔQ_nbr)。
            if (other_ic < nCells) {
                flow_float dqn[5];
                block_dplur::load_block_vec(other_ic, dq_old_0, dq_old_1, dq_old_2, dq_old_3, dq_old_4, dqn);
                #pragma unroll
                for (int i = 0; i < 5; ++i) dqn[i] *= face_area;
                block_dplur::multiply_add_5x5_vec(k_off, dqn, nbr);
            }
        }

        #pragma unroll
        for (int i = 0; i < 5; ++i) b[i] += nbr[i];

        // 軸対称ソースヤコビアン (物理ブロック D0 へ float で加算、既存 block と同式)。
        if (isAxisymmetric == 1) {
            const flow_float A_pl = A_planar[ic];
            const flow_float r_eff = max(v / max(A_pl, static_cast<flow_float>(1.0e-30)), static_cast<flow_float>(1.0e-30));
            const flow_float g1 = gamma - static_cast<flow_float>(1.0);
            const flow_float q2 = vx*vx + vy*vy + vz*vz;
            const flow_float mu_total = laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0));
            const flow_float hoop = static_cast<flow_float>(2.0) * mu_total / (density * r_eff);
            D0[2][0] += -A_pl * (static_cast<flow_float>(0.5)*g1*q2 + hoop*vy);
            D0[2][1] +=  A_pl * (g1*vx);
            D0[2][2] +=  A_pl * (g1*vy + hoop);
            D0[2][3] +=  A_pl * (g1*vz);
            D0[2][4] += -A_pl * g1;
        }

        // Γ_c のランク1寄与: g=(1,u,v,w,H_t), r=∂p/∂Q=(χ_eos+κek,-κu,-κv,-κw,κ), γ=(V/Δτ')·(1-β)/(βc²)。
        // CPG では H_t=c²/(γ-1)+ek・χ_eos=0 で従来式に簡約 (ビット不変)。TP は実 H_t(=local_enthalpy) と
        // χ_eos=c²−κh を使う (build_jacobian_split と同じ一般EOS整合)。κ=γ-1。
        const flow_float ek = static_cast<flow_float>(0.5) * velMag * velMag;
        const flow_float gm1 = gamma - static_cast<flow_float>(1.0);
        const flow_float Htot = local_enthalpy;   // 実 Ht[ic] (CPG/TP 統一。CPG も Ht[ic]=ek+c²/(γ-1))
        const flow_float chi_eos = local_sonic*local_sonic - gm1*(local_enthalpy - ek);  // c²−κh (CPG で ≈0)
        flow_float gvec[5] = { static_cast<flow_float>(1.0), vx, vy, vz, Htot };
        const flow_float rvec[5] = { chi_eos + gm1*ek, -gm1*vx, -gm1*vy, -gm1*vz, gm1 };
        const double dbeta = static_cast<double>(beta);
        const double alpha = (1.0 - dbeta) / (dbeta * static_cast<double>(local_sonic) * static_cast<double>(local_sonic));
        const double gam = static_cast<double>(v_over_dtau) * alpha;   // = (V/Δτ')·α

        // D0 を float で 2 RHS 同時に解く: y=D0⁻¹b, z=D0⁻¹g。
        flow_float y[5], z[5];
        const bool ok = block_dplur::solve_5x5_2rhs(D0, b, gvec, y, z);

        // Sherman-Morrison スカラー (悪条件 1/β はここだけ double): x = y - [γ(rᵀy)/(1+γ(rᵀz))] z。
        double ry = 0.0, rz = 0.0;
        #pragma unroll
        for (int i = 0; i < 5; ++i) { ry += static_cast<double>(rvec[i]) * y[i]; rz += static_cast<double>(rvec[i]) * z[i]; }
        const double denom = 1.0 + gam * rz;
        const double sfac = (fabs(denom) > 1.0e-300) ? gam * ry / denom : 0.0;

        flow_float correction[5];
        #pragma unroll
        for (int i = 0; i < 5; ++i)
            correction[i] = ok ? static_cast<flow_float>((static_cast<double>(y[i]) - sfac * static_cast<double>(z[i]))
                                                         * static_cast<double>(implicit_relax))
                               : static_cast<flow_float>(0.0);

        block_dplur::store_block_vec(ic, correction, dq_new_0, dq_new_1, dq_new_2, dq_new_3, dq_new_4);
    }
}

// block DPLUR の sweep 間バッファ入れ替え。ドライバ側から各 sweep 後に明示的に呼ぶ
// （旧実装は wrapper 内部で暗黙に swap していたが、古典 DPLUR では制御フローを明示化する）。
void swapBlockImplicitCorrectionBuffers(variables& var)
{
    std::swap(var.c_d["dq_block_old_0"], var.c_d["dq_block_new_0"]);
    std::swap(var.c_d["dq_block_old_1"], var.c_d["dq_block_new_1"]);
    std::swap(var.c_d["dq_block_old_2"], var.c_d["dq_block_new_2"]);
    std::swap(var.c_d["dq_block_old_3"], var.c_d["dq_block_new_3"]);
    std::swap(var.c_d["dq_block_old_4"], var.c_d["dq_block_new_4"]);
}

// scalar 対角版 (blockDPLUR==0) の sweep 間バッファ入れ替え。block 版と同様にドライバ側から呼ぶ。
void swapScalarImplicitCorrectionBuffers(variables& var)
{
    std::swap(var.c_d["dq_ro_old"],   var.c_d["dq_ro_new"]);
    std::swap(var.c_d["dq_roUx_old"], var.c_d["dq_roUx_new"]);
    std::swap(var.c_d["dq_roUy_old"], var.c_d["dq_roUy_new"]);
    std::swap(var.c_d["dq_roUz_old"], var.c_d["dq_roUz_new"]);
    std::swap(var.c_d["dq_roe_old"],  var.c_d["dq_roe_new"]);
}

//TODO __global__ void runge_kutta_dual_explicit_d
//TODO // see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
//TODO // N: previous outer step , M: previous inner loop
//TODO ( 
//TODO  geom_int dt ,
//TODO 
//TODO  // mesh structure
//TODO  geom_int nCells_all , geom_int nCells,
//TODO  geom_float* vol ,
//TODO 
//TODO  // variables
//TODO  flow_float* ro  ,
//TODO  flow_float* roUx  ,
//TODO  flow_float* roUy  ,
//TODO  flow_float* roUz  ,
//TODO  flow_float* roe  ,
//TODO 
//TODO  flow_float* roN ,
//TODO  flow_float* roUxN ,
//TODO  flow_float* roUyN ,
//TODO  flow_float* roUzN ,
//TODO  flow_float* roeN ,
//TODO 
//TODO  flow_float* roM ,
//TODO  flow_float* roUxM ,
//TODO  flow_float* roUyM ,
//TODO  flow_float* roUzM ,
//TODO  flow_float* roeM ,
//TODO 
//TODO  flow_float* res_ro,
//TODO  flow_float* res_roUx,
//TODO  flow_float* res_roUy,
//TODO  flow_float* res_roUz,
//TODO  flow_float* res_roe,
//TODO 
//TODO  flow_float* res_ro_dual ,
//TODO  flow_float* res_roUx_dual ,
//TODO  flow_float* res_roUy_dual ,
//TODO  flow_float* res_roUz_dual ,
//TODO  flow_float* res_roe_dual 
//TODO 
//TODO )
//TODO {
//TODO     geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
//TODO 
//TODO     geom_float v = vol[ic];
//TODO 
//TODO     if (ic < nCells) {
//TODO         // N: previous outer step , M: previous inner loop
//TODO         res_ro_dual[ic]   = -(ro[ic]-roN[ic])*v/dt     + res_ro[ic];
//TODO         res_roUx_dual[ic] = -(roUx[ic]-roUxN[ic])*v/dt + res_roUx[ic];
//TODO         res_roUy_dual[ic] = -(roUy[ic]-roUyN[ic])*v/dt + res_roUy[ic];
//TODO         res_roUz_dual[ic] = -(roUz[ic]-roUzN[ic])*v/dt + res_roUz[ic];
//TODO         res_roe_dual[ic]  = -(roe[ic]-roeN[ic])*v/dt   + res_roe[ic];
//TODO     }
//TODO     __syncthreads();
//TODO }

void timeIntegration_d_wrapper(int loop , solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    // 診断 near-axis 安定化係数を env から 1 度だけ device へ設定 (既定 0 = 不変)。
    static bool s_axisAlphaInit = false;
    if (!s_axisAlphaInit) {
        float a = 0.0f;
        if (const char* e = getenv("FORGE_AXIS_DIAG_ALPHA")) a = static_cast<float>(atof(e));
        cudaMemcpyToSymbol(g_axisDiagAlpha, &a, sizeof(float));
        s_axisAlphaInit = true;
    }
    if (cfg.timeIntegration == 4) { // 4th order runge kutta
        runge_kutta_exp_4th_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            loop, 
            cfg.coef_DT_4thRunge[loop],
            cfg.coef_Res_4thRunge[loop],
            cfg.dt ,
            var.c_d["dt_local"],

            // mesh structure
            msh.nCells_all , msh.nCells ,
            var.c_d["volume"],

            // basic variables
            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] ,
            var.c_d["res_ro_m"], var.c_d["res_roUx_m"], var.c_d["res_roUy_m"], var.c_d["res_roUz_m"] , var.c_d["res_roe_m"] 
        ) ;

    } else if (cfg.timeIntegration == 1 or cfg.timeIntegration == 3) { // explicit
        runge_kutta_exp_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            loop,
            cfg.coef_N[loop],
            cfg.coef_M[loop],
            cfg.coef_Res[loop],
            cfg.dt , 
            var.c_d["dt_local"],

            // mesh structure
            msh.nCells_all , msh.nCells ,
            var.c_d["volume"],

            // basic variables
            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] 
        ) ;
    } else if (cfg.timeIntegration == 11) { // implicit defect-correction with diagonal Jacobian approximation
        if (cfg.blockDPLUR == 1) {
            // レジスタ過多のため専用の小さい block サイズで起動（__launch_bounds__ と整合）。
            const int block_threads = BLOCK_DPLUR_THREADS;
            const int block_grid = (msh.nCells_all + block_threads - 1) / block_threads;
            if (cfg.lowMachPrecond == 2) {
              // Phase 4: 完全 Γ⁻¹A 前処理の倍精度カーネル (dt_local は前処理 Δτ' に拡大済)。
              implicit_defect_correction_block_precond_d<<<block_grid , block_threads>>>(
                loop, cfg.dt, var.c_d["dt_local"], cfg.implicitRelax, var.c_d["gamma"], cfg.precondEps,
                msh.nCells_all, msh.nCells, var.c_d["volume"],
                msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d,
                var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
                var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
                cfg.visc, var.c_d["vis_turb"], var.c_d["sonic"],
                var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["Ht"],
                var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"],
                var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
                var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"],
                cfg.isAxisymmetric,
                (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"],
                cfg.unsteadyDiagCoef
              );
            } else {
            // implicitSolvePrecision: 0=float (既定・高速), 1=double (軸対称近軸の根治, 遅い)。
            // 同じテンプレートカーネルを ST=float/double で起動。引数は共通 (FORGE_BDPLUR_ARGS)。
            #define FORGE_BDPLUR_ARGS \
                loop, cfg.dt, var.c_d["dt_local"], cfg.implicitRelax, var.c_d["gamma"], \
                (cfg.thermalMethod == 2 ? 1 : 0), \
                msh.nCells_all, msh.nCells, var.c_d["volume"], \
                msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d, \
                var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"], \
                var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"], \
                var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"], \
                cfg.visc, var.c_d["vis_turb"], var.c_d["sonic"], \
                var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["Ht"], \
                var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"], \
                var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"], \
                var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"], \
                var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"], \
                var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"], \
                var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"], \
                var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"], \
                var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"], \
                var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"], \
                cfg.isAxisymmetric, (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"], cfg.unsteadyDiagCoef, \
                nullptr,  /* axis_flag: in-Jacobian roUy decouple は corner を直さず (corner は多方程式不良) 既定無効 */ \
                ((cfg.discretization == "node" && cfg.nodeWallDirichlet == 1) ? msh.wall_flag_d : nullptr)  /* wall_flag: 壁運動量3行 decouple */
            if (cfg.implicitSolvePrecision == 1)
                implicit_defect_correction_block_d<double><<<block_grid , block_threads>>>(FORGE_BDPLUR_ARGS);
            else
                implicit_defect_correction_block_d<float><<<block_grid , block_threads>>>(FORGE_BDPLUR_ARGS);
            #undef FORGE_BDPLUR_ARGS
            }
        } else {
            implicit_defect_correction_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
                loop,
                cfg.dt,
                var.c_d["dt_local"],
                cfg.implicitRelax,
                var.c_d["gamma"],
                msh.nCells_all,
                msh.nCells,
                var.c_d["volume"],
                msh.map_plane_cells_d,
                msh.map_cell_planes_index_d,
                msh.map_cell_planes_d,
                var.c_d["ccx"],
                var.c_d["ccy"],
                var.c_d["ccz"],
                var.p_d["sx"],
                var.p_d["sy"],
                var.p_d["sz"],
                var.p_d["ss"],
                var.c_d["ro"],
                var.c_d["roUx"],
                var.c_d["roUy"],
                var.c_d["roUz"],
                var.c_d["roe"],
                var.c_d["roN"],
                var.c_d["roUxN"],
                var.c_d["roUyN"],
                var.c_d["roUzN"],
                var.c_d["roeN"],
                cfg.visc,
                var.c_d["vis_turb"],
                var.c_d["sonic"],
                var.c_d["Ux"],
                var.c_d["Uy"],
                var.c_d["Uz"],
                var.c_d["res_ro"],
                var.c_d["res_roUx"],
                var.c_d["res_roUy"],
                var.c_d["res_roUz"],
                var.c_d["res_roe"],
                var.c_d["dq_ro_old"],
                var.c_d["dq_roUx_old"],
                var.c_d["dq_roUy_old"],
                var.c_d["dq_roUz_old"],
                var.c_d["dq_roe_old"],
                var.c_d["dq_ro_new"],
                var.c_d["dq_roUx_new"],
                var.c_d["dq_roUy_new"],
                var.c_d["dq_roUz_new"],
                var.c_d["dq_roe_new"],
                cfg.isAxisymmetric,
                (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"],
                cfg.unsteadyDiagCoef
            );
        }
        // 古典 DPLUR: buffer swap と Q への commit はドライバ側 (main.cpp blockDPLURSolve /
        // applyBlockImplicitCorrection) で明示的に行う。ここでは sweep カーネルの起動のみ。
//TODO    } else if (cfg.timeIntegration == 10) { // implicit (m-time stepping & explicit scheme)
//TODO        runge_kutta_dual_explicit_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
//TODO            cfg.dt , 
//TODO
//TODO            // mesh structure
//TODO            msh.nCells_all , msh.nCells ,
//TODO            var.c_d["volume"],
//TODO
//TODO            // basic variables
//TODO            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
//TODO            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
//TODO            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
//TODO            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] ,
//TODO            var.c_d["res_ro_m"], var.c_d["res_roUx_m"], var.c_d["res_roUy_m"], var.c_d["res_roUz_m"] , var.c_d["res_roe_m"] 
//TODO        ) ;
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

}