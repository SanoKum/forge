#include "axisymmetricSource_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"


// 軸対称 (B 流儀) — 半径方向運動量の圧力・粘性ソース項を residual に加算する。
//
// 円柱座標で書いた半径運動量保存則を r 倍してから (r,z) 平面で積分すると、
// 圧力勾配項のうち r 倍では落ちないテルムとして
//
//     R_{rho u_r} = ∫∫ (p / r) · r dr dz = p_cell · A_planar
//
// が残る。粘性は $-\tau_{\theta\theta}$ を同じ面積で積分した source として加える
// (詳しくは docs/axisymmetric/theory.md)。
//
// res_roUy には -∑(F·S) が累積されており、時間進行は
//   roUy^{n+1} = roUy^n + coef_DT * res_roUy * dt / V
// で書かれる (timeIntegration_d.cu)。よってソース項は res_roUy に
// "+= (p - tau_theta_theta) * A_planar" として加算する。
__global__ void axisymmetricSource_d
(
    geom_int nCells,
    flow_float* P,
    flow_float* A_planar,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float* axisym_divU,
    flow_float* axisym_uy_over_r,
    flow_float* res_roUy,
    geom_int* axis_flag   // node-centered で軸上 CV を示す (nullptr 可)。軸上は半径方向ソースを課さない。
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        // 軸上 CV (R=0) は対称条件で半径方向の正味ソースが 0。node-centered では CV が軸に乗り、
        // P*A_planar/(ro*r_eff) が極端に stiff になって偽の u_r を駆動するため、ソースを課さない。
        if (axis_flag != nullptr && axis_flag[ic] == 1) return;
        const flow_float mu_total = vis_lam[ic] + vis_turb[ic];
        const flow_float tau_theta_theta =
            (flow_float)2.0 * mu_total * axisym_uy_over_r[ic]
            - (flow_float)(2.0 / 3.0) * mu_total * axisym_divU[ic];
        res_roUy[ic] += (P[ic] - tau_theta_theta) * A_planar[ic];
    }
}

// 軸対称ソース項の前計算。半径方向ひずみ u_r/r と発散 divU = ∂xUx + ∂yUy + u_r/r を
// セル毎に求める。これらは turbulent_viscosity (S_θθ) / ransSource / axisymmetricSource が
// 消費するため、それらより前に実行する必要がある。
// 有効半径 r_eff は revolved 体積 (per-radian, V = r·A_planar) を planar 面積で割って得る。
__global__ void axisymmetricGeomTerms_d
(
    geom_int nCells,
    flow_float* Uy,
    flow_float* volume,
    flow_float* A_planar,
    flow_float* dUxdx,
    flow_float* dUydy,
    flow_float* axisym_uy_over_r,
    flow_float* axisym_divU
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float area = A_planar[ic];
        const flow_float r_eff = (area > (flow_float)0.0) ? volume[ic] / area : (flow_float)0.0;

        if (r_eff > (flow_float)0.0) {
            axisym_uy_over_r[ic] = Uy[ic] / r_eff;
        } else {
            axisym_uy_over_r[ic] = (flow_float)0.0;
        }

        axisym_divU[ic] = dUxdx[ic] + dUydy[ic] + axisym_uy_over_r[ic];
    }
}


void axisymmetricSource_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isAxisymmetric != 1) return;

    // 軸 (R≤eps) でソース OFF。ソースは r 重みされない唯一の項で source/volume=P/r→∞ (r→0) が発散源。
    // node-centered 軸対称でのみ作用 (axis_flag = ノード R≈0)。SU2 の y<EPS ソース OFF に相当。
    geom_int* axis_flag = (cfg.discretization == "node" && cfg.isAxisymmetric == 1) ? msh.axis_flag_d : nullptr;
    axisymmetricSource_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["P"],
        var.c_d["A_planar"],
        var.c_d["vis_lam"],
        var.c_d["vis_turb"],
        var.c_d["axisym_divU"],
        var.c_d["axisym_uy_over_r"],
        var.c_d["res_roUy"],
        axis_flag
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

// node-centered 軸対称: 軸上 CV (R=0) で半径方向運動量を 0 に (軸対称の対称条件)。
// エネルギー整合のため、除去する半径方向運動エネルギー 0.5*roUy^2/ro を roe からも引く
// (これを怠ると内部エネルギー=圧力が過大になり発散する)。
__global__ void enforceAxisSymmetry_d
(
    geom_int nCells, geom_int* axis_flag,
    flow_float* ro, flow_float* roUy, flow_float* roe, flow_float* Uy
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells && axis_flag[ic] == 1) {
        const flow_float r = ro[ic];
        if (r > (flow_float)0.0) {
            roe[ic] -= (flow_float)0.5 * roUy[ic] * roUy[ic] / r; // 半径方向 KE を除去
        }
        roUy[ic] = (flow_float)0.0;
        Uy[ic]   = (flow_float)0.0;
    }
}

// SU2 流の軸対称対称面処理 (MARKER_SYM): 軸上 CV (R=0) で半径方向運動量の「残差」を 0 にする。
// res_roUy=0 にすると roUy は更新されず初期値 (=0, 対称条件) を保つ。状態やエネルギーをいじらず、
// flux+source を全て積んだ後に残差を射影するだけなので陰解法とも整合し発散しない (状態強制版・source
// 抑制版は発散したのに対し、これは SU2 と同じ「残差から法線運動量成分を除去」する方式)。
__global__ void zeroAxisRadialResidual_d
(
    geom_int nCells, geom_int* axis_flag, flow_float* res_roUy
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells && axis_flag[ic] == 1) {
        res_roUy[ic] = (flow_float)0.0;
    }
}

void zeroAxisRadialResidual_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isAxisymmetric != 1 || cfg.discretization != "node" || msh.axis_flag_d == nullptr) return;
    zeroAxisRadialResidual_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.axis_flag_d, var.c_d["res_roUy"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void enforceAxisSymmetry_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    // node-centered 軸対称でのみ作用 (軸上ノードが CV 中心になるのは node モード)。
    if (cfg.isAxisymmetric != 1 || cfg.discretization != "node" || msh.axis_flag_d == nullptr) return;

    enforceAxisSymmetry_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.axis_flag_d,
        var.c_d["ro"], var.c_d["roUy"], var.c_d["roe"], var.c_d["Uy"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

// node-centered 壁 Dirichlet (state 初期化): 壁ノードで速度を厳密に 0 にする。除去する運動エネルギー
// 0.5*(roUx²+roUy²+roUz²)/ro を roe からも引き、内部エネルギー=圧力が過大にならないようにする
// (enforceAxisSymmetry_d の 3 成分版)。IC 確定後に一度だけ呼ぶ。
__global__ void enforceWallNoSlip_d
(
    geom_int nCells, geom_int* wall_flag,
    flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz, flow_float* roe,
    flow_float* Ux, flow_float* Uy, flow_float* Uz
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells && wall_flag[ic] == 1) {
        const flow_float r = ro[ic];
        if (r > (flow_float)0.0) {
            roe[ic] -= (flow_float)0.5 * (roUx[ic]*roUx[ic] + roUy[ic]*roUy[ic] + roUz[ic]*roUz[ic]) / r;
        }
        roUx[ic] = (flow_float)0.0; roUy[ic] = (flow_float)0.0; roUz[ic] = (flow_float)0.0;
        Ux[ic]   = (flow_float)0.0; Uy[ic]   = (flow_float)0.0; Uz[ic]   = (flow_float)0.0;
    }
}

// node-centered 壁 Dirichlet (残差射影): 壁ノードで運動量「残差」を 0 にし、速度を 0 に保つ。
// 状態やエネルギーをいじらず flux+source 積算後に残差を射影するだけなので block-DPLUR と整合する
// (zeroAxisRadialResidual_d の 3 成分版)。assembleResidual 末尾 (軸射影の後) で毎反復呼ぶ。
__global__ void zeroWallMomentumResidual_d
(
    geom_int nCells, geom_int* wall_flag,
    flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells && wall_flag[ic] == 1) {
        res_roUx[ic] = (flow_float)0.0;
        res_roUy[ic] = (flow_float)0.0;
        res_roUz[ic] = (flow_float)0.0;
    }
}

void enforceWallNoSlip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || cfg.nodeWallDirichlet == 0 || msh.wall_flag_d == nullptr) return;
    enforceWallNoSlip_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.wall_flag_d,
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
        var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void zeroWallMomentumResidual_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || cfg.nodeWallDirichlet == 0 || msh.wall_flag_d == nullptr) return;
    zeroWallMomentumResidual_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.wall_flag_d,
        var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void axisymmetricGeomTerms_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isAxisymmetric != 1) return;

    axisymmetricGeomTerms_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["Uy"],
        var.c_d["volume"],
        var.c_d["A_planar"],
        var.c_d["dUxdx"],
        var.c_d["dUydy"],
        var.c_d["axisym_uy_over_r"],
        var.c_d["axisym_divU"]
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}
