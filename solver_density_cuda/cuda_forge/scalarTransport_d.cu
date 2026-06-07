#include "scalarTransport_d.cuh"

#include <array>

namespace {

struct ScalarTransportDesc {
    flow_float* phi;
    flow_float* rho_phi;
    flow_float* rho_phi_N;
    flow_float* rho_phi_M;
    flow_float* res_rho_phi;
    flow_float* res_rho_phi_m;
    flow_float* src_jac;       // SST 消散ヤコビアン対角（k:β*ω, ω:2βω）。陽解法 RK の point-implicit 減衰に使用
    flow_float floor;          // realizability 下限（ρk≥0, ρω>1e-20）。陰解法 applySSTPointImplicit と整合
    flow_float sigma;
};

__global__ void scalar_advection_first_order_d(
    geom_int nCells,
    geom_int nNormalHaloPlanes,
    geom_int* normal_halo_planes,
    geom_int* plane_cells,
    flow_float* phi,
    flow_float* massflux,
    flow_float* res_rho_phi)
{
    geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;

    if (ih < nNormalHaloPlanes) {
        const geom_int ip = normal_halo_planes[ih];
        const geom_int ic0 = plane_cells[2 * ip + 0];
        const geom_int ic1 = plane_cells[2 * ip + 1];

        const flow_float mdot = massflux[ip];

        const flow_float phi_upwind = (mdot >= 0.0) ? phi[ic0] : phi[ic1];
        const flow_float flux = mdot * phi_upwind;

        if (ic0 < nCells) {
            atomicAdd(&res_rho_phi[ic0], -flux);
        }
        if (ic1 < nCells) {
            atomicAdd(&res_rho_phi[ic1], flux);
        }
    }
}

__global__ void scalar_diffusion_first_order_d(
    geom_int nCells,
    geom_int nNormalHaloPlanes,
    geom_int* normal_halo_planes,
    geom_int* plane_cells,
    geom_float* ccx,
    geom_float* ccy,
    geom_float* ccz,
    geom_float* fx,
    geom_float* sx,
    geom_float* sy,
    geom_float* sz,
    geom_float* ss,
    flow_float* phi,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float sigma,
    flow_float* res_rho_phi)
{
    geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;

    if (ih < nNormalHaloPlanes) {
        const geom_int ip = normal_halo_planes[ih];
        const geom_int ic0 = plane_cells[2 * ip + 0];
        const geom_int ic1 = plane_cells[2 * ip + 1];

        const geom_float f = fx[ip];
        const geom_float sxx = sx[ip];
        const geom_float syy = sy[ip];
        const geom_float szz = sz[ip];
        const geom_float sss = ss[ip];

        const flow_float dcc_x = ccx[ic1] - ccx[ic0];
        const flow_float dcc_y = ccy[ic1] - ccy[ic0];
        const flow_float dcc_z = ccz[ic1] - ccz[ic0];
        const flow_float dcc = sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z);

        const flow_float denom = dcc_x * sxx + dcc_y * syy + dcc_z * szz;
        const flow_float safe_denom = (abs(denom) < 1.0e-12) ? ((denom >= 0.0) ? 1.0e-12 : -1.0e-12) : denom;
        const flow_float delta = dcc * sss * sss / safe_denom;

        const flow_float mu0 = vis_lam[ic0] + sigma * max(vis_turb[ic0], static_cast<flow_float>(0.0));
        const flow_float mu1 = vis_lam[ic1] + sigma * max(vis_turb[ic1], static_cast<flow_float>(0.0));
        const flow_float mu_face = f * mu0 + (1.0 - f) * mu1;

        const flow_float dphi = phi[ic1] - phi[ic0];
        const flow_float flux = mu_face * (dphi / dcc) * delta;

        if (ic0 < nCells) {
            atomicAdd(&res_rho_phi[ic0], flux);
        }
        if (ic1 < nCells) {
            atomicAdd(&res_rho_phi[ic1], -flux);
        }
    }
}

__global__ void runge_kutta_exp_scalar_4th_d(
    int loop,
    flow_float coef_DT,
    flow_float coef_Res,
    flow_float* dt_local,
    geom_int nCells,
    geom_float* vol,
    flow_float* rho_phi,
    flow_float* rho_phi_N,
    flow_float* rho_phi_M,
    flow_float* res_rho_phi,
    flow_float* res_rho_phi_m)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        const flow_float dt_l = dt_local[ic];
        const geom_float v = vol[ic];

        if (loop == 0) {
            res_rho_phi_m[ic] = 0.0;
        }

        res_rho_phi_m[ic] += coef_Res * res_rho_phi[ic] * dt_l / v;

        if (loop < 3) {
            rho_phi[ic] = rho_phi_N[ic] + coef_DT * res_rho_phi[ic] * dt_l / v;
        } else {
            res_rho_phi[ic] = res_rho_phi_m[ic] * v / dt_l;
            rho_phi[ic] = rho_phi_N[ic] + res_rho_phi_m[ic];
        }
    }
}

__global__ void runge_kutta_exp_scalar_d(
    flow_float coef_N,
    flow_float coef_M,
    flow_float coef_Res,
    flow_float* dt_local,
    geom_int nCells,
    geom_float* vol,
    flow_float* rho_phi,
    flow_float* rho_phi_N,
    flow_float* rho_phi_M,
    flow_float* res_rho_phi,
    flow_float* src_jac,
    flow_float floor)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        const flow_float dt_l = dt_local[ic];
        const geom_float v = vol[ic];
        // SST 消散項の stiff 性を point-implicit で減衰（src_jac=∂D/∂(ρφ)≥0、乱流なしでは 0 で従来と一致）。
        // 減衰係数は block 陰解法対角 D_φ=V/Δτ+V·src_jac に Δτ/V を掛けた形と整合。理論は docs/time_integration。
        const flow_float fac = static_cast<flow_float>(1.0) + coef_Res * dt_l * src_jac[ic];
        const flow_float updated = coef_N * rho_phi_N[ic] + coef_M * rho_phi_M[ic]
                    + (coef_Res * res_rho_phi[ic] * dt_l / v) / fac;
        // realizability 下限（ρk≥0, ρω>1e-20）。applySSTPointImplicit と整合。乱流なしでは floor 以下にならず無影響。
        rho_phi[ic] = max(updated, floor);
    }
}

inline bool scalarTransportEnabled(const solverConfig& cfg)
{
    return cfg.LESorRANS == 2 && cfg.RANSmodel == 1;
}

__global__ void calc_scalar_gradient_face_d(
    geom_int nPlanes,
    geom_int* plane_cells,
    geom_float* fx,
    geom_float* sx,
    geom_float* sy,
    geom_float* sz,
    flow_float* k,
    flow_float* omega,
    flow_float* dKdx,
    flow_float* dKdy,
    flow_float* dKdz,
    flow_float* dOmegadx,
    flow_float* dOmegady,
    flow_float* dOmegadz)
{
    geom_int ip = blockDim.x * blockIdx.x + threadIdx.x;
    if (ip >= nPlanes) return;

    const geom_int ic0 = plane_cells[2 * ip + 0];
    const geom_int ic1 = plane_cells[2 * ip + 1];

    const geom_float f   = fx[ip];
    const flow_float kf  = f * k[ic0]     + (1.0 - f) * k[ic1];
    const flow_float wf  = f * omega[ic0] + (1.0 - f) * omega[ic1];
    const geom_float sxx = sx[ip];
    const geom_float syy = sy[ip];
    const geom_float szz = sz[ip];

    atomicAdd(&dKdx[ic0],     sxx * kf);
    atomicAdd(&dKdy[ic0],     syy * kf);
    atomicAdd(&dKdz[ic0],     szz * kf);
    atomicAdd(&dOmegadx[ic0], sxx * wf);
    atomicAdd(&dOmegady[ic0], syy * wf);
    atomicAdd(&dOmegadz[ic0], szz * wf);

    atomicAdd(&dKdx[ic1],     -sxx * kf);
    atomicAdd(&dKdy[ic1],     -syy * kf);
    atomicAdd(&dKdz[ic1],     -szz * kf);
    atomicAdd(&dOmegadx[ic1], -sxx * wf);
    atomicAdd(&dOmegady[ic1], -syy * wf);
    atomicAdd(&dOmegadz[ic1], -szz * wf);
}

__global__ void calc_scalar_gradient_div_vol_d(
    geom_int nCells,
    geom_float* vol,
    flow_float* dKdx,
    flow_float* dKdy,
    flow_float* dKdz,
    flow_float* dOmegadx,
    flow_float* dOmegady,
    flow_float* dOmegadz)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    const geom_float v = vol[ic];
    dKdx[ic]     /= v;
    dKdy[ic]     /= v;
    dKdz[ic]     /= v;
    dOmegadx[ic] /= v;
    dOmegady[ic] /= v;
    dOmegadz[ic] /= v;
}

std::array<ScalarTransportDesc, 2> buildScalarDescs(variables& var)
{
    return {{
        {var.c_d["k"], var.c_d["roK"], var.c_d["roKN"], var.c_d["roKM"], var.c_d["res_roK"], var.c_d["res_roK_m"], var.c_d["src_jac_k"], static_cast<flow_float>(0.0), static_cast<flow_float>(0.85)},
        {var.c_d["omega"], var.c_d["roOmega"], var.c_d["roOmegaN"], var.c_d["roOmegaM"], var.c_d["res_roOmega"], var.c_d["res_roOmega_m"], var.c_d["src_jac_omega"], static_cast<flow_float>(1.0e-20), static_cast<flow_float>(0.5)}
    }};
}

}

void scalarTransport_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roK"], 0.0, msh.nCells * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roOmega"], 0.0, msh.nCells * sizeof(flow_float)));

    if (!scalarTransportEnabled(cfg)) {
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        return;
    }

    dim3 dimGrid_normal_halo = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
    const auto scalar_descs = buildScalarDescs(var);

    for (const auto& desc : scalar_descs) {
        scalar_advection_first_order_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>>(
            msh.nCells,
            msh.nNormal_halo_Planes,
            msh.normal_halo_planes_d,
            msh.map_plane_cells_d,
            desc.phi,
            var.p_d["massflux"],
            desc.res_rho_phi);

        if (cfg.scalarDiffusion == 1) {
            scalar_diffusion_first_order_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>>(
                msh.nCells,
                msh.nNormal_halo_Planes,
                msh.normal_halo_planes_d,
                msh.map_plane_cells_d,
                var.c_d["ccx"],
                var.c_d["ccy"],
                var.c_d["ccz"],
                var.p_d["fx"],
                var.p_d["sx"],
                var.p_d["sy"],
                var.p_d["sz"],
                var.p_d["ss"],
                desc.phi,
                var.c_d["vis_lam"],
                var.c_d["vis_turb"],
                desc.sigma,
                desc.res_rho_phi);
        }
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void scalarTimeIntegration_d_wrapper(int loop , solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!scalarTransportEnabled(cfg)) {
        return;
    }

    const auto scalar_descs = buildScalarDescs(var);

    for (const auto& desc : scalar_descs) {
        if (cfg.timeIntegration == 4) {
            runge_kutta_exp_scalar_4th_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
                loop,
                cfg.coef_DT_4thRunge[loop],
                cfg.coef_Res_4thRunge[loop],
                var.c_d["dt_local"],
                msh.nCells,
                var.c_d["volume"],
                desc.rho_phi,
                desc.rho_phi_N,
                desc.rho_phi_M,
                desc.res_rho_phi,
                desc.res_rho_phi_m);
        } else if (cfg.timeIntegration == 1 || cfg.timeIntegration == 3) {
            runge_kutta_exp_scalar_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
                cfg.coef_N[loop],
                cfg.coef_M[loop],
                cfg.coef_Res[loop],
                var.c_d["dt_local"],
                msh.nCells,
                var.c_d["volume"],
                desc.rho_phi,
                desc.rho_phi_N,
                desc.rho_phi_M,
                desc.res_rho_phi,
                desc.src_jac,
                desc.floor);
        }
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void calcScalarGradient_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!scalarTransportEnabled(cfg)) {
        return;
    }

    flow_float* grad_vol = (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"];
    flow_float* grad_sx  = (cfg.isAxisymmetric == 1) ? var.p_d["sx_planar"] : var.p_d["sx"];
    flow_float* grad_sy  = (cfg.isAxisymmetric == 1) ? var.p_d["sy_planar"] : var.p_d["sy"];
    flow_float* grad_sz  = (cfg.isAxisymmetric == 1) ? var.p_d["sz_planar"] : var.p_d["sz"];

    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dKdx"],     0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dKdy"],     0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dKdz"],     0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dOmegadx"], 0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dOmegady"], 0, msh.nCells_all * sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dOmegadz"], 0, msh.nCells_all * sizeof(flow_float)));

    calc_scalar_gradient_face_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
        msh.nPlanes,
        msh.map_plane_cells_d,
        var.p_d["fx"],
        grad_sx, grad_sy, grad_sz,
        var.c_d["k"],
        var.c_d["omega"],
        var.c_d["dKdx"], var.c_d["dKdy"], var.c_d["dKdz"],
        var.c_d["dOmegadx"], var.c_d["dOmegady"], var.c_d["dOmegadz"]);

    calc_scalar_gradient_div_vol_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        grad_vol,
        var.c_d["dKdx"], var.c_d["dKdy"], var.c_d["dKdz"],
        var.c_d["dOmegadx"], var.c_d["dOmegady"], var.c_d["dOmegadz"]);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}