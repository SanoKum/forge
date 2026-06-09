#include "speciesTransport_d.cuh"

#include "scalarTransport_d.cuh"
#include "thermo_d.cuh"

#include <string>
#include <vector>

namespace {

// device の roY / res_roY / transport_diag ポインタ配列 (flow_float*[nSpecies])。
// speciesInit_d で 1 度だけ構築 (kinetic 拡散カーネルが化学種ループするため)。
flow_float** g_roY_dev = nullptr;
flow_float** g_resroY_dev = nullptr;
flow_float** g_transdiag_dev = nullptr;
int          g_nSpecies = 0;

constexpr flow_float kSmall = static_cast<flow_float>(1.0e-30);

inline bool speciesEnabled(const variables& var)
{
    return var.nSpeciesRegistered >= 2;
}

// 化学種 s 用のスカラ輸送記述子を構築する (RANS buildScalarDescs と同形)。
// floor=0 (Y_s>=0)、sigma=0 (M2 は移流のみ。拡散は M4)。
ScalarTransportDesc buildSpeciesDesc(variables& var, int s)
{
    const std::string i = std::to_string(s);
    return ScalarTransportDesc{
        var.c_d["Y"+i], var.c_d["roY"+i], var.c_d["roY"+i+"N"], var.c_d["roY"+i+"M"],
        var.c_d["res_roY"+i], var.c_d["res_roY"+i+"_m"],
        var.c_d["src_jac_Y"+i], var.c_d["transport_diag_Y"+i],
        static_cast<flow_float>(0.0), static_cast<flow_float>(0.0),
        0  // 汎用拡散は使わない (化学種 Fick 拡散は species_diffusion_d で別途)
    };
}

__global__ void species_primitive_d(
    geom_int nCells_all,
    flow_float* ro,
    flow_float* roY,
    flow_float* Y)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells_all) {
        Y[ic] = roY[ic] / max(ro[ic], kSmall);
    }
}

// Neumann (zero-gradient) ghost 充填: roY[ig]=roY[ic], Y[ig]=Y[ic]。
__global__ void species_neumann_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* roY,
    flow_float* Y)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];
        roY[ig] = roY[ic];
        Y[ig]   = Y[ic];
    }
}

// 実現可能性 + 再正規化: 各 ρY_s>=0 にクランプ後、Σ_s ρY_s = ρ となるよう再スケール (ΣY_s=1)。
__global__ void species_renormalize_d(
    geom_int nCells,
    int nSpecies,
    flow_float** roY,
    flow_float* ro)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        double sum = 0.0;
        for (int s = 0; s < nSpecies; s++) {
            flow_float v = roY[s][ic];
            if (v < 0.0) v = 0.0;
            roY[s][ic] = v;
            sum += (double)v;
        }
        const double factor = (double)ro[ic] / (sum > (double)kSmall ? sum : (double)kSmall);
        for (int s = 0; s < nSpecies; s++) {
            roY[s][ic] = (flow_float)((double)roY[s][ic] * factor);
        }
    }
}

// M4: 化学種 Fick 拡散 + 質量保存補正 (ΣJ=0) + エンタルピー拡散のエネルギー結合。
//   J_s = ρ D_s ∇Y_s (over-relaxed 法線), 補正 J_s* = J_s - Y_s Σ_k J_k (Σ J_s*=0),
//   res_roY_s += J_s*, res_roe += Σ_s h_s(T) J_s*。
//   D_s: diffMethod==1 で混合平均 (Chapman-Enskog), ==0 で定数 Schmidt D=μ/(ρ Sc)。
__global__ void species_diffusion_d(
    geom_int nCells,
    geom_int nNormalHaloPlanes,
    geom_int* normal_halo_planes,
    geom_int* plane_cells,
    geom_float* ccx, geom_float* ccy, geom_float* ccz,
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    const SpeciesThermo* sp, int nSpecies,
    flow_float** roY, flow_float** res_roY, flow_float** transport_diag,
    flow_float* ro, flow_float* T, flow_float* P, flow_float* vis_lam, flow_float* vis_turb,
    flow_float* res_roe,
    int diffMethod, flow_float Sc, flow_float Sc_t)
{
    geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;

    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2 * ip + 0];
    const geom_int ic1 = plane_cells[2 * ip + 1];

    const geom_float f   = fx[ip];
    const geom_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];

    const flow_float dccx = ccx[ic1] - ccx[ic0];
    const flow_float dccy = ccy[ic1] - ccy[ic0];
    const flow_float dccz = ccz[ic1] - ccz[ic0];
    const flow_float dcc  = sqrt(dccx*dccx + dccy*dccy + dccz*dccz);
    const flow_float denom = dccx*sxx + dccy*syy + dccz*szz;
    const flow_float Dsafe = (fabs(denom) < 1.0e-30) ? ((denom>=0)?1.0e-30:-1.0e-30) : denom;
    const flow_float delta = dcc * sss * sss / Dsafe;       // over-relaxed 法線

    const double ro0 = (double)max(ro[ic0], (flow_float)1.0e-30);
    const double ro1 = (double)max(ro[ic1], (flow_float)1.0e-30);
    const double ro_face = (double)f*ro0 + (1.0-(double)f)*ro1;
    const double T_face  = (double)f*(double)T[ic0] + (1.0-(double)f)*(double)T[ic1];
    const double P_face  = (double)f*(double)P[ic0] + (1.0-(double)f)*(double)P[ic1];

    // 面の質量分率 Yf (正規化) -> モル分率 X
    double Yf[THERMO_MAX_SPECIES], X[THERMO_MAX_SPECIES];
    double ysum = 0.0;
    for (int s=0;s<nSpecies;s++){
        double y = (double)f*((double)roY[s][ic0]/ro0) + (1.0-(double)f)*((double)roY[s][ic1]/ro1);
        if (y<0.0) y=0.0; Yf[s]=y; ysum+=y;
    }
    const double yinv = 1.0/(ysum>1.0e-30?ysum:1.0e-30);
    for (int s=0;s<nSpecies;s++) Yf[s]*=yinv;
    thermo_X_from_Y(sp, nSpecies, Yf, X);

    const double mu_face = (double)f*(double)vis_lam[ic0] + (1.0-(double)f)*(double)vis_lam[ic1];
    const double mut_face = (double)f*(double)vis_turb[ic0] + (1.0-(double)f)*(double)vis_turb[ic1];
    // 乱流化学種拡散 D_t = μ_t/(ρ Sc_t) は全種共通で加える。
    const double Dt = (mut_face > 0.0) ? mut_face/(ro_face*(double)Sc_t) : 0.0;

    // 各化学種の非補正 Fick flux J_s と Σ
    double Js[THERMO_MAX_SPECIES];
    double sumJ = 0.0;
    for (int s=0;s<nSpecies;s++){
        const double Ys0 = (double)roY[s][ic0]/ro0;
        const double Ys1 = (double)roY[s][ic1]/ro1;
        double D;
        if (diffMethod == 1) D = thermo_Dmix_species(sp, nSpecies, X, s, T_face, P_face);
        else                 D = mu_face/(ro_face*(double)Sc);
        D += Dt;  // 層流 (Fick/Sc) + 乱流 (μ_t/Sc_t)
        Js[s] = ro_face * D * ((Ys1 - Ys0)/(double)dcc) * (double)delta;
        sumJ += Js[s];
        // point-implicit 拡散対角 (各セル ρ で正規化)
        const double diag = ro_face * D * fabs((double)delta) / (double)max(dcc,(flow_float)1.0e-30);
        if (ic0 < nCells) atomicAdd(&transport_diag[s][ic0], (flow_float)(diag/ro0));
        if (ic1 < nCells) atomicAdd(&transport_diag[s][ic1], (flow_float)(diag/ro1));
    }

    // 補正 J_s* = J_s - Y_s Σ_k J_k (Σ J_s* = 0) と エネルギー結合 Σ h_s J_s*
    double q = 0.0;
    for (int s=0;s<nSpecies;s++){
        const double Jc = Js[s] - Yf[s]*sumJ;
        if (ic0 < nCells) atomicAdd(&res_roY[s][ic0],  (flow_float)Jc);
        if (ic1 < nCells) atomicAdd(&res_roY[s][ic1], -(flow_float)Jc);
        const double hs = thermo_h_mass(sp[s], T_face);   // NASA 絶対エンタルピー [J/kg]
        q += hs * Jc;
    }
    if (ic0 < nCells) atomicAdd(&res_roe[ic0],  (flow_float)q);
    if (ic1 < nCells) atomicAdd(&res_roe[ic1], -(flow_float)q);
}

}  // namespace

void speciesInit_d(solverConfig& cfg, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) {
        g_roY_dev = nullptr;
        g_nSpecies = 0;
        return;
    }
    g_nSpecies = var.nSpeciesRegistered;

    std::vector<flow_float*> hroY(g_nSpecies), hres(g_nSpecies), htd(g_nSpecies);
    for (int s = 0; s < g_nSpecies; s++) {
        const std::string i = std::to_string(s);
        hroY[s] = var.c_d["roY"+i];
        hres[s] = var.c_d["res_roY"+i];
        htd[s]  = var.c_d["transport_diag_Y"+i];
    }
    const size_t pbytes = g_nSpecies*sizeof(flow_float*);
    gpuErrchk( cudaMalloc((void**)&g_roY_dev, pbytes) );
    gpuErrchk( cudaMalloc((void**)&g_resroY_dev, pbytes) );
    gpuErrchk( cudaMalloc((void**)&g_transdiag_dev, pbytes) );
    gpuErrchk( cudaMemcpy(g_roY_dev, hroY.data(), pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(g_resroY_dev, hres.data(), pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(g_transdiag_dev, htd.data(), pbytes, cudaMemcpyHostToDevice) );

    std::cout << "speciesInit_d: built device roY/res/diag[] for nSpecies=" << g_nSpecies << "\n";
}

flow_float** species_roY_device_ptr()
{
    return g_roY_dev;
}

void speciesPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        species_primitive_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            msh.nCells_all,
            var.c_d["ro"],
            var.c_d["roY"+i],
            var.c_d["Y"+i]);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var)
{
    (void)msh;
    if (!speciesEnabled(var)) return;
    if (bc.iPlanes.empty()) return;

    const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        // M2: 全境界種別を Neumann で扱う (slip 閉領域・内部 contact の検証に十分)。
        species_neumann_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
            nb,
            bc.map_bplane_cell_d,
            bc.map_bplane_cell_ghst_d,
            var.c_d["roY"+i],
            var.c_d["Y"+i]);
    }
}

void applySpeciesBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (auto& bc : msh.bconds) {
        speciesBoundary_d_wrapper(cfg, cuda_cfg, bc, msh, var);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roY"+i], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["transport_diag_Y"+i], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["src_jac_Y"+i], 0, msh.nCells * sizeof(flow_float)));
    }

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const ScalarTransportDesc desc = buildSpeciesDesc(var, s);
        scalarTransportResidual_d(cfg, cuda_cfg, msh, var, desc);
    }

    // M4: 粘性ケースのみ Fick 拡散 + ΣJ=0 補正 + エンタルピー拡散 (res_roe へ加算)。
    if (cfg.viscMethod != 0 && g_roY_dev != nullptr) {
        dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
        species_diffusion_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
            msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
            var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
            thermo_species_device_ptr(), g_nSpecies,
            g_roY_dev, g_resroY_dev, g_transdiag_dev,
            var.c_d["ro"], var.c_d["T"], var.c_d["P"], var.c_d["vis_lam"], var.c_d["vis_turb"],
            var.c_d["res_roe"],
            cfg.speciesDiffusionMethod, cfg.Sc, cfg.Sc_t);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const ScalarTransportDesc desc = buildSpeciesDesc(var, s);
        scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesRenormalize_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var) || g_roY_dev == nullptr) return;

    species_renormalize_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        g_nSpecies,
        g_roY_dev,
        var.c_d["ro"]);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

void speciesUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!speciesEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"N"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"M"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
    }
}

void speciesUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!speciesEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"M"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
    }
}
