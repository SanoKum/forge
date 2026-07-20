#include "nodeWallDirichlet_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/thermo_d.cuh"
#include "cuda_forge/speciesTransport_d.cuh"   // species_roY_device_ptr()
#include "cuda_forge/wmlesWallModel_d.cuh"     // wmlesActiveForBcond / wmlesNodeIsothermalActive (モデル層ゲート)

// node 壁強境界の「状態層」実装。設計は nodeWallDirichlet_d.cuh 冒頭コメントと
// methods/boundary.md「node 等温壁の壁ノード温度ピン」を参照。
// (enforceWallNoSlip / 残差射影は axisymmetricSource_d.cu から、温度ピンは
//  wmlesWallModel_d.cu から移設・共通化した。演算列は移設前と同一。2026-07-20)

namespace {

// 壁ノード状態を u=0 へ確定 (KE を roe から除去)
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

// 壁ノードの Dirichlet 行残差射影: 状態やエネルギー積算をいじらず flux+source 積算後に
// 残差を射影するだけなので block-DPLUR と整合する。どの行をゼロにするかは wrapper が決める
// (運動量 3 行は常時 / ω は SST / roe は 等温ピンまたは WMLES 等温 [Qw_Wall マーカ])。
__global__ void zeroWallDirichletResiduals_d
(
    geom_int nCells, geom_int* wall_flag,
    flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz,
    // SST: 壁ノードで ω は Dirichlet ピン (rans_wall_scalar_boundary_d) なので残差を 0 に射影する。
    // Dirichlet ノードの残差は BC 強制であり物理的不均衡でない → rms_roOmega の汚染 (収束判定の誤検出) を防ぐ。
    // nullptr で無効 (非 SST)。k はノイマンなので res_roK は触らない。
    flow_float* res_roOmega,
    // WMLES 等温壁 (node): 壁ノード温度は温度ピンで Dirichlet されるため res_roe も 0 に射影する。
    // 対象ノードの識別は Qw_Wall マーカ (>-0.5 = 等温 WMLES 壁ノード)。nullptr で無効。
    // 素の等温壁 (非 WMLES) の res_roe は zeroNodeIsothermalEnergyResidual (bcond 単位) が担う。
    flow_float* Qw_Wall, flow_float* res_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells && wall_flag[ic] == 1) {
        res_roUx[ic] = (flow_float)0.0;
        res_roUy[ic] = (flow_float)0.0;
        res_roUz[ic] = (flow_float)0.0;
        if (res_roOmega != nullptr) res_roOmega[ic] = (flow_float)0.0;
        if (Qw_Wall != nullptr && Qw_Wall[ic] > (flow_float)-0.5) res_roe[ic] = (flow_float)0.0;
    }
}

// 等温壁ノードの温度状態ピン: T, roe (=ρ(e(Tw)+ek)), P (=ρR Tw), sonic を整合させる。
// 運動量は enforceWallNoSlip / nodeWallDirichlet が担う (nodeWallDirichlet=1 なら ek=0)。
__global__ void pin_wall_node_temperature_d(
    geom_int nb,
    geom_int* bplane_cell,
    int thermalMethod, flow_float ga, flow_float cp_const,
    const SpeciesThermo* sp, flow_float* const* roY, int nSpecies,
    flow_float* Tsb,
    flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz,
    flow_float* T, flow_float* P, flow_float* roe, flow_float* sonic)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib >= nb) return;
    const geom_int ic = bplane_cell[ib];
    const flow_float Tw = Tsb[ib];
    const flow_float ek = static_cast<flow_float>(0.5)
                        * (Ux[ic]*Ux[ic] + Uy[ic]*Uy[ic] + Uz[ic]*Uz[ic]);
    const GasStateAtT g = thermo_state_at_T(thermalMethod, ga, cp_const, sp, roY, nSpecies, ic, Tw);
    T[ic]     = Tw;
    roe[ic]   = ro[ic] * (g.e + ek);
    P[ic]     = ro[ic] * g.R * Tw;
    sonic[ic] = sqrt(g.gamma * g.R * Tw);
}

// res_roe を bcond の壁ノードで 0 化 (Dirichlet ノードの残差は BC 強制であり物理不均衡でない)
__global__ void zero_res_roe_bplane_d(geom_int nb, geom_int* bplane_cell, flow_float* res_roe)
{
    const geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;
    if (ib >= nb) return;
    res_roe[bplane_cell[ib]] = static_cast<flow_float>(0.0);
}

} // namespace

void enforceWallNoSlip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || cfg.nodeWallDirichlet == 0 || msh.wall_flag_d == nullptr) return;
    enforceWallNoSlip_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.wall_flag_d,
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
        var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void zeroWallDirichletResiduals_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || cfg.nodeWallDirichlet == 0 || msh.wall_flag_d == nullptr) return;
    const bool sst = (cfg.LESorRANS == 2 && cfg.RANSmodel == 1);
    const bool wmlesIso = wmlesNodeIsothermalActive(cfg, msh);
    zeroWallDirichletResiduals_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.wall_flag_d,
        var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"],
        sst ? var.c_d["res_roOmega"] : nullptr,
        wmlesIso ? var.c_d["Qw_Wall"] : nullptr,
        wmlesIso ? var.c_d["res_roe"] : nullptr
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

    // 素の node 等温壁 (非 WMLES): 壁ノード T ピン (applyNodeIsothermalWallPin) と対で res_roe を 0 化
    zeroNodeIsothermalEnergyResidual(cfg , cuda_cfg , msh , var);
}

bool nodeIsothermalPinActive(const solverConfig& cfg, const mesh& msh)
{
    if (cfg.discretization != "node" || cfg.nodeWallDirichlet == 0 || msh.wall_flag_d == nullptr) return false;
    for (const auto& bc : msh.bconds)
        if (bc.bcondKind == "wall_isothermal" && !wmlesActiveForBcond(cfg, bc)) return true;
    return false;
}

void pinWallNodeTemperature_bcond(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , variables& var)
{
    if (bc.iPlanes.empty()) return;
    pin_wall_node_temperature_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
        static_cast<geom_int>(bc.iPlanes.size()),
        bc.map_bplane_cell_d,
        cfg.thermalMethod, cfg.gamma, cfg.cp,
        thermo_species_device_ptr(), species_roY_device_ptr(), cfg.nSpecies,
        bc.bvar_d["Ts"],
        var.c_d["ro"], var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"],
        var.c_d["T"], var.c_d["P"], var.c_d["roe"], var.c_d["sonic"]);
}

// 素の node 等温壁 (非 WMLES): 壁ノード温度を BC 値にピンする (nodeWallDirichlet の熱版)。
// 旧来はエネルギーだけ弱形式で、壁ノード T が壁 CV 平均に緩み第 1 スペーシングの勾配が歪んだ
// (case/24 純伝導検証 2026-07-20)。WMLES 等温壁は applyWmlesWallModel が同じ
// pinWallNodeTemperature_bcond を呼ぶ。
void applyNodeIsothermalWallPin(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!nodeIsothermalPinActive(cfg, msh)) return;
    for (auto& bc : msh.bconds) {
        if (bc.bcondKind != "wall_isothermal") continue;
        if (wmlesActiveForBcond(cfg, bc)) continue;   // WMLES 側で同一ピン適用済み
        pinWallNodeTemperature_bcond(cfg, cuda_cfg, bc, var);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void zeroNodeIsothermalEnergyResidual(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!nodeIsothermalPinActive(cfg, msh)) return;
    for (auto& bc : msh.bconds) {
        if (bc.bcondKind != "wall_isothermal") continue;
        if (wmlesActiveForBcond(cfg, bc)) continue;
        if (bc.iPlanes.empty()) continue;
        zero_res_roe_bplane_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
            static_cast<geom_int>(bc.iPlanes.size()),
            bc.map_bplane_cell_d,
            var.c_d["res_roe"]);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
