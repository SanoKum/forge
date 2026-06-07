#include "dependentVariables_d.cuh"

__global__ void updateVariablesOuter_d
( 
 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* roK  ,
 flow_float* roOmega  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,
 flow_float* roKN ,
 flow_float* roOmegaN ,
 
 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM 

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells_all) {

        roN[ic]   = ro[ic];
        roUxN[ic] = roUx[ic];
        roUyN[ic] = roUy[ic];
        roUzN[ic] = roUz[ic];
        roeN[ic]  = roe[ic];
        roKN[ic] = roK[ic];
        roOmegaN[ic] = roOmega[ic];

        roM[ic]   = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic]  = roe[ic];
        roKM[ic] = roK[ic];
        roOmegaM[ic] = roOmega[ic];

    }
}


void updateVariablesOuter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    updateVariablesOuter_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"]  , var.c_d["roe"]  , var.c_d["roK"] , var.c_d["roOmega"] ,
        var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"] , var.c_d["roeN"] , var.c_d["roKN"] , var.c_d["roOmegaN"] ,
        var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"] , var.c_d["roeM"] , var.c_d["roKM"] , var.c_d["roOmegaM"] 
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

__global__ void updateVariablesInner_d
( 
 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* roK  ,
 flow_float* roOmega  ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM 

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells_all) {

        roM[ic]   = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic]  = roe[ic];
        roKM[ic] = roK[ic];
        roOmegaM[ic] = roOmega[ic];
    }
}


void updateVariablesInner_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    updateVariablesInner_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"]  , var.c_d["roe"]  , var.c_d["roK"] , var.c_d["roOmega"] ,
        var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"] , var.c_d["roeM"] , var.c_d["roKM"] , var.c_d["roOmegaM"] 
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

__global__ void applyScalarImplicitCorrection_d
(
 geom_int nCells_all , geom_int nCells,

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

 flow_float* dq_ro,
 flow_float* dq_roUx,
 flow_float* dq_roUy,
 flow_float* dq_roUz,
 flow_float* dq_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        ro[ic] = roN[ic] + dq_ro[ic];
        roUx[ic] = roUxN[ic] + dq_roUx[ic];
        roUy[ic] = roUyN[ic] + dq_roUy[ic];
        roUz[ic] = roUzN[ic] + dq_roUz[ic];
        roe[ic] = roeN[ic] + dq_roe[ic];
    }
}

void applyScalarImplicitCorrection_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    applyScalarImplicitCorrection_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells_all , msh.nCells,
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
        var.c_d["dq_ro_old"],
        var.c_d["dq_roUx_old"],
        var.c_d["dq_roUy_old"],
        var.c_d["dq_roUz_old"],
        var.c_d["dq_roe_old"]
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

__global__ void applyBlockImplicitCorrection_d
(
 geom_int nCells_all , geom_int nCells,

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

 flow_float* dq_block_0,
 flow_float* dq_block_1,
 flow_float* dq_block_2,
 flow_float* dq_block_3,
 flow_float* dq_block_4
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        ro[ic]   = roN[ic]   + dq_block_0[ic];
        roUx[ic] = roUxN[ic] + dq_block_1[ic];
        roUy[ic] = roUyN[ic] + dq_block_2[ic];
        roUz[ic] = roUzN[ic] + dq_block_3[ic];
        roe[ic]  = roeN[ic]  + dq_block_4[ic];
    }
}

void applyBlockImplicitCorrection_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    // 古典 DPLUR の sweep ループ後、最終補正は dq_block_old に入っている（最後の swap 後）。
    applyBlockImplicitCorrection_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells_all , msh.nCells,
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
        var.c_d["dq_block_old_0"],
        var.c_d["dq_block_old_1"],
        var.c_d["dq_block_old_2"],
        var.c_d["dq_block_old_3"],
        var.c_d["dq_block_old_4"]
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

// SST (k-ω) の segregated point-implicit 更新。
// 平均流 block DPLUR の commit 後に呼び、凍結していたスカラーを陰的に更新する。
// 消散項 + 輸送項(移流+拡散) を point-implicit 化: D_φ = V/Δτ + V·(∂D/∂φ) + transport_diag、
// δ(ρφ) = relax·res_roφ / D_φ。生産・近傍 ΔQ は res_roφ に含む lagged 扱い。ρφ = ρφ_baseline(N) + δ。
__global__ void applySSTPointImplicit_d
(
 geom_int nCells,
 geom_float* vol,
 flow_float* dt_local,
 flow_float implicit_relax,

 flow_float* roK,
 flow_float* roOmega,
 flow_float* roKN,
 flow_float* roOmegaN,

 flow_float* res_roK,
 flow_float* res_roOmega,
 flow_float* src_jac_k,
 flow_float* src_jac_omega,
 flow_float* transport_diag_k,
 flow_float* transport_diag_omega,
 flow_float unsteady_diag
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float v    = vol[ic];
        const flow_float dt_l = max(dt_local[ic], static_cast<flow_float>(1.0e-30));
        // D_φ = V/Δτ + V·src_jac（消散）+ transport_diag（移流+拡散、既に [m³/s]）+ V·unsteady_diag。
        // transport_diag は壁近傍の陽的 k/ω 輸送 stiff 性を陰化し安定 cfl_pseudo を一桁以上引き上げる。
        const flow_float Dk = v / dt_l + v * src_jac_k[ic]     + transport_diag_k[ic]     + v * unsteady_diag;
        const flow_float Dw = v / dt_l + v * src_jac_omega[ic] + transport_diag_omega[ic] + v * unsteady_diag;

        const flow_float dk = implicit_relax * res_roK[ic]     / max(Dk, static_cast<flow_float>(1.0e-30));
        const flow_float dw = implicit_relax * res_roOmega[ic] / max(Dw, static_cast<flow_float>(1.0e-30));

        // realizability: ρk ≥ 0, ρω > 0。dual-time でも正しいよう現在反復値への in-place 加算
        // （定常では roKN==roK のため roKN+dk と等価）。dual-time の roKN/roKNN は BDF 残差項側で使用。
        roK[ic]     = max(roK[ic]     + dk, static_cast<flow_float>(0.0));
        roOmega[ic] = max(roOmega[ic] + dw, static_cast<flow_float>(1.0e-20));
    }
}

void applySSTPointImplicit_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    applySSTPointImplicit_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["volume"],
        var.c_d["dt_local"],
        cfg.implicitRelax,
        var.c_d["roK"],
        var.c_d["roOmega"],
        var.c_d["roKN"],
        var.c_d["roOmegaN"],
        var.c_d["res_roK"],
        var.c_d["res_roOmega"],
        var.c_d["src_jac_k"],
        var.c_d["src_jac_omega"],
        var.c_d["transport_diag_k"],
        var.c_d["transport_diag_omega"],
        cfg.unsteadyDiagCoef
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

// dual-time: 物理時間 BDF 項を残差に加える（assembleResidual の末尾で呼ぶ）。
// res は -R_spatial（陽解法 Q=Q_N+res·dt/V と整合）。dual-time の残差 R* = R_spatial + (V/Δt)(a Q - b Q^n + c Q^{n-1})
// を 0 にするので res* = res - (V/Δt)(a Q - b Q^n + c Q^{n-1}) とする（roN=Q^n, roNN=Q^{n-1}）。
// (a,b,c): BDF1=(1,1,0), BDF2=(1.5,2,0.5)。include_scalar==1 で k/ω も処理。
__global__ void addUnsteadyTimeTerm_d
(
 geom_int nCells,
 geom_float* vol,
 flow_float dt_phys,
 flow_float a, flow_float b, flow_float c,

 flow_float* ro,   flow_float* roUx,   flow_float* roUy,   flow_float* roUz,   flow_float* roe,
 flow_float* roN,  flow_float* roUxN,  flow_float* roUyN,  flow_float* roUzN,  flow_float* roeN,
 flow_float* roNN, flow_float* roUxNN, flow_float* roUyNN, flow_float* roUzNN, flow_float* roeNN,
 flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe,

 int include_scalar,
 flow_float* roK,   flow_float* roOmega,
 flow_float* roKN,  flow_float* roOmegaN,
 flow_float* roKNN, flow_float* roOmegaNN,
 flow_float* res_roK, flow_float* res_roOmega
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float coef = vol[ic] / max(dt_phys, static_cast<flow_float>(1.0e-30));
        res_ro[ic]   -= coef * (a*ro[ic]   - b*roN[ic]   + c*roNN[ic]);
        res_roUx[ic] -= coef * (a*roUx[ic] - b*roUxN[ic] + c*roUxNN[ic]);
        res_roUy[ic] -= coef * (a*roUy[ic] - b*roUyN[ic] + c*roUyNN[ic]);
        res_roUz[ic] -= coef * (a*roUz[ic] - b*roUzN[ic] + c*roUzNN[ic]);
        res_roe[ic]  -= coef * (a*roe[ic]  - b*roeN[ic]  + c*roeNN[ic]);
        if (include_scalar) {
            res_roK[ic]     -= coef * (a*roK[ic]     - b*roKN[ic]     + c*roKNN[ic]);
            res_roOmega[ic] -= coef * (a*roOmega[ic] - b*roOmegaN[ic] + c*roOmegaNN[ic]);
        }
    }
}

void addUnsteadyTimeTerm_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var,
                                   flow_float a, flow_float b, flow_float c, int include_scalar)
{
    addUnsteadyTimeTerm_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["volume"],
        cfg.dt,
        a, b, c,
        var.c_d["ro"],   var.c_d["roUx"],   var.c_d["roUy"],   var.c_d["roUz"],   var.c_d["roe"],
        var.c_d["roN"],  var.c_d["roUxN"],  var.c_d["roUyN"],  var.c_d["roUzN"],  var.c_d["roeN"],
        var.c_d["roNN"], var.c_d["roUxNN"], var.c_d["roUyNN"], var.c_d["roUzNN"], var.c_d["roeNN"],
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"],
        include_scalar,
        var.c_d["roK"],   var.c_d["roOmega"],
        var.c_d["roKN"],  var.c_d["roOmegaN"],
        var.c_d["roKNN"], var.c_d["roOmegaNN"],
        var.c_d["res_roK"], var.c_d["res_roOmega"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

// dual-time: 擬似時間サブ反復の commit。roN=Q^n は BDF 基準として固定するため、
// 定常の Q=roN+dq ではなく現在反復値への in-place 加算 Q += dq を行う。
__global__ void applyBlockImplicitCorrectionInPlace_d
(
 geom_int nCells,
 flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz, flow_float* roe,
 flow_float* dq0, flow_float* dq1, flow_float* dq2, flow_float* dq3, flow_float* dq4
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        ro[ic]   += dq0[ic];
        roUx[ic] += dq1[ic];
        roUy[ic] += dq2[ic];
        roUz[ic] += dq3[ic];
        roe[ic]  += dq4[ic];
    }
}

void applyBlockImplicitCorrectionInPlace_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    applyBlockImplicitCorrectionInPlace_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
        var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"],
        var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

// dual-time の物理時間レベルシフト: roNN ← roN, roN ← ro（mean flow 5 + scalar 2）。
__global__ void shiftDualTimeLevels_d
(
 geom_int nCells_all,
 flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz, flow_float* roe,
 flow_float* roN, flow_float* roUxN, flow_float* roUyN, flow_float* roUzN, flow_float* roeN,
 flow_float* roNN, flow_float* roUxNN, flow_float* roUyNN, flow_float* roUzNN, flow_float* roeNN,
 flow_float* roK, flow_float* roOmega,
 flow_float* roKN, flow_float* roOmegaN,
 flow_float* roKNN, flow_float* roOmegaNN
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells_all) {
        roNN[ic]=roN[ic]; roUxNN[ic]=roUxN[ic]; roUyNN[ic]=roUyN[ic]; roUzNN[ic]=roUzN[ic]; roeNN[ic]=roeN[ic];
        roN[ic]=ro[ic];   roUxN[ic]=roUx[ic];   roUyN[ic]=roUy[ic];   roUzN[ic]=roUz[ic];   roeN[ic]=roe[ic];
        roKNN[ic]=roKN[ic]; roOmegaNN[ic]=roOmegaN[ic];
        roKN[ic]=roK[ic];   roOmegaN[ic]=roOmega[ic];
    }
}

void shiftDualTimeLevels_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    shiftDualTimeLevels_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells_all,
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
        var.c_d["roN"], var.c_d["roUxN"], var.c_d["roUyN"], var.c_d["roUzN"], var.c_d["roeN"],
        var.c_d["roNN"], var.c_d["roUxNN"], var.c_d["roUyNN"], var.c_d["roUzNN"], var.c_d["roeNN"],
        var.c_d["roK"], var.c_d["roOmega"],
        var.c_d["roKN"], var.c_d["roOmegaN"],
        var.c_d["roKNN"], var.c_d["roOmegaNN"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}