#include "dependentVariables_d.cuh"
#include "thermo_d.cuh"

// 温度反転のクランプ範囲 (NASA-9 の有効域より広めに取り, 範囲外は外挿)
#define DEPVAR_TMIN 50.0
#define DEPVAR_TMAX 6000.0

__global__ void dependentVariables_d
(
 // gas properties
 int thermalMethod ,
 flow_float gamma , flow_float cp ,

 // thermally-perfect (thermalMethod==2) 用化学種データ
 const SpeciesThermo* sp , int nSpecies , flow_float** roY ,

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

 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* k   ,
 flow_float* omega,
 flow_float* T   ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 // thermally-perfect 用に毎ステップ更新する per-cell 物性
 flow_float* gam_array ,
 flow_float* cp_array
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    flow_float ek;
    flow_float intE;

    flow_float ro_temp;
    flow_float T_temp;
    flow_float P_temp;

    if (ic < nCells_all) {
        // 密度・圧力フロア: 膨張領域で非物理的な ro→0, P→0 が生じても速度爆発を防ぐ。
        // 1e-6f は実質ゼロと同じ（大気圧比 1e-11）なので物理的に意味のある値にする。
        // ro_min = 1e-4 kg/m³ (大気密度の 1/12000)
        // P_min  = 1.0 Pa    (大気圧の 1/100000)
        ro_temp = max(ro[ic], (flow_float)1.0e-4f);

        Ux[ic] = roUx[ic]/ro_temp;
        Uy[ic] = roUy[ic]/ro_temp;
        Uz[ic] = roUz[ic]/ro_temp;

        ek = 0.5*(Ux[ic]*Ux[ic] +Uy[ic]*Uy[ic] +Uz[ic]*Uz[ic]);
        intE =(roe[ic]/ro_temp -ek);

        if (thermalMethod == 2) {
            // ---- 多成分 thermally-perfect gas (NASA-9) ----
            // 内部計算は全て double。組成 Y を構築 (nSpecies==1 は Y={1})。
            double Y[THERMO_MAX_SPECIES];
            if (nSpecies <= 1 || roY == nullptr) {
                Y[0] = 1.0;
            } else {
                double ysum = 0.0;
                for (int s=0;s<nSpecies;s++){
                    double y = (double)roY[s][ic]/(double)ro_temp;
                    if (y < 0.0) y = 0.0;
                    Y[s] = y; ysum += y;
                }
                double inv = 1.0/(ysum > 1.0e-30 ? ysum : 1.0e-30);
                for (int s=0;s<nSpecies;s++) Y[s] *= inv;
            }

            const double e_in = (double)intE;
            const double Tg   = ((double)T[ic] > DEPVAR_TMIN) ? (double)T[ic] : 300.0; // warm start
            const double Tnew = thermo_T_from_e(sp, nSpecies, Y, e_in, Tg, DEPVAR_TMIN, DEPVAR_TMAX);

            const double Rmix  = thermo_R_mix (sp, nSpecies, Y);
            const double cpmix = thermo_cp_mix(sp, nSpecies, Y, Tnew);
            const double hmix  = thermo_h_mix (sp, nSpecies, Y, Tnew);
            const double cvmix = cpmix - Rmix;
            const double gmix  = cpmix / (cvmix > 1.0e-6 ? cvmix : 1.0e-6);

            double Pnew = (double)ro_temp * Rmix * Tnew;
            if (Pnew < 1.0) Pnew = 1.0;

            T[ic]         = (flow_float)Tnew;
            P[ic]         = (flow_float)Pnew;
            ro[ic]        = ro_temp;
            // roe を (floor 済 ro, 反転 T) と整合させて再構成
            roe[ic]       = (flow_float)((double)ro_temp * ((hmix - Rmix*Tnew) + (double)ek));
            Ht[ic]        = (flow_float)(hmix + (double)ek);
            sonic[ic]     = (flow_float)sqrt(gmix * Rmix * Tnew);
            gam_array[ic] = (flow_float)gmix;
            cp_array[ic]  = (flow_float)cpmix;
        } else {
            // ---- calorically perfect gas (従来経路, ビット不変) ----
            T_temp = max(intE/(cp/gamma), (flow_float)1.0e-4f);
            P_temp = max((gamma-1.0)*(roe[ic]-ro_temp*ek),(flow_float)1.0f);

            T[ic] = T_temp;
            P[ic] = P_temp;

            ro[ic] = ro_temp;
            roe[ic] = P_temp/(gamma-1.0) + ro_temp*ek;

            Ht[ic] = roe[ic]/ro_temp + P_temp/ro_temp;

            sonic[ic] = sqrt(gamma*P_temp/ro_temp);
        }

        k[ic] = max(roK[ic]/ro_temp, static_cast<flow_float>(0.0));
        omega[ic] = max(roOmega[ic]/ro_temp, static_cast<flow_float>(0.0));
    }
}


void dependentVariables_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    dependentVariables_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
        // gas properties
        cfg.thermalMethod ,
        cfg.gamma , cfg.cp ,

        // thermally-perfect 用化学種データ (M1: 単成分なので roY=nullptr)
        thermo_species_device_ptr() , cfg.nSpecies , nullptr ,

        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"], var.c_d["roUy"] , var.c_d["roUz"], var.c_d["roe"] ,
        var.c_d["roK"] , var.c_d["roOmega"],
        var.c_d["P"]   , var.c_d["Ht"]  , var.c_d["sonic"], var.c_d["k"], var.c_d["omega"], var.c_d["T"],
        var.c_d["Ux"]  , var.c_d["Uy"]  , var.c_d["Uz"] ,

        var.c_d["gamma"] , var.c_d["cp"]
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}
