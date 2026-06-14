#include "dependentVariables_d.cuh"
#include "thermo_d.cuh"
#include "speciesTransport_d.cuh"  // species_roY_device_ptr() (多成分 thermalMethod==2)
#include "condensationTransport_d.cuh"  // cond_rog_device_ptr() (二相 EOS)
#include "condensationEOS_d.cuh"        // cond_T_from_e_onetemp (一温度 二相 温度反転)

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

 // 非平衡凝縮 (一温度 二相 EOS)。condensation==0 のとき従来経路 (ビット不変)。
 int condensation , int nCondSpecies , flow_float** rog ,
 int condGasSpecies , int condModel ,   // carrier+condensible: 凝縮気相種 index / モデル(0:N2,1:H2O)

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
 flow_float* cp_array ,
 flow_float* Rmix_array   // M6: 混合比気体定数 R[ic] (SLAU 面エンタルピーが毎面再計算するのを回避)
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

            // 非平衡凝縮 (一温度 二相 EOS): 総液相質量分率 g を集計。condensation==0 で g=0 (従来経路)。
            double g_liq = 0.0;
            const bool carrier = (condensation == 1 && condGasSpecies >= 0);
            if (condensation == 1 && rog != nullptr) {
                for (int s = 0; s < nCondSpecies; ++s) {
                    double gs = (double)rog[s][ic] / (double)ro_temp;
                    if (gs > 0.0) g_liq += gs;
                }
                // realizability: carrier は g≤Y_凝縮種 (蒸気以上は凝縮しない)、pure は g≤0.99。
                if (carrier && roY != nullptr) {
                    double Yw = (double)roY[condGasSpecies][ic]/(double)ro_temp;
                    if (g_liq > Yw) g_liq = (Yw > 0.0 ? Yw : 0.0);
                } else if (g_liq > 0.99) g_liq = 0.99;
                if (g_liq < 0.0) g_liq = 0.0;
            }

            const CondSpeciesProps cprops = (condModel == 1) ? condProps_H2O() : condProps_N2();
            const double Rw = cprops.R;   // 凝縮種の比気体定数 (carrier: 蒸気分圧/EOS に使用)

            // 温度反転。g≈0 は従来 thermo_T_from_e で厳密縮約。
            double Tnew;
            if (g_liq > 1.0e-12 && carrier) {
                Tnew = cond_T_from_e_carrier(sp, nSpecies, Y, e_in, g_liq, Rw, cprops, Tg, DEPVAR_TMIN, DEPVAR_TMAX);
            } else if (g_liq > 1.0e-12) {
                Tnew = cond_T_from_e_onetemp(sp, nSpecies, Y, e_in, g_liq, Tg, DEPVAR_TMIN, DEPVAR_TMAX);
            } else {
                Tnew = thermo_T_from_e(sp, nSpecies, Y, e_in, Tg, DEPVAR_TMIN, DEPVAR_TMAX);
            }

            const double Rmix  = thermo_R_mix (sp, nSpecies, Y);
            double cpmix, hmix;
            thermo_cph_mix(sp, nSpecies, Y, Tnew, &cpmix, &hmix);  // cp,h を 1 スイープ (全蒸気混合)
            const double cvmix = cpmix - Rmix;
            const double gmix  = cpmix / (cvmix > 1.0e-6 ? cvmix : 1.0e-6);

            const double e_v   = hmix - Rmix*Tnew;
            const double Lcond = (g_liq > 1.0e-12) ? cond_latent(cprops, Tnew) : 0.0;
            double e_mix, Pnew, oneMg;
            if (carrier) {
                // carrier+condensible: e_mix=e_全蒸気+g(R_w T-L)、p=ρT(R_mix-g R_w)(凝縮で蒸気モル減)。
                e_mix = e_v + g_liq*(Rw*Tnew - Lcond);
                Pnew  = (double)ro_temp * Tnew * (Rmix - g_liq*Rw);
                oneMg = 1.0;  // 気相質量は別途、p で表現済
            } else {
                // pure-condensible (気相=凝縮種): e_l=e_v+R_vT-L、p=(1-g)ρR T。
                e_mix = e_v + g_liq*Rmix*Tnew - g_liq*Lcond;
                oneMg = 1.0 - g_liq;
                Pnew  = (double)ro_temp * oneMg * Rmix * Tnew;
            }
            if (Pnew < 1.0) Pnew = 1.0;

            T[ic]         = (flow_float)Tnew;
            P[ic]         = (flow_float)Pnew;
            ro[ic]        = ro_temp;
            // roe を (floor 済 ro, 反転 T, 混合内部エネルギー) と整合させて再構成
            roe[ic]       = (flow_float)((double)ro_temp * (e_mix + (double)ek));
            // 総エンタルピー Ht = e_mix + p/ρ + ek (g=0 → hmix+ek と一致)
            Ht[ic]        = (flow_float)(e_mix + (double)Pnew/(double)ro_temp + (double)ek);
            sonic[ic]     = (flow_float)sqrt(gmix * Rmix * Tnew);  // 気相 frozen 音速 (loose coupling 近似)
            gam_array[ic] = (flow_float)gmix;
            cp_array[ic]  = (flow_float)cpmix;
            Rmix_array[ic]= (flow_float)Rmix;
        } else {
            // ---- calorically perfect gas ----
            // 非平衡凝縮 (一温度 二相 EOS, CPG): 総液相質量分率 g を集計。condensation==0 で g=0 (従来経路)。
            double g_liq = 0.0;
            if (condensation == 1 && rog != nullptr) {
                for (int s = 0; s < nCondSpecies; ++s) {
                    double gs = (double)rog[s][ic] / (double)ro_temp;
                    if (gs > 0.0) g_liq += gs;
                }
                if (g_liq > 0.99) g_liq = 0.99;   // realizability
                if (g_liq < 0.0)  g_liq = 0.0;
            }

            const double cv = (double)cp/(double)gamma;            // cv = cp/γ
            const double Rgas = ((double)gamma - 1.0)*cv;          // R = cp - cv = (γ-1)cv

            if (g_liq > 1.0e-12) {
                // 二相: e = cv T - g L(T) = intE を Newton で反転、p=(1-g)ρRT。
                const double e_in = (double)intE;
                const double Tguess = (double)max(intE/(cp/gamma), (flow_float)1.0e-4f);
                const double Tn = cond_T_from_e_cpg(e_in, g_liq, cv, Rgas, Tguess);
                const double L = n2_latent(Tn);
                const double e_mix = (cv + g_liq*Rgas)*Tn - g_liq*L;   // = e_in
                const double oneMg = 1.0 - g_liq;
                double Pn = oneMg*(double)ro_temp*Rgas*Tn;
                if (Pn < 1.0) Pn = 1.0;

                T[ic]   = (flow_float)Tn;
                P[ic]   = (flow_float)Pn;
                ro[ic]  = ro_temp;
                roe[ic] = (flow_float)((double)ro_temp*(e_mix + (double)ek));
                Ht[ic]  = (flow_float)(e_mix + Pn/(double)ro_temp + (double)ek);
                sonic[ic] = (flow_float)sqrt((double)gamma*Rgas*Tn); // 気相 frozen 音速 (loose coupling)
                Rmix_array[ic] = (flow_float)Rgas;
            } else {
                // 単相 CPG (従来経路, ビット不変)
                T_temp = max(intE/(cp/gamma), (flow_float)1.0e-4f);
                P_temp = max((gamma-1.0)*(roe[ic]-ro_temp*ek),(flow_float)1.0f);

                T[ic] = T_temp;
                P[ic] = P_temp;

                ro[ic] = ro_temp;
                roe[ic] = P_temp/(gamma-1.0) + ro_temp*ek;

                Ht[ic] = roe[ic]/ro_temp + P_temp/ro_temp;

                sonic[ic] = sqrt(gamma*P_temp/ro_temp);
                // CPG の混合比気体定数 R = cp - cv = (γ-1)cp/γ (定数。SLAU 単成分経路は未使用だが整合のため埋める)
                Rmix_array[ic] = (gamma-1.0)*cp/gamma;
            }
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

        // thermally-perfect 用化学種データ。多成分 (M2, nSpecies>=2) では device roY 配列を渡し、
        // 単成分のときは nullptr (混合則は Y={1} に縮退)。
        thermo_species_device_ptr() , cfg.nSpecies , species_roY_device_ptr() ,

        // 非平衡凝縮 (二相 EOS)。condensation==0 で rog=nullptr/g=0 → 従来経路ビット不変。
        cfg.condensation , var.nCondSpeciesRegistered , cond_rog_device_ptr() ,
        cfg.condGasSpecies , cfg.condModel ,

        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"], var.c_d["roUy"] , var.c_d["roUz"], var.c_d["roe"] ,
        var.c_d["roK"] , var.c_d["roOmega"],
        var.c_d["P"]   , var.c_d["Ht"]  , var.c_d["sonic"], var.c_d["k"], var.c_d["omega"], var.c_d["T"],
        var.c_d["Ux"]  , var.c_d["Uy"]  , var.c_d["Uz"] ,

        var.c_d["gamma"] , var.c_d["cp"] , var.c_d["Rmix"]
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}
