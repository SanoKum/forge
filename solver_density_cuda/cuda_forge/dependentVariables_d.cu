#include "dependentVariables_d.cuh"

__global__ void dependentVariables_d
( 
 // gas properties
 flow_float gamma , flow_float cp , 

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
 flow_float* Uz  
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
        T_temp = max(intE/(cp/gamma), (flow_float)1.0e-4f);
        P_temp = max((gamma-1.0)*(roe[ic]-ro_temp*ek),(flow_float)1.0f);

        T[ic] = T_temp;
        P[ic] = P_temp;

        ro[ic] = ro_temp;
        roe[ic] = P_temp/(gamma-1.0) + ro_temp*ek;

        Ht[ic] = roe[ic]/ro_temp + P_temp/ro_temp;
        k[ic] = max(roK[ic]/ro_temp, static_cast<flow_float>(0.0));
        omega[ic] = max(roOmega[ic]/ro_temp, static_cast<flow_float>(0.0));

        sonic[ic] = sqrt(gamma*P_temp/ro_temp);

    }
}


void dependentVariables_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    dependentVariables_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // gas properties
        cfg.gamma , cfg.cp , 

        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"], var.c_d["roUy"] , var.c_d["roUz"], var.c_d["roe"] ,
        var.c_d["roK"] , var.c_d["roOmega"],
        var.c_d["P"]   , var.c_d["Ht"]  , var.c_d["sonic"], var.c_d["k"], var.c_d["omega"], var.c_d["T"], 
        var.c_d["Ux"]  , var.c_d["Uy"]  , var.c_d["Uz"] 
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}