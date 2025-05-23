#include <iostream>
#include <vector>
#include <stdio.h>                                                                                       
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "input/solverConfig.hpp"
#include "input/setInitial.hpp"

#include "mesh/gmshReader.hpp"
#include "output/output.hpp"
#include "boundaryCond.hpp"

#include "setStructualVariables.hpp"

#include "gradient.hpp"

#include "dependentVariables.hpp"
//#include "solvePoisson_amgx.hpp"

#include "convectiveFlux.hpp"
#include "implicitCorrection.hpp"
//#include "timeIntegration.hpp"
#include "update.hpp"

#include "common/stringUtil.hpp"
#include "common/vectorUtil.hpp"

//#include "setDT.hpp"

// cuda
#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/cudaConfig.cuh"

#include "cuda_forge/calcGradient_d.cuh"
#include "cuda_forge/convectiveFlux_d.cuh"
#include "cuda_forge/viscousFlux_d.cuh"
#include "cuda_forge/updateCenterVelocity_d.cuh"
#include "cuda_forge/interpVelocity_c2p_d.cuh"
#include "cuda_forge/timeIntegration_d.cuh"
#include "cuda_forge/limiter_d.cuh"
#include "cuda_forge/ducrosSensor_d.cuh"
#include "cuda_forge/turbulent_viscosity_d.cuh"

#include "cuda_forge/fluct_variables_d.cuh"
#include "cuda_forge/gasProperties_d.cuh"

#include "probe/point_probes.cuh"
#include "cuda_forge/setDT_d.cuh"

#include <cuda_runtime.h>


#define CHECK_LAST_CUDA_ERROR() checkLast(__FILE__, __LINE__)
void checkLast(const char* const file, const int line)
{
    cudaError_t err{cudaGetLastError()};
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << std::endl;
        // We don't exit when we encounter CUDA errors in this example.
        std::exit(EXIT_FAILURE);
    }
}

int main(void) {
    clock_t start = clock();

    cout << "Read Solver Config \n";
    solverConfig cfg; 
    cfg.read("solverConfig.yaml");

    cout << "Read Mesh \n";
    mesh msh;

    if (cfg.meshFormat == "hdf5") {
        msh.readMesh(cfg.meshFileName);
    } else {
        cerr << "Error unknown mesh format: " << cfg.meshFormat << endl;
        return 1;
    }

    cout << "Init Matrix (but not used now) \n";
    matrix mat_ns = matrix();
    mat_ns.initMatrix(msh);


    cout << "Read Boundary Conditions \n";
    readBcondConfig(cfg , msh.bconds);

    variables var = variables();
    var.allocVariables(cfg.gpu , msh);

    cout << "Read Initial Values \n";
    var.readValueHDF5(cfg.valueFileName , msh);

    cout << "Set mesh connection map for cuda \n";
    cudaConfig cuda_cfg = cudaConfig(msh);
    msh.setMeshMap_d();

    msh.setPeriodicPartner(); 

    var.setStructuralVariables(cfg , cuda_cfg , msh);

    dependentVariables(cfg , cuda_cfg , msh , var, mat_ns); 

    gasProperties_d_wrapper(cfg , cuda_cfg , msh , var);


    // fluctuating velocity inlet
    fluct_variables fluct = fluct_variables();
    fluct.allocVariables();
    fluct.set_fluctVelocity(cfg , cuda_cfg , msh , var );


    applyBconds(cfg , cuda_cfg , msh , var, mat_ns , fluct);

    calcGradient_d_wrapper(cfg , cuda_cfg , msh , var);

    //copyBcondsGradient(cfg , cuda_cfg , msh , var, mat_ns);

    updateVariablesOuter(cfg , cuda_cfg , msh , var , mat_ns);

    setDT_d_wrapper(cfg , cuda_cfg , msh , var);

    point_probes pprobes;
    pprobes.init(cfg , cuda_cfg , msh);


    outputH5_XDMF(cfg , msh, var, 0);

    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep <cfg.nStep ; iStep++) {

        cout << "----------------------------\n";
        cout << "Step : " << iStep << "  Time : " << cfg.totalTime << "\n";

        for (int iloop = 0 ; iloop <cfg.nStage ; iloop++) {
            updateVariablesInner(cfg , cuda_cfg ,msh , var , mat_ns);

            cout << "       Stage : " << iloop+1 << "\n" ;
            dependentVariables(cfg , cuda_cfg , msh , var, mat_ns);

            gasProperties_d_wrapper(cfg , cuda_cfg , msh , var);

            //applyBconds(cfg , cuda_cfg , msh , var, mat_ns);
            applyBconds(cfg , cuda_cfg , msh , var, mat_ns , fluct);
            
            calcGradient_d_wrapper(cfg , cuda_cfg , msh , var);

            limiter_d_wrapper(cfg , cuda_cfg , msh , var);

            ducrosSensor_d_wrapper(cfg , cuda_cfg , msh , var);

            turbulent_viscosity_d_wrapper(cfg , cuda_cfg , msh , var);

            convectiveFlux_d_wrapper(cfg , cuda_cfg, msh , var, mat_ns);

            viscousFlux_d_wrapper(cfg , cuda_cfg, msh , var, mat_ns);
            
            implicitCorrection(iloop, cfg , msh , var, mat_ns);

            timeIntegration_d_wrapper(iloop, cfg , cuda_cfg , msh , var);
        }

        updateVariablesOuter(cfg , cuda_cfg , msh , var , mat_ns);

        outputH5_XDMF(cfg , msh, var, iStep+1);
        outputBconds_H5_XDMF(cfg , msh, var, iStep+1);
        pprobes.outputProbes(cfg , cuda_cfg , msh , var , iStep+1);

        setDT_d_wrapper(cfg , cuda_cfg, msh , var);

        cfg.totalTime += cfg.dt;
    }

    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time = %.3f s\n", time); 

	return 0;
}