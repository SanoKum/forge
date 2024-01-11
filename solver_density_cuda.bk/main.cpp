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
#include "timeIntegration.hpp"
#include "update.hpp"

#include "common/stringUtil.hpp"
#include "common/vectorUtil.hpp"

#include "setDT.hpp"

// cuda
#include "cuda_nagare/cudaWrapper.cuh"
#include "cuda_nagare/cudaConfig.cuh"

#include "cuda_nagare/calcGradient_d.cuh"
#include "cuda_nagare/calcConvDiff_d.cuh"
#include "cuda_nagare/updateCenterVelocity_d.cuh"
#include "cuda_nagare/interpVelocity_c2p_d.cuh"

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
    } else if (cfg.meshFormat == "gmsh") {
        gmshReader gmsh = gmshReader(cfg.meshFileName);
        msh = gmsh.getMesh();
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

    cout << "Set Initial Values \n";
    setInitial(cfg , msh , var);

    cout << "Set mesh connection map for cuda \n";
    cudaConfig cuda_cfg = cudaConfig(msh);
    msh.setMeshMap_d();

    msh.setPeriodicPartner();

    setStructualVariables(cfg , msh , var);

    dependentVariables(cfg , msh , var, mat_ns);

    applyBconds(cfg , cuda_cfg , msh , var, mat_ns);

    updateVariablesOuter(cfg , msh , var , mat_ns);

    setDT(cfg , msh , var);

    outputH5_XDMF(cfg , msh, var, 0);

    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep <cfg.nStep ; iStep++) {

        cout << "----------------------------\n";
        cout << "Step : " << iStep << "\n";

        for (int iloop = 0 ; iloop <cfg.nLoop ; iloop++) {
            updateVariablesInner(cfg , msh , var , mat_ns);

            cout << "  Inner loop : " << iloop+1 << "\n" ;
            dependentVariables(cfg , msh , var, mat_ns);

            applyBconds(cfg , cuda_cfg , msh , var, mat_ns);
            
            gradientGauss(cfg , msh , var);

            convectiveFlux(iloop, cfg , msh , var, mat_ns);
            
            implicitCorrection(iloop, cfg , msh , var, mat_ns);

            timeIntegration(iloop, cfg , msh , var, mat_ns);
        }

        updateVariablesOuter(cfg , msh , var , mat_ns);

        outputH5_XDMF(cfg , msh, var, iStep+1);

        setDT(cfg , msh , var);
    }

    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time = %.3f s\n", time); 

	return 0;
}