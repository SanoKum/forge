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

//#include "solvePoisson_amgcl.hpp"
#include "solvePoisson_amgx.hpp"
#include "solveNavierStokes.hpp"
#include "update.hpp"

#include "common/stringUtil.hpp"
#include "common/vectorUtil.hpp"

#include "calcCFL.hpp"

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

    cout << "Init Matrix \n";
    matrix mat_poi = matrix();
    mat_poi.initMatrix(msh);

    matrix mat_ns = matrix();
    mat_ns.initMatrix(msh);

    AMGX_SAFE_CALL(AMGX_initialize());


    cout << "Read Boundary Conditions \n";
    readAndSetBcondConfig(msh.bconds);

    cout << "Set Initial Values \n";
    variables v = variables(msh);
    setInitial(cfg , msh , v);

    updateVariablesForNextLoop(cfg , msh , v , mat_poi);
    v.copyVariables_cell_plane_H2D_all();
    CHECK_LAST_CUDA_ERROR();

    cout << "Set mesh map \n";
    // prepare for cuda
    cudaConfig cuda_cfg = cudaConfig(msh);
    msh.setMeshMap_d();
    CHECK_LAST_CUDA_ERROR();

    cout << "Set structual variables \n";
    v.setStructualVariables_d(cuda_cfg , msh );
    CHECK_LAST_CUDA_ERROR();

    setStructualVariables(cfg , msh , v);
    CHECK_LAST_CUDA_ERROR();

    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep <cfg.nStep ; iStep++) {

        cout << " Step =" << iStep << "\n";
        setBcondsValue(cfg , cuda_cfg , msh , v , mat_poi);

        if      (cfg.gpu == 0) { calcGradient(cfg , msh , v ); } 
        else if (cfg.gpu == 1) { calcGradient_d_wrapper(cfg , cuda_cfg , msh , v ); }

        if      (cfg.gpu == 0) { 
            solveNavierStokes(cfg , msh, v , mat_ns);

        } else if (cfg.gpu == 1) { 
            // solve navier stokes 
            calcConvDiff_d_wrapper(cfg , cuda_cfg , msh , v , mat_ns);
            updateCenterVelocity_d_wrapper(cfg , cuda_cfg , msh , v , mat_ns);
            // rhi-chow interpolation
            interpVelocity_c2p_d_wrapper(cfg , cuda_cfg , msh , v , mat_ns);
        }

        if (cfg.gpu == 1) v.copyVariables_cell_plane_D2H_all();

        // set Matrix for poisson equation
        setMatrixPoisson(cfg , msh , v , mat_poi);

        solvePoisson(cfg , msh , v , mat_poi);

        correctPresVel(cfg , msh , v);

        // *** output hdf5 ***
        outputH5_XDMF(cfg , msh, v, iStep);

        // *** update TN, diffT, etc.***
        updateVariablesForNextLoop(cfg , msh , v , mat_poi);

        calcCFL(cfg , msh , v);

        if (cfg.gpu == 1 || cfg.gpu == -1) v.copyVariables_cell_plane_H2D_all();

        cudaThreadSynchronize();
    }

    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time = %.3f s\n", time); 

	return 0;
}