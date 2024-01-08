#include <iostream>
#include <vector>
#include <stdio.h>                                                                                       
#include <fstream>
#include <string>
#include <sstream>

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

#include "solvePoisson_amgcl.hpp"
#include "solveNavierStokes.hpp"
#include "update.hpp"

#include "common/stringUtil.hpp"
#include "common/vectorUtil.hpp"

#include "calcCFL.hpp"


int main(void) {
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

    cout << "Read Boundary Conditions \n";
    readAndSetBcondConfig(msh.bconds);

    cout << "Set Initial Values \n";
    cout << "a\n";
    variables var = variables(msh);
    cout << "b\n";
    setInitial(cfg , msh , var);
    cout << "c\n";

    setStructualVariables(cfg , msh , var);

    updateVariablesForNextLoop(cfg , msh , var , mat_poi);

    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep <cfg.nStep ; iStep++) {

        cout << " Step =" << iStep << "\n";
        setBcondsValue(cfg , msh , var , mat_poi);

        calcGradient(cfg , msh , var);

        solveNavierStokes(cfg , msh , var , mat_ns);

        setMatrixPoisson(cfg , msh , var , mat_poi);

        //solvePoisson(cfg , msh , var , mat_poi);

        if (cfg.gpu == 0) {
            solvePoisson(cfg , msh , var , mat_poi);
        } else if (cfg.gpu == 1) {
            callAmgclCuda(cfg, msh, var, mat_poi );
        } else {
            cerr << "unknown gpu option : " << cfg.gpu << endl;
            return 1;
        }

        correctPresVel(cfg , msh , var);

        // *** output hdf5 ***
        outputH5_XDMF(cfg , msh, var, iStep);

        // *** update TN, diffT, etc.***
        updateVariablesForNextLoop(cfg , msh , var , mat_poi);

        calcCFL(cfg , msh , var);
    }
    
	return 0;
}