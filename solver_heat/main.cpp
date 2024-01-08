#include <iostream>
#include <vector>
#include <stdio.h>                                                                                       
#include <fstream>
#include <string>
#include <sstream>

#include "flowFormat.hpp"
#include "mesh.hpp"
#include "variables.hpp"

#include "readConfig.hpp"

#include "gmshReader.hpp"
#include "output.hpp"
#include "boundaryCond.hpp"

#include "diffusion.hpp"
#include "gradient.hpp"
#include "timeIntegration.hpp"
#include "update.hpp"

#include "solvePoisson_amgcl.hpp"

int main(void) {
    // *** read solver config ***
    solverConfig cfg = solverConfig();

    // *** read mesh ***
    gmshReader gmsh = gmshReader(cfg.meshFileName);
    //mesh msh = gmsh.getMesh();

    // *** init matrix ***
    //matrix mat_poi = matrix();
    //mat_poi.initPoisson(msh);

    // *** read boundary conditon config ***
    readAndSetBcondConfig(gmsh.bconds);

    // *** set initial value ***
    cout << "Set Initial Values\n";
    variables var = variables(gmsh);
    for (geom_int i = 0 ; i<gmsh.nCells ; i++)
    {
        var.c["T"][i] = 300.0;
    }

    gmsh.writeGeometryInput("inputGeometry" , var);


//    updateVariablesForNextLoop(cfg , msh , var , mat_poi);
//
//    cout << "Start Calculation\n";
//
//    for (int iStep = 0 ; iStep <cfg.nStep ; iStep++)
//    {
//        cout << " Step =" << iStep << "\n";
//        // *** boundary conditions ***
//        setBcondsValue(cfg , msh , var , mat_poi);
//
//        // *** calc gradient ***
//        calcGradient(cfg , msh , var);
//
//        // *** calc ***
//        calcDiffusion(cfg , msh , var , mat_poi);
//
//        solvePoisson(cfg , msh , var , mat_poi);
//
//        // *** time integration ***
//        timeIntegration(cfg , msh , var);
//
//        // *** update TN, diffT, etc.***
//        updateVariablesForNextLoop(cfg , msh , var , mat_poi);
//
//        // *** output hdf5 ***
//        outputH5_XDMF(cfg , msh, var, iStep);
//    }

    
	return 0;
}