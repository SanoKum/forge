#include <iostream>
#include <vector>

#include "input/solverConfig.hpp"
#include "input/setInitial.hpp"

#include "mesh/mesh.hpp"
#include "mesh/gmshReader.hpp"
#include "boundaryCond.hpp"

using namespace std;

int main(int argc , char *argv[]) 
{
    if (argc != 3) 
    {
        cerr << "usage: convertGmshToNagare gmshFileName inputMeshName \n"; 
    }

    cout << "-------------------------- \n";
    cout << "*** Read Solver Config *** \n";
    cout << "-------------------------- \n";
    solverConfig cfg = solverConfig();
    string fname = "solverConfig.yaml";
    cfg.read(fname);

    cout << "----------------- \n";
    cout << "*** Read Mesh *** \n";
    cout << "----------------- \n";
    gmshReader gmsh = gmshReader(argv[1]);

    cout << "-------------------------------- \n";
    cout << "*** Read Boundary Conditions *** \n";
    cout << "-------------------------------- \n";
    readBcondConfig(cfg , gmsh.bconds);

    cout << "-------------------------- \n";
    cout << "*** Set Initial Values *** \n";
    cout << "-------------------------- \n";
    variables var = variables();
    var.allocVariables(cfg.gpu , gmsh);
    setInitial(cfg , gmsh , var);

    cout << "------------------------ \n";
    cout << "*** Write Input HDF5 *** \n";
    cout << "------------------------ \n";
    gmsh.writeInputH5(argv[2] , var);

    return 0;
}