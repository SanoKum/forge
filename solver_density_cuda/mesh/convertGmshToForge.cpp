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

    cout << "Read Solver Config \n";
    solverConfig cfg = solverConfig();
    string fname = "solverConfig.yaml";
    cfg.read(fname);

    cout << "Read Mesh \n";
    gmshReader gmsh = gmshReader(argv[1]);

    cout << "Read Boundary Conditions \n";
    readBcondConfig(cfg , gmsh.bconds);

    cout << "Set Initial Values \n";
    variables var = variables();
    var.allocVariables(cfg.gpu , gmsh);
    setInitial(cfg , gmsh , var);

    cout << "Write Input HDF5 \n";
    gmsh.writeInputH5(argv[2] , var);

    return 0;
}