#pragma once

#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

void setMatrixPoisson(solverConfig& , mesh& , variables& , matrix& );

void solvePoisson(solverConfig& , mesh& , variables& , matrix& );

void correctPresVel(solverConfig& , mesh& , variables& );

void correctPresVel_BF(solverConfig& , mesh& , variables& );

void callAmgclCuda(solverConfig& , mesh& , variables& , matrix& );
