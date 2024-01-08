#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"

void solveNavierStokes(solverConfig& , mesh& , variables& , matrix& );
