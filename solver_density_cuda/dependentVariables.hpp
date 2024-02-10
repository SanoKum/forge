#pragma once

#include <cmath>

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"
#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/dependentVariables_d.cuh"

void dependentVariables(solverConfig& , cudaConfig& , mesh& , variables& , matrix& );
