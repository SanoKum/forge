#pragma once

#include <algorithm>
#include <cmath>
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"
#include "cuda_forge/setDT_d.cuh"

void setDT(solverConfig&, cudaConfig&, mesh& , variables& );
