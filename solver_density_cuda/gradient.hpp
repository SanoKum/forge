#pragma once

#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "cuda_forge/calcGradient_d.cuh"

void gradientGauss(solverConfig& , cudaConfig& , mesh& , variables& ); 
