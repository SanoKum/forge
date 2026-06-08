#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

void ransBoundary_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var);