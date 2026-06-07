#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

void scalarTransport_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
void scalarTimeIntegration_d_wrapper(int loop , solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
void calcScalarGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
