#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__device__ void calcDeltaIJ(geom_float pcx , geom_float pcy, geom_float pcz, 
                            flow_float dudx, flow_float dudy,flow_float dudz,
                            flow_float delu_max, flow_float delu_min ,
                            flow_float& delta, flow_float& deltap, flow_float& deltam);

__device__ flow_float venkata_limiter(flow_float delp, flow_float delm , flow_float volume);


void limiter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);