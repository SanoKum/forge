#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void dependentVariables_d
( 
 // gas properties
 geom_int gamma , geom_int cp , 

 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* T   ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  

);

void dependentVariables_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
