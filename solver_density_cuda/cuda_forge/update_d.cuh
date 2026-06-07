#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void updateVariablesOuter_d
( 
 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* roK  ,
 flow_float* roOmega  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,
 flow_float* roKN ,
 flow_float* roOmegaN ,
 
 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM 
);

void updateVariablesOuter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

__global__ void updateVariablesInner_d
( 
 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* roK  ,
 flow_float* roOmega  ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM 
);

void updateVariablesInner_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

__global__ void applyScalarImplicitCorrection_d
(
 geom_int nCells_all , geom_int nCells,

 flow_float* ro,
 flow_float* roUx,
 flow_float* roUy,
 flow_float* roUz,
 flow_float* roe,

 flow_float* roN,
 flow_float* roUxN,
 flow_float* roUyN,
 flow_float* roUzN,
 flow_float* roeN,

 flow_float* dq_ro,
 flow_float* dq_roUx,
 flow_float* dq_roUy,
 flow_float* dq_roUz,
 flow_float* dq_roe
);

void applyScalarImplicitCorrection_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
