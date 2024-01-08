#pragma once

#include "cuda_nagare/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"
__global__ void updateCenterVelocity_d
( 
 flow_float dt,
 // mesh structure
 geom_int nCells,
 geom_float* vol , 

 // variables
 flow_float* convx , flow_float* convy , flow_float* convz,
 flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* Ux  , flow_float* UxN ,
 flow_float* Uy  , flow_float* UyN ,
 flow_float* Uz  , flow_float* UzN ,
 flow_float* ro  ,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz
);

void updateCenterVelocity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns);
