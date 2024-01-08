#pragma once

#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

//#include "cuda_nagare/cudaWrapper.cuh"
#include "cuda_nagare/cudaConfig.cuh"

__global__ void calcStructualVariables_d
//void calcStructualVariables_d
(
 geom_int , geom_int,
 geom_int* ,  
 geom_float* ,  geom_float* ,  geom_float* ,  geom_float* ,  
 geom_float* ,  geom_float* ,  geom_float* ,  
 geom_float* ,  geom_float* ,  geom_float* ,  
 geom_float* ,  geom_float*
);

__global__ void calcStructualVariables_bp_d
//void calcStructualVariables_d
(
 geom_int , 
 geom_int* ,  
 geom_float* ,  geom_float* ,  geom_float* ,  geom_float* ,  
 geom_float* ,  geom_float* ,  geom_float* ,  
 geom_float* ,  geom_float* ,  geom_float* ,  
 geom_float* ,  geom_float*
);


void calcStructualVariables_d_wrapper( cudaConfig& , mesh& ,  variables& );

