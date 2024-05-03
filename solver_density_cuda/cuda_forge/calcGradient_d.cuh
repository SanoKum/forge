#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void calcGradient_1_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* nplane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* Ht  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,
 flow_float* divU   

);


__global__ void calcGradient_2_d
( 
 // mesh structure
 geom_int nCells,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* Ht  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,
 flow_float* divU   
);

__device__ flow_float calcDeltaIJ(geom_float pcx , geom_float pcy, geom_float pcz, 
                                  flow_float dudx, flow_float dudy,flow_float dudz,
                                  flow_float delu_max, flow_float delu_min );

__device__ flow_float venkata_limiter(flow_float x , flow_float volume);


void calcGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
