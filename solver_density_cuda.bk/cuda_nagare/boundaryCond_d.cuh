#pragma once

#include "cuda_nagare/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void slip_d
(
 geom_int nb,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P
);

void slip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);

__global__ void wall_d 
( 
 geom_int nb,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,
 // bvar
 flow_float* Uxb ,  
 flow_float* Uyb ,  
 flow_float* Uzb   
);

void wall_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);

__global__ void outlet_statPress_d 
( 
 geom_int nb,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,
 // bvar
 flow_float* Uxb ,  
 flow_float* Uyb ,  
 flow_float* Uzb   
);

void outlet_statPress_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);

__global__ void inlet_uniformVelocity_statPress_d 
( 
 geom_int nb,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,
 // bvar
 flow_float* Uxb ,  
 flow_float* Uyb ,  
 flow_float* Uzb ,  
 flow_float* Tb ,  
 flow_float* rob 
);

void inlet_uniformVelocity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);
