#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void slip_d
(
 // gas properties
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic  

);

void slip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);

__global__ void wall_d 
( 
 // gas properties
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 geom_float* fx  , 

 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic  
);

void wall_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);

__global__ void outlet_statPress_d 
( 
 // gas properties
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic  ,

 // bvar
 flow_float* Psb ,
 flow_float* Ptb ,
 flow_float* Ttb 
);

void outlet_statPress_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);

__global__ void inlet_uniformVelocity_statPress_d 
( 
 // gas 
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,

 // bvar
 flow_float* rob ,  
 flow_float* Uxb ,  
 flow_float* Uyb ,  
 flow_float* Uzb , 
 flow_float* Psb  
);

void inlet_uniformVelocity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);


__global__ void inlet_Pressure_d 
( 
 // gas 
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* T   ,

 // bvar
 flow_float* Ptb ,
 flow_float* Ttb 
);

void inlet_Pressure_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);



__global__ void inlet_Pressure_dir_d 
( 
 // gas 
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* T   ,

 // bvar
 flow_float* Ptb ,
 flow_float* Ttb ,
 flow_float* Uxb ,
 flow_float* Uyb ,
 flow_float* Uzb 
);

void inlet_Pressure_dir_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);



__global__ void outflow_d
(
 // gas properties
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,

 // bvar
 flow_float* Ptb ,
 flow_float* Ttb 

);

void outflow_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);


__global__ 
void copyBcondsGradient_d 
( 
 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 // variables
 flow_float* dTdx,
 flow_float* dTdy,
 flow_float* dTdz,

 flow_float* dHtdx,
 flow_float* dHtdy,
 flow_float* dHtdz,

 flow_float* drodx,
 flow_float* drody,
 flow_float* drodz,

 flow_float* dUxdx,
 flow_float* dUydx,
 flow_float* dUzdx,

 flow_float* dUxdy,
 flow_float* dUydy,
 flow_float* dUzdy,
 
 flow_float* dUxdz,
 flow_float* dUydz,
 flow_float* dUzdz,
 
 flow_float* dPdx,
 flow_float* dPdy,
 flow_float* dPdz
);

void copyBcondsGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p);