#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void interpVelocity_c2p_d
( 
 flow_float dt,
 // mesh structure
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* Ux ,
 flow_float* Uy ,
 flow_float* Uz ,
 flow_float* ro ,
 flow_float* P  ,
 
 flow_float* divU_vol ,
 
 flow_float* USp,

 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz
);

__global__ void interpVelocity_c2p_outlet_d
( 
 flow_float dt,
 // mesh structure
 geom_int nb, 
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* Ux ,
 flow_float* Uy ,
 flow_float* Uz ,
 flow_float* ro ,
 flow_float* P  ,
 
 flow_float* divU_vol ,
 
 flow_float* USp,
 flow_float* Pp,

 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz
);

__global__ void interpVelocity_c2p_inlet_d
( 
 flow_float dt,
 // mesh structure
 geom_int nb, 
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* Uxp,
 flow_float* Uyp,
 flow_float* Uzp,
 
 flow_float* divU_vol ,
 flow_float* USp
);

void interpVelocity_c2p_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns);
