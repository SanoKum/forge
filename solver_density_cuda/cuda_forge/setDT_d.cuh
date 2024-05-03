#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include <cuda_runtime.h>
#include <cublas_v2.h>

//#include <thrust/extream.h>
//#include <thrust/device_ptr.h>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
//#include "cuda_forge/derived_atomic_functions.cuh"
//#include "cuda_forge/atomic.cuh"

#include "input/solverConfig.hpp"
#include "variables.hpp"

//#include "atomic.cuh"

__global__ void setCFL_pln_d
( 
 flow_float dt,
 flow_float dt_pseudo,
 flow_float visc,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* cfl ,
 flow_float* cfl_pseudo ,
 flow_float* sonic,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 //plane variables
 flow_float* cfl_pln ,
 flow_float* cfl_pseudo_pln
);

__global__ void setCFL_cell_d
( 
 int dtControl, flow_float cfl_target,
 flow_float dt,
 flow_float dt_pseudo,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes,geom_int* cell_planes_index, geom_int* cell_planes,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* cfl ,
 flow_float* cfl_pseudo ,
 flow_float* sonic,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 flow_float* cfl_pln ,
 flow_float* cfl_pseudo_pln 
);

void setDT_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
