#pragma once 

#include <cuda_runtime.h>
#include <vector_types.h>
#include "mesh/mesh.hpp"
#include <cmath>

struct cudaConfig
{
public:
    int blocksize = 512;
    int blocksize_small = 256;
    dim3 dimBlock;
    dim3 dimGrid_cell;
    dim3 dimGrid_normalcell;
    dim3 dimGrid_plane;
    dim3 dimGrid_nplane;
    dim3 dimGrid_bplane;

    dim3 dimBlock_small;
    dim3 dimGrid_normalcell_small;

    cudaConfig(mesh&);
};

inline cudaConfig::cudaConfig(mesh& msh)
{
    this->dimBlock       = dim3(blocksize);
    this->dimBlock_small = dim3(blocksize_small);
    this->dimGrid_cell   = dim3(ceil( msh.nCells_all / (flow_float)blocksize));
    this->dimGrid_normalcell = dim3(ceil( msh.nCells / (flow_float)blocksize));
    this->dimGrid_normalcell_small = dim3(ceil( msh.nCells / (flow_float)blocksize_small));
    this->dimGrid_plane  = dim3(ceil( msh.nPlanes/ (flow_float)blocksize));
    this->dimGrid_nplane = dim3(ceil( msh.nNormalPlanes/ (flow_float)blocksize));
    this->dimGrid_bplane = dim3(ceil( msh.nBPlanes/ (flow_float)blocksize));
};
