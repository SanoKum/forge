#pragma once 

#include <cuda_runtime.h>
#include <vector_types.h>
#include "mesh/mesh.hpp"
#include <cmath>

struct cudaConfig
{
public:
    int blocksize = 512;
    dim3 dimBlock;
    dim3 dimGrid_cell;
    dim3 dimGrid_plane;
    dim3 dimGrid_nplane;
    dim3 dimGrid_bplane;

    cudaConfig(mesh&);
};

inline cudaConfig::cudaConfig(mesh& msh)
{
    this->dimBlock       = dim3(blocksize);
    this->dimGrid_cell   = dim3(ceil( msh.nCells / (float)blocksize));
    this->dimGrid_plane  = dim3(ceil( msh.nPlanes/ (float)blocksize));
    this->dimGrid_nplane = dim3(ceil( msh.nNormalPlanes/ (float)blocksize));
    this->dimGrid_bplane = dim3(ceil( msh.nBPlanes/ (float)blocksize));
};
