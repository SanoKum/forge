#include "cuda_nagare/cudaWrapper.cuh"

#include "variables.hpp"

#include <iostream>
#include <vector>
#include <list>

namespace cudaWrapper {
    // Malloc
    void cudaMalloc_wrapper(flow_float** var_d , geom_int size)
    {
        CHECK_CUDA_ERROR( cudaMalloc(var_d , size*sizeof(flow_float)) );
    };

    void cudaMalloc_wrapper(geom_int** var_d , geom_int size)
    {
        CHECK_CUDA_ERROR( cudaMalloc(var_d , size*sizeof(geom_int)) );
    };

    // Memcpy
    //void cudaMemcpy_vectorToDevice_wrapper(std::vector<flow_float>& vec , flow_float* var_d)
    //{
    //    CHECK_CUDA_ERROR( cudaMemcpy(var_d, vec.data() , vec.size()*sizeof(flow_float), cudaMemcpyHostToDevice) );
    //};

    //void cudaMemcpy_deviceToVector_wrapper(flow_float* var_d , std::vector<flow_float>& vec )
    //{
    //    CHECK_CUDA_ERROR( cudaMemcpy(vec.data() , var_d, vec.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    //};

    void cudaMemcpy_H2D_wrapper(flow_float* vec , flow_float* var_d , geom_int numEle)
    {
        CHECK_CUDA_ERROR( cudaMemcpy(var_d, vec , (size_t)(numEle*sizeof(flow_float)), cudaMemcpyHostToDevice) );
    };

    void cudaMemcpy_D2H_wrapper(flow_float* var_d , flow_float* vec , geom_int numEle)
    {
        CHECK_CUDA_ERROR( cudaMemcpy(vec, var_d, numEle*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    };


    // free
    void cudaFree_wrapper(flow_float* var_d)
    {
        CHECK_CUDA_ERROR( cudaFree(var_d) );
    };

    void cudaFree_wrapper(geom_int* var_d)
    {
        CHECK_CUDA_ERROR( cudaFree(var_d) );
    };

    //void copyVariables_cell_plane_H2D(variables& var)
    //{
    //    for (auto& name : var.cellValNames)
    //    {
    //        cudaMemcpy_vectorToDevice_wrapper(var.c[name] , var.c_d[name]);
    //    }
    //    for (auto& name : var.planeValNames)
    //    {
    //        cudaMemcpy_vectorToDevice_wrapper(var.p[name] , var.p_d[name]);
    //    }
    //}
};
