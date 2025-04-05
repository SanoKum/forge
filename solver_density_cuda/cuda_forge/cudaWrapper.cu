#include "cuda_forge/cudaWrapper.cuh"

#include "variables.hpp"

#include <iostream>
#include <vector>
#include <list>

namespace cudaWrapper {
    // Malloc
    void cudaMalloc_wrapper(flow_float** var_d , geom_int size)
    {
        gpuErrchk( cudaMalloc(var_d , size*sizeof(flow_float)) );
    };

    void cudaMalloc_wrapper(geom_int** var_d , geom_int size)
    {
        gpuErrchk( cudaMalloc(var_d , size*sizeof(geom_int)) );
    };

    void cudaMemcpy_H2D_wrapper(flow_float* vec , flow_float* var_d , geom_int numEle)
    {
        gpuErrchk( cudaMemcpy(var_d, vec , (size_t)(numEle*sizeof(flow_float)), cudaMemcpyHostToDevice) );
    };

    void cudaMemcpy_H2D_wrapper(geom_int* vec , geom_int* var_d , geom_int numEle)
    {
        gpuErrchk( cudaMemcpy(var_d, vec , (size_t)(numEle*sizeof(geom_int)), cudaMemcpyHostToDevice) );
    };


    void cudaMemcpy_D2H_wrapper(flow_float* var_d , flow_float* vec , geom_int numEle)
    {
        gpuErrchk( cudaMemcpy(vec, var_d, (size_t)(numEle*sizeof(flow_float)), cudaMemcpyDeviceToHost) );
    };

    void cudaMemcpy_D2H_wrapper(geom_int* var_d , geom_int* vec , geom_int numEle)
    {
        gpuErrchk( cudaMemcpy(vec, var_d, (size_t)(numEle*sizeof(geom_int)), cudaMemcpyDeviceToHost) );
    };


    // free
//    void cudaFree_wrapper(flow_float* var_d)
//    {
//        int is_dev_ptr = is_device_pointer((const void *) var_d);
//        if (is_dev_ptr == 1) {
//            gpuErrchk( cudaFree(var_d) );
//        }
//    };
//

    void cudaFree_wrapper(int* var_d)
    {
        int is_dev_ptr = is_device_pointer((const void *) var_d);
        if (is_dev_ptr == 1) {
            gpuErrchk( cudaFree(var_d) );
        }
    };

    void cudaFree_wrapper(long* var_d)
    {
        int is_dev_ptr = is_device_pointer((const void *) var_d);
        if (is_dev_ptr == 1) {
            gpuErrchk( cudaFree(var_d) );
        }
    };

    void cudaFree_wrapper(float* var_d)
    {
        int is_dev_ptr = is_device_pointer((const void *) var_d);
        if (is_dev_ptr == 1) {
            gpuErrchk( cudaFree(var_d) );
        }
    };

    void cudaFree_wrapper(double* var_d)
    {
        int is_dev_ptr = is_device_pointer((const void *) var_d);
        if (is_dev_ptr == 1) {
            gpuErrchk( cudaFree(var_d) );
        }
    };

    int is_device_pointer(const void *ptr)
    {
        int is_device_ptr = 0;
        cudaPointerAttributes attributes;
    
        gpuErrchk(cudaPointerGetAttributes(&attributes, ptr));
    
        if(attributes.devicePointer != NULL)
        {
          is_device_ptr = 1;
        }
    
        return is_device_ptr;
    }
};
