#pragma once

#include <iostream>
#include <vector>
#include <list>

#include <cuda_runtime.h>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"

#include "variables.hpp"

//TODO: use template???

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define CHECK_CUDA_ERROR(val) check((val), #val, __FILE__, __LINE__)
template <typename T>
void check(T err, const char* const func, const char* const file,
           const int line)
{
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
        // We don't exit when we encounter CUDA errors in this example.
        // std::exit(EXIT_FAILURE);
    }
}



namespace cudaWrapper{
   void cudaMalloc_wrapper(flow_float** , geom_int );
   void cudaMalloc_wrapper(geom_int** , geom_int );

   void cudaMemcpy_H2D_wrapper(flow_float* , flow_float* , geom_int );
   void cudaMemcpy_H2D_wrapper(geom_int*   , geom_int* , geom_int );
   void cudaMemcpy_D2H_wrapper(flow_float* , flow_float* , geom_int );
   void cudaMemcpy_D2H_wrapper(geom_int*   , geom_int* , geom_int );

   //void cudaFree_wrapper(geom_int* );
   //void cudaFree_wrapper(flow_float* );
   void cudaFree_wrapper(int* );
   void cudaFree_wrapper(long* );
   void cudaFree_wrapper(float* );
   void cudaFree_wrapper(double* );

   int is_device_pointer(const void *ptr);

   //void copyVariables_cell_plane_H2D(variables& );
}
