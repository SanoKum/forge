#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <cstdlib>

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

// per-kernel の cudaDeviceSynchronize は正しさには不要 (同一 default stream のカーネルは順序保証され、
// host が読む値は cudaMemcpy/thrust が intrinsic に同期する)。よって既定では同期を行わず (peek 化)、
// host を待たせず後続カーネルを queue させて per-step の stall を消す。失うのは実行エラーの「即時・行単位の
// 検知」だけ (次の同期点で表面化)。デバッグ時は FORGE_KERNEL_SYNC=1 で従来どおり毎カーネル同期に戻せる。
// 起動エラーの検知は呼び出し側の gpuErrchk(cudaPeekAtLastError()) が引き続き担う。
inline bool forgeKernelSyncEnabled()
{
    static const bool enabled = [](){
        const char* e = std::getenv("FORGE_KERNEL_SYNC");
        return (e != nullptr) && (std::atoi(e) != 0);   // 既定 false = peek 化 (同期しない)
    }();
    return enabled;
}
// 従来の `gpuErrchk( cudaDeviceSynchronize() );` を置換する。既定では no-op、FORGE_KERNEL_SYNC=1 で同期。
#define gpuErrchkKernelSync() do { if (forgeKernelSyncEnabled()) { gpuErrchk(cudaDeviceSynchronize()); } } while (0)

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
