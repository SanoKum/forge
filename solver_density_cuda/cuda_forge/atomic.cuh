#pragma once

#if defined(__CUDA_ARCH__) 


#if __COMPUTE_CAPABILITY__ < 600

static inline __device__ double atomicAdd(double* address, double val)
{
  unsigned long long int* address_as_ull =
                            (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val +
                           __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);

  return __longlong_as_double(old);
}
#endif

#endif

static inline __device__ double2 atomicAdd(double2 *addr, double2 val){
  double2 old = *addr;
  // This is a necessary evil to avoid conflicts between the atomicAdd
  // declared in CUDA 8.0+ headers which are visible for host
  // compilation, which cause a conflict when compiled on clang-cuda.
  // As a result we do not support any architecture without native
  // double precision atomics on clang-cuda.
#if defined(__CUDA_ARCH__) || CUDA_VERSION >= 8000
  old.x = atomicAdd((double*)addr, val.x);
  old.y = atomicAdd((double*)addr + 1, val.y);
#endif
  return old;
}

static inline __device__ float2 atomicAdd(float2 *addr, float2 val){
  float2 old = *addr;
  old.x = atomicAdd((float*)addr, val.x);
  old.y = atomicAdd((float*)addr + 1, val.y);
  return old;
}

static inline __device__ int2 atomicAdd(int2 *addr, int2 val){
  int2 old = *addr;
  old.x = atomicAdd((int*)addr, val.x);
  old.y = atomicAdd((int*)addr + 1, val.y);
  return old;
}

union uint32_short2 { unsigned int i; short2 s; };

static inline __device__ short2 atomicAdd(short2 *addr, short2 val){
  uint32_short2 old, assumed, incremented;
  old.s = *addr;
  do {
    assumed.s = old.s;
    incremented.s = make_short2(val.x + assumed.s.x, val.y + assumed.s.y);
    old.i =  atomicCAS((unsigned int*)addr, assumed.i, incremented.i);
  } while ( assumed.i != old.i );

  return old.s;
}

union uint32_char2 { unsigned short i; char2 s; };

static inline __device__ char2 atomicAdd(char2 *addr, char2 val){
  uint32_char2 old, assumed, incremented;
  old.s = *addr;
  do {
    assumed.s = old.s;
    incremented.s = make_char2(val.x + assumed.s.x, val.y + assumed.s.y);
    old.i =  atomicCAS((unsigned int*)addr, assumed.i, incremented.i);
  } while ( assumed.i != old.i );

  return old.s;
}

static inline __device__ float atomicMax(float *addr, float val){
  unsigned int old = __float_as_uint(*addr), assumed;
  do {
    assumed = old;
    if (__uint_as_float(old) >= val) break;

    old = atomicCAS((unsigned int*)addr,
           assumed,
           __float_as_uint(val));
  } while ( assumed != old );

  return __uint_as_float(old);
}