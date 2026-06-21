#include "residualMonitor_d.cuh"
#include "cudaWrapper.cuh"

#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/logical.h>
#include <thrust/transform_reduce.h>

#include <cmath>
#include <string>
#include <vector>

namespace {

struct SquareValue {
    __host__ __device__
    flow_float operator()(const flow_float value) const
    {
        return value * value;
    }
};

struct IsNonFinite {
    __host__ __device__
    bool operator()(const flow_float value) const
    {
        return !isfinite(value);
    }
};

}

std::array<flow_float, 5> gatherVariableRms_d(
    mesh& msh,
    variables& var,
    const std::array<const char*, 5>& variable_names)
{
    std::array<flow_float, 5> rms{};
    const flow_float cell_count = static_cast<flow_float>(std::max<geom_int>(msh.nCells, 1));

    for (std::size_t i = 0; i < variable_names.size(); ++i) {
        thrust::device_ptr<flow_float> d_ptr = thrust::device_pointer_cast(var.c_d.at(variable_names[i]));
        const flow_float sum_sq = thrust::transform_reduce(
            d_ptr,
            d_ptr + msh.nCells,
            SquareValue{},
            static_cast<flow_float>(0.0),
            thrust::plus<flow_float>()
        );
        rms[i] = std::sqrt(sum_sq / cell_count);
    }

    return rms;
}

std::vector<flow_float> gatherVariableRms_d(
    mesh& msh,
    variables& var,
    const std::vector<std::string>& variable_names)
{
    std::vector<flow_float> rms(variable_names.size(), static_cast<flow_float>(0.0));
    const flow_float cell_count = static_cast<flow_float>(std::max<geom_int>(msh.nCells, 1));

    for (std::size_t i = 0; i < variable_names.size(); ++i) {
        thrust::device_ptr<flow_float> d_ptr = thrust::device_pointer_cast(var.c_d.at(variable_names[i]));
        const flow_float sum_sq = thrust::transform_reduce(
            d_ptr,
            d_ptr + msh.nCells,
            SquareValue{},
            static_cast<flow_float>(0.0),
            thrust::plus<flow_float>()
        );
        rms[i] = std::sqrt(sum_sq / cell_count);
    }

    return rms;
}

bool hasNonFiniteCellValue_d(
    mesh& msh,
    variables& var,
    const std::vector<std::string>& variable_names,
    std::string& offending_name)
{
    offending_name.clear();
    for (const std::string& name : variable_names) {
        auto it = var.c_d.find(name);
        if (it == var.c_d.end()) {
            continue;  // 未確保の変数 (turbulence off の roK 等) はスキップ
        }
        thrust::device_ptr<flow_float> d_ptr = thrust::device_pointer_cast(it->second);
        const bool bad = thrust::any_of(d_ptr, d_ptr + msh.nCells, IsNonFinite{});
        if (bad) {
            offending_name = name;
            return true;
        }
    }
    return false;
}

std::array<flow_float, 5> gatherResidualRms_d(mesh& msh, variables& var)
{
    const std::array<const char*, 5> variable_names = {
        "res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe"
    };

    return gatherVariableRms_d(msh, var, variable_names);
}

// ---- device 常駐の残差 RMS / NaN reducer (per-step host 同期を避ける) ----
namespace {

constexpr int kRmsBlock = 256;

static inline int reduceGrid(geom_int nCells)
{
    int grid = static_cast<int>((nCells + kRmsBlock - 1) / kRmsBlock);
    if (grid < 1) grid = 1;
    if (grid > 1024) grid = 1024;   // atomic 競合を抑えるため grid-stride で上限を設ける
    return grid;
}

// 全変数の sum-of-squares を 1 カーネルで取る (block 縮約 + ブロックごと 1 atomicAdd, double 累算)。
__global__ void residualSumSq_kernel(
    int nVar, geom_int nCells,
    const flow_float* const* vars, double* sumsq)
{
    __shared__ double sdata[kRmsBlock];
    const int tid = threadIdx.x;
    const geom_int gstride = static_cast<geom_int>(gridDim.x) * blockDim.x;

    for (int v = 0; v < nVar; ++v) {
        const flow_float* arr = vars[v];
        double local = 0.0;
        for (geom_int i = static_cast<geom_int>(blockIdx.x) * blockDim.x + tid;
             i < nCells; i += gstride) {
            const double x = static_cast<double>(arr[i]);
            local += x * x;
        }
        sdata[tid] = local;
        __syncthreads();
        for (int s = blockDim.x >> 1; s > 0; s >>= 1) {
            if (tid < s) sdata[tid] += sdata[tid + s];
            __syncthreads();
        }
        if (tid == 0) atomicAdd(&sumsq[v], sdata[0]);
        __syncthreads();
    }
}

__global__ void residualFinalize_kernel(
    int nVar, geom_int nCells, const double* sumsq, flow_float* out)
{
    const int v = blockIdx.x * blockDim.x + threadIdx.x;
    if (v < nVar) {
        const double denom = static_cast<double>(nCells > 0 ? nCells : 1);
        out[v] = static_cast<flow_float>(sqrt(sumsq[v] / denom));
    }
}

__global__ void nonFinite_kernel(
    int nVar, geom_int nCells, const flow_float* const* vars, int* flag)
{
    const geom_int gstride = static_cast<geom_int>(gridDim.x) * blockDim.x;
    for (int v = 0; v < nVar; ++v) {
        const flow_float* arr = vars[v];
        for (geom_int i = static_cast<geom_int>(blockIdx.x) * blockDim.x + threadIdx.x;
             i < nCells; i += gstride) {
            if (!isfinite(arr[i])) { atomicOr(flag, 1); return; }
        }
    }
}

} // namespace

DeviceResidualReducer makeDeviceResidualReducer(
    mesh& msh, variables& var,
    const std::vector<std::string>& names,
    std::vector<std::string>& resolved_names)
{
    resolved_names.clear();
    std::vector<const flow_float*> hptrs;
    for (const std::string& n : names) {
        auto it = var.c_d.find(n);
        if (it == var.c_d.end()) {
            continue;  // 未確保の変数 (turbulence off の res_roK 等) はスキップ
        }
        hptrs.push_back(it->second);
        resolved_names.push_back(n);
    }

    DeviceResidualReducer r;
    r.nVar = static_cast<int>(hptrs.size());
    r.nCells = msh.nCells;
    if (r.nVar == 0) {
        return r;
    }

    gpuErrchk( cudaMalloc((void**)&r.d_ptrs, r.nVar * sizeof(flow_float*)) );
    gpuErrchk( cudaMemcpy((void*)r.d_ptrs, hptrs.data(),
                          r.nVar * sizeof(flow_float*), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMalloc((void**)&r.d_sumsq, r.nVar * sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&r.d_flag, sizeof(int)) );
    return r;
}

void freeDeviceResidualReducer(DeviceResidualReducer& r)
{
    if (r.d_ptrs)  { cudaFree((void*)r.d_ptrs); r.d_ptrs = nullptr; }
    if (r.d_sumsq) { cudaFree(r.d_sumsq); r.d_sumsq = nullptr; }
    if (r.d_flag)  { cudaFree(r.d_flag);  r.d_flag = nullptr; }
    r.nVar = 0;
}

flow_float* allocDeviceRmsBuffer(int capacity, int nVar)
{
    const size_t n = static_cast<size_t>(capacity) * static_cast<size_t>(nVar);
    if (n == 0) {
        return nullptr;
    }
    flow_float* p = nullptr;
    gpuErrchk( cudaMalloc((void**)&p, n * sizeof(flow_float)) );
    return p;
}

void freeDeviceRmsBuffer(flow_float* d_buf)
{
    if (d_buf) cudaFree(d_buf);
}

void reduceResidualToSlot(const DeviceResidualReducer& r, flow_float* d_buf, int slot)
{
    if (r.nVar == 0 || d_buf == nullptr) {
        return;
    }
    gpuErrchk( cudaMemsetAsync(r.d_sumsq, 0, r.nVar * sizeof(double)) );
    residualSumSq_kernel<<<reduceGrid(r.nCells), kRmsBlock>>>(
        r.nVar, r.nCells, r.d_ptrs, r.d_sumsq);
    gpuErrchk( cudaPeekAtLastError() );

    const int fblock = 64;
    const int fgrid = (r.nVar + fblock - 1) / fblock;
    residualFinalize_kernel<<<fgrid, fblock>>>(
        r.nVar, r.nCells, r.d_sumsq, d_buf + static_cast<size_t>(slot) * r.nVar);
    gpuErrchk( cudaPeekAtLastError() );
    // 意図的に同期しない (host 律速回避)。実体化は downloadRmsBuffer (flush 時) で行う。
}

void downloadRmsBuffer(const flow_float* d_buf, flow_float* host_out, int count, int nVar)
{
    if (count <= 0 || nVar <= 0 || d_buf == nullptr) {
        return;
    }
    gpuErrchk( cudaMemcpy(host_out, d_buf,
                          static_cast<size_t>(count) * nVar * sizeof(flow_float),
                          cudaMemcpyDeviceToHost) );
}

void scanNonFiniteToFlag(const DeviceResidualReducer& r)
{
    if (r.nVar == 0) {
        return;
    }
    gpuErrchk( cudaMemsetAsync(r.d_flag, 0, sizeof(int)) );
    nonFinite_kernel<<<reduceGrid(r.nCells), kRmsBlock>>>(
        r.nVar, r.nCells, r.d_ptrs, r.d_flag);
    gpuErrchk( cudaPeekAtLastError() );
}

int downloadNonFiniteFlag(const DeviceResidualReducer& r)
{
    if (r.nVar == 0) {
        return 0;
    }
    int flag = 0;
    gpuErrchk( cudaMemcpy(&flag, r.d_flag, sizeof(int), cudaMemcpyDeviceToHost) );
    return flag;
}