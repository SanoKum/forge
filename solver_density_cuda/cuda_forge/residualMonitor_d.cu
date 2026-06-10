#include "residualMonitor_d.cuh"

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