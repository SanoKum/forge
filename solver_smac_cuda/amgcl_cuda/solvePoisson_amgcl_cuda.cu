#include <vector>
#include <iostream>

#include <amgcl/backend/cuda.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>

#include <amgcl/io/mm.hpp>
#include <amgcl/profiler.hpp>

#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>

#include "amgcl_cuda/solvePoisson_amgcl_cuda.h"

using namespace std;


void solvePoisson_amgcl_cuda( std::vector<ptrdiff_t>& ptr , std::vector<ptrdiff_t>& col,
                              std::vector<flow_float>& val, std::vector<flow_float>& rhs,
                              std::vector<flow_float>& result )
{
    int device;
    cudaDeviceProp prop;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);
    std::cout << prop.name << std::endl;

    // The profiler:
    amgcl::profiler<> prof("poisson3Db");

    //ptrdiff_t rows, cols;
    ptrdiff_t rows = rhs.size();


    auto A = std::tie(rows, ptr, col, val);

    typedef amgcl::backend::cuda<flow_float> Backend;

    typedef amgcl::make_solver<
    amgcl::amg<
        Backend,
        amgcl::runtime::coarsening::wrapper,
        amgcl::runtime::relaxation::wrapper
        >,
    amgcl::runtime::solver::wrapper<Backend>
    > Solver;

    boost::property_tree::ptree prm;

    prm.put("solver.type", "bicgstab");
    prm.put("solver.tol", 1e-8);
    prm.put("solver.maxiter", 100);
    prm.put("precond.coarsening.type", "smoothed_aggregation");
    prm.put("precond.relax.type", "spai0");

    Backend::params bprm;
    cusparseCreate(&bprm.cusparse_handle);

    prof.tic("setup");
    Solver solve(A, prm, bprm);
    prof.toc("setup");

    // Show the mini-report on the constructed solver:
    std::cout << solve << std::endl;

    int iters;
    double error;

    thrust::device_vector<flow_float> f(rhs);
    thrust::device_vector<flow_float> x(rows, 0.0);

    prof.tic("solve");
    std::tie(iters, error) = solve(f, x);
    prof.toc("solve");

    // Output the number of iterations, the relative error,
    // and the profiling data:
    std::cout << "Iters: " << iters << std::endl
              << "Error: " << error << std::endl
              << prof << std::endl;
    
    for (geom_int i=0 ; i<rows ; i++)
    {
        result[i] = x[i];
    }
}

