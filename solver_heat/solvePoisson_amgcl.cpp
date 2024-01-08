#include <vector>                                                                                        
#include <iostream>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>

#include <amgcl/io/mm.hpp>
#include <amgcl/profiler.hpp>

#include "mesh.hpp"
#include "variables.hpp"
#include "readConfig.hpp"

int solvePoisson(const solverConfig& cfg ,
    mesh& msh ,
    variables& var ,
    matrix& mat_p 
    )
{
    std::vector<geom_int>    ptr;
    std::vector<geom_int>    col;
    std::vector<flow_float> val;
    std::vector<flow_float> rhs;


    geom_int nCells = msh.nCells;

    ptr.clear(); ptr.reserve(nCells + 1); ptr.push_back(0);
    col.clear(); col.reserve(nCells * 6); 
    val.clear(); val.reserve(nCells * 6); 

    rhs.resize(nCells);

    geom_int ic0 = 0;

    for (vector<geom_int>& iStr : mat_p.structure)
    {
        col.push_back(ic0);
        val.push_back(mat_p.lhs[ic0][0]);

        for (geom_int& ic1 : iStr)
        {
            col.push_back(ic1);
            val.push_back(mat_p.lhs[ic0][ic1]);
        }

        ptr.push_back(col.size());
        rhs.push_back(mat_p.rhs[ic0]);

        ic0 += 1;
    }

    typedef amgcl::backend::builtin<flow_float> Backend;

    typedef amgcl::make_solver<
        // Use AMG as preconditioner:
        amgcl::amg<
            Backend,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0
            >,
        // And BiCGStab as iterative solver:
        amgcl::solver::bicgstab<Backend>
        > Solver;

    Solver solve( std::tie(nCells, ptr, col, val) );

    std::vector<flow_float> x(nCells, 0.0);
    int    iters;
    flow_float error;
    std::tie(iters, error) = solve(rhs, x);

    return 0;

//    amgcl::profiler<> prof("poisson3Db");
//
//    // We use the tuple of CRS arrays to represent the system matrix.
//    // Note that std::tie creates a tuple of references, so no data is actually
//    // copied here:
//    auto A = std::tie(rows, ptr, col, val);
//
//    // Compose the solver type
//    //   the solver backend:
//    typedef amgcl::backend::builtin<flow_float> SBackend;
//    //   the preconditioner backend:
//#ifdef MIXED_PRECISION
//    typedef amgcl::backend::builtin<float> PBackend;
//#else
//    typedef amgcl::backend::builtin<flow_float> PBackend;
//#endif
//    
//    typedef amgcl::make_solver<
//        amgcl::amg<
//            PBackend,
//            amgcl::coarsening::smoothed_aggregation,
//            amgcl::relaxation::spai0
//            >,
//        amgcl::solver::bicgstab<SBackend>
//        > Solver;
//
//    // Initialize the solver with the system matrix:
//    prof.tic("setup");
//    Solver solve(A);
//    prof.toc("setup");
//
//    // Show the mini-report on the constructed solver:
//    std::cout << solve << std::endl;
//
//    // Solve the system with the zero initial approximation:
//    int iters;
//    double error;
//    std::vector<double> x(rows, 0.0);
//
//    prof.tic("solve");
//    std::tie(iters, error) = solve(A, rhs, x);
//    prof.toc("solve");
//
//    // Output the number of iterations, the relative error,
//    // and the profiling data:
//    std::cout << "Iters: " << iters << std::endl
//              << "Error: " << error << std::endl
//              << prof << std::endl;



}