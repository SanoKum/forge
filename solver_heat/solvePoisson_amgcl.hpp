#ifndef SOLVE_POISSON_AMGCL
#define SOLVE_POISSON_AMGCL

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

void solvePoisson(const solverConfig& cfg ,
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
        geom_int ic_loc = 0;
        for (geom_int& ic1 : iStr)
        {
            col.push_back(ic1);
            val.push_back(mat_p.lhs[ic0][ic_loc]);
            ic_loc += 1;
        }

        ptr.push_back(col.size());
        rhs[ic0] = mat_p.rhs[ic0];

        ic0 += 1;
    }

    typedef amgcl::backend::builtin<flow_float> Backend;

    typedef amgcl::make_solver<
    amgcl::amg<
        Backend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
        >,
    amgcl::solver::bicgstab<Backend>
    > Solver;


    auto A = std::tie(nCells, ptr, col, val);
    Solver solve( A );

    // Show the mini-report on the constructed solver:
    std::cout << solve << std::endl;

    std::vector<flow_float> x(nCells, 0.0);
    int    iters;
    flow_float error;
    std::tie(iters, error) = solve(A, rhs, x);


    std::cout << "Iters: " << iters << std::endl
              << "Error: " << error << std::endl;

    var.c["P"] = x;
}

#endif