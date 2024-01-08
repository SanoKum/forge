#include "solvePoisson_amgx.hpp"

#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

void setMatrixPoisson(solverConfig& cfg , mesh& msh , variables& v , matrix& mat_p )
{
    geom_float ss;
    vector<geom_float> sv(3);
    //vector<geom_float> nv(3); // normal

    vector<geom_float> nv_dia(3); // diagonal
    vector<geom_float> nv_nodia(3); // non-diagonal

    vector<geom_float> c1cent(3);
    vector<geom_float> c2cent(3);
    vector<geom_float> pcent(3);
    geom_int ic0;
    geom_int ic1;
    geom_float dn;


    geom_int ip_loc0;
    geom_int ip_loc1;

    vector<geom_float> dccv(3);
    geom_float dcc;

    geom_float f;

    // normal plane
    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        ic0     = msh.planes[ip].iCells[0];
        ic1     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;
        ss      = msh.planes[ip].surfArea;

        c1cent  = msh.cells[ic0].centCoords;
        c2cent  = msh.cells[ic1].centCoords;

        dccv[0] = c2cent[0] - c1cent[0];
        dccv[1] = c2cent[1] - c1cent[1];
        dccv[2] = c2cent[2] - c1cent[2];
        dcc     = sqrt( pow(dccv[0], 2.0) + pow(dccv[1], 2.0) + pow(dccv[2], 2.0));

        ip_loc0  = mat_p.localPlnOfCell[ip][0];
        ip_loc1  = mat_p.localPlnOfCell[ip][1];

        mat_p.lhs[ic0][0] -= ss/dcc;
        mat_p.lhs[ic0][ip_loc0] = +ss/dcc;

        mat_p.lhs[ic1][0] -= ss/dcc;
        mat_p.lhs[ic1][ip_loc1] = +ss/dcc;
    }

    for (auto& bc : msh.bconds)
    {
        if (bc.valueTypes["P"] == 1)  // dirhicret
        {
            for (auto& ip : bc.iPlanes)
            {
                sv  = msh.planes[ip].surfVect;
                ss  = msh.planes[ip].surfArea;
                ic0 = msh.planes[ip].iCells[0];

                pcent   = msh.planes[ip].centCoords;
                c1cent  = msh.cells[ic0].centCoords;

                dn = ( (pcent[0] - c1cent[0])*sv[0]
                      +(pcent[1] - c1cent[1])*sv[1]
                      +(pcent[2] - c1cent[2])*sv[2] )/ss;

                mat_p.lhs[ic0][0] -= ss/dn;
            }
        } 
    }

    vector<flow_float>& divU_vol = v.c["divU*vol"];
    vector<flow_float>& ro = v.c["ro"];

    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        mat_p.rhs[ic] = divU_vol[ic]*ro[ic]/cfg.dt;
    }

    int neumann_flg = 1;
    for (auto& bc : msh.bconds)
    {
        cout << "valuetypes[p]=" << bc.valueTypes["P"] << endl;
        if (bc.valueTypes["P"] == 1) {
            neumann_flg = 0;
            cout << "neumann_flg=" << neumann_flg << endl;
        };
    }

    if (neumann_flg == 1) {
        mat_p.lhs[0][0] = 1.0e+30;
        mat_p.rhs[0] = 0.0;
    }
};


void solvePoisson(solverConfig& cfg , mesh& msh , variables& var , matrix& mat_p )
{

    std::vector<geom_int>   ptr;
    std::vector<geom_int>   col;
    std::vector<flow_float> val;
    std::vector<flow_float> rhs;

    geom_int nnz = msh.nCells + msh.nNormalPlanes*2;
    geom_int nCells = msh.nCells;

    ptr.clear(); ptr.reserve(nCells + 1); ptr.push_back(0);
    col.clear(); col.reserve(nCells * 6); 
    val.clear(); val.reserve(nCells * 6); 

    rhs.resize(nCells);

    cout << "I'm in solvePoisson" << endl;

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

    AMGX_matrix_handle matrix;
    AMGX_vector_handle rh;
    AMGX_vector_handle soln;
    AMGX_resources_handle rsrc;
    AMGX_solver_handle solver;
    AMGX_config_handle config;
    
    AMGX_SAFE_CALL(AMGX_config_create_from_file(&config, "./GMRES_AMG_D2.json"));
    AMGX_SAFE_CALL(AMGX_resources_create_simple(&rsrc, config));
    AMGX_SAFE_CALL(AMGX_matrix_create(&matrix, rsrc, AMGX_mode_dFFI));
    AMGX_SAFE_CALL(AMGX_vector_create(&rh, rsrc, AMGX_mode_dFFI));
    AMGX_SAFE_CALL(AMGX_vector_create(&soln, rsrc, AMGX_mode_dFFI));
    AMGX_SAFE_CALL(AMGX_solver_create(&solver, rsrc, AMGX_mode_dFFI, config));


    AMGX_SAFE_CALL(AMGX_matrix_upload_all(matrix, msh.nCells, nnz, 1, 1, &ptr[0], &col[0], &val[0], 0));
    AMGX_SAFE_CALL(AMGX_vector_upload(rh, msh.nCells, 1, &rhs[0]));
    AMGX_SAFE_CALL(AMGX_vector_set_zero(soln, msh.nCells, 1));

    AMGX_SAFE_CALL(AMGX_write_system(matrix, rh, soln, "identity_system.mtx"));

    AMGX_SAFE_CALL(AMGX_solver_setup(solver, matrix));
    // slight optimization to tell is that soln is all zeros
    AMGX_SAFE_CALL(AMGX_solver_solve_with_0_initial_guess(solver, rh, soln));
    AMGX_SAFE_CALL(AMGX_vector_download(soln, var.c["dP"].data()));


    //std::vector<geom_int>    ptr;
    //std::vector<geom_int>    col;
    //std::vector<flow_float> val;
    //std::vector<flow_float> rhs;

    //geom_int nCells = msh.nCells;

    //ptr.clear(); ptr.reserve(nCells + 1); ptr.push_back(0);
    //col.clear(); col.reserve(nCells * 6); 
    //val.clear(); val.reserve(nCells * 6); 

    //rhs.resize(nCells);

    //geom_int ic0 = 0;

    //for (vector<geom_int>& iStr : mat_p.structure)
    //{
        //geom_int ic_loc = 0;
        //for (geom_int& ic1 : iStr)
        //{
            //col.push_back(ic1);
            //val.push_back(mat_p.lhs[ic0][ic_loc]);
            //ic_loc += 1;
        //}

        //ptr.push_back(col.size());
        //rhs[ic0] = mat_p.rhs[ic0];

        //ic0 += 1;
    //}

    //typedef amgcl::backend::builtin<flow_float> Backend;

    //typedef amgcl::make_solver<
    //amgcl::amg<
        //Backend,
        //amgcl::runtime::coarsening::wrapper,
        //amgcl::runtime::relaxation::wrapper
        //>,
    //amgcl::runtime::solver::wrapper<Backend>
    //> Solver;

    //boost::property_tree::ptree prm;

    //prm.put("solver.type", "bicgstab");
    //prm.put("solver.tol", 1e-8);
    //prm.put("solver.maxiter", 100);
    //prm.put("precond.coarsening.type", "smoothed_aggregation");
    //prm.put("precond.relax.type", "spai0");

    //auto A = std::tie(nCells, ptr, col, val);
    //Solver solve( A , prm );

    //// Show the mini-report on the constructed solver:
    //std::cout << solve << std::endl;

    //std::vector<flow_float> x(nCells, 0.0);
    //int    iters;
    //flow_float error;
    //std::tie(iters, error) = solve(A, rhs, x);


    //std::cout << "Iters: " << iters << std::endl
              //<< "Error: " << error << std::endl;
    
    //var.c["dP"] = x;
}

void correctPresVel(solverConfig& cfg , mesh& msh , variables& v) 
{
    geom_float ss;
    vector<geom_float> sv(3);
    //vector<geom_float> nv(3); // normal

    //vector<geom_float> nv_dia(3); // diagonal
    //vector<geom_float> nv_nodia(3); // non-diagonal

    vector<geom_float> c1cent(3);
    vector<geom_float> c2cent(3);
    geom_int ic0;
    geom_int ic1;
    geom_float dn;


    geom_int ip_loc0;
    geom_int ip_loc1;

    vector<geom_float> dccv(3);
    geom_float dcc;

    vector<geom_float> dc1pv(3);
    geom_float dc1p;

    vector<geom_float> dc2pv(3);
    geom_float dc2p;

    vector<geom_float> pcent(3);

    geom_float f;

    geom_float dccv_dot_sv ;
    geom_float cosT ;

    vector<flow_float>& divU = v.c["divU"];
    vector<flow_float>& divU_vol = v.c["divU*vol"];
    vector<flow_float>& ro = v.c["ro"];
    vector<flow_float>& rop = v.p["ro"];

    vector<flow_float>& P = v.c["P"];
    vector<flow_float>& dP = v.c["dP"];

    vector<flow_float>& dPPdx = v.c["dPPdx"];
    vector<flow_float>& dPPdy = v.c["dPPdy"];
    vector<flow_float>& dPPdz = v.c["dPPdz"];

    vector<flow_float>& Ux = v.c["Ux"];
    vector<flow_float>& Uy = v.c["Uy"];
    vector<flow_float>& Uz = v.c["Uz"];

    vector<flow_float>& dUxdx = v.c["dUxdx"];
    vector<flow_float>& dUxdy = v.c["dUxdy"];
    vector<flow_float>& dUxdz = v.c["dUxdz"];

    vector<flow_float>& dUydx = v.c["dUydx"];
    vector<flow_float>& dUydy = v.c["dUydy"];
    vector<flow_float>& dUydz = v.c["dUydz"];

    vector<flow_float>& dUzdx = v.c["dUzdx"];
    vector<flow_float>& dUzdy = v.c["dUzdy"];
    vector<flow_float>& dUzdz = v.c["dUzdz"];


    vector<flow_float>& USp = v.p["US"];


    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        divU_vol[ic] = 0.0;
        dPPdx[ic] = 0.0;
        dPPdy[ic] = 0.0;
        dPPdz[ic] = 0.0;
    }

    // normal plane
    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        ic0     = msh.planes[ip].iCells[0];
        ic1     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;
        ss      = msh.planes[ip].surfArea;

        c1cent  = msh.cells[ic0].centCoords;
        c2cent  = msh.cells[ic1].centCoords;
        pcent   = msh.planes[ip].centCoords;

        dccv[0] = c2cent[0] - c1cent[0];
        dccv[1] = c2cent[1] - c1cent[1];
        dccv[2] = c2cent[2] - c1cent[2];
        dcc     = sqrt( pow(dccv[0], 2.0) + pow(dccv[1], 2.0) + pow(dccv[2], 2.0));

        dc1pv[0] = pcent[0] - c1cent[0];
        dc1pv[1] = pcent[1] - c1cent[1];
        dc1pv[2] = pcent[2] - c1cent[2];
        dc1p    = sqrt( pow(dc1pv[0], 2.0) + pow(dc1pv[1], 2.0) + pow(dc1pv[2], 2.0));

        dc2pv[0] = pcent[0] - c2cent[0];
        dc2pv[1] = pcent[1] - c2cent[1];
        dc2pv[2] = pcent[2] - c2cent[2];
        dc2p    = sqrt( pow(dc2pv[0], 2.0) + pow(dc2pv[1], 2.0) + pow(dc2pv[2], 2.0));

        dccv_dot_sv = dccv[0]*sv[0] + dccv[1]*sv[1] + dccv[2]*sv[2];
        cosT = dccv_dot_sv/dcc/ss;

        f = dc2p/dcc;

        USp[ip] += -(dP[ic1] - dP[ic0])/dcc*ss*cfg.dt/(f*ro[ic0] +(1.0-f)*ro[ic1]);

        dPPdx[ic0] += sv[0]*(f*dP[ic0] + (1.0-f)*dP[ic1]);
        dPPdy[ic0] += sv[1]*(f*dP[ic0] + (1.0-f)*dP[ic1]);
        dPPdz[ic0] += sv[2]*(f*dP[ic0] + (1.0-f)*dP[ic1]);

        dPPdx[ic1] -= sv[0]*(f*dP[ic0] + (1.0-f)*dP[ic1]);
        dPPdy[ic1] -= sv[1]*(f*dP[ic0] + (1.0-f)*dP[ic1]);
        dPPdz[ic1] -= sv[2]*(f*dP[ic0] + (1.0-f)*dP[ic1]);

        divU_vol[ic0] += USp[ip];
        divU_vol[ic1] -= USp[ip];
    }

    for (auto& bc : msh.bconds)
    {
        if (bc.valueTypes["P"] == 1)  // dirhicret
        {
            for (auto& ip : bc.iPlanes)
            {
                sv  = msh.planes[ip].surfVect;
                ss  = msh.planes[ip].surfArea;
                ic0 = msh.planes[ip].iCells[0];

                pcent   = msh.planes[ip].centCoords;
                c1cent  = msh.cells[ic0].centCoords;

                dn = ( (pcent[0] - c1cent[0])*sv[0]
                      +(pcent[1] - c1cent[1])*sv[1]
                      +(pcent[2] - c1cent[2])*sv[2] )/ss;

                dPPdx[ic0] += sv[0]*dP[ic0];
                dPPdy[ic0] += sv[1]*dP[ic0];
                dPPdz[ic0] += sv[2]*dP[ic0];

                USp[ip] += -(0.0 - dP[ic0])/dn*ss*cfg.dt/(rop[ip]);
            }
        } 

        for (auto& ip : bc.iPlanes)
        {
            ic0 = msh.planes[ip].iCells[0];
            divU_vol[ic0] += USp[ip];
        }
    }

    for (auto& bc : msh.bconds)
    {
        if (bc.valueTypes["P"] != 1)  // not dirhicret
        {
            for (auto& ip : bc.iPlanes)
            {
                sv  = msh.planes[ip].surfVect;
                ss  = msh.planes[ip].surfArea;
                ic0 = msh.planes[ip].iCells[0];

                dPPdx[ic0] += sv[0]*dP[ic0];
                dPPdy[ic0] += sv[1]*dP[ic0];
                dPPdz[ic0] += sv[2]*dP[ic0];
            }
        } 
    }

    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        dPPdx[ic] /= msh.cells[ic].volume;
        dPPdy[ic] /= msh.cells[ic].volume;
        dPPdz[ic] /= msh.cells[ic].volume;

        Ux[ic] += -dPPdx[ic]/ro[ic]*cfg.dt;
        Uy[ic] += -dPPdy[ic]/ro[ic]*cfg.dt;
        Uz[ic] += -dPPdz[ic]/ro[ic]*cfg.dt;

        P[ic] += dP[ic];

        divU[ic] = divU_vol[ic]/msh.cells[ic].volume;
    }
    flow_float max = *max_element(divU.begin(), divU.end());
    std::cout << "max div = " << max << std::endl;
};

void correctPresVel_BF(solverConfig& cfg , mesh& msh , variables& v) 
{
    geom_float ss;
    vector<geom_float> sv(3);
    //vector<geom_float> nv(3); // normal

    //vector<geom_float> nv_dia(3); // diagonal
    //vector<geom_float> nv_nodia(3); // non-diagonal

    vector<geom_float> c1cent(3);
    vector<geom_float> c2cent(3);
    geom_int ic0;
    geom_int ic1;
    geom_float dn;


    geom_int ip_loc0;
    geom_int ip_loc1;

    vector<geom_float> dccv(3);
    geom_float dcc;

    vector<geom_float> dc1pv(3);
    geom_float dc1p;

    vector<geom_float> dc2pv(3);
    geom_float dc2p;

    vector<geom_float> pcent(3);

    geom_float f;

    geom_float dccv_dot_sv ;
    geom_float cosT ;

    vector<flow_float>& divU = v.c["divU"];
    vector<flow_float>& divU_vol = v.c["divU*vol"];
    vector<flow_float>& ro = v.c["ro"];
    vector<flow_float>& rop = v.p["ro"];

    vector<flow_float>& P = v.c["P"];
    vector<flow_float>& dP = v.c["dP"];

    vector<flow_float>& dPPdx = v.c["dPPdx"];
    vector<flow_float>& dPPdy = v.c["dPPdy"];
    vector<flow_float>& dPPdz = v.c["dPPdz"];

    vector<flow_float>& Ux = v.c["Ux"];
    vector<flow_float>& Uy = v.c["Uy"];
    vector<flow_float>& Uz = v.c["Uz"];

    vector<flow_float>& dUxdx = v.c["dUxdx"];
    vector<flow_float>& dUxdy = v.c["dUxdy"];
    vector<flow_float>& dUxdz = v.c["dUxdz"];

    vector<flow_float>& dUydx = v.c["dUydx"];
    vector<flow_float>& dUydy = v.c["dUydy"];
    vector<flow_float>& dUydz = v.c["dUydz"];

    vector<flow_float>& dUzdx = v.c["dUzdx"];
    vector<flow_float>& dUzdy = v.c["dUzdy"];
    vector<flow_float>& dUzdz = v.c["dUzdz"];


    vector<flow_float>& USp = v.p["US"];


    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        divU_vol[ic] = 0.0;
        dPPdx[ic] = 0.0;
        dPPdy[ic] = 0.0;
        dPPdz[ic] = 0.0;
    }

    // normal plane
    for (geom_int ip=0 ; ip<msh.nNormalPlanes ; ip++)
    {
        ic0     = msh.planes[ip].iCells[0];
        ic1     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;
        ss      = msh.planes[ip].surfArea;

        c1cent  = msh.cells[ic0].centCoords;
        c2cent  = msh.cells[ic1].centCoords;
        pcent   = msh.planes[ip].centCoords;

        dccv[0] = c2cent[0] - c1cent[0];
        dccv[1] = c2cent[1] - c1cent[1];
        dccv[2] = c2cent[2] - c1cent[2];
        dcc     = sqrt( pow(dccv[0], 2.0) + pow(dccv[1], 2.0) + pow(dccv[2], 2.0));

        dc1pv[0] = pcent[0] - c1cent[0];
        dc1pv[1] = pcent[1] - c1cent[1];
        dc1pv[2] = pcent[2] - c1cent[2];
        dc1p    = sqrt( pow(dc1pv[0], 2.0) + pow(dc1pv[1], 2.0) + pow(dc1pv[2], 2.0));

        dc2pv[0] = pcent[0] - c2cent[0];
        dc2pv[1] = pcent[1] - c2cent[1];
        dc2pv[2] = pcent[2] - c2cent[2];
        dc2p    = sqrt( pow(dc2pv[0], 2.0) + pow(dc2pv[1], 2.0) + pow(dc2pv[2], 2.0));

        dccv_dot_sv = dccv[0]*sv[0] + dccv[1]*sv[1] + dccv[2]*sv[2];
        cosT = dccv_dot_sv/dcc/ss;

        f = dc2p/dcc;

        USp[ip] += -(dP[ic1] - dP[ic0])/dcc*ss*cfg.dt/(f*ro[ic0] +(1.0-f)*ro[ic1]);

        dPPdx[ic0] += sv[0]*(f*dP[ic0] + (1.0-f)*dP[ic1]);
        dPPdy[ic0] += sv[1]*(f*dP[ic0] + (1.0-f)*dP[ic1]);
        dPPdz[ic0] += sv[2]*(f*dP[ic0] + (1.0-f)*dP[ic1]);

        dPPdx[ic1] -= sv[0]*(f*dP[ic0] + (1.0-f)*dP[ic1]);
        dPPdy[ic1] -= sv[1]*(f*dP[ic0] + (1.0-f)*dP[ic1]);
        dPPdz[ic1] -= sv[2]*(f*dP[ic0] + (1.0-f)*dP[ic1]);

        divU_vol[ic0] += USp[ip];
        divU_vol[ic1] -= USp[ip];
    }

    for (auto& bc : msh.bconds)
    {
        if (bc.valueTypes["P"] == 1)  // dirhicret
        {
            for (auto& ip : bc.iPlanes)
            {
                sv  = msh.planes[ip].surfVect;
                ss  = msh.planes[ip].surfArea;
                ic0 = msh.planes[ip].iCells[0];

                pcent   = msh.planes[ip].centCoords;
                c1cent  = msh.cells[ic0].centCoords;

                dn = ( (pcent[0] - c1cent[0])*sv[0]
                      +(pcent[1] - c1cent[1])*sv[1]
                      +(pcent[2] - c1cent[2])*sv[2] )/ss;

                dPPdx[ic0] += sv[0]*dP[ic0];
                dPPdy[ic0] += sv[1]*dP[ic0];
                dPPdz[ic0] += sv[2]*dP[ic0];

                USp[ip] += -(0.0 - dP[ic0])/dn*ss*cfg.dt/(rop[ip]);
            }
        } 

        for (auto& ip : bc.iPlanes)
        {
            ic0 = msh.planes[ip].iCells[0];
            divU_vol[ic0] += USp[ip];
        }
    }

    for (auto& bc : msh.bconds)
    {
        if (bc.valueTypes["P"] != 1)  // not dirhicret
        {
            for (auto& ip : bc.iPlanes)
            {
                sv  = msh.planes[ip].surfVect;
                ss  = msh.planes[ip].surfArea;
                ic0 = msh.planes[ip].iCells[0];

                dPPdx[ic0] += sv[0]*dP[ic0];
                dPPdy[ic0] += sv[1]*dP[ic0];
                dPPdz[ic0] += sv[2]*dP[ic0];
            }
        } 
    }

    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        dPPdx[ic] /= msh.cells[ic].volume;
        dPPdy[ic] /= msh.cells[ic].volume;
        dPPdz[ic] /= msh.cells[ic].volume;

        Ux[ic] += -dPPdx[ic]/ro[ic]*cfg.dt;
        Uy[ic] += -dPPdy[ic]/ro[ic]*cfg.dt;
        Uz[ic] += -dPPdz[ic]/ro[ic]*cfg.dt;

        P[ic] += dP[ic];

        divU[ic] = divU_vol[ic]/msh.cells[ic].volume;
    }
    flow_float max = *max_element(divU.begin(), divU.end());
    std::cout << "max div = " << max << std::endl;
};

