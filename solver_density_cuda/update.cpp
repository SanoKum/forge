#include "update.hpp"
#include <vector>

using namespace std;

void updateVariablesOuter(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& v , matrix& mat_p)
//void updateVariablesForNextLoop(solverConfig& cfg , mesh& msh , variables& v)
{
    if (cfg.gpu == 1) {
        updateVariablesOuter_d_wrapper(cfg , cuda_cfg , msh , v);
        return;
    }

    vector<flow_float>& ro   = v.c["ro"];
    vector<flow_float>& roUx = v.c["roUx"];
    vector<flow_float>& roUy = v.c["roUy"];
    vector<flow_float>& roUz = v.c["roUz"];
    vector<flow_float>& roe  = v.c["roe"];

    vector<flow_float>& roN  = v.c["roN"];
    vector<flow_float>& roUxN= v.c["roUxN"];
    vector<flow_float>& roUyN= v.c["roUyN"];
    vector<flow_float>& roUzN= v.c["roUzN"];
    vector<flow_float>& roeN = v.c["roeN"];
    vector<flow_float>& roKN = v.c["roKN"];
    vector<flow_float>& roOmegaN = v.c["roOmegaN"];

    vector<flow_float>& roM  = v.c["roM"];
    vector<flow_float>& roUxM= v.c["roUxM"];
    vector<flow_float>& roUyM= v.c["roUyM"];
    vector<flow_float>& roUzM= v.c["roUzM"];
    vector<flow_float>& roeM = v.c["roeM"];
    vector<flow_float>& roKM = v.c["roKM"];
    vector<flow_float>& roOmegaM = v.c["roOmegaM"];

    vector<flow_float>& roK = v.c["roK"];
    vector<flow_float>& roOmega = v.c["roOmega"];

    //vector<flow_float>& convx = v.c["convx"];
    //vector<flow_float>& convy = v.c["convy"];
    //vector<flow_float>& convz = v.c["convz"];

    //vector<flow_float>& diffx = v.c["diffx"];
    //vector<flow_float>& diffy = v.c["diffy"];
    //vector<flow_float>& diffz = v.c["diffz"];

    vector<flow_float>& divU = v.c["divU"];

    for (geom_int ic=0 ; ic<msh.nCells_all; ic++)
    {

        roN[ic] = ro[ic];
        roUxN[ic] = roUx[ic];
        roUyN[ic] = roUy[ic];
        roUzN[ic] = roUz[ic];
        roeN[ic] = roe[ic];
        roKN[ic] = roK[ic];
        roOmegaN[ic] = roOmega[ic];

        roM[ic] = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic] = roe[ic];
        roKM[ic] = roK[ic];
        roOmegaM[ic] = roOmega[ic];

        divU[ic] = 0.0;
    }

}

void updateVariablesInner(solverConfig& cfg , cudaConfig& cuda_cfg ,mesh& msh , variables& v , matrix& mat_p)
//void updateVariablesForNextLoop(solverConfig& cfg , mesh& msh , variables& v)
{
    if (cfg.gpu==1) {
        updateVariablesInner_d_wrapper(cfg , cuda_cfg , msh , v);
        return;
    }


    vector<flow_float>& ro   = v.c["ro"];
    vector<flow_float>& roUx = v.c["roUx"];
    vector<flow_float>& roUy = v.c["roUy"];
    vector<flow_float>& roUz = v.c["roUz"];
    vector<flow_float>& roe  = v.c["roe"];
    vector<flow_float>& roK  = v.c["roK"];
    vector<flow_float>& roOmega  = v.c["roOmega"];

    vector<flow_float>& roM  = v.c["roM"];
    vector<flow_float>& roUxM= v.c["roUxM"];
    vector<flow_float>& roUyM= v.c["roUyM"];
    vector<flow_float>& roUzM= v.c["roUzM"];
    vector<flow_float>& roeM = v.c["roeM"];
    vector<flow_float>& roKM = v.c["roKM"];
    vector<flow_float>& roOmegaM = v.c["roOmegaM"];

    vector<flow_float>& divU = v.c["divU"];

    for (geom_int ic=0 ; ic<msh.nCells_all; ic++)
    {

        roM[ic] = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic] = roe[ic];
        roKM[ic] = roK[ic];
        roOmegaM[ic] = roOmega[ic];

        divU[ic] = 0.0;
    }

}

void applyScalarImplicitCorrection(solverConfig& cfg , cudaConfig& cuda_cfg ,mesh& msh , variables& v , matrix& mat_p)
{
    if (cfg.gpu==1) {
        applyScalarImplicitCorrection_d_wrapper(cfg , cuda_cfg , msh , v);
        return;
    }

    vector<flow_float>& ro   = v.c["ro"];
    vector<flow_float>& roUx = v.c["roUx"];
    vector<flow_float>& roUy = v.c["roUy"];
    vector<flow_float>& roUz = v.c["roUz"];
    vector<flow_float>& roe  = v.c["roe"];

    vector<flow_float>& roN   = v.c["roN"];
    vector<flow_float>& roUxN = v.c["roUxN"];
    vector<flow_float>& roUyN = v.c["roUyN"];
    vector<flow_float>& roUzN = v.c["roUzN"];
    vector<flow_float>& roeN  = v.c["roeN"];

    vector<flow_float>& dq_ro   = v.c["dq_ro_old"];
    vector<flow_float>& dq_roUx = v.c["dq_roUx_old"];
    vector<flow_float>& dq_roUy = v.c["dq_roUy_old"];
    vector<flow_float>& dq_roUz = v.c["dq_roUz_old"];
    vector<flow_float>& dq_roe  = v.c["dq_roe_old"];

    for (geom_int ic=0 ; ic<msh.nCells; ic++)
    {
        ro[ic] = roN[ic] + dq_ro[ic];
        roUx[ic] = roUxN[ic] + dq_roUx[ic];
        roUy[ic] = roUyN[ic] + dq_roUy[ic];
        roUz[ic] = roUzN[ic] + dq_roUz[ic];
        roe[ic] = roeN[ic] + dq_roe[ic];
    }
}
