#pragma once

#include <iostream>
#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "cuda_forge/cudaConfig.cuh"

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>



class variables {
public:
    std::map<std::string, std::vector<flow_float>> c; // host cell variables
    std::map<std::string, std::vector<flow_float>> p; // host plane variables
    std::map<std::string, flow_float*> c_d; // device cell variables
    std::map<std::string, flow_float*> p_d; // device plane variables


    const std::list<std::string> cellValNames = 
    {
        "Ux"   , "Uy"    , "Uz"    , "T"      , "P" , "deltaP" ,
        "k"    , "omega" ,
        "UxN"  , "UyN"   , "UzN"   , "TN"     ,
        "roN"  , "roUxN" , "roUyN" , "roUzN"  , "roeN" , "roKN" , "roOmegaN" ,
        "roNN" , "roUxNN", "roUyNN", "roUzNN" , "roeNN" , "roKNN" , "roOmegaNN" ,
        "ro"   , "roUx"  , "roUy"  , "roUz"   , "roe"  , "roK" , "roOmega" ,
        "dTdx" , "dTdy"  , "dTdz",
        "dUxdx", "dUxdy" , "dUxdz",
        "dUydx", "dUydy" , "dUydz",
        "dUzdx", "dUzdy" , "dUzdz",
        "dPdx" , "dPdy"  , "dPdz",
        "drodx", "drody" , "drodz",
        "dKdx" , "dKdy"  , "dKdz",
        "dOmegadx", "dOmegady", "dOmegadz",

        "divU*vol" , "divU" , "divU_star",

        "ducros" , 

        "limiter_Ux" , 
        "limiter_Uy" , 
        "limiter_Uz" , 
        "limiter_ro" , 
        "limiter_P" , 
        "limiter_T" , 
        "vis_lam" , 
        "vis_turb" , "wall_dist" , 

        "gamma" , "cp" , 

        "cfl"   , "cfl_pseudo", "dt_local",
        "res_ro"    , "res_roUx"   , "res_roUy"   , "res_roUz"   , "res_roe", "res_roK", "res_roOmega",
        "res_ro_m"  , "res_roUx_m" , "res_roUy_m" , "res_roUz_m" , "res_roe_m", "res_roK_m", "res_roOmega_m",
        "dq_ro_old" , "dq_roUx_old", "dq_roUy_old", "dq_roUz_old", "dq_roe_old",
        "dq_ro_new" , "dq_roUx_new", "dq_roUy_new", "dq_roUz_new", "dq_roe_new",
        "dq_block_old_0", "dq_block_old_1", "dq_block_old_2", "dq_block_old_3", "dq_block_old_4",
        "dq_block_new_0", "dq_block_new_1", "dq_block_new_2", "dq_block_new_3", "dq_block_new_4",
        "rhs_block_0", "rhs_block_1", "rhs_block_2", "rhs_block_3", "rhs_block_4",
        "diag_block_00", "diag_block_01", "diag_block_02", "diag_block_03", "diag_block_04",
        "diag_block_10", "diag_block_11", "diag_block_12", "diag_block_13", "diag_block_14",
        "diag_block_20", "diag_block_21", "diag_block_22", "diag_block_23", "diag_block_24",
        "diag_block_30", "diag_block_31", "diag_block_32", "diag_block_33", "diag_block_34",
        "diag_block_40", "diag_block_41", "diag_block_42", "diag_block_43", "diag_block_44",
        "sonic"   , "Ht" , "ros", // s: entropy

        "roM"  , "roUxM", "roUyM", "roUzM", "roeM" , "roKM", "roOmegaM",

        // Mesh Structure
        "volume" , "ccx" , "ccy" , "ccz",
        "A_planar",  // 軸対称: planar 面積 (= V を r で割った値)。非軸対称時は 0
        "axisym_uy_over_r",
        "axisym_divU",
        // SST 陰解法 (segregated point-implicit) 用: 消散項のヤコビアン対角
        //   src_jac_k     = ∂Dk/∂(ρk)   = β* ω
        //   src_jac_omega = ∂Dω/∂(ρω)   = 2 β ω
        "src_jac_k",
        "src_jac_omega"
    };

    const std::list<std::string> planeValNames = 
    {
        "cfl_pln" , "cfl_pseudo_pln",
        //"US" , "T" ,  
        //"USN", "Ux" , "Uy" , "Uz", "ro" , "P",

        // Mesh Structure
        "sx"  , "sy" , "sz" , "ss" ,
        "sx_planar" , "sy_planar" , "sz_planar" , "ss_planar" ,
        "pcx" , "pcy" , "pcz" ,
        "fx"  , 
        "massflux" ,
        "dcc"   // dcc: distance between two cell centers
    };

    const std::list<std::string> output_cellValNames = 
    {
        "ro"    , "Ux"    , "Uy"    , "Uz"  , "T" , "P" , 
        "roUx"  , "roUy"  , "roUz"  , "roe" , "roK" , "roOmega" , "k" , "omega" ,
        "cfl"   , "volume", "sonic" , 
        "dUxdx" , "dUxdy" , "dUxdz" , 
        "dUydx" , "dUydy" , "dUydz" , 
        "dUzdx" , "dUzdy" , "dUzdz" ,
        "drodx" , "drody" , "drodz" ,
        "dPdx"  , "dPdy"  , "dPdz" ,
        "dKdx"  , "dKdy"  , "dKdz" ,
        "dOmegadx", "dOmegady", "dOmegadz",

        "ducros",

        "limiter_ro" , 
        "limiter_Ux" , 
        "limiter_Uy" , 
        "limiter_Uz" , 
        "limiter_P" , 
        //"limiter_Ht" , 

        "dt_local",

        "axisym_divU",

        "vis_lam" , "vis_turb" , "wall_dist"
    };

    //not yet implemented
    const std::list<std::string> read_cellValNames = 
    {
        "ro"  , "roUx"  , "roUy"  , "roUz"  , "roe" , "wall_dist"
    };

    variables();

    //variables(const int& ,  mesh&);
    ~variables();

    void allocVariables(const int &useGPU , mesh& msh);


    void copyVariables_cell_plane_H2D_all();
    void copyVariables_cell_H2D (std::list<std::string>);
    void copyVariables_plane_H2D(std::list<std::string>);

    void copyVariables_cell_plane_D2H_all();
    void copyVariables_cell_D2H (std::list<std::string>);
    void copyVariables_plane_D2H(std::list<std::string>);

    void setStructuralVariables(solverConfig& cfg, cudaConfig& cuda_cfg , mesh& msh);
    void setStructuralVariables_d(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh);

    void readValueHDF5(std::string fname , mesh& msh);
};