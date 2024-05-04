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
        "UxN"  , "UyN"   , "UzN"   , "TN"     ,
        "roN"  , "roUxN" , "roUyN" , "roUzN"  , "roeN" ,
        "roNN" , "roUxNN", "roUyNN", "roUzNN" , "roeNN" ,
        "ro"   , "roUx"  , "roUy"  , "roUz"   , "roe"  ,
        "dTdx" , "dTdy"  , "dTdz",
        "dUxdx", "dUxdy" , "dUxdz",
        "dUydx", "dUydy" , "dUydz",
        "dUzdx", "dUzdy" , "dUzdz",
        "dPdx" , "dPdy"  , "dPdz",
        "drodx", "drody" , "drodz",
        "dHtdx", "dHtdy" , "dHtdz",
        "divU*vol" , "divU" , "divU_star",
        "cfl"   , "cfl_pseudo",
        "res_ro"    , "res_roUx"   , "res_roUy"   , "res_roUz"   , "res_roe",
        "res_ro_m"  , "res_roUx_m" , "res_roUy_m" , "res_roUz_m" , "res_roe_m",
        "sonic"   , "Ht" , "ros", // s: entropy

        "roM"  , "roUxM", "roUyM", "roUzM", "roeM" , 

        // Mesh Structure
        "volume" , "ccx" , "ccy" , "ccz"
    };

    const std::list<std::string> planeValNames = 
    {
        "cfl_pln" , "cfl_pseudo_pln",
        //"US" , "T" ,  
        //"USN", "Ux" , "Uy" , "Uz", "ro" , "P",

        // Mesh Structure
        "sx"  , "sy" , "sz" , "ss" ,
        "pcx" , "pcy" , "pcz" ,
        "fx"  , 
        "dcc"   // dcc: distance between two cell centers
    };

    const std::list<std::string> output_cellValNames = 
    {
        "ro"    , "Ux"    , "Uy"    , "Uz"  , "T" , "P" , 
        "roUx"  , "roUy"  , "roUz"  , "roe" , 
        "cfl"   , "volume", "sonic" , 
        "dUxdx" , "dUxdy" , "dUxdz" , 
        "dUydx" , "dUydy" , "dUydz" , 
        "dUzdx" , "dUzdy" , "dUzdz" ,
        "drodx" , "drody" , "drodz" ,
        "dPdx"  , "dPdy"  , "dPdz" ,
        "dHtdx" , "dHtdy" , "dHtdz" 
    };

    //not yet implemented
    const std::list<std::string> read_cellValNames = 
    {
        "ro"  , "roUx"  , "roUy"  , "roUz"  , "roe" 
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
    void setStructuralVariables_d(cudaConfig& cuda_cfg , mesh& msh);

    void readValueHDF5(std::string fname , mesh& msh);
};
