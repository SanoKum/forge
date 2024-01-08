#pragma once
#include <iostream>
#include <vector>
#include <list>

#include <flowFormat.hpp>
#include <mesh/mesh.hpp>

class variables {
public:
    std::map<std::string, std::vector<flow_float>> c;
    std::map<std::string, std::vector<flow_float>> p;
    //map<string, vector<flow_float>> b;

    const std::list<std::string> cellValNames = 
    {
        "ro"  , "Ux"  , "Uy"  , "Uz"  , "T" , "P" , "deltaP" ,
        "roN" , "UxN" , "UyN" , "UzN" , "TN" ,
        "dTdx", "dTdy", "dTdz",
        "dUxdx", "dUxdy", "dUxdz",
        "dUydx", "dUydy", "dUydz",
        "dUzdx", "dUzdy", "dUzdz",
        "dPdx", "dPdy", "dPdz",
        "divU*vol" , "divU" , "divU_star",
        "convx" , "convy" , "convz" , "convT",
        "diffx" , "diffy" , "diffz" , "diffT",
//        "src" // , "res", 
        "dP", 
        "dPPdx", "dPPdy", "dPPdz",
//        "volume" , "svsum_x", "svsum_y", "svsum_z",
        "cfl"
    };

    const std::list<std::string> planeValNames = 
    {
        "US" , "T" , "fx" , 
        "USN", "Ux" , "Uy" , "Uz", "ro" , "P"
    };

    const std::list<std::string> output_cellValNames = 
    {
        "ro"  , "Ux"  , "Uy"  , "Uz"  , "T" , "P" , 
        "divU" , "divU_star",
        "convx" , "convy" , "convz" , 
        "diffx" , "diffy" , "diffz" , 
        //"dPPdx" , "dPPdy" , "dPPdz",
        "dPPdx" , "dPPdy"   , "dPPdz",
        "dP"  ,
        //"volume" ,
        //"svsum_x", "svsum_y", "svsum_z",
        "cfl"
    };

    //not yet implemented
    const std::list<std::string> read_cellValNames = 
    {
        "ro"  , "Ux"  , "Uy"  , "Uz"  , "T" , "P" ,
        "dP"
    };

    variables(mesh&);
};
