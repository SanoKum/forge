#ifndef VARIABLES_H
#define VARIABLES_H

#include <iostream>
#include <vector>
#include <list>

#include <flowFormat.hpp>
#include <mesh.hpp>

//using namespace std

class variables {
public:
    map<string, vector<flow_float>> c;
    map<string, vector<flow_float>> p;
    map<string, vector<flow_float>> b;

    list<string> cellValNames = 
    {
        "ro"  , "Ux"  , "Uy"  , "Uz"  , "T" , "P" , 
        "roN" , "UxN" , "UyN" , "UzN" , "TN" ,
        "dTdx", "dTdy", "dTdz",
        "divU" ,
        "convx" , "convy" , "convz" , "convT",
        "diffx" , "diffy" , "diffz" , "diffT",
        "src"  , "rhs", 
    };

    list<string> planeValNames = 
    {
        "T" , "fx"
    };

    list<string> bPlaneValNames = 
    {
        "T"
    };


    //not yet implemented
    list<string> output_cellValNames = 
    {
        "ro"  , "Ux"  , "Uy"  , "Uz"  , "T" , "P" , 
        "divU" ,
    };

    //not yet implemented
    list<string> read_cellValNames = 
    {
        "ro"  , "Ux"  , "Uy"  , "Uz"  , "T" , "P" ,
    };

    variables(mesh msh)
    {
        for (auto cellValName : cellValNames)
        {
            this->c[cellValName].resize(msh.nCells);
        }
        for (auto planeValName : planeValNames)
        {
            this->p[planeValName].resize(msh.nPlanes);
        }
        for (auto bPlaneValName : bPlaneValNames)
        {
            this->b[bPlaneValName].resize(msh.nBPlanes);
        }

    }

};

#endif