#ifndef BOUNDARY_COND_H
#define BOUNDARY_COND_H

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>

#include "boundaryCond.hpp"

#include "flowFormat.hpp"
#include "mesh.hpp"
#include "variables.hpp"

using namespace std;

//struct neumannOrDirchletFlag
//{
//    int ro = 0;
//    int T;
//    int P;
//    int u;
//};

struct bcondConfFormat{
    int physID;
    std::string kind;
    std::map<std::string, int> inputInts ;
    std::map<std::string, flow_float> inputFloats;
    std::map<std::string, int> neuDirFlag ;
};

void readAndSetBcondConfig(vector<bcond>& bconds)
{
    std::string bcondConfigFileName  = "bcondConfig.yaml";
    map<int,bcondConfFormat> bcondConfMap;

    try {
        YAML::Node config = YAML::LoadFile(bcondConfigFileName);

        for(auto it : config) {
            bcondConfFormat bcf;

            cout << "bname =" << it.first.as<std::string>() << endl;
            std::string bname = it.first.as<std::string>();

            int physID      = it.second["physID"].as<int>();

            std::map<std::string, int> inputInts_temp;
            for (auto it_int : it.second["ints"])
            {
                std::string key = it_int.first.as<std::string>();
                int val = it_int.second.as<int>();
                inputInts_temp[key] = val;
            }

            std::map<std::string, flow_float> inputFloats_temp;
            for (auto it_float : it.second["floats"])
            {
                std::string key = it_float.first.as<std::string>();
                int val = it_float.second.as<flow_float>();
                inputFloats_temp[key] = val;
            }

            std::string kind = it.second["kind"].as<std::string>();


            bcf.physID = physID;
            bcf.kind = kind;
            bcf.inputInts = inputInts_temp;
            bcf.inputFloats = inputFloats_temp;

            bcondConfMap[bcf.physID] = bcf;

            cout << "physID=" <<  bcf.physID << endl;
            cout << "kind=" << bcf.kind<< endl;

        }

    } catch(const YAML::BadFile& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);

    } catch(const YAML::ParserException& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);
    }

    // set bconds
    for (bcond& bc : bconds)
    {
        cout << "bc.physID=" << bc.physID << endl;
        bcondConfFormat bcf = bcondConfMap[bc.physID];
        cout << "bcf.kind=" << bcf.kind << endl;
        bc.bcondKind = bcf.kind;
        bc.inputInts = bcf.inputInts;
        bc.inputFloats= bcf.inputFloats;
    }
};


void wall_isothermal(solverConfig& cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p )
{
    //IntName   : none
    //FloatName : T 
    geom_int ib_loc = 0;
    for (auto ib : bc.iBPlanes)
    {
        geom_int ip = bc.iPlanes[ib_loc];
        var.b["T"][ib] = bc.inputFloats["T"];
        var.p["T"][ip] = bc.inputFloats["T"];
        ib_loc += 1;
    }

    vector<geom_float> pcent(3);
    vector<geom_float> c1cent(3);
    geom_int ic0;
    geom_float dn;

    geom_float ss;
    vector<geom_float> sv(3);

    flow_float temp;


    // calculate diffusion term
    for (geom_int& ip : bc.iPlanes)
    {
        ic0     = msh.planes[ip].iCells[0];
        sv      = msh.planes[ip].surfVect;
        ss      = msh.planes[ip].surfArea;
        pcent   = msh.planes[ip].centCoords;
    
        c1cent  = msh.cells[ic0].centCoords;
    
        dn = ( (pcent[0] - c1cent[0])*sv[0]
              +(pcent[1] - c1cent[1])*sv[1]
              +(pcent[2] - c1cent[2])*sv[2] )/ss;
    
        var.c["diffT"][ic0] += (var.p["T"][ip] - var.c["T"][ic0])/dn*ss;

        mat_p.lhs[ic0][0] += ss/dn;
        mat_p.rhs[ic0] += var.p["T"][ip]/dn*ss; 
    }

}


void slip(solverConfig& cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p )
{
    geom_int ib_loc = 0;
    for (auto ib : bc.iBPlanes)
    //for (geom_int ib = 0 ; ib<bc.iBPlanes.size() ; ib++)
    {
        geom_int ip    = bc.iPlanes[ib_loc];
        geom_int icell = bc.iCells [ib_loc];
        var.b["T"][ib] = var.c["T"][icell];
        var.p["T"][ip] = var.c["T"][icell];
        ib_loc += 1;
    }
}

void setBcondsValue(solverConfig& cfg , mesh& msh , variables& var , matrix& mat_p)
{
    for (auto bc : msh.bconds)
    {
        if      (bc.bcondKind == "wall_isothermal") { wall_isothermal(cfg , bc , msh , var , mat_p); } 
        else if (bc.bcondKind == "slip") { slip(cfg , bc , msh , var , mat_p); }
        else {
            cerr << "Error: unknown bcondKind " << bc.bcondKind << endl;
            exit(EXIT_FAILURE);
        };
    }
}

#endif

