#include "boundaryCond.hpp"

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>

#include "yaml-cpp/yaml.h"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

using namespace std;

bcondConfFormat::bcondConfFormat(){};

void readAndSetBcondConfig(vector<bcond>& bconds)
{
    map<int,bcondConfFormat> bcondConfMap;

    std::string bcondConfigFileName  = "bcondConfig.yaml";

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

        map<string,int> valueTypes = bcf.valueTypesOfBC[bcf.kind];
        for (auto& vt : valueTypes)
        {
            bc.valueTypes[vt.first] = vt.second;

            if (vt.second != 1  ) continue; // not dirichlet

            bc.bvar[vt.first].resize(bc.iPlanes.size());
            for (flow_float& var : bc.bvar[vt.first])
            {
                var = bc.inputFloats[vt.first];
            }
        }
    }
};

void wall_isothermal(solverConfig& cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p )
{
//    //IntName   : none
//    //FloatName : T 
//    //geom_int ib_loc = 0;
//    //for (auto ib : bc.iBPlanes)
//    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
//    {
//        geom_int ip = bc.iPlanes[ib];
//        var.p["T"][ip] = bc.bvar["T"][ib];
//    }
//
//    vector<geom_float> pcent(3);
//    vector<geom_float> c1cent(3);
//    geom_int ic0;
//    geom_float dn;
//
//    geom_float ss;
//    vector<geom_float> sv(3);
//
//    flow_float temp;
//
//
//    // calculate diffusion term
//    for (geom_int& ip : bc.iPlanes)
//    {
//        ic0     = msh.planes[ip].iCells[0];
//        sv      = msh.planes[ip].surfVect;
//        ss      = msh.planes[ip].surfArea;
//        pcent   = msh.planes[ip].centCoords;
//    
//        c1cent  = msh.cells[ic0].centCoords;
//    
//        dn = ( (pcent[0] - c1cent[0])*sv[0]
//              +(pcent[1] - c1cent[1])*sv[1]
//              +(pcent[2] - c1cent[2])*sv[2] )/ss;
//    
//        var.c["diffT"][ic0] += (var.p["T"][ip] - var.c["T"][ic0])/dn*ss;
//
//        mat_p.lhs[ic0][0] += ss/dn;
//        mat_p.rhs[ic0] += var.p["T"][ip]/dn*ss; 
//    }
//
};

void wall(solverConfig& cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p )
{
    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip = bc.iPlanes[ib];
        var.p["Ux"][ip] = bc.bvar["Ux"][ib];
        var.p["Uy"][ip] = bc.bvar["Uy"][ib];
        var.p["Uz"][ip] = bc.bvar["Uz"][ib];

        geom_int ic0 = bc.iCells [ib];
        var.p["ro"][ip] = var.c["ro"][ic0];
        var.p["P"][ip]  = var.c["P"][ic0];
        var.p["T"][ip]  = var.c["T"][ic0];
    }
};

void inlet_uniformVelocity(solverConfig& cfg , bcond& bc , mesh& msh , variables& v , matrix& mat_p )
{
    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip = bc.iPlanes[ib];

        v.p["Ux"][ip] = bc.bvar["Ux"][ib];
        v.p["Uy"][ip] = bc.bvar["Uy"][ib];
        v.p["Uz"][ip] = bc.bvar["Uz"][ib];

        geom_int ic0 = bc.iCells [ib];
        v.p["ro"][ip] = v.c["ro"][ic0];
        v.p["P"] [ip] = v.c["P"][ic0];
        v.p["T"] [ip] = v.c["T"][ic0];
    }
}

void outlet_statPress(solverConfig& cfg , bcond& bc , mesh& msh , variables& v , matrix& mat_p )
{
    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip = bc.iPlanes[ib];
        v.p["P"] [ip] = bc.bvar["P"][ib];

        geom_int ic0 = bc.iCells [ib];
        v.p["ro"][ip] = v.c["ro"][ic0];
        v.p["T"] [ip] = v.c["T"][ic0];
        v.p["Ux"][ip] = v.c["Ux"][ic0];
        v.p["Uy"][ip] = v.c["Uy"][ic0];
        v.p["Uz"][ip] = v.c["Uz"][ic0];

    }
}

void slip(solverConfig& cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p )
{
    vector<flow_float>& Tp  = var.p["T"];
    vector<flow_float>& Uxp = var.p["Ux"];
    vector<flow_float>& Uyp = var.p["Uy"];
    vector<flow_float>& Uzp = var.p["Uz"];
    vector<flow_float>& rop = var.p["ro"];
    vector<flow_float>& Pp  = var.p["P"];

    vector<flow_float>& T   = var.c["T"];
    vector<flow_float>& Ux  = var.c["Ux"];
    vector<flow_float>& Uy  = var.c["Uy"];
    vector<flow_float>& Uz  = var.c["Uz"];
    vector<flow_float>& ro  = var.c["ro"];
    vector<flow_float>& P   = var.c["P"];


    for (geom_int ib=0; ib<bc.iPlanes.size(); ib++)
    {
        geom_int ip    = bc.iPlanes[ib];
        geom_int icell = bc.iCells [ib];

        flow_float Ux2 ;
        flow_float Uy2;
        flow_float Uz2;
        flow_float Un;
        vector<flow_float> sv = msh.planes[ip].surfVect;
        flow_float ss = msh.planes[ip].surfArea;

        Un =  (sv[0]*Uxp[icell] + sv[1]*Uyp[icell] + sv[2]*Uzp[icell])/ss;

        Tp[ip]  = T[icell];
        Uxp[ip] = Ux[icell] - Un*sv[0]/ss;
        Uyp[ip] = Uy[icell] - Un*sv[1]/ss;
        Uzp[ip] = Uz[icell] - Un*sv[2]/ss;
        rop[ip] = ro[icell];
        Pp[ip]  = P[icell];

    }
}

void setBcondsValue(solverConfig& cfg , mesh& msh , variables& var , matrix& mat_p)
{
    for (auto bc : msh.bconds)
    {
        if      (bc.bcondKind == "wall_isothermal") { wall_isothermal(cfg , bc , msh , var , mat_p); } 
        else if (bc.bcondKind == "slip") { slip(cfg , bc , msh , var , mat_p); }
        else if (bc.bcondKind == "wall") { wall(cfg , bc , msh , var , mat_p); }
        else if (bc.bcondKind == "inlet_uniformVelocity") { inlet_uniformVelocity(cfg , bc , msh , var , mat_p); }
        else if (bc.bcondKind == "outlet_statPress") { outlet_statPress(cfg , bc , msh , var , mat_p); }
        else {
            cerr << "Error: unknown bcondKind " << bc.bcondKind << endl;
            exit(EXIT_FAILURE);
        };
    }
}

