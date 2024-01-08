#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>

#include "flowFormat.hpp"
#include "elementType.hpp"

#include <stdio.h>                                                                                       

using namespace std;

struct node 
{
    vector<geom_int> iCells;
    vector<geom_int> iPlanes;
    vector<geom_float> coords; 

    node() {}
    node(geom_float &x, geom_float &y, geom_float &z)
    {
        this->coords.push_back(x);
        this->coords.push_back(y);
        this->coords.push_back(z);
    }
};

struct plane 
{
    vector<geom_int> iNodes;
    vector<geom_int> iCells;

    vector<geom_float> surfVect;
    geom_float surfArea;
    vector<geom_float> centCoords;
};

struct cell 
{
    vector<geom_int> iNodes;

    vector<geom_int> iPlanes;
    vector<geom_int> iPlanesDir;

    geom_float volume;

    geom_int ieleType;

    vector<geom_float> centCoords;

    cell() {}
    cell(vector<geom_int> &iNodes) 
    {
        this->iNodes = iNodes;
    }
};

struct bcond
{
    geom_int physID; 
    string physName; 

    vector<geom_int> iPlanes;
    vector<geom_int> iBPlanes;
    vector<geom_int> iCells;

    map<string, int> inputInts;
    map<string, flow_float> inputFloats;

    map<string, int> neuDirFlag; // 0:No Grad, 1:specified Flux, 2: Dirclet

    string bcondKind; 

    bcond() {}
    bcond(const geom_int &physID, const vector<geom_int> &iPlanes, 
          const vector<geom_int> &iCells , const vector<geom_int> &iBPlanes)
    {
        this->physID   = physID;
        this->iPlanes  = iPlanes;
        this->iCells   = iCells;
        this->iBPlanes = iBPlanes;
    }

    void bcondInit(const string &bcondKind,
                   const std::map<std::string, int> &inputInts,
                   const std::map<std::string, flow_float> &inputFloats)
    {
        this->bcondKind = bcondKind;
        this->inputInts = inputInts;
        this->inputFloats = inputFloats;
    }
};

struct mesh 
{
public:
    geom_int  nNodes, nPlanes, nCells, nNormalPlanes, nBPlanes, nBconds;

    vector<node> nodes;
    vector<plane> planes;
    vector<cell> cells;
    vector<bcond> bconds;

    mesh(){}
    mesh(geom_int &nNodes,geom_int &nPlanes,geom_int &nCells, geom_int &nNormalPlanes, 
         geom_int &nBPlanes, geom_int &nBconds,
         vector<node> &nodes , vector<plane> &planes , vector<cell> &cells , vector<bcond> &bconds){
        this->nNodes = nNodes;
        this->nPlanes = nPlanes;
        this->nCells = nCells;
        this->nNormalPlanes= nNormalPlanes;
        this->nBPlanes = nBPlanes;
        this->nBconds = nBconds;
        this->nodes = nodes;
        this->planes = planes;
        this->cells = cells;
        this->bconds = bconds;
    }
};

struct matrix
{
private:
    geom_int itemp;
    geom_int ic0;
    geom_int ic1;
    geom_int ip;

    vector<geom_int> cellPlnCounter;

public:
    vector<vector<geom_int>> structure;
    vector<vector<flow_float>> lhs;
    vector<flow_float> rhs;

    vector<vector<geom_int>> localPlnOfCell;

    matrix(){}

    void initPoisson(mesh& msh)
    {
        structure.resize(msh.nCells);
        lhs.resize(msh.nCells);
        rhs.resize(msh.nCells);

        cellPlnCounter.resize(msh.nCells);
        localPlnOfCell.resize(msh.nNormalPlanes);

        for (ic0=0; ic0<msh.nCells; ic0++)
        {
            structure[ic0].push_back(ic0);
        }

        for (auto& i : cellPlnCounter)
        {
            i = 1;
        }

        //for (ic=0; ic<msh.nCells; ic++)
        for (ip=0; ip<msh.nNormalPlanes; ip++)
        {
            ic0 = msh.planes[ip].iCells[0];
            ic1 = msh.planes[ip].iCells[1];
            structure[ic0].push_back(ic1);
            structure[ic1].push_back(ic0);

            localPlnOfCell[ip].push_back(cellPlnCounter[ic0]);
            localPlnOfCell[ip].push_back(cellPlnCounter[ic1]);

            cellPlnCounter[ic0] += 1;
            cellPlnCounter[ic1] += 1;
        }

        for (geom_int ist=0; ist<structure.size(); ist++)
        {
            lhs[ist].resize(structure[ist].size());
        }
    }
};


#endif