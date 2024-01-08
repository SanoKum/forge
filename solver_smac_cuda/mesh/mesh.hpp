#pragma once

#include <vector>
#include <string>

#include "flowFormat.hpp"
#include "elementType.hpp"

struct node 
{
    std::vector<geom_int> iCells;
    std::vector<geom_int> iPlanes;
    std::vector<geom_float> coords; 

    node();
    node(geom_float& , geom_float& , geom_float& );
};

struct plane 
{
    std::vector<geom_int> iNodes;
    std::vector<geom_int> iCells;

    std::vector<geom_float> surfVect;
    geom_float surfArea;
    std::vector<geom_float> centCoords;
};

struct cell 
{
    std::vector<geom_int> iNodes;

    std::vector<geom_int> iPlanes;
    std::vector<geom_int> iPlanesDir;

    geom_float volume;

    geom_int ieleType;

    std::vector<geom_float> centCoords;

    cell();
    cell(std::vector<geom_int>&) ;
};

struct bcond
{
    geom_int physID; 
    std::string physName; 

    std::vector<geom_int> iPlanes;
    std::vector<geom_int> iBPlanes;
    std::vector<geom_int> iCells;

    std::map<std::string, int> inputInts;
    std::map<std::string, flow_float> inputFloats;

    std::string bcondKind; 

    std::map<std::string, std::vector<flow_float>> bvar;
    std::map<std::string, flow_float* > bvar_d; // cuda

    std::map<std::string,int> valueTypes; // {P, 0}, {T, 1}, etc.

    //cuda
    geom_int* map_bplane_plane_d;
    geom_int* map_bplane_cell_d;

    bcond();
    bcond(const geom_int& , const std::vector<geom_int>& , 
          const std::vector<geom_int>& , const std::vector<geom_int>& );
    ~bcond();

    void bcondInit(const std::string& ,
                   const std::map<std::string, int>& ,
                   const std::map<std::string, flow_float>& );
};

struct mesh 
{
public:
    geom_int  nNodes, nPlanes, nCells, nNormalPlanes, nBPlanes, nBconds;

    std::vector<node> nodes;
    std::vector<plane> planes;
    std::vector<cell> cells;
    std::vector<bcond> bconds;

    // cuda
    geom_int* map_nplane_cells_d;

    mesh();
    ~mesh();
    mesh(geom_int& , geom_int& ,geom_int& , geom_int& , 
         geom_int& , geom_int& ,
         std::vector<node>& , std::vector<plane>& , std::vector<cell>& , std::vector<bcond>& );

    void readMesh(std::string); 
    void setMeshMap_d();
};

struct matrix
{
private:
    geom_int itemp;
    geom_int ic0;
    geom_int ic1;
    geom_int ip;

    std::vector<geom_int> cellPlnCounter;

public:
    std::vector<std::vector<geom_int>> structure;
    std::vector<std::vector<flow_float>> lhs;
    std::vector<flow_float> rhs;

    std::vector<std::vector<geom_int>> localPlnOfCell;

    std::vector<flow_float> row_index;
    std::vector<flow_float> col_index;
    std::vector<flow_float> value;

    matrix();

    void initMatrix(mesh& );
};
