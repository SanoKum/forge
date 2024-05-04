#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <string>
//#include <Eigen/Dense>
#include "flowFormat.hpp"
#include "elementType.hpp"

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

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
    std::vector<geom_int> iCells; //chk

    std::vector<geom_float> surfVect;
    geom_float surfArea;
    std::vector<geom_float> centCoords;
};

struct cell 
{
    std::vector<geom_int> iNodes;

    std::vector<geom_int> iPlanes; //chk
    std::vector<geom_int> iPlanesDir;

    geom_float volume; //chk

    geom_int ieleType;

    std::vector<geom_float> centCoords;

    cell();
    cell(std::vector<geom_int>&) ;
};

struct bcond
{
public:
    std::list<std::string> bplaneValNames =  
    {
        "ro"   , "roUx" , "roUy" , "roUz" , "roUz" , "roe",
        "Ux"   , "Uy"   , "Uz"   , "Tt"   , "Pt"   , "Ts" , "Ps",
        "ypls" , "twall_x" , "twall_y" , "twall_z" ,
        "sx"   , "sy"   , "sz"   , "ss"
    };

    std::list<std::string> bplaneIntNames =  
    {
        "partnerPlnID" // for periodic 
    };


    geom_int physID; 
    std::string physName; 

    std::vector<geom_int> iPlanes;
    std::vector<geom_int> iBPlanes;
    std::vector<geom_int> iCells; //chk
    std::vector<geom_int> iCells_ghst; //ghost cell

    int output_preparation_flg = 0;
    std::vector<node> nodes_local;
    std::vector<plane> planes_local;
    std::vector<geom_int> inodes_l2g; 
    std::vector<geom_int> inodes_g2l; 

    std::map<std::string, int> inputInts;
    std::map<std::string, flow_float> inputFloats;

    std::string bcondKind; 

    std::map<std::string,int> valueTypes; // {P, 0}, {T, 1}, etc.

    std::map<std::string, std::vector<flow_float>> bvar;
    std::map<std::string, flow_float* > bvar_d; // cuda

    // interger values (only for periodic now)
    std::map<std::string, std::vector<flow_float>> bint;
    std::map<std::string, flow_float* > bint_d; // cuda

    //cuda
    geom_int* map_bplane_plane_d;
    geom_int* map_bplane_cell_d;
    geom_int* map_bplane_cell_ghst_d;

    bcond();
    bcond(const geom_int& , const std::vector<geom_int>& , 
          const std::vector<geom_int>& , const std::vector<geom_int>& );
    ~bcond();

    void bcondInitVariables(const int &useGPU);
    void copyVariables_bplane_D2H();
    //void set_nodes_local(mesh& msh);
    void output_preparation(std::vector<node>& nodes, std::vector<plane>& planes);
};

struct mesh 
{
public:
    geom_int  nNodes, nPlanes, nCells, nNormalPlanes, nBPlanes, nBconds;
    geom_int  nCells_all, nCells_ghst , nCells_halo;

    std::vector<node> nodes;
    std::vector<plane> planes;
    std::vector<cell> cells;
    std::vector<bcond> bconds;

    // cuda
    //geom_int* map_nplane_cells_d; // normal plane
    geom_int* map_plane_cells_d; // 

    geom_int* map_cell_planes_index_d; 
    geom_int* map_cell_planes_d; 

    mesh();
    ~mesh();
    mesh(geom_int& , geom_int& ,geom_int& , geom_int& , 
         geom_int& , geom_int& ,
         std::vector<node>& , std::vector<plane>& , std::vector<cell>& , std::vector<bcond>& );

    void readMesh(std::string); 

    void setPeriodicPartner();
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
