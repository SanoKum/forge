#include "mesh.hpp"

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

using namespace std;
using namespace HighFive;

node::node() {};
node::node(geom_float &x, geom_float &y, geom_float &z)
{
    this->coords.push_back(x);
    this->coords.push_back(y);
    this->coords.push_back(z);
}


cell::cell() {}
cell::cell(vector<geom_int> &iNodes) 
{
    this->iNodes = iNodes;
}

bcond::bcond() {}
bcond::bcond(const geom_int& physID, const vector<geom_int>& iPlanes, 
             const vector<geom_int>& iCells , const vector<geom_int>& iBPlanes)
{
    this->physID   = physID;
    this->iPlanes  = iPlanes;
    this->iCells   = iCells;
    this->iBPlanes = iBPlanes;
}

void bcond::bcondInit(const string &bcondKind,
                      const std::map<std::string, int> &inputInts,
                      const std::map<std::string, flow_float> &inputFloats)
{
    this->bcondKind = bcondKind;
    this->inputInts = inputInts;
    this->inputFloats = inputFloats;
}

mesh::mesh(){}
mesh::mesh(geom_int& nNodes,geom_int& nPlanes,geom_int& nCells, geom_int& nNormalPlanes, 
    geom_int& nBPlanes, geom_int& nBconds,
    vector<node> &nodes , vector<plane> &planes , vector<cell>& cells , vector<bcond>& bconds)
{
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

void mesh::readMesh(string fname)
{
    File file(fname, File::ReadOnly);

    // read basic 
    Group group = file.getGroup("/MESH");

    Attribute a = group.getAttribute("nNodes");
    this->nNodes = a.read<geom_int>();

    a = group.getAttribute("nCells");
    this->nCells = a.read<geom_int>();

    a = group.getAttribute("nPlanes");
    this->nPlanes = a.read<geom_int>();

    a = group.getAttribute("nNormalPlanes");
    this->nNormalPlanes = a.read<geom_int>();

    a = group.getAttribute("nBPlanes");
    this->nBPlanes= a.read<geom_int>();

    a = group.getAttribute("nBconds");
    this->nBconds = a.read<geom_int>();

    // nodes
    this->nodes.resize(this->nNodes);
    std::vector<geom_float> coord;
    file.getDataSet("/MESH/COORD").read(coord);

    geom_int ii = 0;
    for (geom_int i=0; i<(this->nNodes); i++)
    {
        nodes[i].coords.resize(3);
        nodes[i].coords[0] = coord[ii];
        ii += 1;
        nodes[i].coords[1] = coord[ii];
        ii += 1;
        nodes[i].coords[2] = coord[ii];
        ii += 1;
    }
    
    // planes
    this->planes.resize(this->nPlanes);
    std::vector<geom_int> strct;
    std::vector<geom_float> surfVect;
    std::vector<geom_float> surfArea;
    std::vector<geom_float> centCoords;
    file.getDataSet("/PLANES/STRUCT").read(strct);
    file.getDataSet("/PLANES/surfVect").read(surfVect);
    file.getDataSet("/PLANES/surfArea").read(surfArea);
    file.getDataSet("/PLANES/centCoords").read(centCoords);

    geom_int ipp = 0;
    for (geom_int ip=0; ip<this->nPlanes; ip++)
    {
        geom_int nn = strct[ipp];
        this->planes[ip].iNodes.resize(nn);
        ipp += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->planes[ip].iNodes[in] = strct[ipp];
            ipp += 1;
        }

        nn = strct[ipp];
        this->planes[ip].iCells.resize(nn);
        ipp += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->planes[ip].iCells[in] = strct[ipp];
            ipp += 1;
        }

        this->planes[ip].surfVect.resize(3);
        this->planes[ip].surfVect[0] = surfVect[3*ip + 0];
        this->planes[ip].surfVect[1] = surfVect[3*ip + 1];
        this->planes[ip].surfVect[2] = surfVect[3*ip + 2];

        this->planes[ip].surfArea    = surfArea[ip];

        this->planes[ip].centCoords.resize(3);
        this->planes[ip].centCoords[0] = centCoords[3*ip + 0];
        this->planes[ip].centCoords[1] = centCoords[3*ip + 1];
        this->planes[ip].centCoords[2] = centCoords[3*ip + 2];
    }
 
    // cells
    this->cells.resize(this->nCells);
    std::vector<geom_int> strct2;
    file.getDataSet("/CELLS/STRUCT").read(strct2);

    std::vector<geom_float> volume;
    std::vector<geom_float> centCoords2;
    file.getDataSet("/CELLS/volume").read(volume);
    file.getDataSet("/CELLS/centCoords").read(centCoords2);

    geom_int icc = 0;
    for (geom_int ic=0; ic<this->nCells; ic++)
    {
        geom_int nn = strct2[icc];
        this->cells[ic].iNodes.resize(nn);
        icc += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->cells[ic].iNodes[in] = strct2[icc];
            icc += 1;
        }

        nn = strct2[icc];
        this->cells[ic].iPlanes.resize(nn);
        icc += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->cells[ic].iPlanes[in] = strct2[icc];
            icc += 1;
        }

        nn = strct2[icc];
        this->cells[ic].iPlanesDir.resize(nn);
        icc += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->cells[ic].iPlanesDir[in] = strct2[icc];
            icc += 1;
        }
        this->cells[ic].ieleType = strct2[icc];
        icc += 1;

        this->cells[ic].volume = volume[ic];

        this->cells[ic].centCoords.resize(3);
        this->cells[ic].centCoords[0] = centCoords2[3*ic + 0];
        this->cells[ic].centCoords[1] = centCoords2[3*ic + 1];
        this->cells[ic].centCoords[2] = centCoords2[3*ic + 2];
    }

    // boundary conditions
    Group grp = file.getGroup("/BCONDS");
    geom_int nb = grp.getNumberObjects();
    bconds.resize(nb);

    geom_int ib = 0;
    for (string oname : grp.listObjectNames())
    {
        this->bconds[ib].physID = stoi(oname);
        grp = file.getGroup("/BCONDS/"+oname);

        a = grp.getAttribute("bcondKind");
        this->bconds[ib].bcondKind = a.read<std::string>();

        // iCells
        std::vector<geom_int> iCells;
        grp.getDataSet("iCells").read(iCells);

        this->bconds[ib].iCells.resize(iCells.size());
        for (geom_int ic = 0 ; ic<iCells.size() ; ic++)
        {
            this->bconds[ib].iCells[ic] = iCells[ic];
        }

        // iPlanes
        std::vector<geom_int> iPlanes;
        grp.getDataSet("iPlanes").read(iPlanes);

        this->bconds[ib].iPlanes.resize(iPlanes.size());
        for (geom_int ip = 0 ; ip<iPlanes.size() ; ip++)
        {
            this->bconds[ib].iPlanes[ip] = iPlanes[ip];
        }

        ib += 1;
    }
}

matrix::matrix(){}
void matrix::initMatrix(mesh& msh)
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