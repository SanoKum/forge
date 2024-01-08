#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include "mesh.hpp"
#include "elementType.hpp"
#include "vectorUtil.hpp"

#include <highfive/H5File.hpp>

using HighFive::File;

//void outputH5_XDMF(const mesh &msh , const variables &var)
void outputH5_XDMF(const solverConfig& cfg , const mesh& msh , const variables& var , const int& iStep)
{
    if (iStep%cfg.outStepInterval != 0) return;
    
    ostringstream oss;
    // ------------
    // *** HDF5 *** 
    // ------------
    oss << iStep;
    string fnameH5 = "vol_"+oss.str()+".h5";
    ofstream ofsH5(fnameH5);

    File file(fnameH5, File::ReadWrite | File::Truncate);

    // write mesh structure
    vector<geom_float> COORD;
    for (auto& nod : msh.nodes)
    {
        COORD.push_back(nod.coords[0]);
        COORD.push_back(nod.coords[1]);
        COORD.push_back(nod.coords[2]);
    }

    file.createDataSet("/MESH/COORD",COORD);

    vector<geom_int> CONNE;
    geom_int CONNE_dim = 0;
    geom_int CONNE0;    
    for (auto& cel : msh.cells)
    {
        geom_int nn = elementTypeMap.mapElementFromGmshID[cel.ieleType].nNodes;
        string name = elementTypeMap.mapElementFromGmshID[cel.ieleType].name;

        if (name == "hex") CONNE0 = 9;
        if (name == "prism") CONNE0 = 8;
        if (name == "pyramid") CONNE0 = 7;
        if (name == "tetra") CONNE0 = 6;

        //CONNE.push_back(nn + 1);
        CONNE.push_back(CONNE0);
        CONNE_dim += nn + 1;

        for (auto& nod : cel.iNodes)
        {
            CONNE.push_back(nod);
        }
    }
    file.createDataSet("/MESH/CONNE",CONNE);

    // write variables
    //for (string name : var.output_cellValNames)
    //TODO: ややこしくしているのでシンプルにoutput_cellValNamesで回したい。が、なぜかエラーになる
    for (auto& v : var.c)
    {
        string name = v.first;

        auto itr = std::find(var.output_cellValNames.begin(), var.output_cellValNames.end(), name);
        if (itr == var.output_cellValNames.end()) {
            continue; // notfound
        }

        file.createDataSet("/VALUE/"+name , v.second);
    }

    // ------------
    // *** XDMF ***
    // ------------
    string fnameXDMF = "vol_"+oss.str()+".xmf";
    ofstream ofs(fnameXDMF);

    ofs << "<?xml version='1.0' ?>\n";
    ofs << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n";
    ofs << "<Xdmf>\n";
    ofs << "  <Domain>\n";
    ofs << "    <Grid  GridType='Collection' CollectionType='Spatial' Name='Mixed'>\n";
    ofs << "      <Grid Name='aaa'>\n";
    ofs << "        <Topology Type='Mixed' NumberOfElements='" << msh.cells.size() << "'>\n";
    ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << CONNE_dim << "'>\n";
    ofs << "            vol_" << oss.str() <<".h5:MESH/CONNE\n";
    ofs << "          </DataItem>\n";
    ofs << "        </Topology>\n";

    ofs << "        <Geometry Type='XYZ'>\n";
    ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << msh.nodes.size()*3 << "'>\n";
    ofs << "            vol_"<< oss.str() <<".h5:MESH/COORD\n";
    ofs << "          </DataItem>\n";
    ofs << "        </Geometry>\n";

    for (string name : var.output_cellValNames)
    {
    //for (auto& v : var.c) {
        //string name = v.first;
        ofs << "        <Attribute Name='"  << name << "' Center='Cell' >\n";
        ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << msh.cells.size() << "'>\n";
        ofs << "            vol_"<< oss.str() << ".h5:VALUE/" << name << "\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Attribute>\n";
    }

    ofs << "      </Grid>\n";
    ofs << "    </Grid>\n";
    ofs << "  </Domain>\n";
    ofs << "</Xdmf>\n";
}

#endif
