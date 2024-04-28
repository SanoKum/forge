#include "output.hpp"

#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include "mesh/mesh.hpp"
#include "mesh/elementType.hpp"
#include "common/vectorUtil.hpp"

#include <highfive/H5File.hpp>

using HighFive::File;

void outputH5_XDMF(const solverConfig& cfg , const mesh& msh , variables& var , const int& iStep)
{
    if (iStep%cfg.outStepInterval != 0) return;

    if (cfg.gpu == 1) {
        var.copyVariables_cell_D2H(var.output_cellValNames);
    }

    elementTypeMap eleTypeMap;
    
    ostringstream oss;

    // ------------
    // *** HDF5 *** 
    // ------------
    oss << iStep;
    string fnameH5 = "res_"+oss.str()+".h5";
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
    //for (auto& cel : msh.cells)
    for (geom_int ic = 0 ; ic<msh.nCells ; ic++)
    {
        auto cell = msh.cells[ic];
        geom_int nn = eleTypeMap.mapElementFromGmshID[cell.ieleType].nNodes;
        string name = eleTypeMap.mapElementFromGmshID[cell.ieleType].name;

        if (name == "hex") CONNE0 = 9;
        if (name == "prism") CONNE0 = 8;
        if (name == "pyramid") CONNE0 = 7;
        if (name == "tetra") CONNE0 = 6;

        //CONNE.push_back(nn + 1);
        CONNE.push_back(CONNE0);
        CONNE_dim += nn + 1;

        for (auto& nod : cell.iNodes)
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

        std::vector<flow_float> vtemp;
        vtemp.resize(msh.nCells);
        copy(v.second.begin(), v.second.begin()+msh.nCells, vtemp.begin());

        //file.createDataSet("/VALUE/"+name , v.second);
        file.createDataSet("/VALUE/"+name , vtemp);
    }

    // ------------
    // *** XDMF ***
    // ------------
    string fnameXDMF = "res_"+oss.str()+".xmf";
    ofstream ofs(fnameXDMF);

    ofs << "<?xml version='1.0' ?>\n";
    ofs << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n";
    ofs << "<Xdmf>\n";
    ofs << "  <Domain>\n";
    ofs << "    <Grid  GridType='Collection' CollectionType='Spatial' Name='Mixed'>\n";
    ofs << "      <Grid Name='aaa'>\n";
    ofs << "        <Topology Type='Mixed' NumberOfElements='" << msh.nCells << "'>\n";
    ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << CONNE_dim << "'>\n";
    ofs << "            res_" << oss.str() <<".h5:MESH/CONNE\n";
    ofs << "          </DataItem>\n";
    ofs << "        </Topology>\n";

    ofs << "        <Geometry Type='XYZ'>\n";
    ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << msh.nNodes*3 << "'>\n";
    ofs << "            res_"<< oss.str() <<".h5:MESH/COORD\n";
    ofs << "          </DataItem>\n";
    ofs << "        </Geometry>\n";

    for (string name : var.output_cellValNames)
    {
    //for (auto& v : var.c) {
        //string name = v.first;
        ofs << "        <Attribute Name='"  << name << "' Center='Cell' >\n";
        ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << msh.nCells << "'>\n";
        ofs << "            res_"<< oss.str() << ".h5:VALUE/" << name << "\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Attribute>\n";
    }

    ofs << "      </Grid>\n";
    ofs << "    </Grid>\n";
    ofs << "  </Domain>\n";
    ofs << "</Xdmf>\n";
}

void outputBconds_H5_XDMF(const solverConfig& cfg , mesh& msh , variables& var , const int& iStep)
{
    if (iStep%cfg.outStepInterval != 0) return;

    for (auto& bc : msh.bconds) {
        if (bc.bcondKind  != "wall") continue;

        bc.copyVariables_bplane_D2H();

        elementTypeMap eleTypeMap;
        ostringstream oss;

        // ------------
        // *** HDF5 *** 
        // ------------
        bc.output_preparation(msh.nodes, msh.planes);

        oss << iStep;
        string fnameH5 = "res_wall_"+oss.str()+".h5";
        ofstream ofsH5(fnameH5);

        File file(fnameH5, File::ReadWrite | File::Truncate);

        // write boundary
        vector<geom_float> COORD;
        geom_int inl;

        for (geom_int inl=0; inl<bc.inodes_l2g.size(); inl++) {
            geom_int ing = bc.inodes_l2g[inl];
            COORD.push_back(msh.nodes[ing].coords[0]);
            COORD.push_back(msh.nodes[ing].coords[1]);
            COORD.push_back(msh.nodes[ing].coords[2]);
        }

        file.createDataSet("/MESH/COORD",COORD);

        vector<geom_int> CONNE;
        geom_int CONNE_dim = 0;
        geom_int CONNE0;    

        for (auto& pln : bc.planes_local)
        {
            geom_int nn = pln.iNodes.size();
            if (nn == 3) { // tri
                string name = "triangle";
                CONNE0 = 4;
            } else if (nn == 4 ) {
                string name = "quad";
                CONNE0 = 5;
            }

            CONNE.push_back(CONNE0);
            CONNE_dim += nn + 1;

            for (auto& ing : pln.iNodes)
            {
                inl = bc.inodes_g2l[ing];
                CONNE.push_back(inl);
            }
        }
        file.createDataSet("/MESH/CONNE",CONNE);

        // write variables
        for (auto& v : bc.bvar)
        {
            string name = v.first;

            //auto itr = std::find(var.output_cellValNames.begin(), var.output_cellValNames.end(), name);
            //if (itr == var.output_cellValNames.end()) {
            //    continue; // notfound
            //}

            std::vector<flow_float> vtemp;
            vtemp.resize(bc.planes_local.size());
            copy(v.second.begin(), v.second.begin()+bc.planes_local.size(), vtemp.begin());

            file.createDataSet("/VALUE/"+name , vtemp);
        }

        // ------------------
        // *** Ghost cell ***
        // ------------------

        vector<geom_float> COORD_ghst;

        // out ghost cell 
        for (auto& ighst : bc.iCells_ghst) {
            // write boundary
            COORD_ghst.push_back(msh.cells[ighst].centCoords[0]);
            COORD_ghst.push_back(msh.cells[ighst].centCoords[1]);
            COORD_ghst.push_back(msh.cells[ighst].centCoords[2]);
        }

        file.createDataSet("/MESH_ghst/COORD",COORD_ghst);

        vector<geom_int> CONNE_ghst;
        vector<geom_int> Indexes_ghst;
        geom_int CONNE_dim_ghst = 0;
        geom_int CONNE0_ghst;    

        geom_int ighst_l = 0;
        for (auto& ighst : bc.iCells_ghst) {
            geom_int nn = 1;
            string name = "point";
            CONNE0_ghst = 1;

            CONNE_ghst.push_back(CONNE0_ghst);
            CONNE_dim_ghst += nn + 1;

            CONNE_ghst.push_back(ighst_l);
            Indexes_ghst.push_back(ighst_l);
            ighst_l++;
        }
        file.createDataSet("/MESH_ghst/CONNE",CONNE_ghst);
        file.createDataSet("/MESH_ghst/Indexes",Indexes_ghst);

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
    
            std::vector<flow_float> vtemp;
            vtemp.resize(bc.iCells_ghst.size());
            for (geom_int igl = 0 ; igl<bc.iCells_ghst.size() ; igl++) {
                geom_int ig = bc.iCells_ghst[igl];
                vtemp[igl] = v.second[ig];
            }
    
            //file.createDataSet("/VALUE/"+name , v.second);
            file.createDataSet("/VALUE_ghst/"+name , vtemp);
        }

        // ------------
        // *** XDMF ***
        // ------------
        string fnameXDMF = "res_wall_"+oss.str()+".xmf";
        ofstream ofs(fnameXDMF);

        ofs << "<?xml version='1.0' ?>\n";
        ofs << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n";
        ofs << "<Xdmf>\n";
        ofs << "  <Domain>\n";
        ofs << "    <Grid  GridType='Collection' CollectionType='Spatial' Name='Mixed'>\n";
        ofs << "      <Grid Name='wall'>\n";
        ofs << "        <Topology Type='Mixed' NumberOfElements='" << bc.iPlanes.size() << "'>\n";
        ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << CONNE_dim << "'>\n";
        ofs << "            res_wall_" << oss.str() <<".h5:MESH/CONNE\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Topology>\n";

        ofs << "        <Geometry Type='XYZ'>\n";
        ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << bc.inodes_l2g.size()*3 << "'>\n";
        ofs << "            res_wall_"<< oss.str() <<".h5:MESH/COORD\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Geometry>\n";

        for (string name : bc.bplaneValNames)
        {
        //for (auto& v : var.c) {
            //string name = v.first;
            ofs << "        <Attribute Name='"  << name << "' Center='Cell' >\n";
            ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << bc.iPlanes.size() << "'>\n";
            ofs << "            res_wall_"<< oss.str() << ".h5:VALUE/" << name << "\n";
            ofs << "          </DataItem>\n";
            ofs << "        </Attribute>\n";
        }

        ofs << "      </Grid>\n";

        // ghost cell
        ofs << "      <Grid Name='wall_ghst'>\n";
        ofs << "        <Topology Type='Polyvertex' Dimensions='" << bc.iCells_ghst.size() <<"' NodesPerElement='1'>\n";
        //ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << CONNE_dim_ghst << "'>\n";
        ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << bc.iCells_ghst.size() << "'>\n";
        //ofs << "            res_wall_" << oss.str() <<".h5:MESH_ghst/CONNE\n";
        ofs << "            res_wall_" << oss.str() <<".h5:MESH_ghst/Indexes\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Topology>\n";

        ofs << "        <Geometry Type='XYZ'>\n";
        ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << bc.iCells_ghst.size()*3 << "'>\n";
        ofs << "            res_wall_"<< oss.str() <<".h5:MESH_ghst/COORD\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Geometry>\n";

        for (string name : var.output_cellValNames)
        {
            ofs << "        <Attribute Name='"  << name << "' Center='Node' >\n";
            ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << bc.iCells_ghst.size() << "'>\n";
            ofs << "            res_wall_"<< oss.str() << ".h5:VALUE_ghst/" << name << "\n";
            ofs << "          </DataItem>\n";
            ofs << "        </Attribute>\n";
        }
        ofs << "      </Grid>\n";

        ofs << "    </Grid>\n";
        ofs << "  </Domain>\n";
        ofs << "</Xdmf>\n";
    }
}
