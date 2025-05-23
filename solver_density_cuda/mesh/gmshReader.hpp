#pragma once

//TODO: public と privateの使い分けがわからない。
//TODO: あと、physnameとかを整える

#include <vector>
#include <cmath>
#include <string>

#include <flowFormat.hpp>
#include <mesh/elementType.hpp>
#include <mesh/mesh.hpp>
#include <common/stringUtil.hpp>
#include "common/vectorUtil.hpp"
#include "variables.hpp"

#include <sstream>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

#include <stdio.h>                                                                                       

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

using namespace std;
//using HighFive::File;
using namespace HighFive;

class gmshReader : public mesh
{
public:
    class surfEnt
    {
    public:
        geom_int entTag;
        geom_float minX, minY, minZ;
        geom_float maxX, maxY, maxZ;
        geom_int nPhysTag;
        //vector<geom_int> physTags;
        geom_int physTag;
        //geom_int nBoundCurves;
        //vector<geom_int> curveTags;


        surfEnt(){};
        surfEnt(geom_int tag, geom_float minX, geom_float minY, geom_float minZ,
                              geom_float maxX, geom_float maxY, geom_float maxZ,
                              geom_int nPhysTag, geom_int physTag)
                              //geom_int nBoundCurves, vector<geom_int> curveTags  )
        {
            this->entTag = tag;
            this->minX = minX; this->minY = minY; this->minZ = minZ;
            this->maxX = maxX; this->maxY = maxY; this->maxZ = maxZ;
            this->nPhysTag = nPhysTag;
            this->physTag  = physTag;
            //this->nBoundCurves = nBoundCurves;
            //this->curveTags = curveTags;
        };
    };

    class volumeEnt
    {
    public:
        geom_int entTag;
        geom_float minX, minY, minZ;
        geom_float maxX, maxY, maxZ;
        geom_int nPhysTag;
        geom_int physTag;

        volumeEnt(){};
        volumeEnt(geom_int tag, geom_float minX, geom_float minY, geom_float minZ,
                                geom_float maxX, geom_float maxY, geom_float maxZ,
                                geom_int nPhysTag, geom_int physTag)
        {
            this->entTag = tag;
            this->minX = minX; this->minY = minY; this->minZ = minZ;
            this->maxX = maxX; this->maxY = maxY; this->maxZ = maxZ;
            this->nPhysTag = nPhysTag;
            this->physTag = physTag;
            //this->nBoundSurfaces = nBoundSurfaces;
            //this->surfaceTags = surfaceTags;
        }
    };

    class element
    {
    public:
        geom_int eleTag;
        vector<geom_int> iNodes;

        //Element(geom_int eleTag, vector<geom_int> &iNodes)
        element(geom_int eleTag, vector<geom_int> &iNodes)
        {
            this->eleTag = eleTag;
            this->iNodes = iNodes;
        }
    };

    class elementsOfEntity
    {
    public:
        geom_int dimension;
        geom_int entTag;
        geom_int ieleType;
        geom_int nEleEnt;
        vector<element> elements;

        elementsOfEntity(geom_int dimension, geom_int entTag, geom_int ieleType, geom_int nElements,
                 vector<element> &elements)
        {
            this->dimension = dimension;
            this->entTag = entTag;
            this->ieleType = ieleType;
            this->nEleEnt = nElements;
            this->elements = elements;
            //this->nBoundSurfaces = nBoundSurfaces;
            //this->surfaceTags = surfaceTags;
        }
    };

    list<elementsOfEntity> elements_summary;

    geom_int nSurfEnt;
    geom_int nVolumeEnt;
    map<geom_int,surfEnt> surfEntMap;
    map<geom_int,volumeEnt> volumeEntMap;

    string gmshVersion;
    geom_int nPhysName;
    list<geom_int> physDim;
    list<geom_int> physID;
    list<string> physName;

    geom_int nElements;

    elementTypeMap eleTypeMap;

    mesh getMesh()
    {
        mesh msh = mesh(this->nNodes, this->nPlanes, this->nCells, this->nNormalPlanes, 
                        this->nBPlanes, this->nBconds,
                        this->nodes, this->planes, this->cells, this->bconds);
        return msh;
    }


    gmshReader(const string inputFileName)
    {
        ifstream inputFile(inputFileName);

        string line;
        vector<string> l_str;
        vector<string> lines;

        // *** read Gmsh ***
        this->readMeshFormat(inputFile);
        this->readPhysicalNames(inputFile);
        this->readEntities(inputFile);
        this->readNodes(inputFile);
        this->readElements(inputFile);

        // *** Make Grid ***
        cout << "makeMesh in gmsh Reader\n";
        this->makeMesh(); // nodes & planes & cells & bconds are made.

    }

    void readMeshFormat(ifstream &inputFile)
    {
        string line;
        vector<string> l_str;

        getline(inputFile, line);
        getline(inputFile, line);

        boost::split(l_str, line, boost::is_space());

        this->gmshVersion = l_str[0];

        if (l_str[0] == "4.1" | l_str[0] == "2.2" ) {
            cout << "GMSH version is " << this->gmshVersion << endl;

        } else {
            cerr << "Error: unknown gmsh format" << gmshVersion << endl;                                                                                                     
            exit(EXIT_FAILURE);
        }

        getline(inputFile, line); // end mesh format
    }

    void readPhysicalNames(ifstream &inputFile)
    {
        string line;
        vector<string> l_str;

        getline(inputFile, line);

        boost::split(l_str, line, boost::is_space());

        if (l_str[0] != "$PhysicalNames") {
            cerr << "Error: unknown gmsh format phys" << endl;                                                                                                     
            exit(EXIT_FAILURE);
        }

        getline(inputFile, line);
        this->nPhysName = stoi(line);

        for (geom_int i = 0 ; i < nPhysName ; i++)
        {
            getline(inputFile, line);
            boost::split(l_str, line, boost::is_space());

            geom_int physDim = stoi(l_str[0]);
            geom_int physID  = stoi(l_str[1]);
            string physName  = l_str[2];

            this->physDim.push_back(physDim);
            this->physID.push_back(physID);
            this->physName.push_back(physName);
        }
        getline(inputFile, line);
    }

    void readEntities(ifstream &inputFile)
    {
        string line;
        vector<string> l_str;

        getline(inputFile, line);
        getline(inputFile, line);

        boost::split(l_str, line, boost::is_space());
        this->nSurfEnt   = stoi(l_str[2]);
        this->nVolumeEnt = stoi(l_str[3]);

        // Point Entities : just ignore
        for (geom_int i = 0 ; i<stoi(l_str[0]); i++)
        {
            getline(inputFile, line);
        }

        // Curve Entities : just ignore
        for (geom_int i = 0 ; i<stoi(l_str[1]); i++)
        {
            getline(inputFile, line);
        }

        // Surface Entities 
        for (geom_int i = 0 ; i<this->nSurfEnt ; i++)
        {
            getline(inputFile, line);
            boost::split(l_str, line, boost::is_space());

            surfEnt surfEnt_temp = surfEnt(stoi(l_str[0]), stof(l_str[1]), stof(l_str[2]), stof(l_str[3]),
                                           stof(l_str[4]), stof(l_str[5]), stof(l_str[6]),
                                           stoi(l_str[7]), stoi(l_str[8]));
            if (surfEnt_temp.nPhysTag > 1)
            {
                std::cerr << "Error: unknown gmsh format. number of PhysTag is not 1" << std::endl;                                                                                                     
                std::cerr << "nPhysTag = " << surfEnt_temp.nPhysTag << " " << surfEnt_temp.physTag  << std::endl;                                                                                                     
                std::cerr << "Ent id = " << i   << std::endl;                                                                                                     

                //BOOST_FOREACH (std::string s, l_str)
                //{
                //    std::cout << s << std::endl;
                //}


                exit(EXIT_FAILURE);
            }

            surfEntMap.insert(std::make_pair(surfEnt_temp.entTag, surfEnt_temp));
        };

        // Volume Entities 
        for (geom_int i = 0 ; i<this->nVolumeEnt ; i++)
        {
            getline(inputFile, line);
            boost::split(l_str, line, boost::is_space());

            volumeEnt volumeEnt_temp = volumeEnt(stoi(l_str[0]), stof(l_str[1]), stof(l_str[2]), stof(l_str[3]),
                                                 stof(l_str[4]), stof(l_str[5]), stof(l_str[6]),
                                                 stoi(l_str[7]), stoi(l_str[8]));
            if (volumeEnt_temp.nPhysTag != 1)
            {
                std::cerr << "Error: unknown gmsh format. number of PhysTag is not 1" << std::endl;                                                                                                     
                exit(EXIT_FAILURE);
            }

            //volumeEntList.push_back(volumeEnt_temp);
            volumeEntMap.insert(std::make_pair(volumeEnt_temp.entTag, volumeEnt_temp));
        };

        getline(inputFile, line); // read $EndEntities

    }

    void readNodes(ifstream &inputFile)
    {
        string line;
        vector<string> l_str;

        getline(inputFile, line);
        getline(inputFile, line);

        boost::split(l_str, line, boost::is_space());
        this->nNodes = stoi(l_str[1]);
        cout << "numNodes = " << this->nNodes << endl;

        if (l_str[2]!="1" | l_str[3]!=l_str[1])
        {
            cerr << "Error: unknown gmsh format. Nodes numbers" << endl;                                                                                                     
            exit(EXIT_FAILURE);
        }

        geom_int id_now = 0;

        while (id_now < this->nNodes)
        {
            getline(inputFile, line);
            boost::split(l_str, line, boost::is_space());
            geom_int nn = stoi(l_str[3]);

            id_now += nn;

            for (geom_int i = 0 ; i < nn ; i++)
            {
                getline(inputFile, line);
            }

            for (geom_int i = 0 ; i < nn ; i++)
            {
                getline(inputFile, line);
                boost::split(l_str, line, boost::is_space());
                geom_float x = stof(l_str[0]);
                geom_float y = stof(l_str[1]);
                geom_float z = stof(l_str[2]);

                node node_temp = node(x, y, z);
                this->nodes.push_back(node_temp);
            }

        }

        while (line.substr(0, 9) != "$EndNodes")
        {
            getline(inputFile, line); 
        }
    }

    void readElements(ifstream &inputFile) 
    {
        string line;
        vector<string> l_str;

        cout << "read elements" << endl;

        getline(inputFile, line); // $Elements
        getline(inputFile, line); // numEntities, numElements, iBegin, iEnd

        boost::split(l_str, line, boost::is_space());
        this->nElements= stoi(l_str[1]);
        cout << "numElements = " << this->nElements << endl;

        if (l_str[2]!="1" | l_str[3]!=l_str[1])
        {
            cerr << "Error: unknown gmsh format. Element numbers" << endl;                                                                                                     
            exit(EXIT_FAILURE);
        }

        geom_int iEle = 0;

        //list<elementsOfEntity> elements_summary;

        while (iEle < this->nElements)
        {
            getline(inputFile, line);
            boost::split(l_str, line, boost::is_space());
            geom_int dimension   = stoi(l_str[0]);
            geom_int entTag      = stoi(l_str[1]);
            geom_int eleTypeGmsh = stoi(l_str[2]);
            geom_int nEleEnt     = stoi(l_str[3]);

            boost::split(l_str, line, boost::is_space());

            elementTypeFormat eleType = this->eleTypeMap.mapElementFromGmshID[eleTypeGmsh];

            iEle += nEleEnt;

            vector<element> elements_temp;

            for (geom_int i = 0 ; i < nEleEnt ; i++)
            {
                getline(inputFile, line);
                boost::split(l_str, line, boost::is_space());

                vector<geom_int> nodes_temp;

                for (geom_int j = 0; j<eleType.nNodes; j++){
                    nodes_temp.push_back( stoi(l_str[j+1]) -1 ); //gmsh starts with 1. 
                }

                element ele_temp = element(stoi(l_str[0])-1, nodes_temp); //gmsh starts with 1.

                elements_temp.push_back(ele_temp);
            }

            elementsOfEntity elementsOfEnt = elementsOfEntity(dimension, entTag, eleTypeGmsh, nEleEnt,
                                                              elements_temp );

            this->elements_summary.push_back(elementsOfEnt);
        }
    }

    void makeMesh()
    {
        // *** Make Grid ***
        geom_int nCells_temp = 0;

        // iterate over each Entity.
        for (auto &itEnt : elements_summary)
        {
            if (itEnt.dimension != 3) continue; // skip surface. only volume

            geom_int physTag = volumeEntMap[itEnt.entTag].physTag; // physical tag (fluid, inlet, etc.)
            elementTypeFormat eleType = this->eleTypeMap.mapElementFromGmshID[itEnt.ieleType]; // quad, hex, etc.

            vector<vector<geom_int>> nOrderPlanes = eleType.nodesOrderPlanes;

            // elementを構成するノード: itEnt.elementsに入っている
            // elementを構成するプレーン: ローカルにはnOrderPlanesに順序が入っており、グローバルには要再構成
            // *やること
            // プレーンを一つのリストにまとめていく。同時に、プレーンが接するelementも保存していくぞ。
            // ノードに接するセルもまとめていく（重複したプレーンの削除に活用予定）
            // まずはセルを記録していくとともに、ノードが保有するセルも整理していく。

            // iterate over each elements of a entity.
            //for (auto itEle = itEnt.elements.begin(), e2=itEnt.elements.end(); itEle != e2; ++itEle)

            for (auto& itEle : itEnt.elements)
            {
                cell cell_now;
                //cell_now.nNodes = eleType.nNodes;
                cell_now.iNodes = itEle.iNodes;
                cell_now.ieleType = itEnt.ieleType;
                this->cells.push_back(cell_now);

                nCells_temp += 1;

                // cellsOfNodes に、nodeが保有するcellを記録していく
                for (auto &nodGlo : itEle.iNodes)
                {
                    this->nodes[nodGlo].iCells.push_back(nCells_temp-1);
                }
                
            }
        }

        vector<int> cellCheckedFlag(cells.size());

        for (auto& cCF : cellCheckedFlag) {
            cCF = 0;
        }

        // create & merge planes.
        geom_int nPlanes_temp = 0;
        geom_int icell = 0;

        for (auto& cel : cells)
        {
            elementTypeFormat eleType = this->eleTypeMap.mapElementFromGmshID[cel.ieleType]; // quad, hex, etc.
            vector<vector<geom_int>> nOrderPlanes = eleType.nodesOrderPlanes;

            for (auto& plnLocal : nOrderPlanes)
            {
                vector<geom_int> nodes_temp;
                for (auto& nodLocal : plnLocal)
                {
                    geom_int nodGlobal = cel.iNodes[nodLocal];
                    nodes_temp.push_back(nodGlobal);
                }

                vector<geom_int> cellPot;
                for (auto& nod : nodes_temp)
                {
                    //cout << "nod=" << nod << endl;
                    if (cellPot.size() == 0)
                    {
                        cellPot = this->nodes[nod].iCells;
                        //print1Dvector(cellPot);
                    } else {
                        vector<geom_int> cellPot2;

                        //cout << "icell=" << icell << endl;

                        for (auto& ic : nodes[nod].iCells)
                        {
                            vector<geom_int>:: iterator itr;
                            itr = find(cellPot.begin() , cellPot.end(), ic);

                            if (itr == cellPot.end()) continue;
                            geom_int wanted_index = distance(cellPot.begin(), itr);
                            cellPot2.push_back(cellPot[wanted_index]);
                        }
                        cellPot = cellPot2;
                    }
                }

                if (cellPot.size() == 2) // normal plane
                {
                    geom_int icell1 = icell;
                    geom_int icell2 = cellPot[0] + cellPot[1] - icell1;

                    if (cellCheckedFlag[icell2] == 1) // already checked. ignore.
                    {
                        continue;
                    } else {  // new normal plane
                        plane plane_now;
                        plane_now.iNodes = nodes_temp;
                        plane_now.iCells.push_back(icell);
                        plane_now.iCells.push_back(icell2);

                        this->planes.push_back(plane_now);

                        this->cells[icell].iPlanes.push_back(nPlanes_temp);
                        this->cells[icell].iPlanesDir.push_back(1);

                        this->cells[icell2].iPlanes.push_back(nPlanes_temp);
                        this->cells[icell2].iPlanesDir.push_back(-1);

                        nPlanes_temp += 1;
                    }
                } else if (cellPot.size() == 1) // boundaryPlane
                {
                    if (cellPot[0] != icell) {
                       cerr << "something wrong in cellPot" << endl;
                       exit(EXIT_FAILURE);
                    }

                    plane plane_now;
                    plane_now.iNodes = nodes_temp;
                    plane_now.iCells.push_back(icell);
                    this->planes.push_back(plane_now);

                    this->cells[icell].iPlanes.push_back(nPlanes_temp);
                    this->cells[icell].iPlanesDir.push_back(1);

                    nPlanes_temp += 1;
                } else {
                    cerr << "something wrong in cellPot" << endl;
                    print1Dvector(cellPot);
                    cout << icell << endl;
                    exit(EXIT_FAILURE);
                }
            }

            cellCheckedFlag[icell] = 1;
            icell += 1;
        }

        this->nCells  = nCells_temp;
//ghost>
        // no ghost cells for converting
        this->nCells_all = this->nCells;
//ghost<
        this->nPlanes = nPlanes_temp;

        cout << "numPlanes " << this->nPlanes << endl;

        // set planes around each node.
        geom_int i = 0;
        for (const auto &pl : planes)
        {
            for (const auto &iNode : pl.iNodes)
            {
                nodes[iNode].iPlanes.push_back(i);
            }
            i++;
        }

        cout << "make boundary" << endl;

        map<geom_int,vector<vector<geom_int>>> bcellMap_physTag_nodes;
        map<geom_int,vector<geom_int>> bcellMap_physTag_planes;

        // -------------------------------
        // *** make boundary conditions***
        // -------------------------------
        geom_int nBPlanes_temp = 0;
        // make bcellMap[phystag, nodes]
        for (const auto &eleOfEnt : elements_summary)
        {
            if (eleOfEnt.dimension != 2) continue; // skip volume. only surface

            geom_int physTag = surfEntMap[eleOfEnt.entTag].physTag; // physical tag (fluid, inlet, etc.)

            for (const auto &iEle : eleOfEnt.elements)
            {
                nBPlanes_temp += 1;
                bcellMap_physTag_nodes[physTag].push_back(iEle.iNodes);
            }
        }

        this->nBPlanes = nBPlanes_temp;
        this->nNormalPlanes = this->nPlanes - this->nBPlanes;
        this->nBconds = bcellMap_physTag_nodes.size();


        cout << " nBPlanes = " << this->nBPlanes << endl;
        cout << " nNormalPlanes = " << this->nNormalPlanes << endl;
        cout << " nBconds = " << this->nBconds << endl;

        // make bcellMap[phystag, planes]
        for (const auto &item : bcellMap_physTag_nodes) 
        {
            for (const auto &ns : item.second) 
            {
                //cout << "ns[0]" << ns[0]<< endl;
                vector<geom_int> planeCandidates = nodes[ns[0]].iPlanes;

                bool ifFound = false;
                for (const auto &plnCand : planeCandidates)
                {
                    vector<geom_int> nodesCand = planes[plnCand].iNodes;

                    if (ifEqualComponent(nodesCand, ns)) // found
                    {
                        ifFound = true;
                        if (bcellMap_physTag_planes.count(item.first)) // physTag exists
                        {
                            bcellMap_physTag_planes[item.first].push_back(plnCand);
                        } else { // physTag not exists
                            vector<geom_int> temp = { plnCand };
                            bcellMap_physTag_planes[item.first] = temp;
                        }
                    } else { // not found
                        continue;
                    }
                }
                if (ifFound == false)
                {
                    cerr << "Error: Coudn't find bplane" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        
        cout << "bcellMap Made" << endl;

        cout << "Change planes order. move boundary planes to last." << endl;

        // change planes order. move bplanes to last.
        // Firstly, prepare for making order of plane indices for reference.
        vector<geom_int> planeIndex;
        vector<geom_int> planeBCflg(this->nPlanes);
        vector<geom_int> bplaneIndex;

        fill(planeBCflg.begin(), planeBCflg.end(), 0);

        for (const auto &item : bcellMap_physTag_planes) 
        {
            for (const auto &pln : item.second) 
            {
                planeBCflg[pln] = 1;
                bplaneIndex.push_back(pln);
            }
        }

        geom_int ip = 0;
        for (const auto &bcflg : planeBCflg) 
        {
            if(bcflg == 0) {
                planeIndex.push_back(ip);
            }
            ip++;
        }

        planeIndex.insert(planeIndex.end(), bplaneIndex.begin(), bplaneIndex.end());

        cout << "plane index Made" << endl;

        vector<geom_int> planeIndex_old_to_new(this->nPlanes);

        vector<plane> planes_new(this->nPlanes);
        planes_new = planes;

        ip = 0;
        for (auto &i : planeIndex)
        {
            planes_new[ip].iNodes = this->planes[i].iNodes;
            planes_new[ip].iCells = this->planes[i].iCells;
            planeIndex_old_to_new[i] = ip;
            ip += 1;
        }
        this->planes = planes_new;
        planes_new.clear();

        vector<cell> cells_new(this->nCells);
        cells_new = cells;

        // cell
        i = 0;
        for (auto &cell : cells)
        {
            int j = 0;
            for (auto &pln : cell.iPlanes)
            {
                cells_new[i].iPlanes[j] = planeIndex_old_to_new[pln];
                j += 1;
            }
            i += 1;
        }

        this->cells = cells_new;
        cells_new.clear();

        vector<node> nodes_new(this->nNodes);
        nodes_new = nodes;

        // node
        i = 0;
        for (auto &node : nodes)
        {
            geom_int j = 0;
            for (auto &pln : node.iPlanes)
            {
                nodes_new[i].iPlanes[j] = planeIndex_old_to_new[pln];
                j += 1;
            }
            i += 1;
        }

        this->nodes = nodes_new;
        nodes_new.clear();

        // boundary cells, bcond
        map<geom_int,vector<geom_int>> bcellMap_physTag_planes_new;
        bcellMap_physTag_planes_new = bcellMap_physTag_planes;

        for (const auto &item : bcellMap_physTag_planes) 
        {
            geom_int j = 0;
            for (const auto &pln : item.second) 
            {
                bcellMap_physTag_planes_new[item.first][j] = planeIndex_old_to_new[pln];
                j += 1;
            }
        }


        bcellMap_physTag_planes = bcellMap_physTag_planes_new;
        bcellMap_physTag_planes_new.clear();

        for (auto &item : bcellMap_physTag_planes)
        {
            vector<geom_int> iCells_temp;
            vector<geom_int> iBPlanes_temp;

            for (auto &pln : item.second)
            {
                iCells_temp.push_back(planes[pln].iCells[0]);
                iBPlanes_temp.push_back(pln - this->nNormalPlanes);
            }

            bcond bcond_temp = bcond(item.first, item.second, iCells_temp, iBPlanes_temp);

            bconds.push_back(bcond_temp);
        }

        // iPlane direction for surface of cells
        i = 0;
        for (auto &pln : planes)
        {
            geom_int c0 = pln.iCells[0];
            int wantedIndex = findIndex(cells[c0].iPlanes , i);
            cells[c0].iPlanesDir[wantedIndex] = 1;

            if (pln.iCells.size()>1)
            {
                geom_int c1 = pln.iCells[1];
                wantedIndex = findIndex(cells[c1].iPlanes , i);
                cells[c1].iPlanesDir[wantedIndex] = -1;
            }

            i += 1;
        }

        cout << "calculate surface vector & area & center\n";

        // ------------------------------------------------
        // *** calculate surface vector & area & center ***
        // ------------------------------------------------
        for (auto& pln : planes)
        {
            if (pln.iNodes.size() == 3) // triangle
            {
                geom_int n0 = pln.iNodes[0];
                geom_int n1 = pln.iNodes[1];
                geom_int n2 = pln.iNodes[2];

                geom_float r01x = nodes[n1].coords[0] - nodes[n0].coords[0];
                geom_float r01y = nodes[n1].coords[1] - nodes[n0].coords[1];
                geom_float r01z = nodes[n1].coords[2] - nodes[n0].coords[2];

                geom_float r02x = nodes[n2].coords[0] - nodes[n0].coords[0];
                geom_float r02y = nodes[n2].coords[1] - nodes[n0].coords[1];
                geom_float r02z = nodes[n2].coords[2] - nodes[n0].coords[2];

                pln.surfVect.resize(3);

                pln.surfVect[0] = -0.5*(r01y*r02z -r01z*r02y);
                pln.surfVect[1] = -0.5*(r01z*r02x -r01x*r02z);
                pln.surfVect[2] = -0.5*(r01x*r02y -r01y*r02x);
                pln.surfArea = std::sqrt(  std::pow(pln.surfVect[0] , 2.0) 
                                         + std::pow(pln.surfVect[1] , 2.0)
                                         + std::pow(pln.surfVect[2] , 2.0) );

            } else if (pln.iNodes.size() == 4) { // quad 

                geom_int n0 = pln.iNodes[0];
                geom_int n1 = pln.iNodes[1];
                geom_int n2 = pln.iNodes[2];
                geom_int n3 = pln.iNodes[3];

                geom_float r02x = nodes[n2].coords[0] - nodes[n0].coords[0];
                geom_float r02y = nodes[n2].coords[1] - nodes[n0].coords[1];
                geom_float r02z = nodes[n2].coords[2] - nodes[n0].coords[2];

                geom_float r13x = nodes[n1].coords[0] - nodes[n3].coords[0];
                geom_float r13y = nodes[n1].coords[1] - nodes[n3].coords[1];
                geom_float r13z = nodes[n1].coords[2] - nodes[n3].coords[2];

                pln.surfVect.resize(3);

                pln.surfVect[0] = -0.5*(r02y*r13z -r02z*r13y);
                pln.surfVect[1] = -0.5*(r02z*r13x -r02x*r13z);
                pln.surfVect[2] = -0.5*(r02x*r13y -r02y*r13x);
                pln.surfArea = std::sqrt(  std::pow(pln.surfVect[0] , 2.0) 
                                         + std::pow(pln.surfVect[1] , 2.0)
                                         + std::pow(pln.surfVect[2] , 2.0) );
            }
            // set face center
            pln.centCoords.resize(3);

            pln.centCoords[0] = 0.0;
            pln.centCoords[1] = 0.0;
            pln.centCoords[2] = 0.0;
            for (auto &n : pln.iNodes)
            {
                pln.centCoords[0] += nodes[n].coords[0];
                pln.centCoords[1] += nodes[n].coords[1];
                pln.centCoords[2] += nodes[n].coords[2];
            }
            pln.centCoords[0] = pln.centCoords[0]/pln.iNodes.size();
            pln.centCoords[1] = pln.centCoords[1]/pln.iNodes.size();
            pln.centCoords[2] = pln.centCoords[2]/pln.iNodes.size();
        }

        cout << "calculate volume" << endl;
        // ------------------------
        // *** calculate volume ***
        // ------------------------
        for (auto &icell : cells)
        {
            geom_float volume = 0.0 ;

            elementTypeFormat eleType;
            eleType = this->eleTypeMap.mapElementFromGmshID[icell.ieleType];

            if (eleType.name == "tet") // tetra
            {
                icell.volume = tetraVolume(icell);
            }

            if (eleType.name == "hex") // hexahedral
            {
                geom_int n0 = icell.iNodes[0];
                geom_int n1 = icell.iNodes[1];
                geom_int n2 = icell.iNodes[2];
                geom_int n3 = icell.iNodes[3];
                geom_int n4 = icell.iNodes[4];
                geom_int n5 = icell.iNodes[5];
                geom_int n6 = icell.iNodes[6];
                geom_int n7 = icell.iNodes[7];

                vector<geom_int> iNodes_temp1 = {n0, n1, n4, n2};
                volume += tetraVolume(cell(iNodes_temp1));

                vector<geom_int> iNodes_temp2 = {n2, n3, n7, n0};
                volume += tetraVolume(cell(iNodes_temp2));

                vector<geom_int> iNodes_temp3 = {n0, n4, n7, n2};
                volume += tetraVolume(cell(iNodes_temp3));

                vector<geom_int> iNodes_temp4 = {n2, n6, n7, n5};
                volume += tetraVolume(cell(iNodes_temp4));

                vector<geom_int> iNodes_temp5 = {n1, n2, n4, n5};
                volume += tetraVolume(cell(iNodes_temp5));

                vector<geom_int> iNodes_temp6 = {n2, n4, n7, n5};
                volume += tetraVolume(cell(iNodes_temp6));

                icell.volume = volume;
            }

            if (eleType.name == "prism") 
            {
                geom_int n0 = icell.iNodes[0];
                geom_int n1 = icell.iNodes[1];
                geom_int n2 = icell.iNodes[2];
                geom_int n3 = icell.iNodes[3];
                geom_int n4 = icell.iNodes[4];
                geom_int n5 = icell.iNodes[5];

                vector<geom_int> iNodes_temp1 = {n0, n4, n3, n5};
                volume += tetraVolume(cell(iNodes_temp1));

                vector<geom_int> iNodes_temp2 = {n0, n1, n4, n5};
                volume += tetraVolume(cell(iNodes_temp2));

                vector<geom_int> iNodes_temp3 = {n0, n2, n1, n5};
                volume += tetraVolume(cell(iNodes_temp3));

                icell.volume = volume;
            }

            if (eleType.name == "pyramid") 
            {
                geom_int n0 = icell.iNodes[0];
                geom_int n1 = icell.iNodes[1];
                geom_int n2 = icell.iNodes[2];
                geom_int n3 = icell.iNodes[3];
                geom_int n4 = icell.iNodes[4];

                vector<geom_int> iNodes_temp1 = {n0, n2, n1, n4};
                volume += tetraVolume(cell(iNodes_temp1));

                vector<geom_int> iNodes_temp2 = {n0, n3, n2, n4};
                volume += tetraVolume(cell(iNodes_temp2));

                icell.volume = volume;
            }

            // set volume center
            icell.centCoords.resize(3);

            icell.centCoords[0] = 0.0;
            icell.centCoords[1] = 0.0;
            icell.centCoords[2] = 0.0;
            for (auto &n : icell.iNodes)
            {
                icell.centCoords[0] += nodes[n].coords[0];
                icell.centCoords[1] += nodes[n].coords[1];
                icell.centCoords[2] += nodes[n].coords[2];
            }
            icell.centCoords[0] = icell.centCoords[0]/icell.iNodes.size();
            icell.centCoords[1] = icell.centCoords[1]/icell.iNodes.size();
            icell.centCoords[2] = icell.centCoords[2]/icell.iNodes.size();
        }

    }

    geom_float tetraVolume(const cell &tetra)
    {
        node n0 = nodes[tetra.iNodes[0]];
        node n1 = nodes[tetra.iNodes[1]];
        node n2 = nodes[tetra.iNodes[2]];
        node n3 = nodes[tetra.iNodes[3]];

        geom_float ax = n2.coords[0] - n0.coords[0]; 
        geom_float ay = n2.coords[1] - n0.coords[1];
        geom_float az = n2.coords[2] - n0.coords[2];

        geom_float bx = n1.coords[0] - n0.coords[0]; 
        geom_float by = n1.coords[1] - n0.coords[1];
        geom_float bz = n1.coords[2] - n0.coords[2];

        geom_float cx = n3.coords[0] - n0.coords[0]; 
        geom_float cy = n3.coords[1] - n0.coords[1];
        geom_float cz = n3.coords[2] - n0.coords[2];

        geom_float volume = ((ay*bz -az*by)*cx +(az*bx -ax*bz)*cy +(ax*by -ay*bx)*cz)/6.0;
        volume = abs(volume);

        if (volume < 1.0e-20) cout << "WARNING : Too small volume\n" ;

        return volume;

    }

    void writeInputH5(const string outFileName ,variables var)
    {
        // ------------
        // *** HDF5 *** 
        // ------------
        string baseName;

        if (outFileName.substr(outFileName.size()-3, 3) == ".h5") {
            baseName = outFileName.substr(0, outFileName.size()-3);
            cout << "baseName=" << baseName << endl;
        } else if (outFileName.substr(outFileName.size()-5, 5) == ".hdf5") {
            baseName = outFileName.substr(0, outFileName.size()-5);
            cout << "baseName=" << baseName << endl;
        } else {
            cerr << "write hdf5 name is wrong. " << outFileName << endl;
            exit(EXIT_FAILURE);
        }

        File file(outFileName, File::ReadWrite | File::Truncate);

        // write mesh structure
        vector<geom_float> COORD;
        for (auto& nod : this->nodes)
        {
            COORD.push_back(nod.coords[0]);
            COORD.push_back(nod.coords[1]);
            COORD.push_back(nod.coords[2]);
        }

        Group group = file.createGroup("/MESH");
        Attribute a = group.createAttribute<geom_int>( "nCells", DataSpace::From(this->cells.size()));
        a.write(this->cells.size());

        a = group.createAttribute<geom_int>( "nNodes", DataSpace::From(this->nodes.size()));
        a.write(this->nodes.size());

        a = group.createAttribute<geom_int>( "nPlanes", DataSpace::From(this->planes.size()));
        a.write(this->planes.size());

        a = group.createAttribute<geom_int>( "nNormalPlanes", DataSpace::From(this->nNormalPlanes));
        a.write(this->nNormalPlanes);

        a = group.createAttribute<geom_int>( "nBPlanes", DataSpace::From(this->nBPlanes));
        a.write(this->nBPlanes);

        a = group.createAttribute<geom_int>( "nBconds", DataSpace::From(this->nBconds));
        a.write(this->nBconds);


        file.createDataSet("/MESH/COORD",COORD);

        vector<geom_int> CONNE;
        geom_int CONNE_dim = 0;
        geom_int CONNE0;    
        for (auto& cel : this->cells)
        {
            geom_int nn = this->eleTypeMap.mapElementFromGmshID[cel.ieleType].nNodes;
            string name = this->eleTypeMap.mapElementFromGmshID[cel.ieleType].name;

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

        // write planes
        vector<geom_int> planes_struct;
        vector<geom_float> surfVect;
        vector<geom_float> surfArea;
        vector<geom_float> centCoords;


        for (auto& pln : this->planes)
        {
            // add nodes
            geom_int nn = pln.iNodes.size();
            planes_struct.push_back(nn);
            for (geom_int in=0 ; in < nn; in++)
            {
                planes_struct.push_back(pln.iNodes[in]);
            }

            nn = pln.iCells.size();
            planes_struct.push_back(nn);
            for (geom_int in=0 ; in < nn; in++)
            {
                planes_struct.push_back(pln.iCells[in]);
            }
            surfVect.push_back(pln.surfVect[0]);
            surfVect.push_back(pln.surfVect[1]);
            surfVect.push_back(pln.surfVect[2]);
            surfArea.push_back(pln.surfArea);
            centCoords.push_back(pln.centCoords[0]);
            centCoords.push_back(pln.centCoords[1]);
            centCoords.push_back(pln.centCoords[2]);
        }
        file.createDataSet("/PLANES/STRUCT",planes_struct);
        file.createDataSet("/PLANES/surfVect",surfVect);
        file.createDataSet("/PLANES/surfArea",surfArea);
        file.createDataSet("/PLANES/centCoords",centCoords);

        // write cells
        vector<geom_int> cells_struct;
        vector<geom_float> volume;
        vector<geom_float> centCoords2;

        for (auto& cel: this->cells)
        {
            // add nodes
            geom_int nn = cel.iNodes.size();
            cells_struct.push_back(nn);
            for (geom_int in=0 ; in < nn; in++)
            {
                cells_struct.push_back(cel.iNodes[in]);
            }

            nn = cel.iPlanes.size();
            cells_struct.push_back(nn);
            for (geom_int in=0 ; in < nn; in++)
            {
                cells_struct.push_back(cel.iPlanes[in]);
            }

            cells_struct.push_back(nn);
            for (geom_int in=0 ; in < nn; in++)
            {
                cells_struct.push_back(cel.iPlanesDir[in]);
            }
            cells_struct.push_back(cel.ieleType);

            volume.push_back(cel.volume);
            centCoords2.push_back(cel.centCoords[0]);
            centCoords2.push_back(cel.centCoords[1]);
            centCoords2.push_back(cel.centCoords[2]);
        }
        file.createDataSet("/CELLS/STRUCT",cells_struct);
        file.createDataSet("/CELLS/volume",volume);
        file.createDataSet("/CELLS/centCoords",centCoords2);

        // write initial values
        //for (string name : var.output_cellValNames)
        //TODO: ややこしくしているのでシンプルにoutput_cellValNamesで回したい。が、なぜかエラーになる
        for (auto& v : var.c)
        {
            string name = v.first;

            auto itr = std::find(var.read_cellValNames.begin(), var.read_cellValNames.end(), name);
            if (itr == var.read_cellValNames.end()) {
                continue; // notfound
            }

            file.createDataSet("/VALUE/"+name , v.second);
        }

        // boundary conditions
        for (auto& bc : this->bconds)
        {
            ostringstream oss;
            oss << bc.physID;


            Group group = file.createGroup("/BCONDS/"+oss.str());

            group.createDataSet("/BCONDS/"+oss.str()+"/iPlanes", bc.iPlanes);
            group.createDataSet("/BCONDS/"+oss.str()+"/iBPlanes", bc.iBPlanes);
            group.createDataSet("/BCONDS/"+oss.str()+"/iCells", bc.iCells);

            string ATTRIBUTE_NAME_NOTE("bcondKind");
            string string_list(bc.bcondKind);
            //cout << string_list << endl;
            Attribute a = group.createAttribute<std::string>(ATTRIBUTE_NAME_NOTE, 
                                                             DataSpace::From(string_list));
            a.write(string_list);

            for (auto& bvar : bc.bvar)
            {
                string name = bvar.first;


                file.createDataSet("/BCONDS/"+oss.str()+"/VALUE/"+name , bvar.second);
            }
        }

        // ------------
        // *** XDMF ***
        // ------------
        string fnameXDMF = baseName + ".xmf";
        ofstream ofs(fnameXDMF);

        ofs << "<?xml version='1.0' ?>\n";
        ofs << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n";
        ofs << "<Xdmf>\n";
        ofs << "  <Domain>\n";
        ofs << "    <Grid  GridType='Collection' CollectionType='Spatial' Name='Mixed'>\n";
        ofs << "      <Grid Name='aaa'>\n";
        ofs << "        <Topology Type='Mixed' NumberOfElements='" << this->cells.size() << "'>\n";
        ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << CONNE_dim << "'>\n";
//        ofs << "            vol_" << oss.str() <<".h5:MESH/CONNE\n";
        ofs << "           " + outFileName + ":MESH/CONNE\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Topology>\n";

        ofs << "        <Geometry Type='XYZ'>\n";
        ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << this->nodes.size()*3 << "'>\n";
        //ofs << "            vol_"<< oss.str() <<".h5:MESH/COORD\n";
        ofs << "           " + outFileName + ":MESH/COORD\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Geometry>\n";

        for (string name : var.read_cellValNames)
        {
        //for (auto& v : var.c) {
            //string name = v.first;
            ofs << "        <Attribute Name='"  << name << "' Center='Cell' >\n";
            ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << this->cells.size() << "'>\n";
            ofs << "            " + outFileName + ":VALUE/" << name<< "\n";
            ofs << "          </DataItem>\n";
            ofs << "        </Attribute>\n";
        }

        ofs << "      </Grid>\n";
        ofs << "    </Grid>\n";
        ofs << "  </Domain>\n";
        ofs << "</Xdmf>\n";


    }
};
