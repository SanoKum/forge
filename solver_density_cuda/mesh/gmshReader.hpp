#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <array>
#include <map>

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
    bool is3D;

    // ---- median-dual (node-centered) 前処理の生成物 ----
    // buildMedianDual() で構築。空なら writeInputH5 は /DUAL を書かない (cell モード)。
    // CV = ノード、双対面 = primal エッジ (2D では plane と 1:1)、CV 重心 = ノード座標。
    bool dualBuilt = false;
    bool axisCentroidShift = true;  // true: CV 面積加重重心を cell 中心に (軸 R=0 のゼロ体積回避)。
                                    // false: node 座標 (軸 R=0)。軸ソース OFF と併用前提。convertGmshToForge が cfg から設定。
    std::vector<geom_float> dualVolume;     // [nNodes] 各ノードの双対 CV 体積 (2D は面積×単位厚み)
    std::vector<geom_float> dualCentroid;   // [nNodes*3] 双対 CV の面積加重重心 (FV セル中心)。
                                            // 軸対称で重要: 軸上ノード(R=0)でも CV 重心 R>0 となり
                                            // 回転体積 volume=A_planar*r が非ゼロになる (theory §3.4)。
    std::vector<geom_int>   dualFaceCells;  // [nDualFaces*2] 双対面が繋ぐ 2 ノード {n0,n1}
    std::vector<geom_float> dualFaceVect;   // [nDualFaces*3] 双対面ベクトル (n0->n1 向き)
    std::vector<geom_float> dualFaceArea;   // [nDualFaces]   |dualFaceVect|
    std::vector<geom_float> dualFaceCent;   // [nDualFaces*3] 双対面重心 (面積加重)
    // 境界半割双対面: 境界ノードに割り当てる外向き面ベクトル (bcond ごと)。
    std::vector<geom_int>   dualBnodeId;    // [nBHalf]   境界ノード id (bcond 境界順に連結)
    std::vector<geom_float> dualBnodeVect;  // [nBHalf*3] 外向き半割面ベクトル
    std::vector<geom_float> dualBnodeCent;  // [nBHalf*3] 半割面の面積加重重心 (境界上, R≥0)。
                                            // 軸対称 r 重みを正しくし、入口/出口 BC が軸近傍 corner CV に届くようにする。
    std::vector<geom_int>   dualBcondOffset;// [nBconds+1] bcond ごとの dualBnode* (emit 済=非壁) への CSR オフセット
    std::vector<geom_int>   dualBcondPhysID;// [nBconds]   各 bcond の physID (順序対応)
    // [nBconds][各 bcond が所有する境界ノード] (壁優先所有後)。iCells 設定用。壁所有 bcond は plane を emit
    // しない (ghost 不要・Dirichlet) が、wall_flag 構築のためノード列はここに残し h5 の iCells に書く。
    std::vector<std::vector<geom_int>> dualBcondNodes;

    class lineEnt
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


        lineEnt(){};
        lineEnt(geom_int tag, geom_float minX, geom_float minY, geom_float minZ,
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

    geom_int nLineEnt;
    geom_int nSurfEnt;
    geom_int nVolumeEnt;
    map<geom_int,lineEnt> lineEntMap;
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
        this->is3D = false;
        for (auto d : physDim) {
            if (d == 3) { 
                is3D = true;
                break;
            }
        }
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
        this->nLineEnt   = stoi(l_str[1]);
        this->nSurfEnt   = stoi(l_str[2]);
        this->nVolumeEnt = stoi(l_str[3]);

        // Point Entities : just ignore
        for (geom_int i = 0 ; i<stoi(l_str[0]); i++)
        {
            getline(inputFile, line);
        }

        // Curve Entities 
        //for (geom_int i = 0 ; i<stoi(l_str[1]); i++)
        for (geom_int i = 0 ; i<nLineEnt; i++)
        {
            getline(inputFile, line);
            boost::split(l_str, line, boost::is_space());
            lineEnt lineEnt_temp = lineEnt(stoi(l_str[0]), stof(l_str[1]), stof(l_str[2]), stof(l_str[3]),
                                           stof(l_str[4]), stof(l_str[5]), stof(l_str[6]),
                                           stoi(l_str[7]), stoi(l_str[8]));

            if (lineEnt_temp.nPhysTag > 1)
            {
                std::cerr << "Error: unknown gmsh format. number of PhysTag is not 1" << std::endl;                                                                                                     
                std::cerr << "nPhysTag = " << lineEnt_temp.nPhysTag << " " << lineEnt_temp.physTag  << std::endl;                                                                                                     
                std::cerr << "Ent id = " << i   << std::endl;                                                                                                     

                exit(EXIT_FAILURE);
            }

            lineEntMap.insert(std::make_pair(lineEnt_temp.entTag, lineEnt_temp));

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
        cout << line << endl;
        getline(inputFile, line);
        cout << line << endl;

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

        // 大規模メッシュで vector の倍々再確保によるメモリスパイク (OOM 要因) を避けるため、
        // cells を事前 reserve する (セル数は対象次元の要素数の総和)。
        {
            geom_int nCellEst = 0;
            for (auto &itEnt : elements_summary)
            {
                if (is3D == true  && itEnt.dimension != 3) continue;
                if (is3D == false && itEnt.dimension != 2) continue;
                nCellEst += itEnt.nEleEnt;
            }
            this->cells.reserve(nCellEst);
        }

        // iterate over each Entity.
        for (auto &itEnt : elements_summary)
        {
            if (is3D == true  && itEnt.dimension != 3) continue; // skip surface. only volume
            if (is3D == false && itEnt.dimension != 2) continue; // skip edge. only surface

            geom_int physTag = is3D ? volumeEntMap[itEnt.entTag].physTag
                                    : surfEntMap[itEnt.entTag].physTag; // physical tag (fluid, outer, etc.)
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
                cell_now.regionId = physTag;
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
        // planes も倍々再確保を避けるため事前 reserve (hex で ~3.0/cell, tet で ~2.0/cell;
        // 余裕をみて 3.5*nCells を確保。過不足は正しさに影響せずスパイク回避のみが目的)。
        this->planes.reserve(this->cells.size() * 3 + this->cells.size() / 2);

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
            //if (is3D == true  && eleOfEnt.dimension != 2) continue; // skip volume. only surface
            geom_int physTag;

            if (is3D == true) {
                if (eleOfEnt.dimension !=2) continue;
                physTag = surfEntMap[eleOfEnt.entTag].physTag; // physical tag (fluid, inlet, etc.)
            } else { // is3D == false
                if (eleOfEnt.dimension !=1) continue;
                physTag = lineEntMap[eleOfEnt.entTag].physTag; // physical tag (fluid, inlet, etc.)
                cout << "physTag = " << physTag << endl;
            }

            for (const auto &iEle : eleOfEnt.elements)
            {
                nBPlanes_temp += 1;
                bcellMap_physTag_nodes[physTag].push_back(iEle.iNodes);
            }
        }

        // 境界ノードマップを作り終えたら elements_summary は以降不要。大規模メッシュでは
        // 全要素のノードリストを保持しており大きいので解放する (writeInputH5 まで保持しない)。
        std::list<elementsOfEntity>().swap(this->elements_summary);

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

        // 並べ替え: planeIndex は全 plane index の置換 (各 index がちょうど 1 回)。
        // 旧実装は planes を 2 重・3 重に deep copy しメモリスパイク (OOM) の主因だった。
        // ここでは sub-vector を move して付け替える (この時点で surfVect/centCoords は
        // 未計算=空なので iNodes/iCells のみ移送すれば挙動は同一)。
        vector<plane> planes_new(this->nPlanes);

        ip = 0;
        for (auto &i : planeIndex)
        {
            planes_new[ip].iNodes = std::move(this->planes[i].iNodes);
            planes_new[ip].iCells = std::move(this->planes[i].iCells);
            planeIndex_old_to_new[i] = ip;
            ip += 1;
        }
        this->planes = std::move(planes_new);

        // cell: iPlanes は「並べ替え」ではなく plane index 値の付け替えのみ。
        // 旧実装は cells を丸ごと deep copy していたが不要。in-place で remap する。
        for (auto &cell : this->cells)
        {
            for (auto &pln : cell.iPlanes)
            {
                pln = planeIndex_old_to_new[pln];
            }
        }

        // node: 同上 (iPlanes の値のみ remap)。in-place 化で deep copy を排除。
        for (auto &node : this->nodes)
        {
            for (auto &pln : node.iPlanes)
            {
                pln = planeIndex_old_to_new[pln];
            }
        }

        // boundary cells, bcond: 同上 (各 physTag の plane index を in-place remap)。
        for (auto &item : bcellMap_physTag_planes)
        {
            for (auto &pln : item.second)
            {
                pln = planeIndex_old_to_new[pln];
            }
        }

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

            } else if (pln.iNodes.size() == 2) { // line (2D)

                geom_int n0 = pln.iNodes[0];
                geom_int n1 = pln.iNodes[1];

                geom_float r01x = nodes[n1].coords[0] - nodes[n0].coords[0];
                geom_float r01y = nodes[n1].coords[1] - nodes[n0].coords[1];
                geom_float r01z = nodes[n1].coords[2] - nodes[n0].coords[2];

                pln.surfVect.resize(3);

                // For CCW-wound 2D cells (Gmsh default), edges traverse the cell
                // counter-clockwise; rotating the edge vector by -90 deg gives the
                // outward normal (consistent with the c0->c1 sv convention used by
                // the flux/BC kernels).
                pln.surfVect[0] =+r01y;
                pln.surfVect[1] =-r01x;
                pln.surfVect[2] = 0.0;
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

            if (eleType.name == "quad") 
            {
                geom_int n0 = icell.iNodes[0];
                geom_int n1 = icell.iNodes[1];
                geom_int n2 = icell.iNodes[2];
                geom_int n3 = icell.iNodes[3];

                geom_float r02x = nodes[n2].coords[0] - nodes[n0].coords[0];
                geom_float r02y = nodes[n2].coords[1] - nodes[n0].coords[1];
                geom_float r02z = nodes[n2].coords[2] - nodes[n0].coords[2];

                geom_float r13x = nodes[n1].coords[0] - nodes[n3].coords[0];
                geom_float r13y = nodes[n1].coords[1] - nodes[n3].coords[1];
                geom_float r13z = nodes[n1].coords[2] - nodes[n3].coords[2];


                geom_float sv1 = -0.5*(          -r02z*r13y);
                geom_float sv2 = -0.5*(r02z*r13x           );
                geom_float sv3 = -0.5*(r02x*r13y -r02y*r13x);
                geom_float ss  = std::sqrt(  std::pow(sv1 , 2.0) 
                                           + std::pow(sv2 , 2.0)
                                           + std::pow(sv3 , 2.0) );

                icell.volume = ss;
            }

            if (eleType.name == "triangle") 
            {
                geom_int n0 = icell.iNodes[0];
                geom_int n1 = icell.iNodes[1];
                geom_int n2 = icell.iNodes[2];

                geom_float r01x = nodes[n1].coords[0] - nodes[n0].coords[0];
                geom_float r01y = nodes[n1].coords[1] - nodes[n0].coords[1];
                geom_float r01z = nodes[n1].coords[2] - nodes[n0].coords[2];

                geom_float r02x = nodes[n2].coords[0] - nodes[n0].coords[0];
                geom_float r02y = nodes[n2].coords[1] - nodes[n0].coords[1];
                geom_float r02z = nodes[n2].coords[2] - nodes[n0].coords[2];

                geom_float sv1 = -0.5*(r01y*r02z -r01z*r02y);
                geom_float sv2 = -0.5*(r01z*r02x -r01x*r02z);
                geom_float sv3 = -0.5*(r01x*r02y -r01y*r02x);
                geom_float ss  = std::sqrt(  std::pow(sv1 , 2.0) 
                                           + std::pow(sv2 , 2.0)
                                           + std::pow(sv3 , 2.0) );

                icell.volume = ss;
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

        // 面法線 (surfVect) の向きをそろえる。
        // - 内部面: cells[ic0] -> cells[ic1] 方向 (D = (cc1-cc0)・sv > 0) にする。
        // - 境界面: 外向き (Db = (faceCenter - cc0)・sv > 0) にする。
        //
        // 3D メッシュ (hex/prism/tet) では純粋な幾何判定を用いる。XY-shoelace による
        // CW 判定はセル断面が XY 平面に乗る Z 押し出し以外の 3D セルで破綻し、
        // 境界面・内部面の法線を誤って反転させて偽の勾配を生むため使わない。
        // 2D メッシュは従来どおり CW(shoelace) ベースの判定を維持する (回帰防止)。
        if (is3D)
        {
            int nFixed = 0;
            int nFixedBnd = 0;
            for (geom_int i = 0; i < (geom_int)planes.size(); ++i)
            {
                auto& pln = planes[i];
                geom_int ic0 = pln.iCells[0];
                const auto& cc0 = cells[ic0].centCoords;

                if (pln.iCells.size() == 1)
                {
                    // 境界面: 面重心がセル重心の外側に来る向きにそろえる
                    geom_float fx = 0.0, fy = 0.0, fz = 0.0;
                    for (geom_int nn : pln.iNodes)
                    {
                        fx += nodes[nn].coords[0];
                        fy += nodes[nn].coords[1];
                        fz += nodes[nn].coords[2];
                    }
                    const geom_float inv = 1.0 / (geom_float)pln.iNodes.size();
                    fx *= inv; fy *= inv; fz *= inv;
                    const geom_float Db = (fx - cc0[0])*pln.surfVect[0]
                                        + (fy - cc0[1])*pln.surfVect[1]
                                        + (fz - cc0[2])*pln.surfVect[2];
                    if (Db < 0.0)
                    {
                        pln.surfVect[0] = -pln.surfVect[0];
                        pln.surfVect[1] = -pln.surfVect[1];
                        pln.surfVect[2] = -pln.surfVect[2];
                        nFixedBnd += 1;
                    }
                    continue;
                }

                // 内部面: ic0 -> ic1 方向にそろえる
                const auto& cc1 = cells[pln.iCells[1]].centCoords;
                const geom_float D = (cc1[0]-cc0[0])*pln.surfVect[0]
                                   + (cc1[1]-cc0[1])*pln.surfVect[1]
                                   + (cc1[2]-cc0[2])*pln.surfVect[2];
                if (D < 0.0)
                {
                    pln.surfVect[0] = -pln.surfVect[0];
                    pln.surfVect[1] = -pln.surfVect[1];
                    pln.surfVect[2] = -pln.surfVect[2];
                    nFixed += 1;
                }
            }
            cout << "face orientation fix (3D geometric): internal=" << nFixed
                 << ", boundary=" << nFixedBnd << "\n";
        }
        else
        {
            // --- 2D: 従来の CW(shoelace) ベース判定 ---
            std::vector<bool> isCW(cells.size(), false);
            for (geom_int ic = 0; ic < (geom_int)cells.size(); ++ic)
            {
                const auto& cel = cells[ic];
                geom_int nn = (geom_int)cel.iNodes.size();
                if (nn < 3) continue; // 辺要素はスキップ
                geom_float signedArea = 0.0;
                for (geom_int k = 0; k < nn; ++k)
                {
                    const auto& n0 = nodes[cel.iNodes[k]].coords;
                    const auto& n1 = nodes[cel.iNodes[(k+1)%nn]].coords;
                    signedArea += n0[0]*n1[1] - n1[0]*n0[1];
                }
                isCW[ic] = (signedArea < 0.0);
            }
            int nCW = 0;
            for (bool b : isCW) if (b) nCW++;
            cout << "CW cells detected: " << nCW << "\n";

            int nFixed = 0;
            int nFixedBnd = 0;
            int nSkipped = 0; // CW でないが D<0 の内部面（高スキュー警告用）
            for (geom_int i = 0; i < (geom_int)planes.size(); ++i)
            {
                auto& pln = planes[i];
                geom_int ic0 = pln.iCells[0];

                if (pln.iCells.size() == 1)
                {
                    // 境界面: ic0 が CW なら sv が内向きなので反転する
                    if (isCW[ic0])
                    {
                        pln.surfVect[0] = -pln.surfVect[0];
                        pln.surfVect[1] = -pln.surfVect[1];
                        pln.surfVect[2] = -pln.surfVect[2];
                        nFixedBnd += 1;
                    }
                    continue;
                }

                // 内部面
                geom_int ic1 = pln.iCells[1];
                if (isCW[ic0])
                {
                    pln.surfVect[0] = -pln.surfVect[0];
                    pln.surfVect[1] = -pln.surfVect[1];
                    pln.surfVect[2] = -pln.surfVect[2];
                    nFixed += 1;
                }
                else
                {
                    const auto& cc0 = cells[ic0].centCoords;
                    const auto& cc1 = cells[ic1].centCoords;
                    const geom_float D = (cc1[0]-cc0[0])*pln.surfVect[0]
                                       + (cc1[1]-cc0[1])*pln.surfVect[1]
                                       + (cc1[2]-cc0[2])*pln.surfVect[2];
                    if (D < 0.0) nSkipped += 1;
                }
            }
            cout << "face orientation fix (CW ic0): internal=" << nFixed << ", boundary=" << nFixedBnd << "\n";
            if (nSkipped > 0)
                cout << "WARNING: " << nSkipped << " internal faces have D<0 but ic0 is CCW (skewed mesh)\n";
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

    // -------------------------------------------------------------------------
    // 中点双対 (median-dual) メッシュを primal から構築する (node-centered 用)。
    //
    // 現状 2D (triangle / quad) のみ実装。3D は警告して何もしない (M4 で対応)。
    // theory: docs/discretization/theory.md §3、impl: docs/discretization/implementation.md §2。
    //
    // 生成物 (メンバ): dualVolume / dualFace* / dualBnode* / dualBcond*。
    // 整合チェック (前処理必須, theory §3.3): ① Σ双対体積 == Σprimal体積、
    // ② 各ノードの符号付き双対面ベクトル和 + 境界半割面 = 0 (閉性)、
    // ③ Σ半割面積 == primal 境界面積。逸脱が大きければ exit。
    // -------------------------------------------------------------------------
    void buildMedianDual()
    {
        if (is3D)
        {
            cout << "[buildMedianDual] 3D mesh: median-dual generation not implemented yet "
                    "(node-centered M4). Skipping /DUAL output.\n";
            dualBuilt = false;
            return;
        }

        cout << "[buildMedianDual] building 2D median-dual mesh ...\n";

        const geom_int nN = this->nNodes;
        const geom_int nP = this->nPlanes; // 2D では plane = primal エッジ (一意)

        dualVolume.assign(nN, 0.0);
        dualCentroid.assign(nN * 3, 0.0);
        dualFaceCells.assign(nP * 2, -1);
        dualFaceVect.assign(nP * 3, 0.0);
        dualFaceArea.assign(nP, 0.0);
        dualFaceCent.assign(nP * 3, 0.0);

        // ---- 双対面: primal エッジ ip ごとに 1 枚。隣接セルの (midpoint -> centroid)
        //      区間を集約。法線は n0->n1 向きに統一 (cell-centered 規約と同型)。
        for (geom_int ip = 0; ip < nP; ++ip)
        {
            const geom_int A = this->planes[ip].iNodes[0];
            const geom_int B = this->planes[ip].iNodes[1];
            dualFaceCells[2*ip + 0] = A;
            dualFaceCells[2*ip + 1] = B;

            const geom_float Mx = 0.5*(nodes[A].coords[0] + nodes[B].coords[0]);
            const geom_float My = 0.5*(nodes[A].coords[1] + nodes[B].coords[1]);

            const geom_float eABx = nodes[B].coords[0] - nodes[A].coords[0];
            const geom_float eABy = nodes[B].coords[1] - nodes[A].coords[1];

            geom_float vx = 0.0, vy = 0.0;       // 集約面ベクトル
            geom_float cx = 0.0, cy = 0.0, wsum = 0.0; // 面積加重重心

            for (const geom_int ic : this->planes[ip].iCells)
            {
                const geom_float Gx = cells[ic].centCoords[0];
                const geom_float Gy = cells[ic].centCoords[1];

                // 区間 midpoint(M) -> centroid(G)。単位厚みの面ベクトルは rotate(-90)。
                const geom_float sx = Gx - Mx;
                const geom_float sy = Gy - My;
                geom_float nx = sy;     // rotate(-90): (sx,sy) -> (sy,-sx)
                geom_float ny = -sx;
                // n0->n1 (A->B) 向きにそろえる
                if (nx*eABx + ny*eABy < 0.0) { nx = -nx; ny = -ny; }

                vx += nx; vy += ny;
                const geom_float seglen = std::sqrt(nx*nx + ny*ny);
                cx += 0.5*(Mx + Gx) * seglen;
                cy += 0.5*(My + Gy) * seglen;
                wsum += seglen;
            }

            dualFaceVect[3*ip + 0] = vx;
            dualFaceVect[3*ip + 1] = vy;
            dualFaceVect[3*ip + 2] = 0.0;
            dualFaceArea[ip] = std::sqrt(vx*vx + vy*vy);
            if (wsum > 0.0) {
                dualFaceCent[3*ip + 0] = cx / wsum;
                dualFaceCent[3*ip + 1] = cy / wsum;
            } else {
                dualFaceCent[3*ip + 0] = Mx;
                dualFaceCent[3*ip + 1] = My;
            }
            dualFaceCent[3*ip + 2] = 0.0;
        }

        // ---- 双対体積: セルごとに、各ノードの sub-CV 四角形 (A, M1, G, M2) の面積を集計。
        //      M1,M2 = そのセルで A に接する 2 エッジの中点。
        for (geom_int ic = 0; ic < this->nCells; ++ic)
        {
            const geom_float Gx = cells[ic].centCoords[0];
            const geom_float Gy = cells[ic].centCoords[1];

            for (const geom_int A : cells[ic].iNodes)
            {
                // A に接する 2 エッジの中点を集める
                geom_float mid[2][2];
                int nmid = 0;
                for (const geom_int ip : cells[ic].iPlanes)
                {
                    const geom_int e0 = this->planes[ip].iNodes[0];
                    const geom_int e1 = this->planes[ip].iNodes[1];
                    if (e0 != A && e1 != A) continue;
                    if (nmid < 2) {
                        mid[nmid][0] = 0.5*(nodes[e0].coords[0] + nodes[e1].coords[0]);
                        mid[nmid][1] = 0.5*(nodes[e0].coords[1] + nodes[e1].coords[1]);
                    }
                    nmid++;
                }
                if (nmid != 2) {
                    cerr << "[buildMedianDual] node " << A << " in cell " << ic
                         << " has " << nmid << " incident edges (expected 2).\n";
                    exit(EXIT_FAILURE);
                }

                // 四角形 A -> M1 -> G -> M2 の面積と面積加重重心 (shoelace)。
                const geom_float px[4] = { nodes[A].coords[0], mid[0][0], Gx, mid[1][0] };
                const geom_float py[4] = { nodes[A].coords[1], mid[0][1], Gy, mid[1][1] };
                geom_float area2 = 0.0, cx2 = 0.0, cy2 = 0.0;
                for (int k = 0; k < 4; ++k) {
                    const int kn = (k + 1) % 4;
                    const geom_float cross = px[k]*py[kn] - px[kn]*py[k];
                    area2 += cross;
                    cx2 += (px[k] + px[kn]) * cross;
                    cy2 += (py[k] + py[kn]) * cross;
                }
                const geom_float subArea = 0.5 * std::fabs(area2);
                dualVolume[A] += subArea;
                // ポリゴン重心 = (1/6A)Σ(p_k+p_{k+1})cross。area2 の符号と整合させ subArea 加重。
                if (std::fabs(area2) > 0.0) {
                    dualCentroid[3*A + 0] += subArea * (cx2 / (3.0 * area2));
                    dualCentroid[3*A + 1] += subArea * (cy2 / (3.0 * area2));
                    // z は 2D で node 値 (押し出し厚みの中央)。
                    dualCentroid[3*A + 2] += subArea * nodes[A].coords[2];
                }
            }
        }
        // 面積加重重心を正規化 (= 双対 CV の FV セル中心)。
        for (geom_int in = 0; in < nN; ++in) {
            const geom_float v = dualVolume[in];
            if (v > 0.0) {
                dualCentroid[3*in + 0] /= v;
                dualCentroid[3*in + 1] /= v;
                dualCentroid[3*in + 2] /= v;
            } else {
                dualCentroid[3*in + 0] = nodes[in].coords[0];
                dualCentroid[3*in + 1] = nodes[in].coords[1];
                dualCentroid[3*in + 2] = nodes[in].coords[2];
            }
        }

        // ---- 境界半割双対面: 各 bcond の境界エッジを構成ノードへ midpoint 分配 (各 1/2)。
        //
        // マルチマーカ方式 (SU2 流): 各境界ノードの半割面寄与は、その面が属する bcond (ow=ib) へ集約する。
        // コーナー (壁∩出口 等、複数マーカに属するノード) は各 incident bcond から 1 枚ずつ半割面を得る
        // (壁側半割面 + 出口側半割面)。両者のベクトル総和がコーナーの境界周を閉じるため CV が閉じ、出口/
        // 入口 BC がコーナー CV に正しく届く。旧 owner 優先所有 (wall>inlet>outlet>slip) では非所有側の
        // 半割面が欠落しコーナー CV が未閉だった (出口流出が壁所有コーナーで欠落)。壁ノードの u=0 は別途
        // Dirichlet (wall_flag + enforceWallNoSlip/zeroWallMomentumResidual) で厳密化する。
        // 閉性チェック (bnodeAccum) は幾何量なので所有に依らず全境界半割面を集計する。
        const geom_int nBc = (geom_int)this->bconds.size();
        dualBcondOffset.assign(nBc + 1, 0);
        dualBcondPhysID.assign(nBc, 0);
        dualBcondNodes.assign(nBc, {});
        std::vector<geom_float> bnodeAccum(nN * 3, 0.0); // 全体集計 (閉性チェック用、所有非依存)

        // 所有 bcond ごとに半割面ベクトル+面積加重重心を集計。
        std::vector<std::map<geom_int, std::array<geom_float,3>>> halfByOwner(nBc);
        std::vector<std::map<geom_int, std::array<geom_float,4>>> hcentByOwner(nBc);
        for (geom_int ib = 0; ib < nBc; ++ib)
        {
            for (const geom_int ip : this->bconds[ib].iPlanes)
            {
                const geom_int A = this->planes[ip].iNodes[0];
                const geom_int B = this->planes[ip].iNodes[1];
                const geom_float sv0 = this->planes[ip].surfVect[0];
                const geom_float sv1 = this->planes[ip].surfVect[1];
                const geom_float sv2 = this->planes[ip].surfVect[2];
                const geom_float hx = 0.5*sv0; // 外向き (makeMesh で整向済)
                const geom_float hy = 0.5*sv1;
                const geom_float hz = 0.5*sv2;
                const geom_float w = 0.5*std::sqrt(sv0*sv0 + sv1*sv1 + sv2*sv2); // 各ノードの半割面積
                for (const geom_int N : {A, B}) {
                    const geom_int O = (N == A) ? B : A;
                    // 閉性は全半割面を集計 (幾何、所有に依らない)
                    bnodeAccum[3*N + 0] += hx; bnodeAccum[3*N + 1] += hy; bnodeAccum[3*N + 2] += hz;
                    // 寄与はこの境界面が属する bcond へ (マルチマーカ: コーナーは各 incident bcond に計上)
                    const int ow = ib;
                    auto& h = halfByOwner[ow][N];
                    h[0] += hx; h[1] += hy; h[2] += hz;
                    // ノード N の半割 (N→エッジ中点 M) の重心 = (3N+O)/4
                    auto& c = hcentByOwner[ow][N];
                    c[0] += w * (3.0*nodes[N].coords[0] + nodes[O].coords[0]) * 0.25;
                    c[1] += w * (3.0*nodes[N].coords[1] + nodes[O].coords[1]) * 0.25;
                    c[2] += w * (3.0*nodes[N].coords[2] + nodes[O].coords[2]) * 0.25;
                    c[3] += w;
                }
            }
        }

        // emit: 各 bcond が自分の境界面に属するノードの半割面 plane を emit する (ghost を生成し、勾配/粘性の
        // 境界閉性を保つ)。マルチマーカなのでコーナー (壁∩出口) は壁・出口の双方から半割面 ghost を 1 枚ずつ
        // 得る (壁側 mirror ghost + 出口側 ghost)。両半割面ベクトルの総和でコーナー CV が閉じる。壁ノードの
        // u=0 はさらに Dirichlet (wall_flag + enforceWallNoSlip/zeroWallMomentumResidual) で厳密化する。
        for (geom_int ib = 0; ib < nBc; ++ib) {
            dualBcondPhysID[ib] = this->bconds[ib].physID;
            for (const auto& kv : halfByOwner[ib]) {
                const geom_int nd = kv.first;
                dualBcondNodes[ib].push_back(nd);
                dualBnodeId.push_back(nd);
                dualBnodeVect.push_back(kv.second[0]);
                dualBnodeVect.push_back(kv.second[1]);
                dualBnodeVect.push_back(kv.second[2]);
                const auto& c = hcentByOwner[ib].at(nd);
                const geom_float invw = (c[3] > 0.0) ? 1.0/c[3] : 0.0;
                dualBnodeCent.push_back(c[0]*invw);
                dualBnodeCent.push_back(c[1]*invw);
                dualBnodeCent.push_back(c[2]*invw);
            }
            dualBcondOffset[ib + 1] = (geom_int)dualBnodeId.size();
        }

        // ---- 整合チェック ----
        // 診断は double で集計 (geom_float=float の累積丸めをアルゴリズム誤差から分離するため)。
        // ① 体積総和
        double sumDual = 0.0, sumPrimal = 0.0;
        for (geom_int in = 0; in < nN; ++in) sumDual += dualVolume[in];
        for (geom_int ic = 0; ic < this->nCells; ++ic) sumPrimal += cells[ic].volume;
        const double volErr = std::fabs(sumDual - sumPrimal) / std::max(sumPrimal, 1e-30);
        cout << "[buildMedianDual] volume sum: dual=" << sumDual << " primal=" << sumPrimal
             << " relErr=" << volErr << "\n";

        // 負体積チェック
        geom_int nNeg = 0;
        for (geom_int in = 0; in < nN; ++in) if (dualVolume[in] <= 0.0) nNeg++;
        if (nNeg > 0) cout << "[buildMedianDual] WARNING: " << nNeg << " non-positive dual volumes\n";

        // ② 閉性: 各ノードで Σ(符号付き双対面ベクトル) + 半割面 = 0
        std::vector<double> clos(nN * 3, 0.0);
        for (geom_int ip = 0; ip < nP; ++ip)
        {
            const geom_int A = dualFaceCells[2*ip + 0];
            const geom_int B = dualFaceCells[2*ip + 1];
            for (int d = 0; d < 3; ++d) {
                clos[3*A + d] += dualFaceVect[3*ip + d]; // A=n0: 外向き +
                clos[3*B + d] -= dualFaceVect[3*ip + d]; // B=n1: 内向き -
            }
        }
        for (geom_int in = 0; in < nN; ++in)
            for (int d = 0; d < 3; ++d) clos[3*in + d] += bnodeAccum[3*in + d];

        double maxClos = 0.0, refArea = 0.0;
        for (geom_int ip = 0; ip < nP; ++ip) refArea = std::max(refArea, (double)dualFaceArea[ip]);
        for (geom_int in = 0; in < nN; ++in) {
            const double c = std::sqrt(clos[3*in+0]*clos[3*in+0]
                                     + clos[3*in+1]*clos[3*in+1]
                                     + clos[3*in+2]*clos[3*in+2]);
            maxClos = std::max(maxClos, c);
        }
        cout << "[buildMedianDual] closure max|sum dS| = " << maxClos
             << " (ref face area=" << refArea << ", normalized=" << maxClos/std::max(refArea,1e-30) << ")\n";

        // ③ 境界半割面積 == primal 境界面積 (bcond ごと)
        for (geom_int ib = 0; ib < nBc; ++ib)
        {
            geom_float primalA = 0.0;
            for (const geom_int ip : this->bconds[ib].iPlanes) primalA += this->planes[ip].surfArea;
            geom_float halfA = 0.0;
            for (geom_int k = dualBcondOffset[ib]; k < dualBcondOffset[ib+1]; ++k) {
                halfA += std::sqrt(dualBnodeVect[3*k+0]*dualBnodeVect[3*k+0]
                                 + dualBnodeVect[3*k+1]*dualBnodeVect[3*k+1]
                                 + dualBnodeVect[3*k+2]*dualBnodeVect[3*k+2]);
            }
            // 注: 半割の合計面積は端点で内向き相殺し得るため、面積総和ではなくベクトル総和=primal
            // で評価する方が厳密。ここでは参考値として両者を表示する。
            cout << "[buildMedianDual]   bcond physID=" << dualBcondPhysID[ib]
                 << " primalArea=" << primalA << " sum|halfVect|=" << halfA << "\n";
        }

        // 許容値は geom_float=float の累積丸めを考慮した値 (実バグは O(0.1-1) で明確に検出される)。
        const double closTol = 1e-3;
        const double volTol  = 1e-3;
        if (volErr > volTol) {
            cerr << "[buildMedianDual] ERROR: dual volume sum mismatch (relErr=" << volErr << ")\n";
            exit(EXIT_FAILURE);
        }
        if (maxClos/std::max(refArea,1e-30) > closTol) {
            cerr << "[buildMedianDual] ERROR: dual faces not closed (normalized="
                 << maxClos/std::max(refArea,1e-30) << ")\n";
            exit(EXIT_FAILURE);
        }

        dualBuilt = true;
        cout << "[buildMedianDual] done. nNodes(CV)=" << nN << " nDualFaces=" << nP
             << " nBHalf=" << dualBnodeId.size() << "\n";
    }

    // -------------------------------------------------------------------------
    // buildMedianDual() で計算した双対量で primal の cells/planes/bconds を置き換え、
    // 「双対メッシュを primary mesh として」書き出せる状態にする (node-centered)。
    //
    // これにより solver 側 (readMesh / setMeshMap_d / 各カーネル / 境界条件) は一切変更なしで
    // node-centered を扱える: CV=ノード、内部面=双対面、境界面=半割双対面 (solver の readMesh が
    // 境界面から自動でゴースト CV を生成する。ゴースト経由で既存 BC カーネル・対流フラックスが動く)。
    //
    // 注意 (M1): 境界半割面の重心 pc は「ノード座標 + h·n_out」に置く (h=0.5√(双対体積))。
    // これは solver のゴースト鏡像生成で非退化 (dcc>0) を保証するための便宜。半割面は厳密には
    // 境界上だが、M1 は非粘性 1 次でゴーストが 1 次精度を強制するため pc 位置はフラックスに影響しない
    // (fx/dcc は再構成・粘性で使われるが M1 では未使用)。厳密な境界は弱形式 (M2+) で扱う。
    // 可視化: 体積は退避した primal セル接続 (vizCONNE, Center='Node')、壁は退避した primal 境界面接続
    // (bcond.vizBfaceNodes, output.cpp で線/面を作り Center='Node') で正しく描ける。
    // -------------------------------------------------------------------------
    void replacePrimalWithDual()
    {
        if (!dualBuilt) {
            cerr << "[replacePrimalWithDual] buildMedianDual() must run first (dualBuilt=false).\n";
            exit(EXIT_FAILURE);
        }
        cout << "[replacePrimalWithDual] replacing primal cells/planes/bconds with median-dual ...\n";

        // ---- 可視化用に primal 境界面トポロジを bcond ごとに退避 (一度だけ) ----
        // 双対置換で境界も半割面(1ノード)になり面を作れなくなるため、置換前に primal 境界面の
        // global ノードid列を保存する (msh.vizCONNE の境界版)。壁出力で線/面を作り Center='Node' で描く。
        for (auto& bc : this->bconds) {
            bc.vizBfaceNodes.clear();
            bc.vizBfaceNodes.reserve(bc.iPlanes.size());
            for (const geom_int ip : bc.iPlanes)
                bc.vizBfaceNodes.push_back(this->planes[ip].iNodes);
        }

        const geom_int nN = this->nNodes;
        const geom_int nDualInternal = (geom_int)this->dualFaceArea.size();
        const geom_int nBHalf        = (geom_int)this->dualBnodeId.size();

        // ---- 可視化用に primal セルトポロジ (CONNE) を退避 ----
        // 双対へ置換すると primal topology が失われるため、node モードの出力で
        // 「primal セル + Center='Node'」可視化に使えるよう CONNE を保存する。
        this->vizCONNE.clear();
        this->vizCONNE_dim = 0;
        this->nVizCells = (geom_int)this->cells.size();
        for (const auto& cel : this->cells)
        {
            const geom_int nn = this->eleTypeMap.mapElementFromGmshID[cel.ieleType].nNodes;
            const std::string name = this->eleTypeMap.mapElementFromGmshID[cel.ieleType].name;
            geom_int code = 0;
            if (name == "hex") code = 9;
            else if (name == "prism") code = 8;
            else if (name == "pyramid") code = 7;
            else if (name == "tetra") code = 6;
            else if (name == "quad") code = 5;
            else if (name == "triangle") code = 4;
            this->vizCONNE.push_back(code);
            this->vizCONNE_dim += nn + 1;
            for (const geom_int nod : cel.iNodes) this->vizCONNE.push_back(nod);
        }

        // ---- 新 cells (CV = ノード) ----
        std::vector<cell> newCells(nN);
        for (geom_int i = 0; i < nN; ++i) {
            // CV 重心: axisCentroidShift=true なら双対 CV の面積加重重心 (軸上でも R>0 で回転体積を稼ぐ)、
            // false なら node 座標 (軸上は R=0)。後者は軸ソース OFF と併用前提 (ソースが唯一の r 非重み項)。
            if (this->axisCentroidShift)
                newCells[i].centCoords = { dualCentroid[3*i+0], dualCentroid[3*i+1], dualCentroid[3*i+2] };
            else
                newCells[i].centCoords = nodes[i].coords;
            newCells[i].volume     = dualVolume[i];
            newCells[i].ieleType   = 2;                     // placeholder (output CONNE 用, viz は M4)
            newCells[i].iNodes      = { i };
            newCells[i].regionId    = 0;
            // iPlanes/iPlanesDir は下で充填
        }

        // ---- 新 planes: 内部双対面 [0,nDualInternal) + 境界半割面 [nDualInternal, +nBHalf) ----
        std::vector<plane> newPlanes;
        newPlanes.reserve(nDualInternal + nBHalf);
        for (geom_int f = 0; f < nDualInternal; ++f) {
            const geom_int A = dualFaceCells[2*f + 0];
            const geom_int B = dualFaceCells[2*f + 1];
            plane p;
            p.iNodes = { A, B };
            p.iCells = { A, B };
            p.surfVect = { dualFaceVect[3*f+0], dualFaceVect[3*f+1], dualFaceVect[3*f+2] };
            p.surfArea = dualFaceArea[f];
            p.centCoords = { dualFaceCent[3*f+0], dualFaceCent[3*f+1], dualFaceCent[3*f+2] };
            newPlanes.push_back(std::move(p));
        }
        for (geom_int k = 0; k < nBHalf; ++k) {
            const geom_int nd = dualBnodeId[k];
            const geom_float vx = dualBnodeVect[3*k+0];
            const geom_float vy = dualBnodeVect[3*k+1];
            const geom_float vz = dualBnodeVect[3*k+2];
            const geom_float area = std::sqrt(vx*vx + vy*vy + vz*vz);
            plane p;
            p.iNodes = { nd };
            p.iCells = { nd };  // 1 セルのみ → solver readMesh がゴーストを追加
            p.surfVect = { vx, vy, vz };
            p.surfArea = area;
            // pc = 半割面の真の面積加重重心 (境界上, R≥0)。軸対称 r 重みが正しくなり、入口/出口 BC が
            // 軸近傍 corner CV に届く (旧 node+h·n_out 便宜は pcy≈0/<0 で BC を r 重み消失させていた)。
            p.centCoords = { dualBnodeCent[3*k+0], dualBnodeCent[3*k+1], dualBnodeCent[3*k+2] };
            newPlanes.push_back(std::move(p));
        }

        this->nNormalPlanes = nDualInternal;
        this->nBPlanes      = nBHalf;
        this->nPlanes       = nDualInternal + nBHalf;
        this->nCells        = nN;
        this->nCells_all    = nN; // ゴーストは solver readMesh で付与

        // ---- cells の iPlanes / iPlanesDir を充填 ----
        for (geom_int ip = 0; ip < this->nPlanes; ++ip) {
            const auto& ic = newPlanes[ip].iCells;
            for (size_t j = 0; j < ic.size(); ++j) {
                newCells[ic[j]].iPlanes.push_back(ip);
                newCells[ic[j]].iPlanesDir.push_back(j == 0 ? 1 : -1);
            }
        }

        this->cells  = std::move(newCells);
        this->planes = std::move(newPlanes);

        // ---- bconds を双対境界半割面で更新 (bcondKind / physID / 入力値は保持) ----
        const geom_int nBc = (geom_int)this->bconds.size();
        for (geom_int ib = 0; ib < nBc; ++ib) {
            if (this->bconds[ib].physID != dualBcondPhysID[ib]) {
                cerr << "[replacePrimalWithDual] bcond order mismatch at ib=" << ib
                     << " (physID " << this->bconds[ib].physID << " vs " << dualBcondPhysID[ib] << ")\n";
                exit(EXIT_FAILURE);
            }
            this->bconds[ib].iPlanes.clear();
            this->bconds[ib].iCells.clear();
            this->bconds[ib].iBPlanes.clear();
            this->bconds[ib].iCells_ghst.clear();
            // iCells = この bcond に属するノード全て (壁含む。コーナーは複数 bcond の iCells に重複出現)。
            //          solver はこの壁 iCells から wall_flag を構築し、コーナーも壁ノードとして拾う。
            this->bconds[ib].iCells = dualBcondNodes[ib];
            // iPlanes/iBPlanes = この bcond が emit した境界半割面 (壁含め全 bcond が ghost を持つ。コーナーは
            // 壁側・出口側で別々の半割面 plane を持つため iCells と 1:1 対応する)。
            // 壁 ghost は mirror で u_face=0 を与え勾配/粘性の境界閉性を保つ (u=0 厳密化は Dirichlet 側)。
            for (geom_int k = dualBcondOffset[ib]; k < dualBcondOffset[ib+1]; ++k) {
                const geom_int gp = nDualInternal + k;
                this->bconds[ib].iPlanes.push_back(gp);
                this->bconds[ib].iBPlanes.push_back(gp - this->nNormalPlanes);
            }
        }

        // 双対は primary に昇格したので /DUAL の重複出力は不要にする
        dualBuilt = false;

        cout << "[replacePrimalWithDual] done. nCells(CV)=" << this->nCells
             << " nNormalPlanes=" << this->nNormalPlanes
             << " nBPlanes=" << this->nBPlanes
             << " nPlanes=" << this->nPlanes << "\n";
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
            if (name == "quad") CONNE0 = 5;
            if (name == "triangle") CONNE0 = 4;

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

        vector<geom_int> regionIds;
        for (auto& cel : this->cells) {
            regionIds.push_back(cel.regionId);
        }
        file.createDataSet("/CELLS/regionId", regionIds);

        // ---- node-centered 可視化トポロジ (/VIZMESH) ----
        // replacePrimalWithDual() で退避した primal CONNE。solver の出力が node モードで
        // 「primal セル + Center='Node'」可視化に使う (docs/discretization/implementation.md §3.5)。
        if (!this->vizCONNE.empty())
        {
            Group vgrp = file.createGroup("/VIZMESH");
            Attribute va = vgrp.createAttribute<geom_int>("nVizCells", DataSpace::From(this->nVizCells));
            va.write(this->nVizCells);
            va = vgrp.createAttribute<geom_int>("vizCONNE_dim", DataSpace::From(this->vizCONNE_dim));
            va.write(this->vizCONNE_dim);
            file.createDataSet("/VIZMESH/CONNE", this->vizCONNE);
            cout << "writeInputH5: wrote /VIZMESH (nVizCells=" << this->nVizCells << ")\n";
        }

        // ---- median-dual (node-centered) データ ----
        // buildMedianDual() 済みのときのみ /DUAL を書く。solver は discretization=="node"
        // で /DUAL を読み CV=ノードとして消費する (docs/discretization/implementation.md §1)。
        if (this->dualBuilt)
        {
            Group dgrp = file.createGroup("/DUAL");
            const geom_int nDualFaces = (geom_int)this->dualFaceArea.size();
            const geom_int nBc        = (geom_int)this->dualBcondPhysID.size();

            Attribute da = dgrp.createAttribute<geom_int>("nNodes", DataSpace::From(this->nNodes));
            da.write(this->nNodes);
            da = dgrp.createAttribute<geom_int>("nDualFaces", DataSpace::From(nDualFaces));
            da.write(nDualFaces);
            da = dgrp.createAttribute<geom_int>("nBconds", DataSpace::From(nBc));
            da.write(nBc);

            file.createDataSet("/DUAL/volume",    this->dualVolume);
            file.createDataSet("/DUAL/faceCells", this->dualFaceCells);
            file.createDataSet("/DUAL/faceVect",  this->dualFaceVect);
            file.createDataSet("/DUAL/faceArea",  this->dualFaceArea);
            file.createDataSet("/DUAL/faceCent",  this->dualFaceCent);
            file.createDataSet("/DUAL/bnodeId",   this->dualBnodeId);
            file.createDataSet("/DUAL/bnodeVect", this->dualBnodeVect);
            file.createDataSet("/DUAL/bcondOffset", this->dualBcondOffset);
            file.createDataSet("/DUAL/bcondPhysID", this->dualBcondPhysID);
            cout << "writeInputH5: wrote /DUAL (nDualFaces=" << nDualFaces << ")\n";
        }

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

            // node モード壁可視化: primal 境界面接続 (flat global ids + 面ごとのノード数)。
            if (!bc.vizBfaceNodes.empty()) {
                std::vector<geom_int> bfNodes, bfSizes;
                for (const auto& f : bc.vizBfaceNodes) {
                    bfSizes.push_back((geom_int)f.size());
                    for (const geom_int n : f) bfNodes.push_back(n);
                }
                group.createDataSet("/BCONDS/"+oss.str()+"/vizBfaceNodes", bfNodes);
                group.createDataSet("/BCONDS/"+oss.str()+"/vizBfaceSizes", bfSizes);
            }

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
        // node-centered では /VIZMESH/CONNE (primal トポロジ) + Center='Node' で可視化する。
        const bool nodeViz = !this->vizCONNE.empty();
        const geom_int xdmfNElem = nodeViz ? this->nVizCells : (geom_int)this->cells.size();
        const geom_int xdmfConneDim = nodeViz ? this->vizCONNE_dim : CONNE_dim;
        const string conneRef = nodeViz ? (outFileName + ":VIZMESH/CONNE")
                                        : (outFileName + ":MESH/CONNE");
        const char* centerAttr = nodeViz ? "Node" : "Cell";

        ofs << "      <Grid Name='aaa'>\n";
        ofs << "        <Topology Type='Mixed' NumberOfElements='" << xdmfNElem << "'>\n";
        ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << xdmfConneDim << "'>\n";
        ofs << "           " + conneRef + "\n";
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
            ofs << "        <Attribute Name='"  << name << "' Center='" << centerAttr << "' >\n";
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
