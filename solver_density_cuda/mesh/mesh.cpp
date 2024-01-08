#include "mesh.hpp"
#include "cuda_nagare/cudaWrapper.cuh"

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

bcond::~bcond() 
{
    for (auto& item : this->bvar_d)
    {
        //gpuErrchk( cudaFree(item.second) );
        cudaWrapper::cudaFree_wrapper(item.second);
    }
    cudaWrapper::cudaFree_wrapper(this->map_bplane_plane_d);
    cudaWrapper::cudaFree_wrapper(this->map_bplane_cell_d);
}
bcond::bcond(const geom_int& physID, const vector<geom_int>& iPlanes, 
             const vector<geom_int>& iCells , const vector<geom_int>& iBPlanes)
{
    this->physID   = physID;
    this->iPlanes  = iPlanes;
    this->iCells   = iCells;
    this->iBPlanes = iBPlanes;
}

void bcond::bcondInitVariables(const int &useGPU)
{
    // ------------------------------------------------------
    // *** allocate and init variables for boundary cells ***
    // ------------------------------------------------------
    // float


    for (auto& bValName : this->bplaneValNames)
    {
        if (valueTypes.find(bValName) == valueTypes.end()) { // not neccesally variable for this bc
            continue;

        } else { // needed
            this->bvar[bValName].resize(this->iPlanes.size());

            int type = valueTypes[bValName];
            cout  << "bValName = " <<  bValName <<  " type=" << type << endl;

            if (type == 1) { // read uniform float from yaml
                cout << "read " << bValName << " of " << bcondKind  << " from config";
                for (flow_float& var : this->bvar[bValName])
                {
                    var = this->inputFloats[bValName];
                }
            }

            if (useGPU == 1){
                gpuErrchk( cudaMalloc(&(this->bvar_d[bValName]) , this->iPlanes.size()*sizeof(flow_float)) );
                int type = valueTypes[bValName];

                if (type == 1) { // read uniform float from yaml
                    gpuErrchk( cudaMemcpy(this->bvar_d[bValName] , &(this->bvar[bValName][0]) ,
                                                 (size_t)(this->iPlanes.size()*sizeof(flow_float)), 
                                                 cudaMemcpyHostToDevice) );
                }
            }
       }
    }

    // ----------------------------------------------------------
    // *** allocate and init INT variables for boundary cells ***
    // ----------------------------------------------------------
    for (auto& bIntName : this->bplaneIntNames)
    {

        if (valueTypes.find(bIntName) == valueTypes.end()) { // not found
            continue;

        } else { // found
            this->bint[bIntName].resize(this->iPlanes.size());
            int type = valueTypes[bIntName];

            if (type == 1) { // read uniform int from yaml
                cout << "read uniform" << bIntName << " of " << bcondKind  << " from config\n";
                for (flow_float& var : this->bint[bIntName])
                {
                    var = this->inputFloats[bIntName];
                }
            } else if (type == 10) { // special 
                continue; // set later
            }

            if (useGPU == 1){
                gpuErrchk( cudaMalloc(&(this->bint_d[bIntName]) , this->iPlanes.size()*sizeof(flow_float)) );

                int type = valueTypes[bIntName];

                if (type == 1) { // read uniform int from yaml
                    gpuErrchk( cudaMemcpy(this->bint_d[bIntName] , &(this->bint[bIntName][0]) ,
                                                 (size_t)(this->iPlanes.size()*sizeof(flow_float)), 
                                                 cudaMemcpyHostToDevice) );
                } // TODO: set special treatment
            }
        }

    }
}


mesh::mesh(){}
mesh::~mesh()
{
    cudaWrapper::cudaFree_wrapper(this->map_nplane_cells_d);
}

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

//ghst>
    this->nCells_ghst = this->nBPlanes;
    this->nCells_all  = this->nCells + this->nCells_ghst;
    //this->nCells_all  = this->nCells ;
//ghst<

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
 
//ghst>
    // Cells including ghost cells 
//    this->cells.resize(this->nCells);
    this->cells.resize(this->nCells_all);
//ghst<
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
    geom_int nGhost = 0;
    for (string oname : grp.listObjectNames())
    {
        this->bconds[ib].physID = stoi(oname);
        grp = file.getGroup("/BCONDS/"+oname);

        a = grp.getAttribute("bcondKind");
        this->bconds[ib].bcondKind = a.read<std::string>();
        cout << "in mesh.cpp  physID=" << this->bconds[ib].physID << endl;
        cout << "in mesh.cpp  bcondKind=" << this->bconds[ib].bcondKind << endl;

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

        // iBPlanes
        std::vector<geom_int> iBPlanes;
        grp.getDataSet("iBPlanes").read(iBPlanes);

        this->bconds[ib].iBPlanes.resize(iBPlanes.size());
        for (geom_int ibp = 0 ; ibp<iBPlanes.size() ; ibp++)
        {
            this->bconds[ib].iBPlanes[ibp] = iBPlanes[ibp];
        }

//ghst>
        // add ghost cells
        for (geom_int ibl = 0 ; ibl<iBPlanes.size() ; ibl++)
        {
            geom_int ip = iPlanes[ibl];
            geom_int ic;

            // change plane information
            this->bconds[ib].iCells_ghst.resize(iBPlanes.size());
            this->bconds[ib].iCells_ghst[ibl] = nCells+nGhost;
            if (this->planes[ip].iCells.size() == 1) {
                this->planes[ip].iCells.resize(2);
                this->planes[ip].iCells[1] =nCells+nGhost;

                ic = this->planes[ip].iCells[0];
            } else {
                std::cerr << "Error: something wrong to add ghost cells" << endl;
                exit(3);
            }

            // make ghost cell
            this->cells[nCells+nGhost].iPlanes.resize(1);
            this->cells[nCells+nGhost].iPlanes[0] = ip;

            this->cells[nCells+nGhost].ieleType = this->cells[ic].ieleType;
            this->cells[nCells+nGhost].volume   = this->cells[ic].volume;
            this->cells[nCells+nGhost].centCoords.resize(3);

            geom_float xc = this->cells[ic].centCoords[0];
            geom_float yc = this->cells[ic].centCoords[1];
            geom_float zc = this->cells[ic].centCoords[2];

            geom_float xp = this->planes[ip].centCoords[0];
            geom_float yp = this->planes[ip].centCoords[1];
            geom_float zp = this->planes[ip].centCoords[2];

            geom_float dx = xp - xc;
            geom_float dy = yp - yc;
            geom_float dz = zp - zc;
            
            this->cells[nCells+nGhost].centCoords[0] = xp + dx;
            this->cells[nCells+nGhost].centCoords[1] = yp + dy;
            this->cells[nCells+nGhost].centCoords[2] = zp + dz;
            nGhost++;
        }
//ghst<
        ib++;
    }
}

void mesh::setPeriodicPartner()
{
    std::list<int> checked_bcIDs;

    for (auto& bc0 : this->bconds) {
        if (bc0.bcondKind == "periodic") {
            int bcID         = bc0.physID;
            int bcID_partner = bc0.inputInts["partnerBCID"];

            if (std::find(checked_bcIDs.begin(), checked_bcIDs.end(), bcID) == checked_bcIDs.end()){ // not saved BC
                checked_bcIDs.push_back(bcID);
                checked_bcIDs.push_back(bcID_partner);

                std::vector<geom_int> map_ib1_iplane;

                // find nearest planes of partner
                for (auto& bc1 : this->bconds) {
                    if (bc1.physID != bcID_partner) continue;

                    Eigen::VectorXd XYZ(3);
                    Eigen::MatrixXd partnerXYZ(3, bc1.iPlanes.size());
                    Eigen::MatrixXd::Index index;

                    geom_int ib1_local = 0;
                    for (geom_int& ip1 : bc1.iPlanes) {
                        geom_float x1 = this->planes[ip1].centCoords[0];
                        geom_float y1 = this->planes[ip1].centCoords[1];
                        geom_float z1 = this->planes[ip1].centCoords[2];

                        partnerXYZ(0, ib1_local) = x1;
                        partnerXYZ(1, ib1_local) = y1;
                        partnerXYZ(2, ib1_local) = z1;

                        map_ib1_iplane.push_back(ip1);

                        ib1_local++;
                    }

                    geom_int ib0_local = 0;
                    for (geom_int& ip0 : bc0.iPlanes) {
                        geom_float x0 = this->planes[ip0].centCoords[0];
                        geom_float y0 = this->planes[ip0].centCoords[1];
                        geom_float z0 = this->planes[ip0].centCoords[2];

                        XYZ(0) = x0;
                        XYZ(1) = y0;
                        XYZ(2) = z0;

                        (partnerXYZ.colwise() - XYZ).colwise().squaredNorm().minCoeff(&index);
                        //cout << "Nearest neibour is column " << index << ":" << endl;
                        //cout << partnerXYZ.col(index) << endl;

                        geom_int ip1 = map_ib1_iplane[index];

                        bc0.bint["partnerPlnID"][ib0_local] = ip1;
                        bc1.bint["partnerPlnID"][index]     = ip0;

                        ib0_local++;
                    }
                } 
            }else {
                continue; // already added bc
            }
        }
    }
};


void mesh::setMeshMap_d()
{
    gpuErrchk(cudaMalloc((void **)&(this->map_nplane_cells_d), sizeof(geom_int)*this->nNormalPlanes*2));

    geom_int* pc_h;
    geom_int* bp_h;
    geom_int* bc_h;
    pc_h = (geom_int *)malloc(sizeof(geom_int)*this->nNormalPlanes*2);

    for (geom_int ip=0; ip<this->nNormalPlanes; ip++)
    {
        pc_h[2*ip + 0] = this->planes[ip].iCells[0]; 
        pc_h[2*ip + 1] = this->planes[ip].iCells[1]; 
    }

    gpuErrchk(cudaMemcpy(this->map_nplane_cells_d  , pc_h , 
                     sizeof(geom_int)*(this->nNormalPlanes*2) , cudaMemcpyHostToDevice));

    free(pc_h); 


    for (bcond& bc : this->bconds)
    {
        bp_h = (geom_int *)malloc(sizeof(geom_int)*bc.iPlanes.size());
        bc_h = (geom_int *)malloc(sizeof(geom_int)*bc.iPlanes.size());
        gpuErrchk(cudaMalloc((void **)&(bc.map_bplane_plane_d), sizeof(geom_int)*bc.iPlanes.size()));
        gpuErrchk(cudaMalloc((void **)&(bc.map_bplane_cell_d) , sizeof(geom_int)*bc.iPlanes.size()));

        for (geom_int ibl=0 ; ibl<bc.iPlanes.size() ; ibl++)
        {
            bp_h[ibl] = bc.iPlanes[ibl];
            bc_h[ibl] = bc.iCells[ibl];
        }

        gpuErrchk(cudaMemcpy(bc.map_bplane_plane_d , bp_h , 
                             sizeof(geom_int)*(bc.iPlanes.size()) , cudaMemcpyHostToDevice));

        gpuErrchk(cudaMemcpy(bc.map_bplane_cell_d , bc_h , 
                             sizeof(geom_int)*(bc.iPlanes.size()) , cudaMemcpyHostToDevice));
        free(bp_h); 
        free(bc_h);
    }
};

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

