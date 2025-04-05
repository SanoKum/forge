#include "mesh.hpp"
#include "cuda_forge/cudaWrapper.cuh"

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

            if (type == 1) { // read uniform float from yaml
                cout << "read " << bValName << " of " << bcondKind  << " from config\n";
                for (flow_float& var : this->bvar[bValName])
                {
                    var = this->inputFloats[bValName];
                }
            }

            if (useGPU == 1){
                cout << "alloc " << bValName << " of " << bcondKind  << " from config\n";
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
                for (geom_int& var : this->bint[bIntName])
                {
                    var = this->inputFloats[bIntName];
                }
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

void bcond::copyVariables_bplane_D2H()
{
    // float
    for (auto& bValName : this->bplaneValNames)
    {
        if (valueTypes.find(bValName) == valueTypes.end()) { // not neccesally variable for this bc
            continue;

        } else { // needed
            cudaWrapper::cudaMemcpy_D2H_wrapper(this->bvar_d[bValName], this->bvar[bValName].data(), this->bvar[bValName].size());
        }
    }
}

void bcond::output_preparation(std::vector<node>& nodes, std::vector<plane>& planes)
{
    if (this->output_preparation_flg == 1) return;

    vector<geom_int> nodes_bc_flg;

    nodes_bc_flg.resize(nodes.size());
    this->inodes_g2l.resize(nodes.size());
    std::fill(nodes_bc_flg.begin(), nodes_bc_flg.end(), 0);
    std::fill(inodes_g2l.begin(), inodes_g2l.end(), -1);

    geom_int ipl = 0;
    for (auto& ip : this->iPlanes) {
        this->planes_local.push_back(planes[ip]);

        for (auto& in : planes[ip].iNodes) {
            //cout << "inl=" << inl << endl;
            nodes_bc_flg[in] = 1;
            //this->inodes_g2l[in] = inl;
            //this->inodes_l2g.push_back(in);
            //inl++;
        }
        ipl++;
    }

    geom_int inl = -1;
    for (geom_int ing=0 ; ing<nodes.size(); ing++) {
        if (nodes_bc_flg[ing] == 1){
            inl++;
            this->inodes_g2l[ing] = inl;
            this->inodes_l2g.push_back(ing);
        }
    }

    this->output_preparation_flg = 1;

}

mesh::mesh(){}
mesh::~mesh()
{
    cudaWrapper::cudaFree_wrapper(this->map_plane_cells_d);
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
    cout << "Number of Cells: " << this->nCells << endl;

    a = group.getAttribute("nPlanes");
    this->nPlanes = a.read<geom_int>();
    cout << "Number of Planes: " << this->nPlanes << endl;

    a = group.getAttribute("nNormalPlanes");
    this->nNormalPlanes = a.read<geom_int>();
    cout << "Number of Normal Planes: " << this->nNormalPlanes << endl;

    a = group.getAttribute("nBPlanes");
    this->nBPlanes= a.read<geom_int>();
    cout << "Number of Boundary Planes: " << this->nBPlanes << endl;

    a = group.getAttribute("nBconds");
    this->nBconds = a.read<geom_int>();
    cout << "Number of Boundary Conditions: " << this->nBconds << endl;

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
        cout << "               ip min=" << iPlanes[0] << ", ip max=" << iPlanes[iPlanes.size()-1]<< endl;

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

            geom_float ss = this->planes[ip].surfArea;
            geom_float nx = this->planes[ip].surfVect[0]/ss;
            geom_float ny = this->planes[ip].surfVect[1]/ss;
            geom_float nz = this->planes[ip].surfVect[2]/ss;

            geom_float dnx = (dx*nx +dy*ny + dz*nz)*nx;
            geom_float dny = (dx*nx +dy*ny + dz*nz)*ny;
            geom_float dnz = (dx*nx +dy*ny + dz*nz)*nz;
            
            this->cells[nCells+nGhost].centCoords[0] = xc + 2*dnx;
            this->cells[nCells+nGhost].centCoords[1] = yc + 2*dny;
            this->cells[nCells+nGhost].centCoords[2] = zc + 2*dnz;
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

            flow_float dx;
            flow_float dy;
            flow_float dz;
            flow_float dtheta;

            if (bc0.inputInts["type"] == 0) { // Cartesian
                dx = bc0.inputFloats["dx"];
                dy = bc0.inputFloats["dy"];
                dz = bc0.inputFloats["dz"];

            } else if (bc0.inputInts["type"] == 1) { // rotation
                dtheta = bc0.inputFloats["dtheta"];
            }

            if (std::find(checked_bcIDs.begin(), checked_bcIDs.end(), bcID) == checked_bcIDs.end()){ // not saved BC
                checked_bcIDs.push_back(bcID);
                checked_bcIDs.push_back(bcID_partner);

                std::vector<geom_int> map_ib1_iplane;

                // find nearest planes of partner
                for (auto& bc1 : this->bconds) {
                    if (bc1.physID != bcID_partner) continue;

                    //Eigen::VectorXd XYZ(3);
                    //Eigen::MatrixXd partnerXYZ(3, bc1.iPlanes.size());
                    //Eigen::MatrixXd::Index index;

                    vector<geom_float> XYZ(3);
                    vector<vector<geom_float>> partnerXYZ;
                    geom_float index;

                    partnerXYZ.resize(3);
                    partnerXYZ[0].resize(bc1.iPlanes.size());
                    partnerXYZ[1].resize(bc1.iPlanes.size());
                    partnerXYZ[2].resize(bc1.iPlanes.size());

                    geom_float x1;
                    geom_float y1;
                    geom_float z1;

                    geom_float x1_r1;
                    geom_float y1_r1;
                    geom_float z1_r1;

                    geom_int ib1_local = 0;
                    for (geom_int& ip1 : bc1.iPlanes) {
                        x1 = this->planes[ip1].centCoords[0];
                        y1 = this->planes[ip1].centCoords[1];
                        z1 = this->planes[ip1].centCoords[2];

                        if (bc0.inputInts["type"] == 0) {
                            x1_r1 = x1 - dx;
                            y1_r1 = y1 - dy;
                            z1_r1 = z1 - dz;

                        } else if (bc0.inputInts["type"] == 1) {
                            x1_r1 = x1;
                            y1_r1 = cos(-dtheta)*y1 -sin(-dtheta)*z1;
                            z1_r1 = sin(-dtheta)*y1 +cos(-dtheta)*z1;
                        }

                        partnerXYZ[0][ib1_local] = x1_r1;
                        partnerXYZ[1][ib1_local] = y1_r1;
                        partnerXYZ[2][ib1_local] = z1_r1;

                        map_ib1_iplane.push_back(ip1);

                        ib1_local++;
                    }

                    geom_int ib0_local = 0;
                    for (geom_int& ip0 : bc0.iPlanes) {
                        geom_float x0 = this->planes[ip0].centCoords[0];
                        geom_float y0 = this->planes[ip0].centCoords[1];
                        geom_float z0 = this->planes[ip0].centCoords[2];

                        XYZ[0] = x0;
                        XYZ[1] = y0;
                        XYZ[2] = z0;

                        geom_float dist2 = 1e+30;
                        geom_float dist2_temp;
                        geom_int ib1_local = 0;
                        for (geom_int& ip1 : bc1.iPlanes) {
                            x1 = partnerXYZ[0][ib1_local] ;
                            y1 = partnerXYZ[1][ib1_local] ;
                            z1 = partnerXYZ[2][ib1_local] ;

                            dist2_temp = pow((x1-x0),2) +pow((y1-y0),2) +pow((z1-z0),2);

                            if (dist2>dist2_temp) {
                                dist2 = dist2_temp;
                                index = ib1_local;
                            }

                            ib1_local++;
                        }

                        // don't use eigen because too many warnings
                        //(partnerXYZ.colwise() - XYZ).colwise().squaredNorm().minCoeff(&index);

                        geom_int ip1 = map_ib1_iplane[index];

                        geom_int ic0 = this->planes[ip0].iCells[0];
                        geom_int ic1 = this->planes[ip1].iCells[0];

                        bc0.bint["partnerPlnID"][ib0_local] = ip1;
                        bc1.bint["partnerPlnID"][index]     = ip0;

                        bc0.bint["partnerCellID"][ib0_local] = ic1;
                        bc1.bint["partnerCellID"][index]     = ic0;


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
    gpuErrchk(cudaMalloc((void **)&(this->map_plane_cells_d), sizeof(geom_int)*this->nPlanes*2));

    geom_int n_normal_ghst_planes =  this->nNormalPlanes;

    for (auto& bc : this->bconds)
    {
        if (bc.bcondKind == "periodic") {
            n_normal_ghst_planes += bc.iPlanes.size();
        }
    }

    this->nNormal_ghst_Planes = n_normal_ghst_planes;


    gpuErrchk(cudaMalloc((void **)&(this->normal_ghst_planes_d), sizeof(geom_int)*n_normal_ghst_planes));

    geom_int* normal_ghst_planes;
    normal_ghst_planes = (geom_int *)malloc(sizeof(geom_int)*n_normal_ghst_planes);


    geom_int ip_sum = 0;
    for (geom_int ip=0; ip<this->nNormalPlanes; ip++ ) {
        normal_ghst_planes[ip] = ip_sum;
        ip_sum += 1;

    }

     for (auto& bc : this->bconds)
    {
        if (bc.bcondKind == "periodic") {
            for (auto& ip : bc.iPlanes){
                normal_ghst_planes[ip_sum] = ip;
                ip_sum += 1;
            }
        }
    }

    gpuErrchk(cudaMemcpy(this->normal_ghst_planes_d  , normal_ghst_planes , 
                         sizeof(geom_int)*n_normal_ghst_planes , cudaMemcpyHostToDevice));

    free(normal_ghst_planes); 


    geom_int* pc_h;
    geom_int* bp_h;
    geom_int* bc_h;
    geom_int* bcg_h;
    pc_h = (geom_int *)malloc(sizeof(geom_int)*this->nPlanes*2);

    for (geom_int ip=0; ip<this->nPlanes; ip++)
    {
        pc_h[2*ip + 0] = this->planes[ip].iCells[0]; 
        pc_h[2*ip + 1] = this->planes[ip].iCells[1]; 
        //printf("ip=%d, ic1=%d, ic2=%d\n", ip, pc_h[2*ip + 0], pc_h[2*ip + 1]);
    }

    gpuErrchk(cudaMemcpy(this->map_plane_cells_d  , pc_h , 
                     sizeof(geom_int)*(this->nPlanes*2) , cudaMemcpyHostToDevice));

    free(pc_h); 


    // cell -> planes
    gpuErrchk(cudaMalloc((void **)&(this->map_cell_planes_index_d), sizeof(geom_int)*(this->nCells+1)));

    geom_int* cpi_h;
    cpi_h = (geom_int *)malloc(sizeof(geom_int)*(this->nCells+1));

    cpi_h[0] = 0; 
    for (geom_int ic=0; ic<this->nCells; ic++)
    {
        geom_int nplane_local = this->cells[ic].iPlanes.size();
        cpi_h[ic+1] = cpi_h[ic] + nplane_local; 
    }

    gpuErrchk(cudaMalloc((void **)&(this->map_cell_planes_d), sizeof(geom_int)*cpi_h[this->nCells]));
    geom_int* cp_h;
    cp_h = (geom_int *)malloc(sizeof(geom_int)*cpi_h[this->nCells]);

    geom_int ipln = 0;
    for (geom_int ic=0; ic<this->nCells; ic++)
    {
        for (auto ip :this->cells[ic].iPlanes) {
            cp_h[ipln] = ip; 
            ipln += 1;
        }
    }

    gpuErrchk(cudaMemcpy(this->map_cell_planes_index_d  , cpi_h , 
                         sizeof(geom_int)*(this->nCells+1) , cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(this->map_cell_planes_d  , cp_h , 
                     sizeof(geom_int)*cpi_h[this->nCells] , cudaMemcpyHostToDevice));

    free(cpi_h); 
    free(cp_h); 


    for (bcond& bc : this->bconds)
    {
        bp_h = (geom_int *)malloc(sizeof(geom_int)*bc.iPlanes.size());
        bc_h = (geom_int *)malloc(sizeof(geom_int)*bc.iPlanes.size());
        bcg_h = (geom_int *)malloc(sizeof(geom_int)*bc.iPlanes.size()); // ghost cell
        gpuErrchk(cudaMalloc((void **)&(bc.map_bplane_plane_d), sizeof(geom_int)*bc.iPlanes.size()));
        gpuErrchk(cudaMalloc((void **)&(bc.map_bplane_cell_d) , sizeof(geom_int)*bc.iPlanes.size()));
        gpuErrchk(cudaMalloc((void **)&(bc.map_bplane_cell_ghst_d) , sizeof(geom_int)*bc.iPlanes.size()));

        for (geom_int ibl=0 ; ibl<bc.iPlanes.size() ; ibl++)
        {
            bp_h[ibl]  = bc.iPlanes[ibl];
            bc_h[ibl]  = bc.iCells[ibl];
            bcg_h[ibl] = bc.iCells_ghst[ibl];
        }

        gpuErrchk(cudaMemcpy(bc.map_bplane_plane_d , bp_h , 
                             sizeof(geom_int)*(bc.iPlanes.size()) , cudaMemcpyHostToDevice));

        gpuErrchk(cudaMemcpy(bc.map_bplane_cell_d , bc_h , 
                             sizeof(geom_int)*(bc.iPlanes.size()) , cudaMemcpyHostToDevice));

        gpuErrchk(cudaMemcpy(bc.map_bplane_cell_ghst_d , bcg_h , 
                             sizeof(geom_int)*(bc.iPlanes.size()) , cudaMemcpyHostToDevice));
 
        free(bp_h); 
        free(bc_h);
        free(bcg_h);
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

