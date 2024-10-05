#include <iostream>
#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/calcStructualVariables_d.cuh"



variables::variables() {};

variables::~variables() {
    for (auto& cellValName : cellValNames)
    {
        cudaWrapper::cudaFree_wrapper(this->c_d[cellValName]);
    }

    for (auto& planeValName : planeValNames)
    {
        cudaWrapper::cudaFree_wrapper(this->p_d[planeValName]);
    }
}

void variables::allocVariables(const int &useGPU , mesh& msh)
{
    for (auto& cellValName : cellValNames)
    {
        //this->c[cellValName].resize(msh.nCells);
        this->c[cellValName].resize(msh.nCells_all); // including ghost cells

        if (useGPU == 1) 
        {
            gpuErrchk( cudaMalloc((void**) &(this->c_d[cellValName]), (msh.nCells_all)*sizeof(flow_float)) );
        }

    }
    for (auto& planeValName : planeValNames)
    {
        this->p[planeValName].resize(msh.nPlanes);

        if (useGPU == 1) 
        {
            gpuErrchk( cudaMalloc((void**) &(this->p_d[planeValName]), msh.nPlanes*sizeof(flow_float)) );
        }
    }
}



void variables::copyVariables_cell_plane_H2D_all()
{
    for (auto& name : this->cellValNames)
    {
        cudaWrapper::cudaMemcpy_H2D_wrapper(this->c[name].data() , this->c_d[name], this->c[name].size());
    }
    for (auto& name : this->planeValNames)
    {
        cudaWrapper::cudaMemcpy_H2D_wrapper(this->p[name].data() , this->p_d[name], this->p[name].size());
    }
}

//void variables::copyVariables_cell_H2D(std::string name)
void variables::copyVariables_cell_H2D(std::list<std::string> names)
{
    for (auto& name : names) {
        cudaWrapper::cudaMemcpy_H2D_wrapper(this->c[name].data() , this->c_d[name], this->c[name].size());
    }
}
void variables::copyVariables_plane_H2D(std::list<std::string> names)
{
    for (auto& name : names) {
        cudaWrapper::cudaMemcpy_H2D_wrapper(this->p[name].data() , this->p_d[name], this->c[name].size());
    }
}

void variables::copyVariables_cell_plane_D2H_all()
{
    for (auto& name : this->cellValNames)
    {
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->c_d[name], this->c[name].data() , this->c[name].size());
    }
    for (auto& name : this->planeValNames)
    {
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->p_d[name], this->p[name].data() , this->p[name].size());
    }
}

void variables::copyVariables_cell_D2H(std::list<std::string> names)
{
    for (auto& name : names) {
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->c_d[name], this->c[name].data(), this->c[name].size());
    }
}

void variables::copyVariables_plane_D2H(std::list<std::string> names)
{
    for (auto& name : names) {
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->p_d[name], this->p[name].data(), this->p[name].size());
    }
}

void variables::setStructuralVariables(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh)
{
    if (cfg.gpu==1) {
        variables::setStructuralVariables_d(cuda_cfg, msh);
        return;
    }

    //geom_float ss;
    std::vector<geom_float> sv(3);

    std::vector<geom_float> dccv(3);
    geom_float dcc;

    std::vector<geom_float> dc1pv(3);
    geom_float dc1p;

    std::vector<geom_float> dc2pv(3);
    geom_float dc2p;

    std::vector<geom_float> pcent(3);
    std::vector<geom_float> c1cent(3);
    std::vector<geom_float> c2cent(3);

    geom_float volume;

    geom_int ic1;
    geom_int ic2;

    geom_float f;
    std::vector<flow_float>& fxp  = this->p["fx"]; 
    std::vector<flow_float>& dccp = this->p["dcc"];
    std::vector<flow_float>& vvol = this->c["volume"];
    std::vector<flow_float>& vccx = this->c["ccx"];
    std::vector<flow_float>& vccy = this->c["ccy"];
    std::vector<flow_float>& vccz = this->c["ccz"];

    for (geom_int ip=0 ; ip<msh.nPlanes ; ip++)
    {
        ic1     = msh.planes[ip].iCells[0];
        ic2     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;

        pcent   = msh.planes[ip].centCoords;

        c1cent  = msh.cells[ic1].centCoords;
        c2cent  = msh.cells[ic2].centCoords;

        dccv[0] = c2cent[0] - c1cent[0];
        dccv[1] = c2cent[1] - c1cent[1];
        dccv[2] = c2cent[2] - c1cent[2];
        dcc     = sqrt( pow(dccv[0], 2.0) + pow(dccv[1], 2.0) + pow(dccv[2], 2.0));

        dc2pv[0] = pcent[0] - c2cent[0];
        dc2pv[1] = pcent[1] - c2cent[1];
        dc2pv[2] = pcent[2] - c2cent[2];
        dc2p     = sqrt( pow(dc2pv[0], 2.0) + pow(dc2pv[1], 2.0) + pow(dc2pv[2], 2.0));

        fxp[ip]  = dc2p/dcc;
        dccp[ip] = dcc;
    }

    // cell
    for (geom_int ic=0 ; ic<msh.nCells ; ic++)
    {
        c1cent   = msh.cells[ic].centCoords;
        vvol[ic] = msh.cells[ic].volume;
        vccx[ic] = c1cent[0];
        vccy[ic] = c1cent[1];
        vccz[ic] = c1cent[2];
    }

}

void variables::setStructuralVariables_d(cudaConfig& cuda_cfg , mesh& msh )
{
    geom_float* sx;
    geom_float* sy;
    geom_float* sz;
    geom_float* ss;
    geom_float* pcx;
    geom_float* pcy;
    geom_float* pcz;
    geom_float* ccx;
    geom_float* ccy;
    geom_float* ccz;
    geom_float* volume;

    sx = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sy = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sz = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    ss = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcx = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcy = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcz = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);

    ccx = (geom_float*)malloc(sizeof(geom_float)*(msh.nCells_all));
    ccy = (geom_float*)malloc(sizeof(geom_float)*(msh.nCells_all));
    ccz = (geom_float*)malloc(sizeof(geom_float)*(msh.nCells_all));
    volume = (geom_float*)malloc(sizeof(geom_float)*msh.nCells_all);


    for (geom_int ip=0; ip<msh.nPlanes; ip++)
    {
        sx[ip] = msh.planes[ip].surfVect[0];
        sy[ip] = msh.planes[ip].surfVect[1];
        sz[ip] = msh.planes[ip].surfVect[2];
        ss[ip] = msh.planes[ip].surfArea;
        pcx[ip] = msh.planes[ip].centCoords[0];
        pcy[ip] = msh.planes[ip].centCoords[1];
        pcz[ip] = msh.planes[ip].centCoords[2];
    }

    for (geom_int ic=0; ic<msh.nCells_all; ic++)
    {
        ccx[ic] = msh.cells[ic].centCoords[0];
        ccy[ic] = msh.cells[ic].centCoords[1];
        ccz[ic] = msh.cells[ic].centCoords[2];
        volume[ic] = msh.cells[ic].volume;
    }

    cudaMemcpy(this->p_d["sx"] , sx , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->p_d["sy"] , sy , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->p_d["sz"] , sz , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->p_d["ss"] , ss , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);

    cudaMemcpy(this->p_d["pcx"] , pcx , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->p_d["pcy"] , pcy , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->p_d["pcz"] , pcz , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);

    cudaMemcpy(this->c_d["ccx"] , ccx , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->c_d["ccy"] , ccy , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->c_d["ccz"] , ccz , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->c_d["volume"] , volume , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);

    calcStructualVariables_d_wrapper(cuda_cfg , msh , *this);

    free(sx) ; free(sy) ; free(sz) ; free(ss);
    free(pcx); free(pcy); free(pcz);
    free(ccx); free(ccy); free(ccz);
    free(volume); 

}

void variables::readValueHDF5(std::string fname , mesh& msh)
{
    HighFive::File file(fname, HighFive::File::ReadOnly);

    // read basic 
    HighFive::Group group = file.getGroup("/VALUE");

    // nodes
    std::vector<geom_float> ro;
    file.getDataSet("/VALUE/ro").read(ro);
    std::vector<geom_float> roUx;
    file.getDataSet("/VALUE/roUx").read(roUx);
    std::vector<geom_float> roUy;
    file.getDataSet("/VALUE/roUy").read(roUy);
    std::vector<geom_float> roUz;
    file.getDataSet("/VALUE/roUz").read(roUz);
    std::vector<geom_float> roe;
    file.getDataSet("/VALUE/roe").read(roe);
    std::vector<geom_float> wall_dist;
    file.getDataSet("/VALUE/wall_dist").read(wall_dist);
 
   
    for (geom_int i=0; i<msh.nCells; i++)
    {
        this->c["ro"][i] = ro[i];
        this->c["roUx"][i] = roUx[i];
        this->c["roUy"][i] = roUy[i];
        this->c["roUz"][i] = roUz[i];
        this->c["roe"][i] = roe[i];
        this->c["wall_dist"][i] = wall_dist[i];
    }

    std::list<std::string> names = {"ro", "roUx", "roUy", "roUz", "roe", "wall_dist"};
    this->copyVariables_cell_H2D(names);
}