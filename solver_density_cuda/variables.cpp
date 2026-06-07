#include <iostream>
#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/calcStructualVariables_d.cuh"



variables::variables()
{
    for (const auto& cellValName : cellValNames)
    {
        this->c.emplace(cellValName, std::vector<flow_float>{});
        this->c_d.emplace(cellValName, nullptr);
    }

    for (const auto& planeValName : planeValNames)
    {
        this->p.emplace(planeValName, std::vector<flow_float>{});
        this->p_d.emplace(planeValName, nullptr);
    }
};

variables::~variables() {
    for (auto& cellValName : cellValNames)
    {
        cudaWrapper::cudaFree_wrapper(this->c_d.at(cellValName));
    }

    for (auto& planeValName : planeValNames)
    {
        cudaWrapper::cudaFree_wrapper(this->p_d.at(planeValName));
    }
}

void variables::allocVariables(const int &useGPU , mesh& msh)
{
    for (auto& cellValName : cellValNames)
    {
        //this->c[cellValName].resize(msh.nCells);
        std::vector<flow_float>& cellValues = this->c.at(cellValName);
        cellValues.resize(msh.nCells_all); // including ghost cells

        if (useGPU == 1) 
        {
            gpuErrchk( cudaMalloc((void**) &(this->c_d.at(cellValName)), (msh.nCells_all)*sizeof(flow_float)) );
        }

    }
    for (auto& planeValName : planeValNames)
    {
        std::vector<flow_float>& planeValues = this->p.at(planeValName);
        planeValues.resize(msh.nPlanes);

        if (useGPU == 1) 
        {
            gpuErrchk( cudaMalloc((void**) &(this->p_d.at(planeValName)), msh.nPlanes*sizeof(flow_float)) );
        }
    }
}



void variables::copyVariables_cell_plane_H2D_all()
{
    for (auto& name : this->cellValNames)
    {
        std::vector<flow_float>& cellValues = this->c.at(name);
        cudaWrapper::cudaMemcpy_H2D_wrapper(cellValues.data() , this->c_d.at(name), cellValues.size());
    }
    for (auto& name : this->planeValNames)
    {
        std::vector<flow_float>& planeValues = this->p.at(name);
        cudaWrapper::cudaMemcpy_H2D_wrapper(planeValues.data() , this->p_d.at(name), planeValues.size());
    }
}

//void variables::copyVariables_cell_H2D(std::string name)
void variables::copyVariables_cell_H2D(std::list<std::string> names)
{
    for (auto& name : names) {
        std::vector<flow_float>& cellValues = this->c.at(name);
        cudaWrapper::cudaMemcpy_H2D_wrapper(cellValues.data() , this->c_d.at(name), cellValues.size());
    }
}
void variables::copyVariables_plane_H2D(std::list<std::string> names)
{
    for (auto& name : names) {
        std::vector<flow_float>& planeValues = this->p.at(name);
        cudaWrapper::cudaMemcpy_H2D_wrapper(planeValues.data() , this->p_d.at(name), planeValues.size());
    }
}

void variables::copyVariables_cell_plane_D2H_all()
{
    for (auto& name : this->cellValNames)
    {
        std::vector<flow_float>& cellValues = this->c.at(name);
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->c_d.at(name), cellValues.data() , cellValues.size());
    }
    for (auto& name : this->planeValNames)
    {
        std::vector<flow_float>& planeValues = this->p.at(name);
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->p_d.at(name), planeValues.data() , planeValues.size());
    }
}

void variables::copyVariables_cell_D2H(std::list<std::string> names)
{
    for (auto& name : names) {
        std::vector<flow_float>& cellValues = this->c.at(name);
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->c_d.at(name), cellValues.data(), cellValues.size());
    }
}

void variables::copyVariables_plane_D2H(std::list<std::string> names)
{
    for (auto& name : names) {
        std::vector<flow_float>& planeValues = this->p.at(name);
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->p_d.at(name), planeValues.data(), planeValues.size());
    }
}

void variables::setStructuralVariables(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh)
{
    if (cfg.gpu==1) {
        variables::setStructuralVariables_d(cfg, cuda_cfg, msh);
        return;
    }

    //geom_float ss;
    std::vector<geom_float> sv(3);

    std::vector<geom_float> dccv(3);
    geom_float dcc;

    std::vector<geom_float> dc1pv(3);
    std::vector<geom_float> dc2pv(3);
    geom_float dc2p;

    std::vector<geom_float> pcent(3);
    std::vector<geom_float> c1cent(3);
    std::vector<geom_float> c2cent(3);

    geom_int ic1;
    geom_int ic2;
    std::vector<flow_float>& fxp  = this->p.at("fx");
    std::vector<flow_float>& dccp = this->p.at("dcc");
    std::vector<flow_float>& vvol = this->c.at("volume");
    std::vector<flow_float>& vccx = this->c.at("ccx");
    std::vector<flow_float>& vccy = this->c.at("ccy");
    std::vector<flow_float>& vccz = this->c.at("ccz");

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

void variables::setStructuralVariables_d(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh )
{
    geom_float* sx;
    geom_float* sy;
    geom_float* sz;
    geom_float* ss;
    geom_float* sx_planar;
    geom_float* sy_planar;
    geom_float* sz_planar;
    geom_float* ss_planar;
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
    sx_planar = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sy_planar = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sz_planar = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    ss_planar = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcx = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcy = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcz = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);

    ccx = (geom_float*)malloc(sizeof(geom_float)*(msh.nCells_all));
    ccy = (geom_float*)malloc(sizeof(geom_float)*(msh.nCells_all));
    ccz = (geom_float*)malloc(sizeof(geom_float)*(msh.nCells_all));
    volume = (geom_float*)malloc(sizeof(geom_float)*msh.nCells_all);
    geom_float* A_planar_h = (geom_float*)malloc(sizeof(geom_float)*msh.nCells_all);


    for (geom_int ip=0; ip<msh.nPlanes; ip++)
    {
        sx[ip] = msh.planes[ip].surfVect[0];
        sy[ip] = msh.planes[ip].surfVect[1];
        sz[ip] = msh.planes[ip].surfVect[2];
        ss[ip] = msh.planes[ip].surfArea;
        sx_planar[ip] = sx[ip];
        sy_planar[ip] = sy[ip];
        sz_planar[ip] = sz[ip];
        ss_planar[ip] = ss[ip];
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
        A_planar_h[ic] = 0.0;
    }

    if (cfg.isAxisymmetric == 1) {
        // B 流儀: 幾何量に r 重み付け、半径方向の圧力ソース用に planar 面積を保存。
        // 軸 (r=0) 上の face で S を厳密に 0 にすると、下流の flux/BC カーネルで
        // n = S/|S| = 0/0 = NaN になるため、極小の r フロアを入れて n の方向を
        // 保ちつつ寄与を実質ゼロにする。
        const geom_float r_floor = (geom_float)1.0e-20;
        for (geom_int ip=0; ip<msh.nPlanes; ip++) {
            const geom_float r_face = (pcy[ip] > r_floor) ? pcy[ip] : r_floor;
            sx[ip] *= r_face;
            sy[ip] *= r_face;
            sz[ip] *= r_face;
            ss[ip] *= r_face;
        }
        for (geom_int ic=0; ic<msh.nCells_all; ic++) {
            const geom_float r_cell = (ccy[ic] > r_floor) ? ccy[ic] : r_floor;
            A_planar_h[ic] = volume[ic];
            volume[ic]     = volume[ic] * r_cell;
        }
    }

    flow_float* p_sx = this->p_d.at("sx");
    flow_float* p_sy = this->p_d.at("sy");
    flow_float* p_sz = this->p_d.at("sz");
    flow_float* p_ss = this->p_d.at("ss");
    flow_float* p_sx_planar = this->p_d.at("sx_planar");
    flow_float* p_sy_planar = this->p_d.at("sy_planar");
    flow_float* p_sz_planar = this->p_d.at("sz_planar");
    flow_float* p_ss_planar = this->p_d.at("ss_planar");
    flow_float* p_pcx = this->p_d.at("pcx");
    flow_float* p_pcy = this->p_d.at("pcy");
    flow_float* p_pcz = this->p_d.at("pcz");
    flow_float* c_ccx = this->c_d.at("ccx");
    flow_float* c_ccy = this->c_d.at("ccy");
    flow_float* c_ccz = this->c_d.at("ccz");
    flow_float* c_volume = this->c_d.at("volume");

    cudaMemcpy(p_sx , sx , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sy , sy , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sz , sz , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_ss , ss , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sx_planar , sx_planar , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sy_planar , sy_planar , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sz_planar , sz_planar , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_ss_planar , ss_planar , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);

    cudaMemcpy(p_pcx , pcx , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_pcy , pcy , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_pcz , pcz , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);

    cudaMemcpy(c_ccx , ccx , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(c_ccy , ccy , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(c_ccz , ccz , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(c_volume , volume , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);

    if (cfg.isAxisymmetric == 1) {
        flow_float* c_A_planar = this->c_d.at("A_planar");
        cudaMemcpy(c_A_planar , A_planar_h , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    }

    calcStructualVariables_d_wrapper(cuda_cfg , msh , *this);

    free(sx) ; free(sy) ; free(sz) ; free(ss);
    free(sx_planar) ; free(sy_planar) ; free(sz_planar) ; free(ss_planar);
    free(pcx); free(pcy); free(pcz);
    free(ccx); free(ccy); free(ccz);
    free(volume);
    free(A_planar_h);

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
    std::vector<geom_float> roK;
    std::vector<geom_float> roOmega;
    const bool has_roK = file.exist("/VALUE/roK");
    const bool has_roOmega = file.exist("/VALUE/roOmega");
    if (has_roK) {
        file.getDataSet("/VALUE/roK").read(roK);
    }
    if (has_roOmega) {
        file.getDataSet("/VALUE/roOmega").read(roOmega);
    }
 
   
    std::vector<flow_float>& v_ro = this->c.at("ro");
    std::vector<flow_float>& v_roUx = this->c.at("roUx");
    std::vector<flow_float>& v_roUy = this->c.at("roUy");
    std::vector<flow_float>& v_roUz = this->c.at("roUz");
    std::vector<flow_float>& v_roe = this->c.at("roe");
    std::vector<flow_float>& v_wall_dist = this->c.at("wall_dist");
    std::vector<flow_float>& v_roK = this->c.at("roK");
    std::vector<flow_float>& v_roOmega = this->c.at("roOmega");

    for (geom_int i=0; i<msh.nCells; i++)
    {
        v_ro[i] = ro[i];
        v_roUx[i] = roUx[i];
        v_roUy[i] = roUy[i];
        v_roUz[i] = roUz[i];
        v_roe[i] = roe[i];
        v_wall_dist[i] = wall_dist[i];
        v_roK[i] = has_roK ? roK[i] : static_cast<flow_float>(0.0);
        v_roOmega[i] = has_roOmega ? roOmega[i] : static_cast<flow_float>(0.0);
    }

    std::list<std::string> names = {"ro", "roUx", "roUy", "roUz", "roe", "wall_dist", "roK", "roOmega"};
    this->copyVariables_cell_H2D(names);
}