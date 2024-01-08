#include <iostream>
#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

variables::variables(mesh& msh)
{
    for (auto cellValName : cellValNames)
    {
        this->c[cellValName].resize(msh.nCells);
    }
    for (auto planeValName : planeValNames)
    {
        this->p[planeValName].resize(msh.nPlanes);
    }

    //for (geom_int ic=0 ; ic<msh.nCells ; ic++)
    //{
    //    this->c["volume"][ic] = msh.cells[ic].volume;
    //}

    //for (geom_int ic=0 ; ic<msh.nCells ; ic++)
    //{
    //    this->c["svsum_x"][ic] = 0.0;
    //    this->c["svsum_y"][ic] = 0.0;
    //    this->c["svsum_z"][ic] = 0.0;
    //}

    //for (auto ip=0 ; ip<msh.nNormalPlanes ; ip++) 
    //{
    //    geom_int ic0 = msh.planes[ip].iCells[0];
    //    geom_int ic1 = msh.planes[ip].iCells[1];

    //    this->c["svsum_x"][ic0] += msh.planes[ip].surfVect[0];
    //    this->c["svsum_y"][ic0] += msh.planes[ip].surfVect[1];
    //    this->c["svsum_z"][ic0] += msh.planes[ip].surfVect[2];

    //    this->c["svsum_x"][ic1] -= msh.planes[ip].surfVect[0];
    //    this->c["svsum_y"][ic1] -= msh.planes[ip].surfVect[1];
    //    this->c["svsum_z"][ic1] -= msh.planes[ip].surfVect[2];

    //}

    //for (auto ip=msh.nNormalPlanes ; ip<msh.nPlanes ; ip++) 
    //{
    //    geom_int ic0 = msh.planes[ip].iCells[0];
    //    //geom_int ic1 = msh.planes[ip].iCells[1];

    //    this->c["svsum_x"][ic0] += msh.planes[ip].surfVect[0];
    //    this->c["svsum_y"][ic0] += msh.planes[ip].surfVect[1];
    //    this->c["svsum_z"][ic0] += msh.planes[ip].surfVect[2];

    //    //this->c["svsum_x"][ic1] -= msh.planes[ip].surfVect[0];
    //    //this->c["svsum_y"][ic1] -= msh.planes[ip].surfVect[1];
    //    //this->c["svsum_z"][ic1] -= msh.planes[ip].surfVect[2];

    //}
    

}
