#pragma once

#include <iostream>
#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "cuda_forge/cudaConfig.cuh"
#include "cudaWrapper.cuh"


class fluct_variables {
public:

    // fluctuating velocity
    int fluct_init_flag = 0;
    int knum_max = 200;
    flow_float Lref = 0.51;
    flow_float Re = 5.52e6;

    flow_float Mach = 6.0;
    flow_float ramda = 7.739e-5;
    flow_float IntTurb = 0.5;
    flow_float T0 = 420.0;
    flow_float P0 = 965264.3;

    std::vector<flow_float> ubar;
    std::vector<flow_float> sigma_x;
    std::vector<flow_float> sigma_y;
    std::vector<flow_float> sigma_z;
    std::vector<flow_float> k_vec_x;
    std::vector<flow_float> k_vec_y;
    std::vector<flow_float> k_vec_z;
    std::vector<flow_float> Psi;
 

    flow_float* ubar_d;
    flow_float* sigma_x_d;
    flow_float* sigma_y_d;
    flow_float* sigma_z_d;
    flow_float* k_vec_x_d;
    flow_float* k_vec_y_d;
    flow_float* k_vec_z_d;
 
    flow_float* Psi_d;

    fluct_variables();
    ~fluct_variables();

    void allocVariables();
    void set_fluctVelocity(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var );


};


__global__ void inlet_fluctVelocity_statPress_d 
( 
 // gas 
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,

 // bvar
 flow_float* rob ,
 flow_float* roUxb ,
 flow_float* roUyb ,
 flow_float* roUzb ,
 flow_float* roeb ,
 flow_float* Uxb ,
 flow_float* Uyb ,
 flow_float* Uzb ,
 flow_float* Ux0b ,
 flow_float* Uy0b ,
 flow_float* Uz0b ,
 //flow_float* Ttb ,
 //flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb ,

 flow_float time,
 int knum_max,

 flow_float* ubar,

 flow_float* sigma_x,
 flow_float* sigma_y,
 flow_float* sigma_z,

 flow_float* k_vec_x,
 flow_float* k_vec_y,
 flow_float* k_vec_z,

 flow_float* Psi,

 flow_float* px,
 flow_float* py,
 flow_float* pz

);

void inlet_fluctVelocity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p, fluct_variables& fluct);
