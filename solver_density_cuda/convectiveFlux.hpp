#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>


void calcFluxJacobianMatrix (flow_float ga, flow_float nx, flow_float ny, flow_float nz, 
                             flow_float u , flow_float v , flow_float w,
                             //flow_float p , flow_float ro, flow_float H, flow_float c,
                             flow_float H, flow_float c,
                             Eigen::MatrixXd& A,Eigen::MatrixXd& lambda, 
                             Eigen::MatrixXd& R,Eigen::MatrixXd& Rinv);

void calcRoeAverage (flow_float ga, flow_float roL, flow_float roR, flow_float uL, flow_float uR, 
                     flow_float vL, flow_float vR, 
                     flow_float wL, flow_float wR, 
                     flow_float HL, flow_float HR, flow_float cL, flow_float cR,
                     flow_float& ro_a, flow_float& Ux_a, flow_float& Uy_a, flow_float& Uz_a, flow_float& H_a, flow_float& c_a);

void convectiveFlux(int, solverConfig& , mesh& , variables& , matrix& );


