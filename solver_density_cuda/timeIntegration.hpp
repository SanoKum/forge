#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"
#include "cuda_forge/timeIntegration_d.cuh"
#include <algorithm>
#include <vector>
//#include <Eigen/Dense>
#include <iostream>
#include <cmath>

void timeIntegration(int, solverConfig& , cudaConfig &, mesh& , variables& , matrix& );