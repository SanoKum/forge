#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

void timeIntegration(int, solverConfig& , mesh& , variables& , matrix& );