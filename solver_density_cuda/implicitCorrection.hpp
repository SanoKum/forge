#pragma once

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"

void implicitCorrection(int, solverConfig& , mesh& , variables& , matrix& );