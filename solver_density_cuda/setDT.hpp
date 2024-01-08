#pragma once

#include <algorithm>
#include <cmath>
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"

void setDT(solverConfig& , mesh& , variables& );
