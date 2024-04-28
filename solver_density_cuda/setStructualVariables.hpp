#pragma once

#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "bits/stdc++.h"

//void calcGradient(solverConfig& , mesh& , variables& ); 
void setStructualVariables(solverConfig& , cudaConfig& , mesh& , variables& );
