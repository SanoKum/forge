#pragma once

#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

//void outputH5_XDMF(const mesh &msh , const variables &var)
void outputH5_XDMF(const solverConfig& , const mesh& , const variables& , const int& );


