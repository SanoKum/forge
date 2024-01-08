#pragma once

#include <vector>
#include "flowFormat.hpp"

void solvePoisson_amgcl_cuda( std::vector<ptrdiff_t>& , std::vector<ptrdiff_t>& ,
                              std::vector<flow_float>&, std::vector<flow_float>& ,
                              std::vector<flow_float>& );
