#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"
#include "cuda_forge/update_d.cuh"

void updateVariablesOuter(solverConfig&, cudaConfig&, mesh&, variables&, matrix& );
void updateVariablesInner(solverConfig&, cudaConfig&, mesh&, variables&, matrix& );
