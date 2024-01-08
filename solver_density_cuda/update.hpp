#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"

void updateVariablesOuter(solverConfig&, mesh&, variables&, matrix& );
void updateVariablesInner(solverConfig&, mesh&, variables&, matrix& );
