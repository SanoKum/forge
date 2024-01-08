#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "input/solverConfig.hpp"

void updateVariablesForNextLoop(solverConfig&, mesh&, variables&, matrix& );
