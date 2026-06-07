#pragma once

#include <array>
#include <string>
#include <vector>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

std::array<flow_float, 5> gatherResidualRms_d(mesh& msh, variables& var);
std::array<flow_float, 5> gatherVariableRms_d(
	mesh& msh,
	variables& var,
	const std::array<const char*, 5>& variable_names);
std::vector<flow_float> gatherVariableRms_d(
	mesh& msh,
	variables& var,
	const std::vector<std::string>& variable_names);