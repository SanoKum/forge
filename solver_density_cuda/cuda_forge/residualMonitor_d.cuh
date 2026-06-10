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

// 指定した保存量/派生量のセル配列 (内部セル msh.nCells) に NaN/Inf が 1 つでもあれば true を返し、
// 最初に該当した変数名を offending_name に格納する。detectNaN 診断モードから毎ステップ呼ばれる。
bool hasNonFiniteCellValue_d(
	mesh& msh,
	variables& var,
	const std::vector<std::string>& variable_names,
	std::string& offending_name);