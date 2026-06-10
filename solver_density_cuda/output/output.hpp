#pragma once

#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

//void outputH5_XDMF(const mesh &msh , const variables &var)
void outputH5_XDMF(const solverConfig& , const mesh& , variables& , const int& );

// detectNaN 診断モード用: 出力間隔ガードを無視して現在の解を prefix 付きで強制ダンプする。
void dumpSolutionH5_force(const solverConfig& , const mesh& , variables& , const int& , const std::string& );


void outputBconds_H5_XDMF(const solverConfig& cfg , mesh& msh , variables& var , const int& iStep);
