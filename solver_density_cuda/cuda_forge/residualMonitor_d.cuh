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

// ---- device 常駐の残差 RMS / NaN reducer (per-step host 同期を避ける) ----
// 詳細: plans/accepted/architecture-residual-monitor-async.md
struct DeviceResidualReducer {
	int nVar = 0;
	geom_int nCells = 0;
	const flow_float** d_ptrs = nullptr; // [nVar] 対象変数の device ポインタ配列 (device 上)
	double* d_sumsq = nullptr;           // [nVar] sum-of-squares scratch (device, double 累算)
	int* d_flag = nullptr;               // [1] 非有限フラグ (device, detectNaN 用)
};

// names のうち var.c_d に存在する変数の device ポインタを集めて reducer を構築する (device メモリ確保)。
// 実際に採用した名前を resolved_names に返す (存在しない名前は除外し nVar を詰める)。
DeviceResidualReducer makeDeviceResidualReducer(
	mesh& msh, variables& var,
	const std::vector<std::string>& names,
	std::vector<std::string>& resolved_names);
void freeDeviceResidualReducer(DeviceResidualReducer& r);

// device の rms バッファ (capacity*nVar 要素) を確保 / 解放。
flow_float* allocDeviceRmsBuffer(int capacity, int nVar);
void freeDeviceRmsBuffer(flow_float* d_buf);

// 現在の状態の sqrt(sum_sq/nCells) を d_buf[slot*nVar ..] へ async に書き込む (同期しない)。
void reduceResidualToSlot(const DeviceResidualReducer& r, flow_float* d_buf, int slot);
// d_buf 先頭 count*nVar 要素を host へ一括コピー (同期する)。
void downloadRmsBuffer(const flow_float* d_buf, flow_float* host_out, int count, int nVar);

// reducer の全変数に非有限値があれば d_flag を 1 にする (fused, 同期しない; 事前に flag=0 する)。
void scanNonFiniteToFlag(const DeviceResidualReducer& r);
// d_flag を host へ読み出す (同期する)。非0 なら非有限あり。
int downloadNonFiniteFlag(const DeviceResidualReducer& r);