#pragma once

#include <cuda_runtime.h>

// 非平衡凝縮: 凝縮種ごとの物性相関 (飽和蒸気圧・凝縮相密度・潜熱・表面張力)。
// 物性評価は exp/log で桁が飛ぶため内部は **double** で計算する (flow_float が float でも安全)。
//
// 【相の方針 (重要)】Lin 2014 のモデルは **液相 (過冷却=supercooled liquid)**。非平衡凝縮は
// 三重点以下でも準安定な過冷却液滴を作るため、本実装も**全域で液相フィットを使う** (固相昇華へは
// 切替えない)。r* = 2σ/(ρ_l RT ln(p_v/p_sat)) も液密度 ρ_l を使う。詳細は docs/condensation/theory.md 8 節。
//
// 【低温外挿ガード (重要)】液相フィット (Jacobsen p_sat, 潜熱多項式) は ~40K 以下に外挿すると破綻する
// (p_sat が ~33K で非単調、潜熱が崩落)。有効域は概ね 45–126K。凝縮 ON では潜熱放出で活性域が ~40K 以上に
// 留まるが、過渡/dry セルが低温に落ちても破綻しないよう、物性評価温度を [T_PROP_FLOOR, Tc) にクランプする。
// クランプは「より低温=より低 p_sat」側に倒れる安全側 (過剰核生成を起こさない)。既知の近似として docs に明記。
//
// 【気相 thermo の制約】凝縮ケースの気相は **thermalMethod 0 (熱量的完全気体, cp/γ 一定)** を使うこと。
// NASA-9/CEA は ~200K 未満で無効 (出口 ~27K)。論文も calorically perfect gas。
//
// 種の切替は enum + 係数構造体 CondSpeciesProps で行う。N2 を最初に実装、H2O は係数差し替えで足す。

// 液相フィットの有効下限 (これ未満は外挿破綻するためクランプ)。
#define COND_T_PROP_FLOOR 45.0

enum CondPropModel {
    COND_MODEL_N2  = 0,
    COND_MODEL_H2O = 1,   // Phase 3
};

// 凝縮種の物性パラメータ (device へ value 渡しできる POD)。
struct CondSpeciesProps {
    int    model;   // CondPropModel
    double R;       // 蒸気の気体定数 [J/(kg·K)]
    double cv;      // 蒸気の定積比熱 [J/(kg·K)] (calorically perfect)
    double cp;      // 蒸気の定圧比熱 [J/(kg·K)]
    double Tt;      // 三重点温度 [K] (液/固 切替)
    double Tc;      // 臨界温度 [K]
    double M;       // 分子量 [kg/mol]
};

// N2 既定パラメータ。R=296.8, γ=1.4 → cv=R/(γ-1)=742, cp=γcv=1038.8。
__host__ __device__ inline CondSpeciesProps condProps_N2()
{
    CondSpeciesProps s;
    s.model = COND_MODEL_N2;
    s.R  = 296.8;
    s.cv = 742.0;
    s.cp = 1038.8;
    s.Tt = 63.15;
    s.Tc = 126.192;
    s.M  = 0.0280134;
    return s;
}

// 物性評価温度を液相フィット有効域 [T_PROP_FLOOR, Tc) にクランプ。
__host__ __device__ inline double cond_clamp_Tprop(double T, double Tc)
{
    double Tcl = T;
    if (Tcl < COND_T_PROP_FLOOR) Tcl = COND_T_PROP_FLOOR;
    const double Tmax = Tc - 0.5; // 臨界点直下で発散させない
    if (Tcl > Tmax) Tcl = Tmax;
    return Tcl;
}

// --- N2 飽和蒸気圧 p_sat(T) [Pa] --- 液 (Jacobsen, 式22, atm)。過冷却液として全域使用 (T はクランプ)。
__host__ __device__ inline double n2_psat(double T)
{
    const double Tc = 126.192;
    const double Tcl = cond_clamp_Tprop(T, Tc);
    const double n1 = 8394.409444, n2 = -1890.045259, n3 = -7.282229165;
    const double n4 = 0.01022850966, n5 = 5.556063825e-4, n6 = -5.944544662e-6;
    const double n7 = 2.715433932e-8, n8 = -4.879535904e-11, n9 = 509.5360824;
    const double dTc = Tc - Tcl;
    const double lnP = n1/Tcl + n2 + n3*Tcl + n4*pow(dTc, 1.95)
                     + n5*Tcl*Tcl*Tcl + n6*Tcl*Tcl*Tcl*Tcl + n7*Tcl*Tcl*Tcl*Tcl*Tcl
                     + n8*Tcl*Tcl*Tcl*Tcl*Tcl*Tcl + n9*log(Tcl);
    return exp(lnP) * 101325.0; // atm -> Pa
}

// --- N2 凝縮相 (液) 密度 ρ_l(T) [kg/m³] --- Nowak (式24)。
__host__ __device__ inline double n2_rho_cond(double T)
{
    const double Tc = 126.192, rhoc = 313.3;
    const double Tcl = cond_clamp_Tprop(T, Tc);
    double tau = 1.0 - Tcl/Tc;
    if (tau < 0.0) tau = 0.0;
    const double n1 = 1.48654237, n2 = -0.280476066, n3 = 0.0894143085, n4 = -0.119879866;
    const double lnr = n1*pow(tau, 0.3294) + n2*pow(tau, 4.0/6.0)
                     + n3*pow(tau, 16.0/6.0) + n4*pow(tau, 35.0/6.0);
    return rhoc * exp(lnr);
}

// --- N2 潜熱 (蒸発) L(T) [J/kg] --- 式26 (4次多項式 MJ/kg)。
__host__ __device__ inline double n2_latent(double T)
{
    const double Tcl = cond_clamp_Tprop(T, 126.192);
    const double p1 = -2.137e-8, p2 = 7.18e-6, p3 = -9.142e-4, p4 = 0.05069, p5 = -0.809;
    const double L = p1*Tcl*Tcl*Tcl*Tcl + p2*Tcl*Tcl*Tcl + p3*Tcl*Tcl + p4*Tcl + p5; // MJ/kg
    const double Lj = L * 1.0e6;
    return (Lj > 0.0) ? Lj : 0.0;
}

// --- N2 表面張力 (液) σ(T) [N/m] --- Stansfield (式28, dyn/cm)。
__host__ __device__ inline double n2_sigma(double T)
{
    const double Tc = 126.0, sig0 = 29.06; // dyn/cm
    const double Tcl = cond_clamp_Tprop(T, 126.192);
    double x = 1.0 - Tcl/Tc;
    if (x < 0.0) x = 0.0;
    return sig0 * pow(x, 1.247) * 1.0e-3; // dyn/cm -> N/m
}

// --- 種ディスパッチ (将来 H2O を足す枠) ---
__host__ __device__ inline double cond_psat(const CondSpeciesProps& s, double T)
{
    // 現状 N2 のみ。H2O は Phase 3 で分岐追加。
    return n2_psat(T);
    (void)s;
}
__host__ __device__ inline double cond_rho_cond(const CondSpeciesProps& s, double T)
{
    return n2_rho_cond(T);
    (void)s;
}
__host__ __device__ inline double cond_latent(const CondSpeciesProps& s, double T)
{
    return n2_latent(T);
    (void)s;
}
__host__ __device__ inline double cond_sigma(const CondSpeciesProps& s, double T)
{
    return n2_sigma(T);
    (void)s;
}
