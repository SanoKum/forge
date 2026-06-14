#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <stdio.h>                                                                                       
#include <fstream>
#include <string>
#include <time.h>
#include <limits>
#include <cstdio>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "input/solverConfig.hpp"
#include "input/setInitial.hpp"

#include "mesh/gmshReader.hpp"
#include "output/output.hpp"
#include "boundaryCond.hpp"

#include "setStructualVariables.hpp"

#include "gradient.hpp"

#include "dependentVariables.hpp"
//#include "solvePoisson_amgx.hpp"

#include "convectiveFlux.hpp"
//#include "timeIntegration.hpp"
#include "update.hpp"

#include "common/stringUtil.hpp"
#include "common/vectorUtil.hpp"

//#include "setDT.hpp"

// cuda
#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/cudaConfig.cuh"

#include "cuda_forge/calcGradient_d.cuh"
#include "cuda_forge/convectiveFlux_d.cuh"
#include "cuda_forge/ransTransport_d.cuh"
#include "cuda_forge/ransSource_d.cuh"
#include "cuda_forge/speciesTransport_d.cuh"
#include "cuda_forge/condensationTransport_d.cuh"
#include "cuda_forge/viscousFlux_d.cuh"
#include "cuda_forge/updateCenterVelocity_d.cuh"
#include "cuda_forge/interpVelocity_c2p_d.cuh"
#include "cuda_forge/timeIntegration_d.cuh"
#include "cuda_forge/limiter_d.cuh"
#include "cuda_forge/ducrosSensor_d.cuh"
#include "cuda_forge/axisymmetricSource_d.cuh"
#include "cuda_forge/turbulent_viscosity_d.cuh"
#include "cuda_forge/residualMonitor_d.cuh"

#include "cuda_forge/fluct_variables_d.cuh"
#include "cuda_forge/gasProperties_d.cuh"
#include "cuda_forge/thermo_d.cuh"

#include "probe/point_probes.cuh"
#include "cuda_forge/setDT_d.cuh"

#include <cuda_runtime.h>


#define CHECK_LAST_CUDA_ERROR() checkLast(__FILE__, __LINE__)
void checkLast(const char* const file, const int line)
{
    cudaError_t err{cudaGetLastError()};
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << std::endl;
        // We don't exit when we encounter CUDA errors in this example.
        std::exit(EXIT_FAILURE);
    }
}

namespace {

enum class ProfileSection {
    UpdateInner,
    DependentVariables,
    GasProperties,
    ApplyBconds,
    CalcGradient,
    Limiter,
    DucrosSensor,
    TurbulenceModel,
    ConvectiveFlux,
    AxisymmetricSource,
    ViscousFlux,
    ImplicitCorrection,
    TimeIntegration,
    UpdateOuter,
    WriteOutputs,
    SetDt,
    StepTotal,
    Count
};

constexpr std::size_t kProfileSectionCount = static_cast<std::size_t>(ProfileSection::Count);
enum class ResidualPhase {
    OuterBegin,
    InnerBegin,
    InnerIter,
    OuterEnd,
};

constexpr std::array<const char*, 5> kResidualEquationNames = {
    "ro", "roUx", "roUy", "roUz", "roe"
};

constexpr std::array<const char*, 2> kScalarResidualEquationNames = {
    "roK", "roOmega"
};

bool scalarResidualEnabled(const solverConfig& cfg)
{
    return cfg.LESorRANS == 2 && cfg.RANSmodel == 1;
}

// 非平衡凝縮 (Phase 1): 4 モーメントの輸送が有効か。
bool condensationEnabled(const solverConfig& cfg)
{
    return cfg.condensation == 1 && cfg.nCondSpecies >= 1;
}

// 凝縮モーメントの保存量名 (登録順)。variables::registerCondensation と同じ命名規約で cfg のみから導出。
// 1 凝縮種あたり 4 本: rog_{s}, roQ2_{s}, roQ1_{s}, roQ0_{s}。
std::vector<std::string> condMomentConsNames(const solverConfig& cfg)
{
    std::vector<std::string> names;
    if (!condensationEnabled(cfg)) return names;
    const std::array<const char*, 4> bases = {"g", "Q2", "Q1", "Q0"};
    for (int s = 0; s < cfg.nCondSpecies; ++s) {
        for (const auto* b : bases) {
            names.emplace_back(std::string("ro") + b + "_" + std::to_string(s));
        }
    }
    return names;
}

std::vector<std::string> residualEquationNames(const solverConfig& cfg)
{
    std::vector<std::string> names;
    names.reserve(kResidualEquationNames.size() + kScalarResidualEquationNames.size());

    for (const auto* name : kResidualEquationNames) {
        names.emplace_back(name);
    }

    if (scalarResidualEnabled(cfg)) {
        for (const auto* name : kScalarResidualEquationNames) {
            names.emplace_back(name);
        }
    }

    if (condensationEnabled(cfg)) {
        for (const auto& name : condMomentConsNames(cfg)) {
            names.emplace_back(name);
        }
    }

    return names;
}

struct ResidualSnapshot {
    std::vector<flow_float> rms;
};

struct CorrectionSnapshot {
    std::array<flow_float, kResidualEquationNames.size()> rms{};
};

struct ImplicitDiagSnapshot {
    double mean_pseudo_diag = 0.0;
    double mean_face_diag = 0.0;
    double mean_pseudo_ratio = 0.0;
    double min_pseudo_ratio = 0.0;
    double max_pseudo_ratio = 0.0;
    double mean_face_to_pseudo = 0.0;
    double mean_dt_local = 0.0;
};

class ResidualCsvLogger {
public:
    ResidualCsvLogger(const std::string& file_name, const solverConfig& cfg)
        : residual_names_(residualEquationNames(cfg)),
          stream_(file_name, std::ios::out | std::ios::trunc)
    {
        if (!stream_) {
            throw std::runtime_error("Failed to open residual log file: " + file_name);
        }

        stream_ << "step,inner,phase";
        for (const auto& name : residual_names_) {
            stream_ << ",rms_" << name;
        }
        for (const auto* name : kResidualEquationNames) {
            stream_ << ",rms_dq_" << name;
        }
        stream_ << "\n";

        last_snapshot_.rms.assign(residual_names_.size(), static_cast<flow_float>(0.0));
    }

    void logOuterBegin(
        int step,
        int inner,
        const ResidualSnapshot& snapshot,
        const CorrectionSnapshot& correction_snapshot = {})
    {
        last_snapshot_ = snapshot;
        last_correction_snapshot_ = correction_snapshot;

        writeRow(step, -1, ResidualPhase::OuterBegin, snapshot, correction_snapshot);
        writeRow(step, inner, ResidualPhase::InnerBegin, snapshot, correction_snapshot);
    }

    void logInnerIter(
        int step,
        int inner_index,
        const ResidualSnapshot& snapshot,
        const CorrectionSnapshot& correction_snapshot = {})
    {
        last_snapshot_ = snapshot;
        last_correction_snapshot_ = correction_snapshot;
        writeRow(step, inner_index, ResidualPhase::InnerIter, snapshot, correction_snapshot);
    }

    void logOuterEnd(
        int step,
        const ResidualSnapshot* snapshot = nullptr,
        const CorrectionSnapshot* correction_snapshot = nullptr)
    {
        if (snapshot != nullptr) {
            last_snapshot_ = *snapshot;
        }
        if (correction_snapshot != nullptr) {
            last_correction_snapshot_ = *correction_snapshot;
        }
        writeRow(step, -1, ResidualPhase::OuterEnd, last_snapshot_, last_correction_snapshot_);
    }

private:
    static const char* phaseName(ResidualPhase phase)
    {
        switch (phase) {
            case ResidualPhase::OuterBegin: return "outer_begin";
            case ResidualPhase::InnerBegin: return "inner_begin";
            case ResidualPhase::InnerIter: return "inner_iter";
            case ResidualPhase::OuterEnd: return "outer_end";
        }

        return "unknown";
    }

    void writeRow(
        int step,
        int inner,
        ResidualPhase phase,
        const ResidualSnapshot& snapshot,
        const CorrectionSnapshot& correction_snapshot)
    {
        stream_ << step << ',' << inner << ',' << phaseName(phase);
        for (flow_float value : snapshot.rms) {
            stream_ << ',' << std::setprecision(16) << value;
        }
        for (flow_float value : correction_snapshot.rms) {
            stream_ << ',' << std::setprecision(16) << value;
        }
        stream_ << "\n";
        stream_.flush();
    }

    std::ofstream stream_;
    std::vector<std::string> residual_names_;
    ResidualSnapshot last_snapshot_{};
    CorrectionSnapshot last_correction_snapshot_{};
};

class ImplicitDiagLogger {
public:
    ImplicitDiagLogger()
    {
        const char* path = std::getenv("FORGE_IMPLICIT_DIAG_CSV");
        if (path == nullptr || path[0] == '\0') {
            return;
        }

        stream_.open(path, std::ios::out | std::ios::trunc);
        if (!stream_) {
            throw std::runtime_error("Failed to open implicit diag log file: " + std::string(path));
        }

        enabled_ = true;
        stream_ << "step,inner,mean_pseudo_diag,mean_face_diag,mean_pseudo_ratio,min_pseudo_ratio,max_pseudo_ratio,mean_face_to_pseudo,mean_dt_local\n";
    }

    bool enabled() const
    {
        return enabled_;
    }

    void log(int step, int inner, const ImplicitDiagSnapshot& snapshot)
    {
        if (!enabled_) {
            return;
        }

        stream_ << step << ',' << inner
                << ',' << std::setprecision(16) << snapshot.mean_pseudo_diag
                << ',' << snapshot.mean_face_diag
                << ',' << snapshot.mean_pseudo_ratio
                << ',' << snapshot.min_pseudo_ratio
                << ',' << snapshot.max_pseudo_ratio
                << ',' << snapshot.mean_face_to_pseudo
                << ',' << snapshot.mean_dt_local
                << "\n";
        stream_.flush();
    }

private:
    bool enabled_ = false;
    std::ofstream stream_;
};

ResidualSnapshot gatherResidualSnapshot(solverConfig& cfg, mesh& msh, variables& var)
{
    std::vector<std::string> variable_names = {
        "res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe"
    };
    if (scalarResidualEnabled(cfg)) {
        variable_names.emplace_back("res_roK");
        variable_names.emplace_back("res_roOmega");
    }
    if (condensationEnabled(cfg)) {
        for (const auto& name : condMomentConsNames(cfg)) {
            variable_names.emplace_back("res_" + name);
        }
    }

    ResidualSnapshot snapshot;
    if (cfg.gpu == 1) {
        snapshot.rms = gatherVariableRms_d(msh, var, variable_names);
        return snapshot;
    }

    snapshot.rms.assign(variable_names.size(), static_cast<flow_float>(0.0));

    for (std::size_t i = 0; i < variable_names.size(); ++i) {
        const std::vector<flow_float>& values = var.c.at(variable_names[i]);
        flow_float sum_sq = 0.0;
        for (geom_int ic = 0; ic < msh.nCells; ++ic) {
            const flow_float value = values[ic];
            sum_sq += value * value;
        }
        snapshot.rms[i] = std::sqrt(sum_sq / static_cast<flow_float>(std::max<geom_int>(msh.nCells, 1)));
    }

    return snapshot;
}

// frozen scalar 陰解法（次フェーズ）の inner 補正ログ用。現状の経路では未使用だが温存する。
[[maybe_unused]] CorrectionSnapshot gatherCorrectionSnapshot(solverConfig& cfg, mesh& msh, variables& var)
{
    CorrectionSnapshot snapshot;

    if (cfg.timeIntegration != 11) {
        return snapshot;
    }

    const std::array<const char*, kResidualEquationNames.size()> variable_names = {
        cfg.blockDPLUR == 1 ? "dq_block_old_0" : "dq_ro_old",
        cfg.blockDPLUR == 1 ? "dq_block_old_1" : "dq_roUx_old",
        cfg.blockDPLUR == 1 ? "dq_block_old_2" : "dq_roUy_old",
        cfg.blockDPLUR == 1 ? "dq_block_old_3" : "dq_roUz_old",
        cfg.blockDPLUR == 1 ? "dq_block_old_4" : "dq_roe_old"
    };

    if (cfg.gpu == 1) {
        snapshot.rms = gatherVariableRms_d(msh, var, variable_names);
        return snapshot;
    }

    for (std::size_t i = 0; i < variable_names.size(); ++i) {
        const std::vector<flow_float>& values = var.c.at(variable_names[i]);
        flow_float sum_sq = 0.0;
        for (geom_int ic = 0; ic < msh.nCells; ++ic) {
            const flow_float value = values[ic];
            sum_sq += value * value;
        }
        snapshot.rms[i] = std::sqrt(sum_sq / static_cast<flow_float>(std::max<geom_int>(msh.nCells, 1)));
    }

    return snapshot;
}

ImplicitDiagSnapshot gatherImplicitDiagSnapshot(solverConfig& cfg, mesh& msh, variables& var)
{
    const std::list<std::string> variable_names = {
        "ro", "Ux", "Uy", "Uz", "sonic", "vis_turb", "dt_local"
    };

    if (cfg.gpu == 1) {
        var.copyVariables_cell_D2H(variable_names);
    }

    const std::vector<flow_float>& ro = var.c.at("ro");
    const std::vector<flow_float>& ux = var.c.at("Ux");
    const std::vector<flow_float>& uy = var.c.at("Uy");
    const std::vector<flow_float>& uz = var.c.at("Uz");
    const std::vector<flow_float>& sonic = var.c.at("sonic");
    const std::vector<flow_float>& vis_turb = var.c.at("vis_turb");
    const std::vector<flow_float>& dt_local = var.c.at("dt_local");

    ImplicitDiagSnapshot snapshot;
    snapshot.min_pseudo_ratio = std::numeric_limits<double>::infinity();

    for (geom_int ic = 0; ic < msh.nCells; ++ic) {
        const flow_float volume = static_cast<flow_float>(msh.cells[ic].volume);
        const flow_float dt_l = std::max(dt_local[ic], static_cast<flow_float>(1.0e-30));
        const flow_float density = std::max(ro[ic], static_cast<flow_float>(1.0e-30));
        const flow_float nu_eff = (cfg.visc + std::max(vis_turb[ic], static_cast<flow_float>(0.0))) / density;
        const flow_float local_sonic = std::max(sonic[ic], static_cast<flow_float>(0.0));

        double face_diag_sum = 0.0;
        for (const geom_int ip : msh.cells[ic].iPlanes) {
            const plane& pln = msh.planes[ip];
            const flow_float sx = static_cast<flow_float>(pln.surfVect[0]);
            const flow_float sy = static_cast<flow_float>(pln.surfVect[1]);
            const flow_float sz = static_cast<flow_float>(pln.surfVect[2]);
            const flow_float face_area = std::max(static_cast<flow_float>(pln.surfArea), static_cast<flow_float>(1.0e-30));

            const flow_float advective_radius = std::fabs(
                ux[ic] * sx + uy[ic] * sy + uz[ic] * sz
            ) / face_area + local_sonic;

            const geom_int ic0 = pln.iCells[0];
            const geom_int ic1 = pln.iCells[1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;

            const geom_float dcc_x = msh.cells[other_ic].centCoords[0] - msh.cells[ic].centCoords[0];
            const geom_float dcc_y = msh.cells[other_ic].centCoords[1] - msh.cells[ic].centCoords[1];
            const geom_float dcc_z = msh.cells[other_ic].centCoords[2] - msh.cells[ic].centCoords[2];
            const flow_float dcc = std::max(
                static_cast<flow_float>(std::sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z)),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float dcc_dot_s = std::max(
                static_cast<flow_float>(std::fabs(dcc_x * sx + dcc_y * sy + dcc_z * sz)),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float delta = std::max(
                dcc * face_area * face_area / dcc_dot_s,
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float viscous_radius = static_cast<flow_float>(2.0) * nu_eff / delta;
            const flow_float face_coeff = face_area * (advective_radius + viscous_radius);
            face_diag_sum += static_cast<double>(face_coeff);
        }

        const double pseudo_diag = static_cast<double>(volume / dt_l);
        const double total_diag = std::max(pseudo_diag + face_diag_sum, 1.0e-30);
        const double pseudo_ratio = pseudo_diag / total_diag;
        const double face_to_pseudo = face_diag_sum / std::max(pseudo_diag, 1.0e-30);

        snapshot.mean_pseudo_diag += pseudo_diag;
        snapshot.mean_face_diag += face_diag_sum;
        snapshot.mean_pseudo_ratio += pseudo_ratio;
        snapshot.min_pseudo_ratio = std::min(snapshot.min_pseudo_ratio, pseudo_ratio);
        snapshot.max_pseudo_ratio = std::max(snapshot.max_pseudo_ratio, pseudo_ratio);
        snapshot.mean_face_to_pseudo += face_to_pseudo;
        snapshot.mean_dt_local += static_cast<double>(dt_l);
    }

    const double inv_n_cells = 1.0 / static_cast<double>(std::max<geom_int>(msh.nCells, 1));
    snapshot.mean_pseudo_diag *= inv_n_cells;
    snapshot.mean_face_diag *= inv_n_cells;
    snapshot.mean_pseudo_ratio *= inv_n_cells;
    snapshot.mean_face_to_pseudo *= inv_n_cells;
    snapshot.mean_dt_local *= inv_n_cells;
    if (!std::isfinite(snapshot.min_pseudo_ratio)) {
        snapshot.min_pseudo_ratio = 0.0;
    }

    return snapshot;
}

struct ProfileStat {
    double total_ms = 0.0;
    std::size_t calls = 0;
};

class RuntimeProfiler {
public:
    RuntimeProfiler()
        : enabled_(readFlag("FORGE_PROFILE")),
          verbose_(readFlag("FORGE_PROFILE_VERBOSE"))
    {
        if (enabled_) {
            CHECK_CUDA_ERROR(cudaEventCreate(&start_event_));
            CHECK_CUDA_ERROR(cudaEventCreate(&stop_event_));
        }
    }

    ~RuntimeProfiler()
    {
        if (!enabled_) {
            return;
        }

        cudaEventDestroy(start_event_);
        cudaEventDestroy(stop_event_);
    }

    bool enabled() const
    {
        return enabled_;
    }

    template <typename Func>
    void measureWall(ProfileSection section, Func&& func)
    {
        if (!enabled_) {
            func();
            return;
        }

        const auto start = std::chrono::steady_clock::now();
        func();
        const auto end = std::chrono::steady_clock::now();
        const double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
        add(section, elapsed_ms);
    }

    template <typename Func>
    void measureCuda(ProfileSection section, Func&& func)
    {
        if (!enabled_) {
            func();
            return;
        }

        CHECK_CUDA_ERROR(cudaEventRecord(start_event_));
        func();
        CHECK_CUDA_ERROR(cudaEventRecord(stop_event_));
        CHECK_CUDA_ERROR(cudaEventSynchronize(stop_event_));

        float elapsed_ms = 0.0f;
        CHECK_CUDA_ERROR(cudaEventElapsedTime(&elapsed_ms, start_event_, stop_event_));
        add(section, static_cast<double>(elapsed_ms));
    }

    void printSummary() const
    {
        if (!enabled_) {
            return;
        }

        const ProfileStat& total = stats_[index(ProfileSection::StepTotal)];
        std::cout << "\n=== Runtime Profile Summary ===\n";
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "Profiled steps: " << total.calls << "\n";
        if (total.calls == 0 || total.total_ms <= 0.0) {
            std::cout << "No profiled steps were recorded.\n";
            return;
        }

        for (std::size_t i = 0; i < kProfileSectionCount; ++i) {
            const auto section = static_cast<ProfileSection>(i);
            if (section == ProfileSection::StepTotal) {
                continue;
            }

            const ProfileStat& stat = stats_[i];
            if (stat.calls == 0) {
                continue;
            }

            const double avg_ms = stat.total_ms / static_cast<double>(stat.calls);
            const double share_pct = (stat.total_ms / total.total_ms) * 100.0;
            std::cout << std::setw(20) << std::left << sectionName(section)
                      << " total_ms=" << std::setw(12) << stat.total_ms
                      << " avg_ms=" << std::setw(10) << avg_ms
                      << " share=" << std::setw(8) << share_pct << "%"
                      << " calls=" << stat.calls << "\n";
        }

        if (verbose_) {
            std::cout << "Profiling source: CPU sections use steady_clock, GPU wrapper sections use CUDA events on the default stream.\n";
        }
    }

private:
    static bool readFlag(const char* name)
    {
        const char* value = std::getenv(name);
        if (value == nullptr) {
            return false;
        }

        return std::string(value) != "0";
    }

    static constexpr std::size_t index(ProfileSection section)
    {
        return static_cast<std::size_t>(section);
    }

    void add(ProfileSection section, double elapsed_ms)
    {
        ProfileStat& stat = stats_[index(section)];
        stat.total_ms += elapsed_ms;
        stat.calls += 1;
    }

    static const char* sectionName(ProfileSection section)
    {
        switch (section) {
            case ProfileSection::UpdateInner: return "update_inner";
            case ProfileSection::DependentVariables: return "dependent_vars";
            case ProfileSection::GasProperties: return "gas_properties";
            case ProfileSection::ApplyBconds: return "apply_bconds";
            case ProfileSection::CalcGradient: return "calc_gradient";
            case ProfileSection::Limiter: return "limiter";
            case ProfileSection::DucrosSensor: return "ducros_sensor";
            case ProfileSection::TurbulenceModel: return "turbulence_model";
            case ProfileSection::ConvectiveFlux: return "convective_flux";
            case ProfileSection::AxisymmetricSource: return "axisym_source";
            case ProfileSection::ViscousFlux: return "viscous_flux";
            case ProfileSection::ImplicitCorrection: return "implicit_corr";
            case ProfileSection::TimeIntegration: return "time_integr";
            case ProfileSection::UpdateOuter: return "update_outer";
            case ProfileSection::WriteOutputs: return "write_outputs";
            case ProfileSection::SetDt: return "set_dt";
            case ProfileSection::StepTotal: return "step_total";
            case ProfileSection::Count: return "count";
        }

        return "unknown";
    }

    bool enabled_ = false;
    bool verbose_ = false;
    cudaEvent_t start_event_{};
    cudaEvent_t stop_event_{};
    std::array<ProfileStat, kProfileSectionCount> stats_{};
};

cudaConfig initializeSimulation(
    solverConfig& cfg,
    mesh& msh,
    matrix& mat_ns,
    variables& var,
    fluct_variables& fluct,
    point_probes& pprobes)
{
    cout << "Read Solver Config \n";
    cfg.read("solverConfig.yaml");

    cout << "Init Thermo DB \n";
    thermo_init_db(cfg);   // NASA-9/LJ 化学種 DB を構築し device へアップロード (thermalMethod==2 用)

    cout << "Read Mesh \n";
    if (cfg.meshFormat == "hdf5") {
        msh.readMesh(cfg.meshFileName);
    } else {
        cerr << "Error unknown mesh format: " << cfg.meshFormat << endl;
        std::exit(EXIT_FAILURE);
    }

    cout << "Init Matrix (but not used now) \n";
    mat_ns.initMatrix(msh);

    cout << "Read Boundary Conditions \n";
    readBcondConfig(cfg , msh.bconds);

    // 化学種変数を登録 (allocVariables より前)。nSpecies<=1 では no-op。
    var.registerSpecies(cfg.nSpecies);

    // 非平衡凝縮モーメント変数を登録 (allocVariables より前)。condensation==0 では no-op。
    var.registerCondensation(cfg.nCondSpecies);

    var.allocVariables(cfg.gpu , msh);

    // device roY[] ポインタ配列を構築 (c_d 確保後, dependentVariables より前)。
    speciesInit_d(cfg , var);

    cout << "Read Initial Values \n";
    var.readValueHDF5(cfg.valueFileName , msh);

    cout << "Set mesh connection map for cuda \n";
    cudaConfig cuda_cfg(msh);
    msh.setMeshMap_d();
    msh.setPeriodicPartner();

    var.setStructuralVariables(cfg , cuda_cfg , msh);
    speciesPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);  // Y_s = ρY_s/ρ (roY を読込済)
    condensationPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);  // φ = ρφ/ρ (液相モーメント読込済)
    dependentVariables(cfg , cuda_cfg , msh , var, mat_ns);
    // node-centered 壁 Dirichlet: IC の壁ノード速度を厳密 0 に初期化 (KE を roe から除去)。
    // この後 gasProperties が補正 roe から P/T を再計算する。cell/非 node では no-op。
    enforceWallNoSlip_d_wrapper(cfg , cuda_cfg , msh , var);
    gasProperties_d_wrapper(cfg , cuda_cfg , msh , var);

    fluct.allocVariables();
    fluct.set_fluctVelocity(cfg , cuda_cfg , msh , var);

    applyBconds(cfg , cuda_cfg , msh , var, mat_ns , fluct);
    applyRansScalarBoundaries(cfg , cuda_cfg , msh , var);
    applySpeciesBoundaries(cfg , cuda_cfg , msh , var);
    applyCondensationBoundaries(cfg , cuda_cfg , msh , var);
    calcGradient_d_wrapper(cfg , cuda_cfg , msh , var);
    axisymmetricGeomTerms_d_wrapper(cfg , cuda_cfg , msh , var);
    updateVariablesOuter(cfg , cuda_cfg , msh , var , mat_ns);
    speciesUpdateOuter_d_wrapper(cfg , cuda_cfg , msh , var);  // roY{s}N/M ベースライン
    condensationUpdateOuter_d_wrapper(cfg , cuda_cfg , msh , var);  // 液相モーメント N/M ベースライン
    setDT_d_wrapper(cfg , cuda_cfg , msh , var);

    pprobes.init(cfg , cuda_cfg , msh);

    return cuda_cfg;
}

void writeStepOutputs(
    solverConfig& cfg,
    cudaConfig& cuda_cfg,
    mesh& msh,
    variables& var,
    point_probes& pprobes,
    int iStep)
{
    outputH5_XDMF(cfg , msh, var, iStep);
    outputBconds_H5_XDMF(cfg , msh, var, iStep);
    pprobes.outputProbes(cfg , cuda_cfg , msh , var , iStep);
}

void writeInitialOutputs(
    solverConfig& cfg,
    mesh& msh,
    variables& var)
{
    outputH5_XDMF(cfg , msh, var, 0);
}

// 1 ステップ分の状態を束ねる軽量コンテキスト。旧 advanceOneStep の [&] capture を置換し、
// 自由関数間でデータフローを明示する。
struct StepContext {
    solverConfig&       cfg;
    cudaConfig&         cuda_cfg;
    mesh&               msh;
    matrix&             mat_ns;   // 陰解法では未使用だが非陰解法のシグネチャに必要なため保持
    variables&          var;
    fluct_variables&    fluct;
    point_probes&       pprobes;
    RuntimeProfiler&    profiler;
    ResidualCsvLogger&  residual_logger;
    ImplicitDiagLogger& implicit_diag_logger;
    int                 iStep;
};

// 残差組み立ての単一情報源（旧 assembleCurrentState）。保存量から派生量・境界・勾配・各フラックス・
// ソース項を計算し res_* を確定する。explicit / implicit 双方が呼ぶ。
void assembleResidual(StepContext& s, int stage_index)
{
    (void)stage_index;  // 現状カーネルは stage_index を使わない（dual-time 拡張用に interface 保持）
    s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
        updateVariablesInner(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
    });
    s.profiler.measureCuda(ProfileSection::DependentVariables, [&]() {
        speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // Y_s = ρY_s/ρ (混合則 thermo の前)
        condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // φ = ρφ/ρ (スカラ移流の上流値)
    });
    s.profiler.measureWall(ProfileSection::DependentVariables, [&]() {
        dependentVariables(s.cfg , s.cuda_cfg , s.msh , s.var, s.mat_ns);
    });
    s.profiler.measureCuda(ProfileSection::GasProperties, [&]() {
        gasProperties_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
        applyBconds(s.cfg , s.cuda_cfg , s.msh , s.var, s.mat_ns , s.fluct);
    });
    s.profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
        applyRansScalarBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
        applySpeciesBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
        applyCondensationBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::CalcGradient, [&]() {
        calcGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
        axisymmetricGeomTerms_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::Limiter, [&]() {
        limiter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::DucrosSensor, [&]() {
        ducrosSensor_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        turbulent_viscosity_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::ConvectiveFlux, [&]() {
        convectiveFlux_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, s.mat_ns);
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        ransTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        speciesTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);  // 化学種移流残差
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        condensationTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);  // 液相モーメント移流残差 (Phase 1)
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        ransGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        ransSource_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
        axisymmetricSource_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::ViscousFlux, [&]() {
        viscousFlux_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, s.mat_ns);
    });
    // SU2 流の軸対称対称面 (MARKER_SYM): 軸上 CV の半径方向運動量残差を 0 に射影し roUy=0 を保つ。
    // explicit では ΔroUy=0 になり直接効く (implicit/block-DPLUR では連成 solve が補正を漏らすため Jacobian
    // 整合が別途必要・open issue, docs §7.1)。cell/非軸対称/平面では no-op。
    zeroAxisRadialResidual_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // node-centered 壁 Dirichlet: 壁ノードの運動量残差を 0 に射影し u=0 を保つ (壁ゴースト撤廃の代替)。
    // 軸射影の後に置き、コーナー (壁∩軸はまれだが) でも壁 no-slip を最終確定する。cell/非 node では no-op。
    zeroWallMomentumResidual_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // TODO(dual-time): unsteady のとき addUnsteadyTimeTerm(s) で BDF 物理時間項を res_* と
    // 対角に加える。定常では no-op。本体は次フェーズ。
}

// 残差スナップショットを CSV へ記録（inner_index==0 で outer_begin、それ以外で inner_iter）。
void logResidualSnapshot(StepContext& s, int inner_index)
{
    const ResidualSnapshot residual_snapshot = gatherResidualSnapshot(s.cfg, s.msh, s.var);
    if (inner_index == 0) {
        s.residual_logger.logOuterBegin(s.iStep, inner_index, residual_snapshot);
    } else {
        s.residual_logger.logInnerIter(s.iStep, inner_index, residual_snapshot);
    }
    if (s.implicit_diag_logger.enabled() && s.cfg.timeIntegration == 11) {
        s.implicit_diag_logger.log(s.iStep, inner_index, gatherImplicitDiagSnapshot(s.cfg, s.msh, s.var));
    }
}

// 古典 DPLUR 線形ソルバ。固定残差 res_* に対し Q を更新せず dq_block を nStepInner 回 Jacobi 緩和する。
// 各 sweep 後にバッファを swap し、最終補正は dq_block_old に残る（commit は呼び出し側）。
void blockDPLURSolve(StepContext& s)
{
    // 古典 DPLUR は dq=0 から開始する。前ステップの残留値による近傍参照を避けるため明示ゼロ化。
    // blockDPLUR==1: 5×5 block 版 (dq_block_*)、blockDPLUR==0: scalar 対角版 (dq_ro_* 等)。
    const size_t bytes = static_cast<size_t>(s.msh.nCells_all) * sizeof(flow_float);
    const bool useBlock = (s.cfg.blockDPLUR == 1);
    if (useBlock) {
        cudaMemset(s.var.c_d["dq_block_old_0"], 0, bytes);
        cudaMemset(s.var.c_d["dq_block_old_1"], 0, bytes);
        cudaMemset(s.var.c_d["dq_block_old_2"], 0, bytes);
        cudaMemset(s.var.c_d["dq_block_old_3"], 0, bytes);
        cudaMemset(s.var.c_d["dq_block_old_4"], 0, bytes);
    } else {
        cudaMemset(s.var.c_d["dq_ro_old"],   0, bytes);
        cudaMemset(s.var.c_d["dq_roUx_old"], 0, bytes);
        cudaMemset(s.var.c_d["dq_roUy_old"], 0, bytes);
        cudaMemset(s.var.c_d["dq_roUz_old"], 0, bytes);
        cudaMemset(s.var.c_d["dq_roe_old"],  0, bytes);
    }

    const int nSweep = std::max(1, s.cfg.nStepInner);
    for (int iSweep = 0; iSweep < nSweep; ++iSweep) {
        s.profiler.measureCuda(ProfileSection::TimeIntegration, [&]() {
            timeIntegration_d_wrapper(iSweep, s.cfg , s.cuda_cfg , s.msh , s.var);
        });
        if (useBlock) {
            swapBlockImplicitCorrectionBuffers(s.var);
        } else {
            swapScalarImplicitCorrectionBuffers(s.var);
        }
    }
}

// 1 回の非線形（擬似時間）更新。定常・dual-time 共有の核。
// 残差 1 回構築 → 局所擬似時間 dτ → 古典 DPLUR 線形解 → Q への commit。
void implicitNonlinearUpdate(StepContext& s, int inner_index)
{
    assembleResidual(s, 1);
    logResidualSnapshot(s, inner_index);
    s.profiler.measureCuda(ProfileSection::SetDt, [&]() {
        setDT_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);
    });
    blockDPLURSolve(s);
    s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
        if (s.cfg.blockDPLUR == 1) {
            applyBlockImplicitCorrection(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
        } else {
            applyScalarImplicitCorrection(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
        }
    });
    // 注: 軸上 roUy=0 の commit 後強制は block-DPLUR と非整合で Mach~1000 に発散 (外部状態手術不可)。
    // SU2 流の対称面は Jacobian 内で対称化する必要がある (open issue, docs §7.1)。暫定で無効。
    // enforceAxisSymmetry_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // SST (k-ω) を segregated point-implicit で更新（凍結解除）。残差・消散ヤコビアンは
    // 直前の assembleResidual (ransSource) で確定済み、dt_local は setDT 済み。
    if (scalarResidualEnabled(s.cfg)) {
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            applySSTPointImplicit(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
        });
    }

    // 化学種 (多成分 TP) を segregated point-implicit で更新（陰解法での凍結解除）。
    // res_roY/transport_diag/src_jac は assembleResidual の speciesTransport で確定済み。
    // baseline roY_N=roY を取り、1 回 point-implicit forward-Euler 更新 → 実現可能性再正規化。
    // 各 wrapper は nSpecies<2 で no-op (単成分/CPG は不変)。
    s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
        speciesUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // roY_N = roY_M = roY
        speciesTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
        speciesRenormalize_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);     // Y=roY/ρ (出力/次残差用に同期)
    });

    // 液相モーメント (非平衡凝縮) を segregated point-implicit で更新 (Phase 1 はソース=0 の純移流)。
    // res_/transport_diag は assembleResidual の condensationTransport で確定済み。各 wrapper は
    // condensation==0 で no-op。
    s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
        condensationUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // ro*_N = ro*_M = ro*
        condensationTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
        condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);     // φ=ρφ/ρ (出力/次残差用に同期)
    });
}

// 陽解法 Runge-Kutta（tI 1/3/4）。steady/unsteady 両対応。
void advanceExplicitRK(StepContext& s)
{
    const int iteration_count = s.cfg.perStepIterationCount();
    const char* iteration_label = s.cfg.perStepIterationLabel();

    for (int iloop = 0 ; iloop < iteration_count ; iloop++) {
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            updateVariablesInner(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
            speciesUpdateInner_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // roY{s}M ステージ始点
            condensationUpdateInner_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // 液相モーメント M ステージ始点
        });

        cout << "       " << iteration_label << " : " << iloop+1 << "\n";
        assembleResidual(s, iloop + 1);
        logResidualSnapshot(s, iloop);
        s.profiler.measureCuda(ProfileSection::TimeIntegration, [&]() {
            timeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);
            ransTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);
            speciesTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);
            speciesRenormalize_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // ρY_s>=0, ΣρY_s=ρ
            condensationTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);  // 液相モーメント (Phase 1 ソース=0)
        });
        // 注: explicit 軸対称は軸 CV が step1 で発散するため (recipe 併用でも不変)、enforce は呼んでも
        // 検証できない。explicit の near-axis 安定化は別途要 (open issue)。暫定で無効。
        // enforceAxisSymmetry_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    }

    s.profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
        updateVariablesOuter(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
        speciesUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // roY{s}N/M 次ステップ用
        speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);    // 出力 Y_s を最終 roY_s と同期
        condensationUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // 液相モーメント N/M 次ステップ用
        condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);    // 出力 φ を最終 ρφ と同期
    });
    s.profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
        writeStepOutputs(s.cfg , s.cuda_cfg , s.msh , s.var , s.pprobes , s.iStep+1);
    });
    s.profiler.measureCuda(ProfileSection::SetDt, [&]() {
        setDT_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);
    });
    s.residual_logger.logOuterEnd(s.iStep);
    if (s.cfg.unsteady == 1) {
        s.cfg.totalTime += s.cfg.dt;
    }
}

// 定常 block DPLUR 陰解法。メインループ（nStepOuter）を擬似時間とし、1 ステップ = 1 非線形更新の縮退形。
void advanceImplicitSteady(StepContext& s)
{
    // baseline (roN) は前ステップ末尾 / 初期化の updateVariablesOuter で設定済み（ro == roN）。
    implicitNonlinearUpdate(s, 0);

    s.profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
        updateVariablesOuter(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
    });
    s.profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
        writeStepOutputs(s.cfg , s.cuda_cfg , s.msh , s.var , s.pprobes , s.iStep+1);
    });
    s.profiler.measureCuda(ProfileSection::SetDt, [&]() {
        setDT_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);
    });
    s.residual_logger.logOuterEnd(s.iStep);
}

// 非定常 dual-time 陰解法。1 物理ステップ = 時間レベルシフト → 擬似時間サブ反復（BDF 物理時間項つき
// 非線形更新を nSubIterDualTime 回）→ 物理時間前進。BDF1（初回）/ BDF2（以降, bdfOrder==2 のとき）。
void advanceImplicitDualTime(StepContext& s)
{
    if (s.cfg.dualTime != 1) {
        throw std::runtime_error(
            "Unsteady implicit (timeIntegration==11 && unsteady==1) requires dualTime=1.");
    }
    if (s.cfg.blockDPLUR != 1) {
        throw std::runtime_error(
            "Dual-time implicit currently supports blockDPLUR=1 (5x5 block) only.");
    }
    if (s.cfg.dtControl != 0) {
        // 物理 Δt は固定でなければならない（setDT が cfg.dt を CFL 適応すると BDF 項が壊れる）。
        throw std::runtime_error(
            "Dual-time implicit requires time.deltaT.control=0 (fixed physical dt).");
    }

    // 物理時間レベルシフト: roNN ← roN, roN ← ro（現在の ro = Q^n）。
    s.profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
        shiftDualTimeLevels_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });

    // BDF 係数: 初回ステップ or bdfOrder==1 は BDF1 (1,1,0)、以降 BDF2 (3/2,2,1/2)。
    const bool useBDF2 = (s.cfg.bdfOrder >= 2) && (s.iStep > 0);
    const flow_float a = useBDF2 ? static_cast<flow_float>(1.5) : static_cast<flow_float>(1.0);
    const flow_float b = useBDF2 ? static_cast<flow_float>(2.0) : static_cast<flow_float>(1.0);
    const flow_float c = useBDF2 ? static_cast<flow_float>(0.5) : static_cast<flow_float>(0.0);
    const int include_scalar = scalarResidualEnabled(s.cfg) ? 1 : 0;
    // 対角へ加える物理時間項係数 a/Δt（block/scalar/SST カーネルが cfg 経由で参照）。
    s.cfg.unsteadyDiagCoef = a / std::max(s.cfg.dt, static_cast<flow_float>(1.0e-30));

    const int nSub = std::max(1, s.cfg.nSubIterDualTime);
    for (int m = 0; m < nSub; ++m) {
        assembleResidual(s, 1);
        // 残差に物理時間 BDF 項を加える: res* = res - (V/Δt)(a Q - b Q^n + c Q^{n-1})。
        addUnsteadyTimeTerm_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, a, b, c, include_scalar);
        logResidualSnapshot(s, m);
        s.profiler.measureCuda(ProfileSection::SetDt, [&]() {
            setDT_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);
        });
        blockDPLURSolve(s);
        // dual-time の commit は in-place（roN=Q^n は BDF 基準で固定のため roN+dq は使えない）。
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            applyBlockImplicitCorrectionInPlace_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        });
        if (include_scalar) {
            s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
                applySSTPointImplicit(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
            });
        }
        // 液相モーメント (非平衡凝縮) を segregated point-implicit で更新。condensation==0 で no-op。
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            condensationUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            condensationTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
            condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        });
    }

    s.cfg.unsteadyDiagCoef = 0.0; // 定常側へ影響しないようリセット
    s.cfg.totalTime += s.cfg.dt;

    s.profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
        writeStepOutputs(s.cfg , s.cuda_cfg , s.msh , s.var , s.pprobes , s.iStep+1);
    });
    s.residual_logger.logOuterEnd(s.iStep);
}

// detectNaN==1 のときのみ毎ステップ終端で呼ぶ診断ルーチン。保存量 (ro,roUx,roUy,roUz,roe) と圧力 P
// (RANS 時は roK,roOmega も) を内部セルにわたって検査し、NaN/Inf があれば解を res_nan_<step>.h5 に
// 強制ダンプして即停止する。off のときは一切呼ばれないので通常実行の性能・結果には影響しない。
void checkNonFiniteAndHalt(StepContext& s)
{
    std::vector<std::string> names = {"ro", "roUx", "roUy", "roUz", "roe", "P"};
    if (scalarResidualEnabled(s.cfg)) {
        names.emplace_back("roK");
        names.emplace_back("roOmega");
    }

    std::string offending;
    bool bad = false;
    if (s.cfg.gpu == 1) {
        bad = hasNonFiniteCellValue_d(s.msh, s.var, names, offending);
    } else {
        for (const std::string& name : names) {
            auto it = s.var.c.find(name);
            if (it == s.var.c.end()) continue;
            for (geom_int ic = 0; ic < s.msh.nCells; ic++) {
                if (!std::isfinite(it->second[ic])) { bad = true; offending = name; break; }
            }
            if (bad) break;
        }
    }

    if (!bad) return;

    const int dumpStep = s.iStep + 1;
    std::cerr << "[detectNaN] Non-finite value detected in '" << offending
              << "' at step " << dumpStep
              << ". Dumping solution to res_nan_" << dumpStep
              << ".h5 and halting." << std::endl;
    dumpSolutionH5_force(s.cfg, s.msh, s.var, dumpStep, "res_nan_");
    std::exit(EXIT_FAILURE);
}

void advanceOneStep(
    solverConfig& cfg,
    cudaConfig& cuda_cfg,
    mesh& msh,
    matrix& mat_ns,
    variables& var,
    fluct_variables& fluct,
    point_probes& pprobes,
    RuntimeProfiler& profiler,
    ResidualCsvLogger& residual_logger,
    ImplicitDiagLogger& implicit_diag_logger,
    int iStep)
{
    cout << "----------------------------\n";
    cout << "Step : " << iStep;
    if (cfg.unsteady == 1) {
        cout << "  Time : " << cfg.totalTime;
    }
    cout << "\n";

    StepContext s{cfg, cuda_cfg, msh, mat_ns, var, fluct, pprobes,
                  profiler, residual_logger, implicit_diag_logger, iStep};

    profiler.measureWall(ProfileSection::StepTotal, [&]() {
        if (cfg.isImplicit == 1) {
            if (cfg.unsteady == 1) {
                advanceImplicitDualTime(s);
            } else {
                advanceImplicitSteady(s);
            }
        } else {
            advanceExplicitRK(s);
        }
    });

    // detectNaN 診断モード (既定 off): 保存量+P を検査し NaN/Inf があればダンプして停止する。
    if (cfg.detectNaN == 1) {
        checkNonFiniteAndHalt(s);
    }
}

}

int main(void) {
    clock_t start = clock();
    RuntimeProfiler profiler;

    solverConfig cfg;
    mesh msh;
    matrix mat_ns;
    variables var;
    fluct_variables fluct;
    point_probes pprobes;
    ImplicitDiagLogger implicit_diag_logger;

    cudaConfig cuda_cfg = initializeSimulation(cfg, msh, mat_ns, var, fluct, pprobes);
    ResidualCsvLogger residual_logger("residual_history.csv", cfg);

    writeInitialOutputs(cfg , msh , var);

    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep < cfg.mainLoopCount() ; iStep++) {
        advanceOneStep(cfg , cuda_cfg , msh , mat_ns , var , fluct , pprobes , profiler , residual_logger , implicit_diag_logger , iStep);
    }

    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time = %.3f s\n", time); 
    profiler.printSummary();

	return 0;
}