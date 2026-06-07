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
#include "cuda_forge/scalarTransport_d.cuh"
#include "cuda_forge/ransSource_d.cuh"
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
    TurbulentViscosity,
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

CorrectionSnapshot gatherCorrectionSnapshot(solverConfig& cfg, mesh& msh, variables& var)
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
            case ProfileSection::TurbulentViscosity: return "turb_viscosity";
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

    var.allocVariables(cfg.gpu , msh);

    cout << "Read Initial Values \n";
    var.readValueHDF5(cfg.valueFileName , msh);

    cout << "Set mesh connection map for cuda \n";
    cudaConfig cuda_cfg(msh);
    msh.setMeshMap_d();
    msh.setPeriodicPartner();

    var.setStructuralVariables(cfg , cuda_cfg , msh);
    dependentVariables(cfg , cuda_cfg , msh , var, mat_ns);
    gasProperties_d_wrapper(cfg , cuda_cfg , msh , var);

    fluct.allocVariables();
    fluct.set_fluctVelocity(cfg , cuda_cfg , msh , var);

    applyBconds(cfg , cuda_cfg , msh , var, mat_ns , fluct);
    applyRansScalarBoundaries(cfg , cuda_cfg , msh , var);
    calcGradient_d_wrapper(cfg , cuda_cfg , msh , var);
    axisymmetricGeomTerms_d_wrapper(cfg , cuda_cfg , msh , var);
    updateVariablesOuter(cfg , cuda_cfg , msh , var , mat_ns);
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

    profiler.measureWall(ProfileSection::StepTotal, [&]() {
        const bool steady_implicit = (cfg.isImplicit == 1 && cfg.unsteady == 0);
        const bool frozen_scalar_implicit = (
            cfg.timeIntegration == 11 &&
            cfg.isImplicit == 1 &&
            cfg.unsteady == 0 &&
            cfg.blockDPLUR == 0
        );

        const auto assembleCurrentState = [&](int stage_index) {
            profiler.measureWall(ProfileSection::UpdateInner, [&]() {
                updateVariablesInner(cfg , cuda_cfg , msh , var , mat_ns);
            });
            profiler.measureWall(ProfileSection::DependentVariables, [&]() {
                dependentVariables(cfg , cuda_cfg , msh , var, mat_ns);
            });
            profiler.measureCuda(ProfileSection::GasProperties, [&]() {
                gasProperties_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
                applyBconds(cfg , cuda_cfg , msh , var, mat_ns , fluct);
            });
            profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
                applyRansScalarBoundaries(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::CalcGradient, [&]() {
                calcGradient_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
                axisymmetricGeomTerms_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::Limiter, [&]() {
                limiter_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::DucrosSensor, [&]() {
                ducrosSensor_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::TurbulentViscosity, [&]() {
                turbulent_viscosity_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::ConvectiveFlux, [&]() {
                convectiveFlux_d_wrapper(cfg , cuda_cfg, msh , var, mat_ns);
            });
            profiler.measureCuda(ProfileSection::ConvectiveFlux, [&]() {
                scalarTransport_d_wrapper(cfg , cuda_cfg, msh , var);
            });
            profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
                calcScalarGradient_d_wrapper(cfg , cuda_cfg , msh , var);
                ransSource_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
                axisymmetricSource_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::ViscousFlux, [&]() {
                viscousFlux_d_wrapper(cfg , cuda_cfg, msh , var, mat_ns);
            });
        };

        const auto refreshCommittedState = [&]() {
            profiler.measureWall(ProfileSection::DependentVariables, [&]() {
                dependentVariables(cfg , cuda_cfg , msh , var, mat_ns);
            });
            profiler.measureCuda(ProfileSection::GasProperties, [&]() {
                gasProperties_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
                applyBconds(cfg , cuda_cfg , msh , var, mat_ns , fluct);
            });
            profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
                applyRansScalarBoundaries(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::CalcGradient, [&]() {
                calcGradient_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
                axisymmetricGeomTerms_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::Limiter, [&]() {
                limiter_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::DucrosSensor, [&]() {
                ducrosSensor_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            profiler.measureCuda(ProfileSection::TurbulentViscosity, [&]() {
                turbulent_viscosity_d_wrapper(cfg , cuda_cfg , msh , var);
            });
        };

        const auto logResidualSnapshot = [&](int inner_index) {
            const ResidualSnapshot residual_snapshot = gatherResidualSnapshot(cfg, msh, var);
            if (inner_index == 0) {
                residual_logger.logOuterBegin(iStep, inner_index, residual_snapshot);
            } else {
                residual_logger.logInnerIter(iStep, inner_index, residual_snapshot);
            }
            if (implicit_diag_logger.enabled() && cfg.timeIntegration == 11) {
                implicit_diag_logger.log(iStep, inner_index, gatherImplicitDiagSnapshot(cfg, msh, var));
            }
        };

        if (frozen_scalar_implicit) {
            assembleCurrentState(1);

            const ResidualSnapshot outer_residual_snapshot = gatherResidualSnapshot(cfg, msh, var);
            residual_logger.logOuterBegin(iStep, 0, outer_residual_snapshot);
            if (implicit_diag_logger.enabled() && cfg.timeIntegration == 11) {
                implicit_diag_logger.log(iStep, 0, gatherImplicitDiagSnapshot(cfg, msh, var));
            }

            for (int iInner = 0; iInner < std::max(1, cfg.nStepInner); ++iInner) {
                cout << "       Inner : " << iInner + 1 << "\n";
                profiler.measureCuda(ProfileSection::TimeIntegration, [&]() {
                    timeIntegration_d_wrapper(iInner, cfg , cuda_cfg , msh , var);
                });
                residual_logger.logInnerIter(
                    iStep,
                    iInner + 1,
                    outer_residual_snapshot,
                    gatherCorrectionSnapshot(cfg, msh, var));
                if (implicit_diag_logger.enabled() && cfg.timeIntegration == 11) {
                    implicit_diag_logger.log(iStep, iInner + 1, gatherImplicitDiagSnapshot(cfg, msh, var));
                }
            }

            profiler.measureWall(ProfileSection::UpdateInner, [&]() {
                applyScalarImplicitCorrection(cfg , cuda_cfg , msh , var , mat_ns);
            });
            refreshCommittedState();

            profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
                updateVariablesOuter(cfg , cuda_cfg , msh , var , mat_ns);
            });
            profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
                writeStepOutputs(cfg , cuda_cfg , msh , var , pprobes , iStep+1);
            });
            profiler.measureCuda(ProfileSection::SetDt, [&]() {
                setDT_d_wrapper(cfg , cuda_cfg, msh , var);
            });
            const ResidualSnapshot committed_residual_snapshot = gatherResidualSnapshot(cfg, msh, var);
            const CorrectionSnapshot committed_correction_snapshot = gatherCorrectionSnapshot(cfg, msh, var);
            residual_logger.logOuterEnd(iStep, &committed_residual_snapshot, &committed_correction_snapshot);
            if (cfg.unsteady == 1) {
                cfg.totalTime += cfg.dt;
            }
            return;
        }

        const int iteration_count = cfg.perStepIterationCount();
        const char* iteration_label = cfg.perStepIterationLabel();

        for (int iloop = 0 ; iloop < iteration_count ; iloop++) {
            profiler.measureWall(ProfileSection::UpdateInner, [&]() {
                updateVariablesInner(cfg , cuda_cfg , msh , var , mat_ns);
            });

            cout << "       " << iteration_label << " : " << iloop+1 << "\n";
            assembleCurrentState(iloop + 1);
            logResidualSnapshot(iloop);
            profiler.measureCuda(ProfileSection::TimeIntegration, [&]() {
                timeIntegration_d_wrapper(iloop, cfg , cuda_cfg , msh , var);
                scalarTimeIntegration_d_wrapper(iloop, cfg , cuda_cfg , msh , var);
            });

            if (steady_implicit) {
                const int pseudo_step = iStep * iteration_count + iloop + 1;
                profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
                    writeStepOutputs(cfg , cuda_cfg , msh , var , pprobes , pseudo_step);
                });
                profiler.measureCuda(ProfileSection::SetDt, [&]() {
                    setDT_d_wrapper(cfg , cuda_cfg, msh , var);
                });
            }
        }

        profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
            updateVariablesOuter(cfg , cuda_cfg , msh , var , mat_ns);
        });
        if (!steady_implicit) {
            profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
                writeStepOutputs(cfg , cuda_cfg , msh , var , pprobes , iStep+1);
            });
            profiler.measureCuda(ProfileSection::SetDt, [&]() {
                setDT_d_wrapper(cfg , cuda_cfg, msh , var);
            });
        }
        residual_logger.logOuterEnd(iStep);
        if (cfg.unsteady == 1) {
            cfg.totalTime += cfg.dt;
        }
    });
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