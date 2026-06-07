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
#include <sstream>
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

void dumpAxisymmetricGradientDebug(
    solverConfig& cfg,
    mesh& msh,
    variables& var,
    int step,
    const char* tag);

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
    axisymmetricDiagnostics_d_wrapper(cfg , cuda_cfg , msh , var);
    dumpAxisymmetricGradientDebug(cfg , msh , var , 0 , "init");
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

void dumpAxisymmetricGradientDebug(
    solverConfig& cfg,
    mesh& msh,
    variables& var,
    int step,
    const char* tag)
{
    if (cfg.isAxisymmetric != 1) {
        return;
    }

    bcond* axis_bc = nullptr;
    for (auto& bc : msh.bconds) {
        if (bc.bcondKind == "axis") {
            axis_bc = &bc;
            break;
        }
    }
    if (axis_bc == nullptr || axis_bc->iPlanes.empty()) {
        return;
    }

    var.copyVariables_cell_D2H({
        "ccx", "ccy", "ccz", "volume", "A_planar",
        "ro", "Ux", "Uy", "Uz", "P", "dt_local", "cfl",
        "dUxdx", "dUxdy", "dUxdz",
        "dUydx", "dUydy", "dUydz",
        "dPdx", "dPdy", "dPdz",
        "axisym_r", "axisym_p_over_r", "axisym_uy_over_r", "axisym_divU", "axisym_roUy_source"
    });
    var.copyVariables_plane_D2H({"sx", "sy", "sz", "ss", "pcx", "pcy", "pcz"});
    axis_bc->copyVariables_bplane_D2H();

    char plane_fname[256];
    char cell_fname[256];
    std::snprintf(plane_fname, sizeof(plane_fname), "axis_debug_planes_%s_step%04d.csv", tag, step);
    std::snprintf(cell_fname, sizeof(cell_fname), "axis_debug_cells_%s_step%04d.csv", tag, step);

    std::ofstream plane_out(plane_fname);
    std::ofstream cell_out(cell_fname);

    plane_out
        << "ib,ip,ic,ig,pcx,pcy,pcz,sx,sy,sz,ss,ccy_owner,ccy_ghost,volume_owner,A_planar_owner,"
        << "axisym_r,p_over_r,uy_over_r,axisym_divU,roUy_source,radial_momentum_balance,"
        << "Ux_owner,Uy_owner,Uz_owner,P_owner,ro_owner,Un_owner,"
        << "Ux_ghost,Uy_ghost,Uz_ghost,"
        << "Ux_bface,Uy_bface,Uz_bface,Ps_bface,"
        << "axis_face_dUydy,axis_face_dPdy,total_dUydy,total_dPdy,dt_local,cfl\n";
    cell_out
        << "ic,ccx,ccy,ccz,volume,A_planar,ro,Ux,Uy,Uz,P,"
        << "axisym_r,axisym_p_over_r,axisym_uy_over_r,axisym_divU,axisym_roUy_source,radial_momentum_balance,"
        << "dUxdx,dUxdy,dUxdz,dUydx,dUydy,dUydz,dPdx,dPdy,dPdz,dt_local,cfl\n";

    std::vector<int> dumped_cell(msh.nCells_all, 0);
    for (std::size_t ib = 0; ib < axis_bc->iPlanes.size(); ++ib) {
        const geom_int ip = axis_bc->iPlanes[ib];
        const geom_int ic = axis_bc->iCells[ib];
        const geom_int ig = axis_bc->iCells_ghst[ib];

        const flow_float sx = var.p["sx"][ip];
        const flow_float sy = var.p["sy"][ip];
        const flow_float sz = var.p["sz"][ip];
        const flow_float ss = var.p["ss"][ip];
        const flow_float volume = var.c["volume"][ic];

        const flow_float Ux = var.c["Ux"][ic];
        const flow_float Uy = var.c["Uy"][ic];
        const flow_float Uz = var.c["Uz"][ic];
        const flow_float Un = (Ux * sx + Uy * sy + Uz * sz) / ss;
        const flow_float p_over_r = var.c["axisym_p_over_r"][ic];
        const flow_float radial_momentum_balance = -var.c["dPdy"][ic] + p_over_r;

        const flow_float axis_face_dUydy = sy * axis_bc->bvar["Uy"][ib] / volume;
        const flow_float axis_face_dPdy = sy * axis_bc->bvar["Ps"][ib] / volume;

        plane_out
            << ib << ',' << ip << ',' << ic << ',' << ig << ','
            << var.p["pcx"][ip] << ',' << var.p["pcy"][ip] << ',' << var.p["pcz"][ip] << ','
            << sx << ',' << sy << ',' << sz << ',' << ss << ','
            << var.c["ccy"][ic] << ',' << var.c["ccy"][ig] << ','
            << volume << ',' << var.c["A_planar"][ic] << ','
            << var.c["axisym_r"][ic] << ',' << p_over_r << ',' << var.c["axisym_uy_over_r"][ic] << ','
            << var.c["axisym_divU"][ic] << ',' << var.c["axisym_roUy_source"][ic] << ',' << radial_momentum_balance << ','
            << Ux << ',' << Uy << ',' << Uz << ','
            << var.c["P"][ic] << ',' << var.c["ro"][ic] << ',' << Un << ','
            << var.c["Ux"][ig] << ',' << var.c["Uy"][ig] << ',' << var.c["Uz"][ig] << ','
            << axis_bc->bvar["Ux"][ib] << ',' << axis_bc->bvar["Uy"][ib] << ',' << axis_bc->bvar["Uz"][ib] << ','
            << axis_bc->bvar["Ps"][ib] << ','
            << axis_face_dUydy << ',' << axis_face_dPdy << ','
            << var.c["dUydy"][ic] << ',' << var.c["dPdy"][ic] << ','
            << var.c["dt_local"][ic] << ',' << var.c["cfl"][ic] << '\n';

        if (dumped_cell[ic] == 0) {
            dumped_cell[ic] = 1;
            cell_out
                << ic << ','
                << var.c["ccx"][ic] << ',' << var.c["ccy"][ic] << ',' << var.c["ccz"][ic] << ','
                << var.c["volume"][ic] << ',' << var.c["A_planar"][ic] << ','
                << var.c["ro"][ic] << ',' << var.c["Ux"][ic] << ',' << var.c["Uy"][ic] << ',' << var.c["Uz"][ic] << ',' << var.c["P"][ic] << ','
                << var.c["axisym_r"][ic] << ',' << var.c["axisym_p_over_r"][ic] << ',' << var.c["axisym_uy_over_r"][ic] << ','
                << var.c["axisym_divU"][ic] << ',' << var.c["axisym_roUy_source"][ic] << ',' << radial_momentum_balance << ','
                << var.c["dUxdx"][ic] << ',' << var.c["dUxdy"][ic] << ',' << var.c["dUxdz"][ic] << ','
                << var.c["dUydx"][ic] << ',' << var.c["dUydy"][ic] << ',' << var.c["dUydz"][ic] << ','
                << var.c["dPdx"][ic] << ',' << var.c["dPdy"][ic] << ',' << var.c["dPdz"][ic] << ','
                << var.c["dt_local"][ic] << ',' << var.c["cfl"][ic] << '\n';
        }
    }
}

void dumpAxisymmetricResidualDebug(
    solverConfig& cfg,
    mesh& msh,
    variables& var,
    int step,
    int stage,
    const char* tag)
{
    if (cfg.isAxisymmetric != 1) {
        return;
    }

    bcond* axis_bc = nullptr;
    for (auto& bc : msh.bconds) {
        if (bc.bcondKind == "axis") {
            axis_bc = &bc;
            break;
        }
    }
    if (axis_bc == nullptr || axis_bc->iPlanes.empty()) {
        return;
    }

    var.copyVariables_cell_D2H({
        "ccx", "ccy", "ccz", "volume", "A_planar",
        "ro", "Ux", "Uy", "Uz", "P",
        "dPdy", "dUydy",
        "axisym_r", "axisym_p_over_r", "axisym_uy_over_r", "axisym_divU", "axisym_roUy_source",
        "res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe"
    });

    char cell_fname[256];
    std::snprintf(cell_fname, sizeof(cell_fname), "axis_debug_residual_%s_stage%02d_step%04d.csv", tag, stage, step);
    std::ofstream cell_out(cell_fname);
    cell_out
        << "ic,ccx,ccy,ccz,volume,A_planar,ro,Ux,Uy,Uz,P,"
        << "axisym_r,axisym_p_over_r,axisym_uy_over_r,axisym_divU,axisym_roUy_source,"
        << "dUydy,dPdy,res_ro,res_roUx,res_roUy,res_roUz,res_roe\n";

    std::vector<int> dumped_cell(msh.nCells_all, 0);
    for (std::size_t ib = 0; ib < axis_bc->iPlanes.size(); ++ib) {
        const geom_int ic = axis_bc->iCells[ib];
        if (dumped_cell[ic] != 0) {
            continue;
        }
        dumped_cell[ic] = 1;

        cell_out
            << ic << ','
            << var.c["ccx"][ic] << ',' << var.c["ccy"][ic] << ',' << var.c["ccz"][ic] << ','
            << var.c["volume"][ic] << ',' << var.c["A_planar"][ic] << ','
            << var.c["ro"][ic] << ',' << var.c["Ux"][ic] << ',' << var.c["Uy"][ic] << ',' << var.c["Uz"][ic] << ',' << var.c["P"][ic] << ','
            << var.c["axisym_r"][ic] << ',' << var.c["axisym_p_over_r"][ic] << ',' << var.c["axisym_uy_over_r"][ic] << ','
            << var.c["axisym_divU"][ic] << ',' << var.c["axisym_roUy_source"][ic] << ','
            << var.c["dUydy"][ic] << ',' << var.c["dPdy"][ic] << ','
            << var.c["res_ro"][ic] << ',' << var.c["res_roUx"][ic] << ',' << var.c["res_roUy"][ic] << ',' << var.c["res_roUz"][ic] << ',' << var.c["res_roe"][ic] << '\n';
    }
}

flow_float interp1stUpHost(
    int,
    int,
    flow_float phiC,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float)
{
    return phiC;
}

flow_float interpMuscl2ndHost(
    int,
    int,
    flow_float phiC,
    flow_float,
    flow_float dphidx,
    flow_float dphidy,
    flow_float dphidz,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float cpdx,
    flow_float cpdy,
    flow_float cpdz,
    flow_float,
    flow_float limiter)
{
    return phiC + limiter * (dphidx * cpdx + dphidy * cpdy + dphidz * cpdz);
}

flow_float interpMuscl3rdHost(
    int,
    int,
    flow_float phiC,
    flow_float phiD,
    flow_float dphidx,
    flow_float dphidy,
    flow_float dphidz,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float,
    flow_float cpdx,
    flow_float cpdy,
    flow_float cpdz,
    flow_float,
    flow_float limiter)
{
    const flow_float k = static_cast<flow_float>(1.0 / 3.0);
    return phiC + limiter * (static_cast<flow_float>(0.5) * k * (phiD - phiC)
        + (static_cast<flow_float>(1.0) - k) * (dphidx * cpdx + dphidy * cpdy + dphidz * cpdz));
}

flow_float interpMinmodHost(
    int,
    int,
    flow_float phiC,
    flow_float phiD,
    flow_float dphidxC,
    flow_float dphidyC,
    flow_float dphidzC,
    flow_float dphidxD,
    flow_float dphidyD,
    flow_float dphidzD,
    flow_float dx,
    flow_float dy,
    flow_float dz,
    flow_float cpdx,
    flow_float cpdy,
    flow_float cpdz,
    flow_float f,
    flow_float)
{
    const flow_float dd2dx = (static_cast<flow_float>(2.0) * f - static_cast<flow_float>(1.0)) * dx;
    const flow_float dd2dy = (static_cast<flow_float>(2.0) * f - static_cast<flow_float>(1.0)) * dy;
    const flow_float dd2dz = (static_cast<flow_float>(2.0) * f - static_cast<flow_float>(1.0)) * dz;

    const flow_float phiDD = phiD - (dd2dx * dphidxD + dd2dy * dphidyD + dd2dz * dphidzD);
    const flow_float phiU = phiDD - static_cast<flow_float>(4.0) * (static_cast<flow_float>(1.0) - f)
        * (dx * dphidxC + dy * dphidyC + dz * dphidzC);
    const flow_float denom = phiD - phiC;
    const flow_float r = (std::fabs(denom) > static_cast<flow_float>(1.0e-30))
        ? (phiC - phiU) / denom
        : static_cast<flow_float>(0.0);
    const flow_float limit = std::max(static_cast<flow_float>(0.0), std::min(static_cast<flow_float>(1.0), r));
    return phiC + limit * (dphidxC * cpdx + dphidyC * cpdy + dphidzC * cpdz);
}

flow_float interpDispatchHost(
    int scheme,
    int limit_scheme,
    flow_float phiC,
    flow_float phiD,
    flow_float dphidxC,
    flow_float dphidyC,
    flow_float dphidzC,
    flow_float dphidxD,
    flow_float dphidyD,
    flow_float dphidzD,
    flow_float dx,
    flow_float dy,
    flow_float dz,
    flow_float cpdx,
    flow_float cpdy,
    flow_float cpdz,
    flow_float f,
    flow_float limiter)
{
    if (scheme == 0 || scheme == -1) {
        return interp1stUpHost(scheme, limit_scheme, phiC, phiD,
            dphidxC, dphidyC, dphidzC,
            dphidxD, dphidyD, dphidzD,
            dx, dy, dz, cpdx, cpdy, cpdz, f, limiter);
    }
    if (scheme == 1 && limit_scheme >= 0) {
        return interpMuscl2ndHost(scheme, limit_scheme, phiC, phiD,
            dphidxC, dphidyC, dphidzC,
            dphidxD, dphidyD, dphidzD,
            dx, dy, dz, cpdx, cpdy, cpdz, f, limiter);
    }
    if (scheme == 2 && limit_scheme >= 0) {
        return interpMuscl3rdHost(scheme, limit_scheme, phiC, phiD,
            dphidxC, dphidyC, dphidzC,
            dphidxD, dphidyD, dphidzD,
            dx, dy, dz, cpdx, cpdy, cpdz, f, limiter);
    }
    return interpMinmodHost(scheme, limit_scheme, phiC, phiD,
        dphidxC, dphidyC, dphidzC,
        dphidxD, dphidyD, dphidzD,
        dx, dy, dz, cpdx, cpdy, cpdz, f, limiter);
}

flow_float signSanoHost(flow_float value)
{
    return value > static_cast<flow_float>(0.0)
        ? static_cast<flow_float>(1.0)
        : (value < static_cast<flow_float>(0.0) ? static_cast<flow_float>(-1.0) : static_cast<flow_float>(0.0));
}

flow_float betaPlsSlauHost(flow_float mach)
{
    if (std::fabs(mach) >= static_cast<flow_float>(1.0)) {
        return static_cast<flow_float>(0.25) * (static_cast<flow_float>(2.0) - mach) * std::pow(mach + static_cast<flow_float>(1.0), static_cast<flow_float>(2.0));
    }
    return static_cast<flow_float>(0.5) * (static_cast<flow_float>(1.0) + signSanoHost(mach));
}

flow_float betaMnsSlauHost(flow_float mach)
{
    if (std::fabs(mach) >= static_cast<flow_float>(1.0)) {
        return static_cast<flow_float>(0.25) * (static_cast<flow_float>(2.0) + mach) * std::pow(mach - static_cast<flow_float>(1.0), static_cast<flow_float>(2.0));
    }
    return static_cast<flow_float>(0.5) * (static_cast<flow_float>(1.0) + signSanoHost(-mach));
}

void dumpAxisymmetricConvectiveFaceDebug(
    solverConfig& cfg,
    mesh& msh,
    variables& var,
    int step,
    int stage)
{
    if (cfg.isAxisymmetric != 1 || cfg.solver != "SLAU") {
        return;
    }

    bcond* axis_bc = nullptr;
    for (auto& bc : msh.bconds) {
        if (bc.bcondKind == "axis") {
            axis_bc = &bc;
            break;
        }
    }
    if (axis_bc == nullptr || axis_bc->iPlanes.empty()) {
        return;
    }

    var.copyVariables_cell_D2H({
        "ccx", "ccy", "ccz", "ro", "Ux", "Uy", "Uz", "P", "Ht", "sonic",
        "res_roUy",
        "limiter_ro", "limiter_Ux", "limiter_Uy", "limiter_Uz", "limiter_P",
        "drodx", "drody", "drodz",
        "dUxdx", "dUxdy", "dUxdz",
        "dUydx", "dUydy", "dUydz",
        "dUzdx", "dUzdy", "dUzdz",
        "dPdx", "dPdy", "dPdz"
    });
    var.copyVariables_plane_D2H({"pcx", "pcy", "pcz", "fx", "sx", "sy", "sz", "ss"});

    char cell_fname[256];
    std::snprintf(cell_fname, sizeof(cell_fname), "axis_debug_convective_faces_stage%02d_step%04d.csv", stage, step);
    std::ofstream out(cell_fname);
    out << "ic,ip,other_ic,is_boundary,cell_side,sx,sy,sz,ss,fx,Uy_L,Uy_R,P_L,P_R,Vn_L,Vn_R,mdot,p_tilde,res_roUy_face,cell_contribution,cell_res_roUy_after_convective\n";

    std::vector<int> dumped_cell(msh.nCells_all, 0);
    for (std::size_t ib = 0; ib < axis_bc->iPlanes.size(); ++ib) {
        const geom_int ic = axis_bc->iCells[ib];
        if (dumped_cell[ic] != 0) {
            continue;
        }
        dumped_cell[ic] = 1;

        for (geom_int ip : msh.cells[ic].iPlanes) {
            const geom_int ic0 = msh.planes[ip].iCells[0];
            const geom_int ic1 = msh.planes[ip].iCells[1];
            const bool boundary_face = (ic1 >= msh.nCells);
            int conv_scheme = cfg.convMethod;
            if (boundary_face) {
                conv_scheme = -1;
            }

            const flow_float f = var.p["fx"][ip];
            const flow_float sxx = var.p["sx"][ip];
            const flow_float syy = var.p["sy"][ip];
            const flow_float szz = var.p["sz"][ip];
            const flow_float sss = var.p["ss"][ip];
            if (std::fabs(sss) <= static_cast<flow_float>(1.0e-30)) {
                continue;
            }

            const flow_float ccx_0 = var.c["ccx"][ic0];
            const flow_float ccy_0 = var.c["ccy"][ic0];
            const flow_float ccz_0 = var.c["ccz"][ic0];
            const flow_float ccx_1 = var.c["ccx"][ic1];
            const flow_float ccy_1 = var.c["ccy"][ic1];
            const flow_float ccz_1 = var.c["ccz"][ic1];

            const flow_float dcc_x = ccx_1 - ccx_0;
            const flow_float dcc_y = ccy_1 - ccy_0;
            const flow_float dcc_z = ccz_1 - ccz_0;
            const flow_float dc0p_x = var.p["pcx"][ip] - ccx_0;
            const flow_float dc0p_y = var.p["pcy"][ip] - ccy_0;
            const flow_float dc0p_z = var.p["pcz"][ip] - ccz_0;
            const flow_float dc1p_x = var.p["pcx"][ip] - ccx_1;
            const flow_float dc1p_y = var.p["pcy"][ip] - ccy_1;
            const flow_float dc1p_z = var.p["pcz"][ip] - ccz_1;

            const flow_float ro_L = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["ro"][ic0], var.c["ro"][ic1],
                var.c["drodx"][ic0], var.c["drody"][ic0], var.c["drodz"][ic0],
                var.c["drodx"][ic1], var.c["drody"][ic1], var.c["drodz"][ic1],
                dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, var.c["limiter_ro"][ic0]);
            const flow_float Ux_L = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["Ux"][ic0], var.c["Ux"][ic1],
                var.c["dUxdx"][ic0], var.c["dUxdy"][ic0], var.c["dUxdz"][ic0],
                var.c["dUxdx"][ic1], var.c["dUxdy"][ic1], var.c["dUxdz"][ic1],
                dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, var.c["limiter_Ux"][ic0]);
            const flow_float Uy_L = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["Uy"][ic0], var.c["Uy"][ic1],
                var.c["dUydx"][ic0], var.c["dUydy"][ic0], var.c["dUydz"][ic0],
                var.c["dUydx"][ic1], var.c["dUydy"][ic1], var.c["dUydz"][ic1],
                dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, var.c["limiter_Uy"][ic0]);
            const flow_float Uz_L = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["Uz"][ic0], var.c["Uz"][ic1],
                var.c["dUzdx"][ic0], var.c["dUzdy"][ic0], var.c["dUzdz"][ic0],
                var.c["dUzdx"][ic1], var.c["dUzdy"][ic1], var.c["dUzdz"][ic1],
                dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, var.c["limiter_Uz"][ic0]);
            const flow_float P_L = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["P"][ic0], var.c["P"][ic1],
                var.c["dPdx"][ic0], var.c["dPdy"][ic0], var.c["dPdz"][ic0],
                var.c["dPdx"][ic1], var.c["dPdy"][ic1], var.c["dPdz"][ic1],
                dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, var.c["limiter_P"][ic0]);

            const flow_float ro_R = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["ro"][ic1], var.c["ro"][ic0],
                var.c["drodx"][ic1], var.c["drody"][ic1], var.c["drodz"][ic1],
                var.c["drodx"][ic0], var.c["drody"][ic0], var.c["drodz"][ic0],
                -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, static_cast<flow_float>(1.0) - f, var.c["limiter_ro"][ic1]);
            const flow_float Ux_R = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["Ux"][ic1], var.c["Ux"][ic0],
                var.c["dUxdx"][ic1], var.c["dUxdy"][ic1], var.c["dUxdz"][ic1],
                var.c["dUxdx"][ic0], var.c["dUxdy"][ic0], var.c["dUxdz"][ic0],
                -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, static_cast<flow_float>(1.0) - f, var.c["limiter_Ux"][ic1]);
            const flow_float Uy_R = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["Uy"][ic1], var.c["Uy"][ic0],
                var.c["dUydx"][ic1], var.c["dUydy"][ic1], var.c["dUydz"][ic1],
                var.c["dUydx"][ic0], var.c["dUydy"][ic0], var.c["dUydz"][ic0],
                -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, static_cast<flow_float>(1.0) - f, var.c["limiter_Uy"][ic1]);
            const flow_float Uz_R = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["Uz"][ic1], var.c["Uz"][ic0],
                var.c["dUzdx"][ic1], var.c["dUzdy"][ic1], var.c["dUzdz"][ic1],
                var.c["dUzdx"][ic0], var.c["dUzdy"][ic0], var.c["dUzdz"][ic0],
                -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, static_cast<flow_float>(1.0) - f, var.c["limiter_Uz"][ic1]);
            const flow_float P_R = interpDispatchHost(conv_scheme, cfg.limiter,
                var.c["P"][ic1], var.c["P"][ic0],
                var.c["dPdx"][ic1], var.c["dPdy"][ic1], var.c["dPdz"][ic1],
                var.c["dPdx"][ic0], var.c["dPdy"][ic0], var.c["dPdz"][ic0],
                -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, static_cast<flow_float>(1.0) - f, var.c["limiter_P"][ic1]);

            const flow_float Vn_L = (Ux_L * sxx + Uy_L * syy + Uz_L * szz) / sss;
            const flow_float Vn_R = (Ux_R * sxx + Uy_R * syy + Uz_R * szz) / sss;
            const flow_float c_hat = static_cast<flow_float>(0.5) * (var.c["sonic"][ic0] + var.c["sonic"][ic1]);
            const flow_float M_L = Vn_L / c_hat;
            const flow_float M_R = Vn_R / c_hat;

            flow_float beta_p;
            if (std::fabs(M_L) >= static_cast<flow_float>(1.0)) {
                beta_p = static_cast<flow_float>(0.5) * (M_L + std::fabs(M_L)) / M_L;
            } else {
                beta_p = static_cast<flow_float>(0.25) * std::pow(M_L + static_cast<flow_float>(1.0), static_cast<flow_float>(2.0))
                    * (static_cast<flow_float>(2.0) - M_L);
            }

            flow_float beta_m;
            if (std::fabs(M_R) >= static_cast<flow_float>(1.0)) {
                beta_m = static_cast<flow_float>(0.5) * (M_R - std::fabs(M_R)) / M_R;
            } else {
                beta_m = static_cast<flow_float>(0.25) * std::pow(M_R - static_cast<flow_float>(1.0), static_cast<flow_float>(2.0))
                    * (static_cast<flow_float>(2.0) + M_R);
            }

            const flow_float g = -std::max(std::min(M_L, static_cast<flow_float>(0.0)), static_cast<flow_float>(-1.0))
                * std::min(std::max(M_R, static_cast<flow_float>(0.0)), static_cast<flow_float>(1.0));
            const flow_float Vn_hat_abs = (ro_L * std::fabs(Vn_L) + ro_R * std::fabs(Vn_R)) / (ro_L + ro_R);
            const flow_float Vn_hat_p_abs = (static_cast<flow_float>(1.0) - g) * Vn_hat_abs + g * std::fabs(Vn_L);
            const flow_float Vn_hat_m_abs = (static_cast<flow_float>(1.0) - g) * Vn_hat_abs + g * std::fabs(Vn_R);

            const flow_float velocity2_avg = static_cast<flow_float>(0.5)
                * ((Ux_R * Ux_R + Uy_R * Uy_R + Uz_R * Uz_R) + (Ux_L * Ux_L + Uy_L * Uy_L + Uz_L * Uz_L));
            const flow_float M_hat = std::min(static_cast<flow_float>(1.0), std::sqrt(velocity2_avg) / c_hat);
            const flow_float chi = (static_cast<flow_float>(1.0) - M_hat) * (static_cast<flow_float>(1.0) - M_hat);
            const flow_float pressure_sum = P_L + P_R;
            const flow_float p_tilde = static_cast<flow_float>(0.5) * pressure_sum
                + static_cast<flow_float>(0.5) * (beta_p - beta_m) * (P_L - P_R)
                + (static_cast<flow_float>(1.0) - chi) * (beta_p + beta_m - static_cast<flow_float>(1.0))
                    * static_cast<flow_float>(0.5) * pressure_sum;

            const flow_float mdot = sss * static_cast<flow_float>(0.5)
                * ((ro_L * (Vn_L + Vn_hat_p_abs) + ro_R * (Vn_R - Vn_hat_m_abs))
                    - chi / c_hat * (P_R - P_L));
            const flow_float res_roUy_face = static_cast<flow_float>(0.5) * (mdot + std::fabs(mdot)) * Uy_L
                + static_cast<flow_float>(0.5) * (mdot - std::fabs(mdot)) * Uy_R
                + p_tilde * syy;
            const flow_float cell_contribution = (ic == ic0)
                ? -res_roUy_face
                : (ic == ic1 ? res_roUy_face : static_cast<flow_float>(0.0));

            out << ic << ',' << ip << ',' << (ic == ic0 ? ic1 : ic0) << ','
                << (boundary_face ? 1 : 0) << ',' << (ic == ic0 ? "left" : "right") << ','
                << sxx << ',' << syy << ',' << szz << ',' << sss << ',' << f << ','
                << Uy_L << ',' << Uy_R << ',' << P_L << ',' << P_R << ','
                << Vn_L << ',' << Vn_R << ',' << mdot << ',' << p_tilde << ','
                << res_roUy_face << ',' << cell_contribution << ',' << var.c["res_roUy"][ic] << '\n';
        }
    }
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
                axisymmetricDiagnostics_d_wrapper(cfg , cuda_cfg , msh , var);
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
            if (cfg.isAxisymmetric == 1 && iStep == 0) {
                dumpAxisymmetricResidualDebug(cfg, msh, var, iStep + 1, stage_index, "after_convective");
                if (stage_index == 1) {
                    dumpAxisymmetricConvectiveFaceDebug(cfg, msh, var, iStep + 1, stage_index);
                }
            }
            profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
                axisymmetricSource_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            if (cfg.isAxisymmetric == 1 && iStep == 0) {
                dumpAxisymmetricResidualDebug(cfg, msh, var, iStep + 1, stage_index, "after_axis_source");
            }
            profiler.measureCuda(ProfileSection::ViscousFlux, [&]() {
                viscousFlux_d_wrapper(cfg , cuda_cfg, msh , var, mat_ns);
            });
            if (cfg.isAxisymmetric == 1 && iStep == 0) {
                dumpAxisymmetricResidualDebug(cfg, msh, var, iStep + 1, stage_index, "after_viscous");
            }
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
                axisymmetricDiagnostics_d_wrapper(cfg , cuda_cfg , msh , var);
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
        if (cfg.isAxisymmetric == 1 && iStep == 0) {
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
                axisymmetricDiagnostics_d_wrapper(cfg , cuda_cfg , msh , var);
            });
            dumpAxisymmetricGradientDebug(cfg , msh , var , 1 , "poststep");
        }
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