#include "boundaryCond.hpp"

#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/boundaryCond_d.cuh"
#include "cuda_forge/ransScalarBoundary_d.cuh"
#include "cuda_forge/fluct_variables_d.cuh"

using namespace std;

bcondConfFormat::bcondConfFormat(){};

void readBcondConfig(solverConfig& cfg , vector<bcond>& bconds)
{
    map<int,bcondConfFormat> bcondConfMap;
    const auto isInletKind = [](const std::string& kind) {
        return kind.rfind("inlet_", 0) == 0;
    };

    std::string bcondConfigFileName  = "bcondConfig.yaml";

    try {
        YAML::Node config = YAML::LoadFile(bcondConfigFileName);

        for(auto bcInYaml : config) {
            bcondConfFormat bcf;

            std::string bname = bcInYaml.first.as<std::string>();
            cout << "bname = " << bname << endl;

            // akirameta
            //if (bcInYaml.find("physID") == bcInYaml.end()) {
            //    // not found
            //    cerr << "Error: couldn't find physID of " << bname << " in bcondConfig.yaml\n";
            //    exit();
            //}

            int physID = bcInYaml.second["physID"].as<int>();

            std::map<std::string, int> inputInts_temp;
            for (auto bcInYaml_int : bcInYaml.second["ints"])
            {
                std::string key = bcInYaml_int.first.as<std::string>();
                int val = bcInYaml_int.second.as<int>();
                inputInts_temp[key] = val;
            }

            std::map<std::string, flow_float> inputFloats_temp;
            for (auto bcInYaml_float : bcInYaml.second["floats"])
            {
                std::string key = bcInYaml_float.first.as<std::string>();
                flow_float val = bcInYaml_float.second.as<flow_float>();
                inputFloats_temp[key] = val;
            }

            std::string kind = bcInYaml.second["kind"].as<std::string>();

            int outputHDFflg_temp = bcInYaml.second["outputHDFflg"].as<int>();

            bcf.physID = physID;
            bcf.physName = bname;
            bcf.kind = kind;
            bcf.inputInts = inputInts_temp;
            bcf.inputFloats = inputFloats_temp;
            bcf.outputHDFflg = outputHDFflg_temp;

            if (cfg.LESorRANS == 2 && isInletKind(kind)) {
                const bool hasK = inputFloats_temp.find("k") != inputFloats_temp.end();
                const bool hasOmega = inputFloats_temp.find("omega") != inputFloats_temp.end();
                if (!hasK || !hasOmega) {
                    cerr << "Error: inlet boundary '" << bname << "' of kind '" << kind
                         << "' must define floats 'k' and 'omega' when LESorRANS=2." << endl;
                    exit(EXIT_FAILURE);
                }
            }

            bcondConfMap[bcf.physID] = bcf;

            cout << "physID=" <<  bcf.physID << endl;
            cout << "kind=" << bcf.kind<< endl;
            cout << "outputHDF=" << bcf.outputHDFflg<< endl;
        }

    } catch(const YAML::BadFile& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);

    } catch(const YAML::ParserException& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);
    }

    // set bconds from cofig file (YAML file)
    for (bcond& bc : bconds)
    {
        cout << "bc.physID=" << bc.physID << endl;

        const auto confIt = bcondConfMap.find(bc.physID);
        if (confIt == bcondConfMap.end()) {
            // not found
            cerr << "Error: couldn't find physID " << bc.physID << " in bcondConfig.yaml\n";
            exit(EXIT_FAILURE);
        }
        const bcondConfFormat& bcf = confIt->second;

        cout << "bcf.kind=" << bcf.kind << endl;

        bc.bcondKind    = bcf.kind;
        bc.physName     = bcf.physName;
        bc.inputInts    = bcf.inputInts;
        bc.inputFloats  = bcf.inputFloats;
        bc.outputHDFflg = bcf.outputHDFflg;

        // copy value types to bcond
        const auto valueTypesIt = bcf.valueTypesOfBC.find(bcf.kind);
        if (valueTypesIt == bcf.valueTypesOfBC.end()) {
            cerr << "Error: unsupported boundary condition kind " << bcf.kind << " for physID " << bc.physID << "\n";
            exit(EXIT_FAILURE);
        }

        const map<string,int>& valueTypes = valueTypesIt->second;
        for (auto& vt : valueTypes)
        {
            bc.valueTypes[vt.first] = vt.second;
        }

        bc.bcondInitVariables(cfg.gpu); // allocate and set boundary variables
    }
};

void applyBconds(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p , fluct_variables& fluct)
{
    if (cfg.gpu == 0) {
        cerr << "Error: gpu=0 is not supported in applyBconds" << endl;
        exit(EXIT_FAILURE);
    } else if (cfg.gpu ==1) { // gpu
        for (auto& bc : msh.bconds)
        {
            //if      (bc.bcondKind == "wall_isothermal") { wall_isothermal_d(cfg , bc , msh , var , mat_p); } 
            if      (bc.bcondKind == "slip") slip_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
            else if (bc.bcondKind == "axis") slip_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
            else if (bc.bcondKind == "wall") wall_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
            else if (bc.bcondKind == "wall_isothermal") wall_isothermal_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
            else if (bc.bcondKind == "inlet_uniformVelocity") { inlet_uniformVelocity_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "inlet_fluctVelocity") { inlet_fluctVelocity_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p , fluct); }
            else if (bc.bcondKind == "outlet_statPress") { outlet_statPress_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "inlet_Pressure") { inlet_Pressure_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "inlet_Pressure_dir") { inlet_Pressure_dir_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "outflow") { outflow_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
            else if (bc.bcondKind == "periodic") { periodic_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); }
        }
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
    } else {
        cerr << "Error: unsupported gpu flag " << cfg.gpu << endl;
        exit(EXIT_FAILURE);
    }
}

void applyRansScalarBoundaries(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!(cfg.gpu == 1 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1)) {
        return;
    }

    for (auto& bc : msh.bconds)
    {
        ransScalarBoundary_d_wrapper(cfg , cuda_cfg , bc , msh , var);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

//void copyBcondsGradient(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p)
//{
//    if (cfg.gpu == 0) { // cpu
//        cerr << "Error: can't use with cpu" << endl;
//        exit(EXIT_FAILURE);
//
//    } else if (cfg.gpu ==1) { // gpu
//        for (auto& bc : msh.bconds)
//        {
//            copyBcondsGradient_d_wrapper(cfg , cuda_cfg , bc , msh , var , mat_p); 
//        }
//        gpuErrchk( cudaPeekAtLastError() );
//        gpuErrchk( cudaDeviceSynchronize() );
//
//    }
//}
//