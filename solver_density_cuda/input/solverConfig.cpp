#include "input/solverConfig.hpp"


solverConfig::solverConfig(){};

// getValidatedValue 関数の定義
template <typename T>
T getValidatedValue(const YAML::Node& node, const std::string& key, const std::string& parent = "") {
    try {
        if (!node[key].IsDefined()) {
            throw std::runtime_error("Key '" + key + "' is missing" + (parent.empty() ? "" : " in '" + parent + "'."));
        }
        T value = node[key].as<T>();
        std::cout << "Loaded value for '" << key << "'" << (parent.empty() ? "" : " in '" + parent + "'") << ": " << value << std::endl;
        return value;
    } catch (const YAML::BadConversion& e) {
        throw std::runtime_error("Key '" + key + "' in '" + parent + "' has an incorrect type.");
    }
}


void solverConfig::read(std::string fname)
{
//    try {
//
//
//        this->solConfigFileName = fname;
//        YAML::Node config = YAML::LoadFile(this->solConfigFileName);
//
//        
//        std::string meshFormat    = config["mesh"]["meshFormat"].as<std::string>();
//        std::string meshFileName  = config["mesh"]["meshFileName"].as<std::string>();
//        std::string valueFileName = config["mesh"]["valueFileName"].as<std::string>();
//        std::string solver = config["solver"].as<std::string>();
//
//        this->gpu = config["gpu"].as<int>();
//
//        int unsteady = config["time"]["unsteady"].as<int>();
//        int dualTime = config["time"]["dualTime"].as<int>();
//
//        auto last = config["time"]["last"];
//            int endTimeControl = last["control"].as<int>();
//            int nStep = last["nStep"].as<int>();
//            int endTime = last["time"].as<flow_float>();
//
//        auto deltaT = config["time"]["deltaT"];
//            int dtControl = deltaT["control"].as<int>();
//            flow_float dt = deltaT["dt"].as<flow_float>();
//            flow_float cfl = deltaT["cfl"].as<flow_float>();
//            flow_float cfl_pseudo = deltaT["cfl_pseudo"].as<flow_float>();
//            flow_float dt_max= deltaT["dt_max"].as<flow_float>();
//            flow_float dt_min= deltaT["dt_min"].as<flow_float>();
//
//        int outStepInterval = config["time"]["outStepInterval"].as<int>();
//        int outStepStart    = config["time"]["outStepStart"].as<int>();
//        int timeIntegration = config["time"]["timeIntegration"].as<int>();
//        int nInnerLoop      = config["time"]["nInnerLoop"].as<int>();
//
//        auto space = config["space"];
//            int convMethod = space["convMethod"].as<flow_float>();
//            int limiter    = space["limiter"].as<int>();
//
//        auto turb = config["turbulence"];
//            int LESorRANS = turb["LESorRANS"].as<int>();
//            int LESmodel  = turb["LESmodel"].as<int>();
//
//        auto physProp = config["physProp"];
//            int isCompressible = physProp["isCompressible"].as<int>();
//            int thermalMethod = physProp["thermalMethod"].as<int>();
//            int viscMethod = physProp["viscMethod"].as<int>();
//            printf("viscMethod=%d\n", viscMethod);
//            flow_float ro = physProp["ro"].as<flow_float>();
//            flow_float visc = physProp["visc"].as<flow_float>();
//            flow_float thermCond = physProp["thermCond"].as<flow_float>();
//            flow_float cp = physProp["cp"].as<flow_float>();
//            flow_float gamma = physProp["gamma"].as<flow_float>();
//
//
//        this->initial = config["initial"].as<std::string>();
//
//        std::cout << "Mesh Name : " << meshFileName << '\n';
//        std::cout << "Mesh Name size : " << meshFileName.size() << '\n';
//        std::cout << "Value File Name : " << valueFileName << '\n';
//
//        std::cout << "Solver Name : " << solver << '\n';
//        if        (endTimeControl == 0) { std::cout << "End Step : " << nStep << '\n';
//        } else if (endTimeControl == 1) { std::cout << "End Time : " << endTime << '\n';
//        } else {
//            std::cerr << "something wrong in solver yaml" << std::endl;
//            exit(EXIT_FAILURE);
//        }
//
//        if        (dtControl == 0) { std::cout << "delta Time : " << dt << '\n'; 
//        } else if (dtControl == 1) { std::cout << "CFL : " << cfl << '\n';
//        } else {
//            std::cerr << "something wrong in solver yaml" << std::endl;
//            exit(EXIT_FAILURE);
//        }
//
//        std::cout << "convection method : " << convMethod<< '\n';
//        std::cout << "limiter : " << limiter<< '\n';
//
//        this->meshFormat   = meshFormat;
//        this->meshFileName = meshFileName;
//        this->valueFileName = valueFileName;
//        this->solver = solver;
//
//        this->endTimeControl = endTimeControl;
//        this->nStep = nStep;
//        this->outStepInterval = outStepInterval;
//        this->outStepStart    = outStepStart;
//
//        this->dtControl = dtControl;
//        this->dt = dt;
//        this->cfl = cfl;
//        this->cfl_pseudo = cfl_pseudo;
//        this->dt_max = dt_max;
//        this->dt_min = dt_min;
//        this->unsteady = unsteady;
//        this->dualTime = dualTime;
//        this->timeIntegration = timeIntegration;
//        this->nInnerLoop = nInnerLoop;
//        this->cfl_pseudo = cfl_pseudo;
//
//        initTimeIntegrationScheme(timeIntegration);
//
//        this->convMethod = convMethod;
//        this->limiter = limiter;
//
//
//        this->LESorRANS = LESorRANS;
//        this->LESmodel  = LESmodel;
//
//
//        this->isCompressible = isCompressible ;
//
//        this->thermalMethod = thermalMethod;
//        this->viscMethod = viscMethod;
//        this->ro = ro;
//        this->visc = visc;
//        this->thermCond = thermCond;
//        this->cp = cp;
//        this->gamma = gamma;
//
//        this->initial = initial;
//
//    } catch (const YAML::BadFile& e) {
//        std::cerr << e.msg << std::endl;
//        exit(EXIT_FAILURE);
//
//    } catch (const YAML::ParserException& e) {
//        std::cerr << e.msg << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//
    try {
        //YAML::Node config = YAML::LoadFile("config.yaml");
        this->solConfigFileName = fname;
        YAML::Node config = YAML::LoadFile(this->solConfigFileName);

        // mesh関連
        this->meshFormat = getValidatedValue<std::string>(config["mesh"], "meshFormat", "mesh");
        this->meshFileName = getValidatedValue<std::string>(config["mesh"], "meshFileName", "mesh");
        this->valueFileName = getValidatedValue<std::string>(config["mesh"], "valueFileName", "mesh");

        // solver関連
        this->solver = getValidatedValue<std::string>(config, "solver");

        // GPU設定
        this->gpu = getValidatedValue<int>(config, "gpu");

        // time関連
        this->unsteady = getValidatedValue<int>(config["time"], "unsteady", "time");
        this->dualTime = getValidatedValue<int>(config["time"], "dualTime", "time");

        auto last = config["time"]["last"];
        this->endTimeControl = getValidatedValue<int>(last, "control", "time.last");
        this->nStep = getValidatedValue<int>(last, "nStep", "time.last");

        auto deltaT = config["time"]["deltaT"];
        this->dtControl = getValidatedValue<int>(deltaT, "control", "time.deltaT");
        this->dt = getValidatedValue<double>(deltaT, "dt", "time.deltaT");
        this->cfl = getValidatedValue<double>(deltaT, "cfl", "time.deltaT");
        this->cfl_pseudo = getValidatedValue<double>(deltaT, "cfl_pseudo", "time.deltaT");
        this->dt_max = getValidatedValue<double>(deltaT, "dt_max", "time.deltaT");
        this->dt_min = getValidatedValue<double>(deltaT, "dt_min", "time.deltaT");

        this->outStepInterval = getValidatedValue<int>(config["time"], "outStepInterval", "time");
        this->outStepStart = getValidatedValue<int>(config["time"], "outStepStart", "time");
        this->timeIntegration = getValidatedValue<int>(config["time"], "timeIntegration", "time");
        this->nInnerLoop = getValidatedValue<int>(config["time"], "nInnerLoop", "time");

        // 空間設定
        auto space = config["space"];
        this->convMethod = getValidatedValue<int>(space, "convMethod", "space");
        this->limiter = getValidatedValue<int>(space, "limiter", "space");

        // turbulence model
        auto turb = config["turbulence"];
        this->LESorRANS = getValidatedValue<int>(turb, "LESorRANS", "turbulence");
        this->LESmodel = getValidatedValue<int>(turb, "LESmodel", "turbulence");

        // 物理プロパティ
        auto physProp = config["physProp"];
        this->isCompressible = getValidatedValue<int>(physProp, "isCompressible", "physProp");
        this->thermalMethod = getValidatedValue<int>(physProp, "thermalMethod", "physProp");
        this->viscMethod = getValidatedValue<int>(physProp, "viscMethod", "physProp");
        this->ro = getValidatedValue<double>(physProp, "ro", "physProp");
        this->visc = getValidatedValue<double>(physProp, "visc", "physProp");
        this->thermCond = getValidatedValue<double>(physProp, "thermCond", "physProp");
        this->cp = getValidatedValue<double>(physProp, "cp", "physProp");
        //double gamma = getValidatedValue<double>(physProp, "gamma", "physProp");
        this->gamma = getValidatedValue<double>(physProp, "gamma", "physProp");

        // 初期条件
        this->initial = getValidatedValue<std::string>(config, "initial");

        // デバッグ用出力
        std::cout << "Mesh Format: " << meshFormat << '\n';
        std::cout << "Mesh File Name: " << meshFileName << '\n';
        std::cout << "Solver: " << solver << '\n';

    } catch (const std::runtime_error& e) {
        std::cerr << "Configuration Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    } catch (const YAML::Exception& e) {
        std::cerr << "YAML Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    initTimeIntegrationScheme(timeIntegration);
}

void solverConfig::initTimeIntegrationScheme(int timeIntegration){

    int ns;
    switch (timeIntegration) {
        case 1: // 1st Euler Explicit
            ns = 1;
            this->isImplicit = 0;
            this->nStage = ns;
            this->coef_N.resize(ns);
            this->coef_M.resize(ns);
            this->coef_Res.resize(ns);

            coef_N[0]   = 1.0;
            coef_M[0]   = 0.0;
            coef_Res[0] = 1.0;
            break;
        case 3: // 3rd TVD Runge-Kutta Explicit
            ns = 3;
            this->isImplicit = 0;
            this->nStage = ns;
            this->coef_N.resize(ns);
            this->coef_M.resize(ns);
            this->coef_Res.resize(ns);

            coef_N[0]   = 1.0 ; coef_N[1]   = 0.75 ; coef_N[2]   = 1.0/3.0;
            coef_M[0]   = 0.0 ; coef_M[1]   = 0.25 ; coef_M[2]   = 2.0/3.0;
            coef_Res[0] = 1.0 ; coef_Res[1] = 0.25 ; coef_Res[2] = 2.0/3.0;
            break;

        case 4: // 4th Runge-Kutta Explicit
            ns = 4;
            this->isImplicit = 0;
            this->nStage = ns;
            this->coef_DT_4thRunge.resize(ns);
            this->coef_Res_4thRunge.resize(ns);

            this->coef_DT_4thRunge[0] = 0.5;
            this->coef_DT_4thRunge[1] = 0.5;
            this->coef_DT_4thRunge[2] = 1.0;
            this->coef_DT_4thRunge[3] = -1e+30;

            this->coef_Res_4thRunge[0] = 1.0/6.0;
            this->coef_Res_4thRunge[1] = 1.0/3.0;
            this->coef_Res_4thRunge[2] = 1.0/3.0;
            this->coef_Res_4thRunge[3] = 1.0/6.0;

            break;

    }

}
