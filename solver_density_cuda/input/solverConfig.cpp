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
        std::cout << "'" << key << "'" << (parent.empty() ? "" : " in '" + parent + "'") << ": " << value << std::endl;
        return value;
    } catch (const YAML::BadConversion& e) {
        throw std::runtime_error("Key '" + key + "' in '" + parent + "' has an incorrect type.");
    }
}

template <typename T>
T getOptionalValidatedValue(const YAML::Node& node, const std::string& key, const T& default_value, const std::string& parent = "") {
    if (!node[key].IsDefined()) {
        return default_value;
    }

    try {
        T value = node[key].as<T>();
        std::cout << "'" << key << "'" << (parent.empty() ? "" : " in '" + parent + "'") << ": " << value << std::endl;
        return value;
    } catch (const YAML::BadConversion&) {
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

        // 離散化レイアウト (任意, 既定 "cell")。"cell": cell-centered、"node": node-centered
        // (中点双対 median-dual)。docs/discretization/ 参照。
        if (config["mesh"]["discretization"]) {
            this->discretization = config["mesh"]["discretization"].as<std::string>();
            if (this->discretization != "cell" && this->discretization != "node") {
                throw std::runtime_error("Key 'discretization' in 'mesh' must be 'cell' or 'node'.");
            }
        } else {
            this->discretization = "cell";
        }
        if (config["mesh"]["bndFirstOrder"]) {
            this->bndFirstOrder = config["mesh"]["bndFirstOrder"].as<int>();
        }
        if (config["mesh"]["axisCentroidShift"]) {
            this->axisCentroidShift = config["mesh"]["axisCentroidShift"].as<int>();
        }
        if (config["mesh"]["nodeWallDirichlet"]) {
            this->nodeWallDirichlet = config["mesh"]["nodeWallDirichlet"].as<int>();
        }
        if (config["mesh"]["gradLSQ"]) {
            this->gradLSQ = config["mesh"]["gradLSQ"].as<int>();
        }
        if (config["mesh"]["nodeMidpointFx"]) {
            this->nodeMidpointFx = config["mesh"]["nodeMidpointFx"].as<int>();
        }
        if (config["mesh"]["nodeWallViscGradFlux"]) {
            this->nodeWallViscGradFlux = config["mesh"]["nodeWallViscGradFlux"].as<int>();
        }
        if (config["mesh"]["nodeWallDistFloorCoef"]) {
            this->nodeWallDistFloorCoef = config["mesh"]["nodeWallDistFloorCoef"].as<flow_float>();
        }

        // solver関連
        this->solver = getValidatedValue<std::string>(config, "solver");

        // GPU設定
        this->gpu = getValidatedValue<int>(config, "gpu");

        // time関連
        this->unsteady = getValidatedValue<int>(config["time"], "unsteady", "time");
        this->dualTime = getValidatedValue<int>(config["time"], "dualTime", "time");

        auto last = config["time"]["last"];
        this->endTimeControl = getValidatedValue<int>(last, "control", "time.last");
        if (last["nStep"].IsDefined()) {
            throw std::runtime_error(
                "Key 'nStep' in 'time.last' is no longer supported. Rename it to 'nStepOuter'."
            );
        }
        this->nStepOuter = getValidatedValue<int>(last, "nStepOuter", "time.last");

        auto deltaT = config["time"]["deltaT"];
        this->dtControl = getValidatedValue<int>(deltaT, "control", "time.deltaT");
        this->dt = getValidatedValue<double>(deltaT, "dt", "time.deltaT");
        this->cfl = getValidatedValue<double>(deltaT, "cfl", "time.deltaT");
        this->cfl_pseudo = getValidatedValue<double>(deltaT, "cfl_pseudo", "time.deltaT");
        this->implicitRelax = getOptionalValidatedValue<double>(deltaT, "implicitRelax", 1.0, "time.deltaT");
        {
            double raw = getOptionalValidatedValue<double>(deltaT, "implicitRelaxSST", -1.0, "time.deltaT");
            this->implicitRelaxSST = (raw < 0.0) ? this->implicitRelax : (flow_float)raw;
        }
        this->blockDPLUR = getOptionalValidatedValue<int>(deltaT, "blockDPLUR", 0, "time.deltaT");
        // 多成分 TP 陰解法の化学種更新方式: 既定 0 (従来 segregated 点陰的・ビット不変)。
        // 1 で緩和整合 scalar-DPLUR (流れ block と同一緩和。plan thermophysics-species-implicit-coupling.md)。
        this->speciesImplicitCoupling = getOptionalValidatedValue<int>(deltaT, "speciesImplicitCoupling", 0, "time.deltaT");
        // 多成分 face 整合再構成: 既定 0 (mixed-order・ビット不変)。1 で Y を ρ と同じ再構成し thermo/species 流束整合。
        this->speciesFaceReconstruction = getOptionalValidatedValue<int>(deltaT, "speciesFaceReconstruction", 0, "time.deltaT");
        // 軸対称 near-axis 安定化係数 β_axis: 既定 0 (不変)。擬似時間スペクトル半径に λ_axis=β(|u_r|+c)A_planar を加える。
        this->axisTimestepBeta = getOptionalValidatedValue<flow_float>(deltaT, "axisTimestepBeta", 0.0, "time.deltaT");
        // block-DPLUR 線形 solve の内部精度: 既定 0 (float・従来高速)。1 で double 化 (軸対称近軸の根治用)。
        this->implicitSolvePrecision = getOptionalValidatedValue<int>(deltaT, "implicitSolvePrecision", 0, "time.deltaT");
        // NaN 検知診断モード: 既定 0 で従来挙動 (検査なし・ビット不変)。1 で毎ステップ終端に検査。
        this->detectNaN = getOptionalValidatedValue<int>(deltaT, "detectNaN", 0, "time.deltaT");
        // 低マッハ前処理 (Weiss-Smith): 既定 0 で従来挙動 (ビット不変)。
        this->lowMachPrecond = getOptionalValidatedValue<int>(deltaT, "lowMachPrecond", 0, "time.deltaT");
        this->precondEps = getOptionalValidatedValue<double>(deltaT, "precondEps", 0.15, "time.deltaT");
        // Thornber 型低マッハ再構成補正: 既定 0 で従来挙動 (ビット不変)。lowMachPrecond と直交・併用可。
        this->lowMachThornber = getOptionalValidatedValue<int>(deltaT, "lowMachThornber", 0, "time.deltaT");
        this->dt_max = getValidatedValue<double>(deltaT, "dt_max", "time.deltaT");
        this->dt_min = getValidatedValue<double>(deltaT, "dt_min", "time.deltaT");

        this->outStepInterval = getValidatedValue<int>(config["time"], "outStepInterval", "time");
        this->outStepStart = getValidatedValue<int>(config["time"], "outStepStart", "time");
        this->timeIntegration = getValidatedValue<int>(config["time"], "timeIntegration", "time");
        if (config["time"]["nInnerLoop"].IsDefined()) {
            throw std::runtime_error(
                "Key 'nInnerLoop' in 'time' is no longer supported. Rename it to 'nStepInner'."
            );
        }
        if (config["time"]["dualTime_InnerLoop"].IsDefined()) {
            throw std::runtime_error(
                "Key 'dualTime_InnerLoop' in 'time' is no longer supported. Rename it to 'nStepInner'."
            );
        }
        this->nStepInner = getOptionalValidatedValue<int>(config["time"], "nStepInner", 1, "time");
        // dual-time (非定常陰解法): 1 物理ステップあたりの擬似時間サブ反復数と BDF 次数。
        this->nSubIterDualTime = getOptionalValidatedValue<int>(config["time"], "nSubIterDualTime", 20, "time");
        this->bdfOrder = getOptionalValidatedValue<int>(config["time"], "bdfOrder", 2, "time");

        // 空間設定
        auto space = config["space"];
        this->convMethod = getValidatedValue<int>(space, "convMethod", "space");
        this->limiter = getValidatedValue<int>(space, "limiter", "space");
        if (this->limiter != 0 && this->limiter != 1 && this->limiter != 2 && this->limiter != -1) {
            throw std::runtime_error("Key 'limiter' in 'space' must be one of 0, 1, 2, or -1.");
        }
        // free-stream 保存用の基準静圧 (既定 0.0 = 従来挙動・ビット不変)
        this->pRef = getOptionalValidatedValue<double>(space, "pRef", 0.0, "space");

        // turbulence model
        auto turb = config["turbulence"];
        this->LESorRANS = getValidatedValue<int>(turb, "LESorRANS", "turbulence");
        this->LESmodel = getValidatedValue<int>(turb, "LESmodel", "turbulence");
        this->RANSmodel = getOptionalValidatedValue<int>(turb, "RANSmodel", 0, "turbulence");
        this->scalarDiffusion = getOptionalValidatedValue<int>(turb, "scalarDiffusion", 1, "turbulence");
        this->dilatationCorrection = getOptionalValidatedValue<int>(turb, "dilatationCorrection", 2, "turbulence");

        if (this->dilatationCorrection < 0 || this->dilatationCorrection > 2) {
            throw std::runtime_error("Key 'dilatationCorrection' in 'turbulence' must be one of 0, 1, or 2.");
        }

        this->katoLaunder = getOptionalValidatedValue<int>(turb, "katoLaunder", 0, "turbulence");

        if (this->katoLaunder < 0 || this->katoLaunder > 1) {
            throw std::runtime_error("Key 'katoLaunder' in 'turbulence' must be 0 or 1.");
        }

        // 乱流プラントル数 Pr_t (乱流熱伝導 k_t = cp*mu_t/Pr_t)。既定 0.9。
        this->turbulentPrandtl = getOptionalValidatedValue<flow_float>(turb, "turbulentPrandtl", 0.85, "turbulence");
        if (this->turbulentPrandtl <= 0.0) {
            throw std::runtime_error("Key 'turbulentPrandtl' in 'turbulence' must be positive.");
        }

        if (this->LESorRANS < 0 || this->LESorRANS > 2) {
            throw std::runtime_error("Key 'LESorRANS' in 'turbulence' must be one of 0, 1, or 2.");
        }
        if (this->LESorRANS == 1 && this->LESmodel <= 0) {
            throw std::runtime_error("Key 'LESmodel' in 'turbulence' must be positive when LESorRANS == 1.");
        }
        if (this->LESorRANS == 2 && this->RANSmodel <= 0) {
            throw std::runtime_error("Key 'RANSmodel' in 'turbulence' must be set when LESorRANS == 2.");
        }
        if (this->RANSmodel < 0 || this->RANSmodel > 1) {
            throw std::runtime_error("Key 'RANSmodel' in 'turbulence' must be one of 0 or 1.");
        }

        // 物理プロパティ
        auto physProp = config["physProp"];
        this->isCompressible = getValidatedValue<int>(physProp, "isCompressible", "physProp");
        if (physProp["isAxisymmetric"]) {
            this->isAxisymmetric = physProp["isAxisymmetric"].as<int>();
        } else {
            this->isAxisymmetric = 0;
        }
        this->thermalMethod = getValidatedValue<int>(physProp, "thermalMethod", "physProp");
        this->viscMethod = getValidatedValue<int>(physProp, "viscMethod", "physProp");
        this->ro = getValidatedValue<double>(physProp, "ro", "physProp");
        this->visc = getValidatedValue<double>(physProp, "visc", "physProp");
        this->thermCond = getValidatedValue<double>(physProp, "thermCond", "physProp");
        this->cp = getValidatedValue<double>(physProp, "cp", "physProp");
        //double gamma = getValidatedValue<double>(physProp, "gamma", "physProp");
        this->gamma = getValidatedValue<double>(physProp, "gamma", "physProp");

        // 多成分 thermally-perfect gas 設定 (任意, thermalMethod==2 で使用)
        if (physProp["species"]) {
            this->speciesNames.clear();
            for (const auto& sn : physProp["species"]) this->speciesNames.push_back(sn.as<std::string>());
            this->nSpecies = static_cast<int>(this->speciesNames.size());
        }
        if (physProp["speciesDBFile"])          this->speciesDBFile = physProp["speciesDBFile"].as<std::string>();
        if (physProp["speciesDiffusionMethod"]) this->speciesDiffusionMethod = physProp["speciesDiffusionMethod"].as<int>();
        if (physProp["thermoHrefTemp"])         this->thermoHrefTemp = physProp["thermoHrefTemp"].as<double>();
        if (physProp["Sc"])                     this->Sc = physProp["Sc"].as<double>();
        if (physProp["Sc_t"])                   this->Sc_t = physProp["Sc_t"].as<double>();
        // 乱流シュミット数は turbulence.turbulentSchmidt でも設定可 (turbulentPrandtl と同じ場所)。
        // 後方互換: physProp.Sc_t も有効。両方あれば turbulence 側を優先する。
        if (config["turbulence"] && config["turbulence"]["turbulentSchmidt"]) {
            this->Sc_t = config["turbulence"]["turbulentSchmidt"].as<flow_float>();
            std::cout << "'turbulentSchmidt' in 'turbulence': " << this->Sc_t << std::endl;
        }
        if (this->Sc_t <= 0.0) {
            throw std::runtime_error("Turbulent Schmidt number (Sc_t / turbulentSchmidt) must be positive.");
        }
        if (this->thermalMethod == 2 && this->speciesNames.empty()) {
            // 既定: 単成分 N2 (CEA 検証用)
            this->speciesNames.push_back("N2");
            this->nSpecies = 1;
        }

        // 非平衡凝縮 (任意セクション)。docs/condensation/ 参照。
        // Phase 1 は受動スカラー輸送のみ (核生成/成長係数はまだ読まない)。
        if (config["condensation"]) {
            auto cond = config["condensation"];
            this->condensation = getOptionalValidatedValue<int>(cond, "condensation", 0, "condensation");
            this->nCondSpecies = getOptionalValidatedValue<int>(cond, "nCondSpecies", 0, "condensation");
            this->condModel = getOptionalValidatedValue<int>(cond, "condModel", 0, "condensation");
            this->condGasSpecies = getOptionalValidatedValue<int>(cond, "condGasSpecies", -1, "condensation");
            this->condKantrowitz = getOptionalValidatedValue<int>(cond, "condKantrowitz", 0, "condensation");
            this->condGrowthModel = getOptionalValidatedValue<int>(cond, "condGrowthModel", 0, "condensation");
            this->condGyarmathyC = getOptionalValidatedValue<double>(cond, "condGyarmathyC", 3.18, "condensation");
            this->condTwoTemp = getOptionalValidatedValue<int>(cond, "condTwoTemp", 0, "condensation");
        }
        if (this->condensation != 0 && this->condensation != 1) {
            throw std::runtime_error("Key 'condensation' in 'condensation' must be 0 or 1.");
        }
        if (this->condensation == 1 && this->nCondSpecies < 1) {
            throw std::runtime_error("Key 'nCondSpecies' in 'condensation' must be >= 1 when condensation == 1.");
        }
        if (this->condensation == 0) {
            this->nCondSpecies = 0; // off のときは登録しない
        }

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

        case 11: // implicit pseudo-time defect-correction (block DPLUR)
            this->isImplicit = 1;
            this->nStage = 1;
            this->coef_N.clear();
            this->coef_M.clear();
            this->coef_Res.clear();
            this->coef_DT_4thRunge.clear();
            this->coef_Res_4thRunge.clear();
            // blockDPLUR==1: 5×5 block (LU-SGS A±)、blockDPLUR==0: scalar 対角 (スペクトル半径) 版。
            // どちらも古典 DPLUR 制御フロー (残差固定 + nStepInner sweep + 単一 commit) に対応。
            if (this->blockDPLUR != 0 && this->blockDPLUR != 1) {
                throw std::runtime_error(
                    "timeIntegration==11: blockDPLUR must be 0 (scalar-diagonal) or 1 (5x5 block).");
            }
            break;

        default:
            throw std::runtime_error("Unsupported timeIntegration: " + std::to_string(timeIntegration));

    }

    // implicitSolvePrecision (block-DPLUR 線形 solve の double 化) のサポート範囲チェック。
    // 未対応の組み合わせはフラグが黙って無視され「効いていない」事故になるため、明示エラーで停止する。
    // 実装済みは block 経路 (blockDPLUR==1) かつ非 precond (lowMachPrecond 0/1) の implicit (timeIntegration==11) のみ。
    // 詳細: .github/plans/precision-mixed-axisym.md §13。
    if (this->implicitSolvePrecision != 0 && this->implicitSolvePrecision != 1) {
        throw std::runtime_error(
            "implicitSolvePrecision must be 0 (float) or 1 (double).");
    }
    if (this->implicitSolvePrecision == 1) {
        if (timeIntegration != 11) {
            throw std::runtime_error(
                "implicitSolvePrecision==1 (double solve) requires timeIntegration==11 "
                "(implicit block DPLUR). It has no effect for explicit schemes.");
        }
        if (this->blockDPLUR != 1) {
            throw std::runtime_error(
                "implicitSolvePrecision==1 (double solve) is only supported with blockDPLUR==1 "
                "(5x5 block DPLUR); scalar-diagonal DPLUR (blockDPLUR==0) is not supported.");
        }
        if (this->lowMachPrecond == 2) {
            throw std::runtime_error(
                "implicitSolvePrecision==1 (double solve) is not supported with lowMachPrecond==2 "
                "(full preconditioned block DPLUR). Use lowMachPrecond 0 or 1.");
        }
    }

}

int solverConfig::mainLoopCount() const
{
    return std::max(1, this->nStepOuter);
}

int solverConfig::perStepIterationCount() const
{
    if (this->isImplicit == 1) {
        return std::max(1, this->nStepInner);
    }

    return std::max(1, this->nStage);
}

const char* solverConfig::perStepIterationLabel() const
{
    if (this->isImplicit == 1) {
        return "Inner";
    }

    return "Stage";
}
